"""NEMO domain, coordinate, time, and wet-mask data structures.

The rest of the assimilation code talks to these compact classes rather than
raw NetCDF files. To couple pyPDAF to another complex model, provide an
equivalent domain object with local dimensions, coordinates, a wet/active-point
mask, and the time/domain attributes needed by the writer.
"""
import dataclasses
import datetime
import functools
import math

import mpi4py.MPI
import numpy as np
import numpy.typing as npt

import log
import io_nemo
import parallel

@dataclasses.dataclass(frozen=True)
class LocalDomain:
    """The class stores NEMO grid and domain information for a single process.

    Attributes
    ----------
    ni_p : int
        The number of grid points in the i-direction for the local sub-domain
    nj_p : int
        The number of grid points in the j-direction for the local sub-domain
    nk_p : int
        The number of grid points in the z-direction for the local sub-domain
    i0 : int
        The starting index in the i-direction for the local sub-domain
    j0 : int
        The starting index in the j-direction for the local sub-domain
    halo0 : int
        The size of the halo region at the start of the local sub-domain
    halo1 : int
        The size of the halo region at the end of the local sub-domain
    """
    ni_p: int
    nj_p: int
    nk_p: int
    i0: int
    j0: int
    halo0: npt.NDArray
    halo1: npt.NDArray

    @functools.cached_property
    def dim_2d_p(self) -> int:
        """The number of grid points for local 2D state."""
        return self.ni_p * self.nj_p

    @functools.cached_property
    def dim_3d_p(self) -> int:
        """The number of grid points for local 3D state."""
        return self.ni_p * self.nj_p * self.nk_p

    def local_slices(self) -> dict[str, slice]:
        """Return the slices for the local sub-domain in the global domain."""
        return {
            "x": slice(self.i0 - 1, self.i0 - 1 + self.ni_p),
            "y": slice(self.j0 - 1, self.j0 - 1 + self.nj_p),
        }

@dataclasses.dataclass(frozen=True)
class GridCoordinates:
    """The class stores NEMO grid coordinates.

    Attributes
    ----------
    nav_lev : npt.NDArray
        The vertical grid levels
    nav_lon : npt.NDArray
        The longitudes of the grid points
    nav_lat : npt.NDArray
        The latitudes of the grid points
    glamt : npt.NDArray
        The longitudes of the grid points on t-points
    glamu : npt.NDArray
        The longitudes of the grid points on u-points
    glamv : npt.NDArray
        The longitudes of the grid points on v-points
    gphit : npt.NDArray
        The latitudes of the grid points on t-points
    gphiu : npt.NDArray
        The latitudes of the grid points on u-points
    gphiv : npt.NDArray
        The latitudes of the grid points on v-points
    """
    nav_lev: npt.NDArray
    nav_lon: npt.NDArray
    nav_lat: npt.NDArray
    glamt: npt.NDArray
    glamu: npt.NDArray
    glamv: npt.NDArray
    gphit: npt.NDArray
    gphiu: npt.NDArray
    gphiv: npt.NDArray


@dataclasses.dataclass
class WetGrid:
    """The class stores NEMO wet grid information.

    Attributes
    ----------
    tmask : npt.NDArray[np.bool_]
        The wet mask for the grid points.
        True indicates a wet point, False indicates a land point.
    wet_pts : npt.NDArray[np.int64]
        - 0 Global y index,
        - 1 Global x index,
        - 2 number wet layers at given latlon
        - 3 index in 2d grid box
        - 4 starting index in all wet points for vertical column
        - 5 local y index in subdomain
        - 6 local x index in subdomain
    idx_wet_2d : npt.NDArray[np.int64]
        The indices of the wet grid points in the 2D grid
    idx_nwet_3d : npt.NDArray[np.int64]
        The indices of the wet grid points in the 3D grid
    nlev_wet_2d : npt.NDArray[np.int64]
        The number of wet layers for each grid point in the 2D grid
    """
    tmask: npt.NDArray[np.bool_]
    wet_pts: npt.NDArray[np.int64]
    idx_wet_2d: npt.NDArray[np.int64]
    idx_nwet_3d: npt.NDArray[np.int64]
    nlev_wet_2d: npt.NDArray[np.int64]

    def __init__(self, domvars: dict[str, np.ndarray], local_dom: LocalDomain):
        self.tmask = self.get_tmask(domvars, local_dom.nk_p)
        self.wet_pts = self.get_wet_pts(local_dom)

    def get_tmask(self, domvars:dict[str, np.ndarray], nk_p:int) -> npt.NDArray[np.bool_]:
        """Return the wet mask with dimensions ``(nk, ny, nx)``."""
        # create the wet mask
        k_top = domvars["top_level"]
        k_bot = domvars["bottom_level"]
        k = np.arange(nk_p)[:, None, None]
        return (k_top[None, :, :] != 0) & (k >= k_top[None, :, :] - 1) & (k < k_bot[None, :, :])

    def get_wet_pts(self, local_dom: LocalDomain):
        """Return the wet grid points information."""
        if self.nwet2d > 0:
            j, i = np.nonzero(self.tmask[0])
            wet_pts = np.zeros((7, self.nwet2d), dtype=int)
            wet_pts[0] = j + local_dom.j0
            wet_pts[1] = i + local_dom.i0
            wet_pts[2] = self.tmask[:, j, i].sum(axis=0)
            wet_pts[3] = np.ravel_multi_index((j, i), (local_dom.nj_p, local_dom.ni_p))
            wet_pts[4, 0] = 0
            wet_pts[4, 1:] = np.cumsum(wet_pts[2, :-1])
            wet_pts[5] = j
            wet_pts[6] = i
        else:
            wet_pts = np.zeros((7, 0), dtype=int)
        return wet_pts

    @functools.cached_property
    def nwet2d(self) -> int:
        """Number of surface wet grid points"""
        return int(np.count_nonzero(self.tmask[0]))

    @functools.cached_property
    def nwet3d(self) -> int:
        """The number of wet grid points in the 3D grid."""
        return int(np.count_nonzero(self.tmask))

    @functools.cached_property
    def sdim3d(self) -> int:
        """The 3D dimension of the state vector."""
        return self.nwet3d

    @functools.cached_property
    def sdim2d(self) -> int:
        """The 2D dimension of the state vector."""
        return self.nwet2d


@dataclasses.dataclass(frozen=True)
class SpatialTime:
    """The class stores the dimension of the NEMO domain.

    Attributes
    ----------
    time_counter : float
        The time counter from the restart file
    ndastp : float
        The NEMO time string specified in the NEMO namelist namrun
    jpiglo : int
        The global NEMO grid dimension in the i-direction
    jpjglo : int
        The global NEMO grid dimension in the j-direction
    jpk : int
        The global NEMO grid dimension in the k-direction
    numcat : int
        The number of categories in the sea ice model
    nn_time0 : int
        The initial time of day in hhmm
    """
    time_counter: npt.NDArray
    ndastp: float
    jpiglo: int
    jpjglo: int
    jpk: int
    numcat: int
    nn_time0: int


class NemoDomain:
    """The class stores NEMO grid and domain information.

    The object is almost always created by the input reader from NEMO files
    using the class method ``from_input``.

    One can adapt this object for other models by creating equivalent domain,
    coordinate, active-cell, and time information. The key contract is that
    ``wet_grid`` identifies which local grid cells are part of the state vector.
    """

    def __init__(self, local_dom: LocalDomain, coords: GridCoordinates,
                 spatial_time: SpatialTime, wet_grid: WetGrid,
                 pe: parallel.Parallel) -> None:
        self.local_dom = local_dom
        self.coords = coords
        self.spatial_time = spatial_time
        self.wet_grid = wet_grid
        self.pe = pe
        self.bkg_ens = np.zeros((1))

        if self.wet_grid.nwet2d == 0:
            msg = f"No valid local domains, PE={self.pe.mype_filter:3d}"
            log.warning(msg, ranks=self.pe.mype_filter)
        self._log_grid_counts()

    @classmethod
    def from_input(cls, numcat:int, reader:io_nemo.Reader,
                   pe:parallel.Parallel) -> "NemoDomain":
        """
        Create a NemoDomain instance from restart and NEMO domain data .

        Parameters
        ----------
        numcat : int
            The number of categories in the sea ice model.
            It can be an arbitrary value if sea ice model is not used, or
            sea ice state vector is not included in the state vector.
        reader : io_pdaf.Reader
            The reader for accessing input data
        pe : parallel.Parallelisation
            The parallelisation object

        Returns
        -------
        NemoDomain
            The created NemoDomain instance
        """
        # reading restart file for process domain information
        attnames = ["DOMAIN_size_local", "DOMAIN_position_first",
                    "DOMAIN_halo_size_start", "DOMAIN_halo_size_end",
                    "DOMAIN_number_total"]
        vnames = ["ndastp", "time_counter", "nav_lev", "nav_lon", "nav_lat",
                  "ntime", 'adatrj']
        rst_vars, rst_attrs = reader.read_local_domain_coords(vnames, attnames)
        assert rst_attrs['DOMAIN_number_total'] == pe.npes_filter, \
                (f"Number of domains in restart file ({rst_attrs['DOMAIN_number_total']}) "
                f"does not match number of processes ({pe.npes_filter})")
        vnames = ["glamt", "glamu", "glamv", "gphit", "gphiu", "gphiv",
                  "bottom_level", "top_level"]

        local_dom = LocalDomain(
            ni_p=int(rst_attrs["DOMAIN_size_local"][0]),
            nj_p=int(rst_attrs["DOMAIN_size_local"][1]),
            nk_p=len(rst_vars['nav_lev']),
            i0=int(rst_attrs["DOMAIN_position_first"][0]),
            j0=int(rst_attrs["DOMAIN_position_first"][1]),
            halo0=rst_attrs["DOMAIN_halo_size_start"],
            halo1=rst_attrs["DOMAIN_halo_size_end"]
        )
        # read global domain config file for grid coordinates and wet grid information
        dims, domvars = reader.read_global_domain(vnames, local_dom.local_slices())
        coords = GridCoordinates(
            nav_lev=rst_vars["nav_lev"],
            nav_lon=rst_vars["nav_lon"],
            nav_lat=rst_vars["nav_lat"],
            glamt=domvars["glamt"],
            glamu=domvars["glamu"],
            glamv=domvars["glamv"],
            gphit=domvars["gphit"],
            gphiu=domvars["gphiu"],
            gphiv=domvars["gphiv"]
        )

        # set up the global domain and time information
        hour, minute = divmod(int(np.rint(rst_vars['ntime'])), 100)
        result = datetime.datetime(2000, 1, 1, hour, minute) + \
                datetime.timedelta(seconds=(float(rst_vars['adatrj']) -
                                            math.floor(float(rst_vars['adatrj']))
                                            ) * 86_400)
        spatial_time = SpatialTime(
            time_counter=rst_vars["time_counter"],
            ndastp=float(rst_vars["ndastp"]),
            jpiglo=dims['x'],
            jpjglo=dims['y'],
            jpk=dims['z'],
            numcat=numcat,
            nn_time0=result.hour * 100 + result.minute
        )

        # compute wet grid information for the local domain
        wet_grid = WetGrid(domvars, local_dom)

        return cls(local_dom, coords, spatial_time, wet_grid, pe)

    def _log_grid_counts(self):
        if self.pe.npes_filter > 1:
            nwet2d_g = self.pe.comm_filter.reduce(self.wet_grid.nwet2d, op=mpi4py.MPI.SUM, root=0)
            nwet3d_g = self.pe.comm_filter.reduce(self.wet_grid.nwet3d, op=mpi4py.MPI.SUM, root=0)
            if self.pe.mype_filter == 0:
                msg = f"Number of global wet surface points   {nwet2d_g:11d}"
                log.info(msg, ranks=0, screen=1)
                msg = f"Number of global 3D wet points        {nwet3d_g:11d}"
                log.info(msg, ranks=0, screen=1)
                msg = ("Number of global 2D wet points * nlayers        "
                       f"{nwet2d_g * self.local_dom.nk_p:11d}")
                log.info(msg, ranks=0, screen=1)

        msg = (f"PE {self.pe.mype_filter:4d}  Number of wet surface points   "
               f"{self.wet_grid.nwet2d:11d}")
        log.info(msg, screen=2)
        msg = (f"PE {self.pe.mype_filter:4d}  Number of 3D wet points        "
               f"{self.wet_grid.nwet3d:11d}")
        log.info(msg, screen=2)
        msg = (f"PE {self.pe.mype_filter:4d}  2D wet points * nlayers        "
               f"{self.wet_grid.nwet2d * self.local_dom.nk_p:11d}")
        log.info(msg, screen=2)
