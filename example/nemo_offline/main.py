"""Entry point for the offline pyPDAF/NEMO analysis workflow."""
import configparser

import log
import parallel
import da_sys

def main():
    """Run one offline PDAF analysis cycle from ``config.ini``.

    The script expects NEMO restart/domain files and ensemble directories to
    match the paths configured in ``config.ini``. It initialises PDAF/MPI,
    reads the forecast ensemble, assimilates configured observations, writes
    NEMO ASM files, and finalises MPI.
    """
    config = configparser.ConfigParser()
    config.read('config.ini')
    screen = config['pdaf'].getint('screen')
    pe = parallel.Parallel(dim_ens=config['pdaf'].getint('dim_ens'),
                         screen=screen)
    log.set_printing_process(pe, ranks=0, screen=screen)
    pdaf_sys = da_sys.DAsystem(pe, config)
    pdaf_sys.init_pdaf(screen)
    pdaf_sys.assimilate()
    pdaf_sys.finalise()
    pe.finalize_parallel()

if __name__ == '__main__':
    main()
