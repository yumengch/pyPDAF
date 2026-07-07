module pdafomi_c_type
use PDAFomi_obs_f, only: obs_f
use PDAFomi_obs_l, only: obs_l
implicit none

type(obs_f), allocatable, target :: thisobs(:)
type(obs_l), allocatable, target :: thisobs_l(:)
integer :: n_obs_omi

!$OMP THREADPRIVATE(thisobs_l)

end module pdafomi_c_type