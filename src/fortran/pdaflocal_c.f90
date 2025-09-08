MODULE pdaflocal_c
use iso_c_binding, only: c_int, c_double, c_bool
use PDAF
use pdaf_c_cb_interface

implicit none

contains
   SUBROUTINE c__PDAFlocal_set_indices(dim_l, map) bind(c)
      ! Dimension of local state vector
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Index array for mapping between local and global state vector
      INTEGER(c_int), DIMENSION(dim_l), INTENT(in) :: map


      call PDAFlocal_set_indices(dim_l, map)

   END SUBROUTINE c__PDAFlocal_set_indices

   SUBROUTINE c__PDAFlocal_set_increment_weights(dim_l, weights) bind(c)
      ! Dimension of local state vector
      INTEGER(c_int), INTENT(in) :: dim_l
      ! Weights array
      REAL(c_double), DIMENSION(dim_l), INTENT(in) :: weights


      call PDAFlocal_set_increment_weights(dim_l, weights)

   END SUBROUTINE c__PDAFlocal_set_increment_weights

   SUBROUTINE c__PDAFlocal_clear_increment_weights() bind(c)
      call PDAFlocal_clear_increment_weights()

   END SUBROUTINE c__PDAFlocal_clear_increment_weights
END MODULE pdaflocal_c
