      module vrthov_mod

      implicit none

      real*4,  allocatable, target :: cb_vrthov(:)
      integer, allocatable, target :: cb_vrthov_i(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nvrthov, nvrthov_i
c        ic.f (mfi_gets/mfi_getv: w2(lrbs), lrbs=20*lx1*ly1*lz1*lelt)
         integer sz_w2
c        ic.f (mfi_gets/mfi_getv: vi(2+lrbs_loc,lelt))
         integer sz_vi

         sz_w2 = 20*lx1*ly1*lz1*lelt
         sz_vi = (2+20*lx1*ly1*lz1)*lelt

         nvrthov   = sz_w2
         nvrthov_i = sz_vi

         allocate(cb_vrthov(nvrthov), stat=ierr)
         allocate(cb_vrthov_i(nvrthov_i), stat=ierr)

      end subroutine init
      end module vrthov_mod
