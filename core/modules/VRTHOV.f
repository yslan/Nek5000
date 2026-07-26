      module vrthov_mod

      implicit none

      real*4,  allocatable, target :: cb_vrthov(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nvrthov
c        ic.f (mfi_gets/mfi_getv: w2(lrbs), lrbs=20*lx1*ly1*lz1*lelt)
         integer sz_w2

         sz_w2 = 20*lx1*ly1*lz1*lelt

         nvrthov = sz_w2

         allocate(cb_vrthov(nvrthov), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vrthov$',ierr)

c        cb_vrthov_i (the pre-fusion CR receive tuple `vi`) is gone: the fused
c        reader receives into wk (/scrns/), not a dedicated integer buffer.

      end subroutine init
      end module vrthov_mod
