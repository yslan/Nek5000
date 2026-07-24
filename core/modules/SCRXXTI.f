      module scrxxti_mod

      implicit none

      integer, allocatable, target :: cb_scrxxti(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, n
c        subs1.f (set_up_h1_crs_strs: ia,ja(ldim*ldim*lcr*lcr*lelv));
c        dominates navier8.f (set_up_h1_crs: ia,ja(lcr,lcr,lelv))
         integer sz
         integer lxc, lcr

         lxc = 2
         lcr = lxc**ldim

         sz = 2*ldim*ldim*lcr*lcr*lelv

         n = sz

         allocate(cb_scrxxti(n), stat=ierr)

      end subroutine init
      end module scrxxti_mod
