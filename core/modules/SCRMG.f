      module scrmg_mod

      implicit none

      real, allocatable, target :: cb_scrmg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrmg
         integer sz_lelt4   ! hmholtz.f/hsmg.f/navier6.f (4 planes)
         integer sz_avg     ! prepost.f (tax: LX1*LY1*LELT*LDIMT)

         sz_lelt4 = 4*lx1*ly1*lz1*lelt
         sz_avg   = lx1*ly1*lelt*ldimt

         nscrmg = max(sz_lelt4, sz_avg)

         allocate(cb_scrmg(nscrmg), stat=ierr)

      end subroutine init
      end module scrmg_mod
