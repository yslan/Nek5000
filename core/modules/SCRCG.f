      module scrcg_mod

      implicit none

      real, allocatable, target :: cb_scrcg(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrcg

c        hmholtz.f (hmh_flex_cg: r,z,p,w, each LX1*LY1*LZ1*LELT);
c        dominates every other SCRCG site in the source
         nscrcg = 4*lx1*ly1*lz1*lelt

         allocate(cb_scrcg(nscrcg), stat=ierr)

      end subroutine init
      end module scrcg_mod
