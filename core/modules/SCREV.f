      module screv_mod

      implicit none

      real, allocatable, target :: cb_screv(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrev

c        dominant family: two LX1*LY1*LZ1*LELT/LELV "planes"
c        (drive2.f SII/SIII, eigsolv.f X1/Y1, subs1.f dpc/p1,
c        subs2.f ei2/ei3, navier1.f H1/H2 & RKZ1/RKZ2,
c        navier6.f x(2*ltotd)); dominates every other SCREV site
         nscrev = 2*lx1*ly1*lz1*lelt

         allocate(cb_screv(nscrev), stat=ierr)

      end subroutine init
      end module screv_mod
