      module scrvh_mod

      implicit none

      real, allocatable, target :: cb_scrvh(:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, nscrvh
         integer sz_lelt2   ! conduct.f/navier4.f (h1(lt),h2(lt))
         integer sz_lelv3   ! hsmg.f (mg_set_h1/h2: h1,h2,h2inv)
         integer sz_lelv4   ! convect.f (char_conv: bmsk,bdwt,bmst,u1)

         sz_lelt2 = 2*lx1*ly1*lz1*lelt
         sz_lelv3 = 3*lx1*ly1*lz1*lelv
         sz_lelv4 = 4*lx1*ly1*lz1*lelv

         nscrvh = max(sz_lelt2, sz_lelv3, sz_lelv4)

         allocate(cb_scrvh(nscrvh), stat=ierr)

      end subroutine init
      end module scrvh_mod
