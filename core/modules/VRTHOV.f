      module vrthov_mod
      include 'SIZE'            ! = use size_mod: lelt for the /mfi_hs/ bound
      implicit none

      real*4,  allocatable, target :: cb_vrthov(:)

c     mfi restart: crystal handshake work arrays, shared by
c     mfi_redist_plan / mfi_redist_round_cr / mfi_redist_round_rma.
c     All arrays bound by lhs_mx=lelt:
c       send-side: nb<=lelt elems/batch;
c       recv-side (it): #rows n1 <= recvcnt <= nelt_hr0 <= lelt (indep. of np)
c       actual usage << lelt (depends on #dests/#sources per rank)
      integer, parameter :: lhs_mx = lelt
      integer kv,ord,ioff,dstlist,cnt,boff,it,ndest
      common /mfi_hs/ kv(2,lhs_mx),ord(lhs_mx),ioff(lhs_mx+1),
     $               dstlist(lhs_mx),cnt(lhs_mx),boff(lhs_mx),
     $               it(3,lhs_mx),ndest

      contains

      subroutine init
         implicit none

         integer ierr, nvrthov
c        ic.f (mfi_gets/mfi_getv: w2(lrbs), lrbs=20*lx1*ly1*lz1*lelt)
         integer sz_w2

         sz_w2 = 20*lx1*ly1*lz1*lelt

         nvrthov = sz_w2

         allocate(cb_vrthov(nvrthov), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vrthov$',ierr)

      end subroutine init
      end module vrthov_mod
