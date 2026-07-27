      module vrthov_mod
      include 'SIZE'            ! = use size_mod: lelt for the /mfi_hs/ bound
      implicit none

      real*4,  allocatable, target :: cb_vrthov(:)

c     mfi restart: crystal fan-in/fan-out handshake index (CSR), shared by
c     mfi_redist_plan / mfi_redist_round / mfi_redist_round_rma. Kept a static
c     common (small, constant in problem size ~56 KB/rank), consolidated into
c     this restart-scratch module instead of a separate MFI_HS.f.
      integer, parameter :: lbrst_max = 1024 ! send-side bound: batch elems
c     handshake recv bound: a dest gets <= its whole local field per batch, so
c     #contributing sources <= recvcnt <= nelt_hr0 <= lelt (independent of np).
      integer, parameter :: lhs_mx = lelt
      integer kv,ord,ioff,dstlist,cnt,boff,it,ndest
      common /mfi_hs/ kv(2,lbrst_max),ord(lbrst_max),ioff(lbrst_max+1),
     $               dstlist(lbrst_max),cnt(lbrst_max),boff(lbrst_max),
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

c        cb_vrthov_i (the pre-fusion CR receive tuple `vi`) is gone: the fused
c        reader receives into wk (/scrns/), not a dedicated integer buffer.

      end subroutine init
      end module vrthov_mod
