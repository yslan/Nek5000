      module vrthov_mod
      include 'SIZE'            ! = use size_mod: lelt for the /mfi_hs/ bound
      implicit none

c     restart buffers (dedicated, problem-size independent):
c       cb_vrthov = w2 : file-order read buffer for one batch
c       cb_wk     = wk : CR/RMA redistribution buffer (also the RMA window)
c     both sized to (at least) lrst_targ words (= lrst_mb MB) in reserve()
      real*4,  allocatable, target :: cb_vrthov(:)
      real*4,  allocatable, target :: cb_wk(:)

c     per-element floor (worst single element, real*4 words): FP64 + mapab
c     source-order margin (nxr<=lx1+6). one element must fit either buffer.
      integer, parameter :: lelem_mx = 2*(lx1+6)*(ly1+6)*(lz1+6)

c     fixed byte budget per buffer (knob lrst_mb, in MB); lrst_targ = words
      integer, parameter :: lrst_mb   = 16
      integer, parameter :: lrst_targ = lrst_mb*1024*1024/4

c     lidst: RMA idstage (round_rma) length = max read-batch; mfi caps
c            lbrst <= lidst ('d' check). Sibling of lhs_mx (both = lelt).
      integer, parameter :: lidst = lelt

c     redist handshake work arrays, bound by lhs_mx=lelt:
c       send-side: nb<=lelt elems/batch;
c       recv-side (it): #rows n1 <= recvcnt <= nelt_hr0 <= lelt (indep. of np)
c       actual usage << lelt (depends on #dests/#sources per rank)
      integer, parameter :: lhs_mx = lelt
      integer kv,ord,ioff,dstlist,cnt,boff,it,ndest
      common /mfi_hs/ kv(2,lhs_mx),ord(lhs_mx),ioff(lhs_mx+1),
     $               dstlist(lhs_mx),cnt(lhs_mx),boff(lhs_mx),
     $               it(3,lhs_mx),ndest

      contains

c     vrthov_reserve: ensure restart buffers have capacity >= max(lrst_targ,need)
c     grown only if `need` (one elem) exceeds the budget (rare, high order).
c     MUST NOT run while cb_wk backs a live window.
      subroutine vrthov_reserve(need)
         implicit none
         integer, optional :: need
         integer ns, ierr

         ns = lrst_targ                    ! default = 16MB budget
         if (present(need)) then
            if (need.gt.ns) ns = need      ! increase only if needed
         endif
         if (mod(ns,2).ne.0) ns = ns + 1   ! keep even (2*lwk)

         if (allocated(cb_vrthov)) then
            if (size(cb_vrthov).ge.ns) return ! big enough: no realloc
            deallocate(cb_vrthov)             ! grow -> both together
         endif
         if (allocated(cb_wk)) deallocate(cb_wk)

         allocate(cb_vrthov(ns), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_vrthov$',ns)
         allocate(cb_wk(ns), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_wk$',ns)
      end subroutine vrthov_reserve

c     escape hatch: release both buffers. free any MPI window over cb_wk first.
      subroutine vrthov_free
         implicit none
         if (allocated(cb_vrthov)) deallocate(cb_vrthov)
         if (allocated(cb_wk))     deallocate(cb_wk)
      end subroutine vrthov_free
      end module vrthov_mod
