      module mfi_hs_mod
      include 'SIZE'
c
c     mfi restart: crystal fan-in/fan-out handshake index (CSR) shared by
c     mfi_redist_plan / mfi_redist_round / mfi_redist_round_rma.
c     Style-1 (kept static common, wrapped in a module): small and constant in
c     problem size (~56 KB/rank), so no allocatable backing / init.
c
      parameter(lbrst_max=1024) ! send-side bound: batch elems (nb<=lbrst)
c     handshake recv bound: a dest gets <= its whole local field per batch, so
c     #contributing sources <= recvcnt <= nelt_hr0 <= lelt (independent of np).
      parameter(lhs_mx=lelt)
      integer kv,ord,ioff,dstlist,cnt,boff,it,ndest
      common /mfi_hs/ kv(2,lbrst_max),ord(lbrst_max),ioff(lbrst_max+1),
     $               dstlist(lbrst_max),cnt(lbrst_max),boff(lbrst_max),
     $               it(3,lhs_mx),ndest
      end module mfi_hs_mod
