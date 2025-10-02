c-----------------------------------------------------------------------
      subroutine set_up_hmg_crs
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP' ! ifield
      include 'REFINEMG'

      if (nhref.eq.0) return
      if (nhref.ne.1) call exitti('hmg only support one level$',nhref)

      ncut = hrefcuts(1)
      nblk = ncut**ldim
      call lim_chk(ncut,lcut,'ncut ','lcut ',' h_ref_mg ')
      if (ncut.le.1) call exitti('hmg invalid ncut$',ncut)
      if (nio.eq.0) write(*,*)'hmg crs setup',ncut,nblk

      call hmg_set_interp_mat(ncut,lxc,hmg_pc,hmg_pt)

      call refine_cbc_r2o(hmg_CBCo,hmg_nelv_o,ncut,ifield)

      call set_up_hmg_crs_matrix

      return
      end
c-----------------------------------------------------------------------
      subroutine hmg_set_interp_mat(ncut,nxc,pc,pt)
      implicit none
      include 'SIZE'

      integer nxc,ncut
      real pc(nxc*nxc,ncut) ! Interpolation matrix
      real pt(nxc*nxc,ncut)

      integer i,k
      real zh(nxc), dr, r0

      real z_out, wk, wk2
      common /qtmp0/ z_out(lx1),wk(lx1*lx1),wk2(lx1*lx1)

      integer ncut_save
      save ncut_save
      data ncut_save /0/

      if (ncut.le.1)
     $  call exitti('invalid ncut in hmg_set_interp_mat$',ncut)

      call lim_chk(nxc,lx1,'nxc  ','lx1  ',' hmg_intp ')

      if (ncut.ne.ncut_save) then

        call zwgll(zh,wk,nxc)

        dr = 2./ncut
        do k=1,ncut
          r0 = -1. + (k-1)*dr
          do i=1,nxc
            z_out(i) = r0 + dr*(zh(i)+1)/2.
          enddo
          call interp_mat(pc(1,k),z_out,nxc,zh,nxc,wk,wk2)
          call transpose (pt(1,k),nxc,pc(1,k),nxc)
        enddo

        ncut_save = ncut

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hmg_interp_o2r(ur,uo)
c     uo -> ur
      implicit none
      include 'SIZE'
      include 'REFINEMG'

      real ur(lxyzc,nelv)
      real uo(lxyzc,hmg_nelv_o)

      integer kcut,ic,jc,kc,e,el,er

      real wk
      common /qtmp0/ wk(lxc*lxc)

      kcut=ncut
      if (ldim.eq.2) kcut=1

      call rzero(ur,lxyzc*nelv)

      do e = hmg_nelv_o,1,-1

        el = 0
        do kc=1,kcut
        do jc=1,ncut
        do ic=1,ncut

          el = el + 1
          er = nblk*(e-1) + el

          call tensr3(ur(1,er),lxc,uo(1,e),lxc
     $               ,hmg_pc(1,ic),hmg_pt(1,jc),hmg_pt(1,kc),wk)

        enddo
        enddo
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine hmg_interp_r2o(uo,ur)
c     ur -> uo
      implicit none
      include 'SIZE'
      include 'REFINEMG'

      real ur(lxyzc,nelv)
      real uo(lxyzc,hmg_nelv_o)
      real u1(lxyzc)

      integer kcut,ic,jc,kc,e,el,er

      real wk
      common /qtmp0/ wk(lxc*lxc)

      kcut=ncut
      if (ldim.eq.2) kcut=1

      call rzero(uo,lxyzc*hmg_nelv_o)

      do e = hmg_nelv_o,1,-1

        el = 0
        do kc=1,kcut
        do jc=1,ncut
        do ic=1,ncut

          el = el + 1
          er = nblk*(e-1) + el

          call tensr3(u1,lxc,ur(1,er),lxc
     $               ,hmg_pt(1,ic),hmg_pc(1,jc),hmg_pc(1,kc),wk)

          call add2(uo(1,e),u1,lxyzc)

        enddo
        enddo
        enddo

      enddo

      call cmult(uo,1.0/nblk,lxyzc*hmg_nelv_o)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_hmg_crs_matrix
      implicit none
      include 'SIZE'
      include 'GEOM' ! ifvcor
      include 'REFINEMG'
      include 'INPUT'
      include 'PARALLEL' ! xxth
      include 'TSTEP' ! ifield

      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer null_space

      character*3 cb
      integer ia,ja
      common /scrhmgxxti/ ia(lxyzc*lxyzc*lelv), ja(lxyzc*lxyzc*lelv)
      real a
      common /scrhmgxxtr/ a(lxyzc*lxyzc*lelv)

      real h1,h2,w1,w2
      common /scrhmghx/ h1(lx1*ly1*lz1*lelv),h2(lx1*ly1*lz1*lelv)
      common /scrhmgxx/ w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelv)

      real z
      integer*8 ngv
      character*132 amgfile_c
      character*1   fname1(132)
      equivalence  (fname1,amgfile_c)

      integer nxc,nzc,ncr,ntot,nz,nfaces,n,na,la
      integer ie,ierr,iface,isolver
      integer lamgn,ltrunc,iglmax,iglmax_ms
      real*8 t0, dnekclock

      t0 = dnekclock()

      nfaces = 2*ldim
      nxc  = 2
      nzc  = 1
      if (if3d) nzc = nxc
      ncr  = nxc**ldim
      ntot = ncr * hmg_nelv_o

      na = ncr*ncr*hmg_nelv_o
      la = lxyzc*lxyzc*lelv
      call lim_chk(na,la,'a    ','27*lt',' hmg_crsm ')

      if(nio.eq.0) write(6,*) 'setup hmg coarse grid, nx_crs='
     $                        ,nxc,hmg_nelv_o,ifield

c
c     Set SEM_to_GLOB
c
      call set_vert(hmg_se_to_gcrs,ngv
     $             ,nxc,hmg_nelv_o,hmg_vertex_o,.true.)

c     Set mask
      z=0
      call rone(hmg_crs_mask,ntot)
      call rone(hmg_crs_cmlt,ntot)
      if (ifield.eq.1) then
         do ie=1,hmg_nelv_o
         do iface=1,nfaces
            cb=hmg_cbco(iface,ie)
            if (cb.eq.'o  '  .or.  cb.eq.'on '  .or.
     $          cb.eq.'O  '  .or.  cb.eq.'ON '  .or.  cb.eq.'MM '  .or.
     $          cb.eq.'mm '  .or.  cb.eq.'ms '  .or.  cb.eq.'MS ')
     $           call facev(hmg_crs_mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
         enddo
         enddo
      elseif (ifield.eq.ifldmhd) then   ! no ifmhd ?avo?
         do ie=1,hmg_nelv_o
         do iface=1,nfaces
            cb=hmg_cbco(iface,ie)
            if (cb.eq.'ndd'  .or.  cb.eq.'dnd'  .or.  cb.eq.'ddn')
     $          call facev(hmg_crs_mask,ie,iface,z,nxc,nxc,nzc)
         enddo
         enddo
      endif

c     Set global index of dirichlet nodes to zero; xxt will ignore them

      call fgslib_gs_setup(hmg_crs_gsh,hmg_se_to_gcrs,ntot,nekcomm,mp)
      call fgslib_gs_op   (hmg_crs_gsh,hmg_crs_mask,1,2,0)  !  "*"
      call fgslib_gs_op   (hmg_crs_gsh,hmg_crs_cmlt,1,1,0)  !  "+"
c     call fgslib_gs_free (hmg_crs_gsh)
      call set_jl_crs_mask(ntot,hmg_crs_mask,hmg_se_to_gcrs)

      call invcol1(hmg_crs_cmlt,ntot)

c     Setup local SEM-based Neumann operators (for now, just full...)

c      if (param(51).eq.1) then     ! old coarse grid
c         nxyz1=lx1*ly1*lz1
c         lda = 27*nxyz1*lelt
c         ldw =  7*nxyz1*lelt
c         call get_local_crs(a,lda,nxc,h1,h2,w,ldw)
c      else
c        NOTE: a(),h1,...,w2() must all be large enough
         n = lx1*ly1*lz1*nelv
         call rone (h1,n)
         call rzero(h2,n)
         call hmg_get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)
c      endif

      call set_mat_ij(ia,ja,ncr,hmg_nelv_o) ! (ia,ja) are 0-base
      null_space=0
      if (ifield.eq.1) then
         if (ifvcor)  null_space=1
      elseif (ifield.eq.ifldmhd) then
         if (ifbcor)  null_space=1
      endif

      nz=ncr*ncr*hmg_nelv_o
      isolver = param(40)

      call blank(fname1,132)
      lamgn = ltrunc(amgfile,len(amgfile))
      call chcopy(fname1,amgfile,lamgn)
      call chcopy(fname1(lamgn+1),char(0),1)

      ierr = 0
      call crs_setup(xxth(ifield),isolver,nekcomm,mp,ntot,
     $     hmg_se_to_gcrs,nz,ia,ja,a, null_space, crs_param,
     $     amgfile_c,ierr)
      ierr = iglmax(ierr,1)
      if (ifneknek) ierr = iglmax_ms(ierr,1)
      if (ierr.eq.1) then
         call exitt
      endif

      t0 = dnekclock()-t0
      if (nio.eq.0) then
         write(6,*) 'done :: setup hmg coarse grid ',t0, ' sec'
         write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hmg_get_local_crs_galerkin(a,ncl,nxc,h1,h2,w1,w2)
      implicit none
      include 'SIZE'
      include 'REFINEMG'

      real    a(ncl,ncl,1),h1(1),h2(1)
      real    w1(lx1*ly1*lz1,nelv),w2(lx1*ly1*lz1,nelv)

      real b
      common /ctmp1z2/ b(lx1*ly1*lz1,lblk,8)

      integer e,i,j,ec,er,ic,jc,kc,kcut,isd,imsh,ncl,nxc,nxyz
      real vlsc2

      call lim_chk(nblk,lblk,'nblk ','lblk ',' hmg_crsl ')

      kcut = ncut
      if (ldim.eq.2) kcut = 1

      isd  = 1
      imsh = 1

      nxyz = lx1*ly1*lz1

      call rzero(a,ncl*ncl*hmg_nelv_o)

      do j=1,ncl
         ec = 0
         do kc = 1,kcut
         do jc = 1,ncut
         do ic = 1,ncut
            ec = ec + 1
            call hmg_gen_crs_basis(b(1,ec,j),ncut,ic,jc,kc,j)
         enddo
         enddo
         enddo
      enddo

      do j = 1,ncl

         do e = 1,hmg_nelv_o
            ec = 0
            do kc = 1,kcut
            do jc = 1,ncut
            do ic = 1,ncut
               ec = ec + 1
               er = (e-1)*nblk + ec
               call copy(w1(1,er),b(1,ec,j),nxyz)
            enddo
            enddo
            enddo
         enddo

         call axhelm (w2,w1,h1,h2,imsh,isd)        ! A^e * bj

         do e = 1,hmg_nelv_o
            ec = 0
            do kc = 1,kcut
            do jc = 1,ncut
            do ic = 1,ncut
               ec = ec + 1
               er = (e-1)*nblk + ec
               do i = 1,ncl
                  a(i,j,e) = a(i,j,e) + vlsc2(b(1,ec,i),w2(1,er),nxyz)  ! bi^T * A^e * bj
               enddo
            enddo
            enddo
            enddo
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine hmg_gen_crs_basis(b,ncut,ic,jc,kc,iv) ! bi- tri-linear
      implicit none
      include 'SIZE'

      real b(lx1,ly1,lz1)
      real zh(lx1)
      real zr(lx1),zs(lx1),zt(lx1)

      integer i,p,q,r,ic,jc,kc,iv,ncut

      ! get (zr,zs,zt) in [0,1]^3
      call zwgll(zh,zs,lx1)
      do i=1,lx1
         zh(i) = 0.5*(1.0+zh(i)) / (1.0*ncut) ! [0, 1/ncut]
         zr(i) = (ic-1) / (1.0*ncut) + zh(i)
         zs(i) = (jc-1) / (1.0*ncut) + zh(i)
         zt(i) = (kc-1) / (1.0*ncut) + zh(i)
      enddo

      ! get linear basis in each direction
      do i=1,lx1
         zr(i) = 1.0 - zr(i)
         zs(i) = 1.0 - zs(i)
         zt(i) = 1.0 - zt(i)
      enddo
      if (mod(iv,2).eq.0) then
         do i=1,lx1
            zr(i) = 1.0 - zr(i)
         enddo
      endif
      if (iv.eq.3.or.iv.eq.4.or.iv.eq.7.or.iv.eq.8) then
         do i=1,lx1
            zs(i) = 1.0 - zs(i)
         enddo
      endif
      if (iv.gt.4) then
         do i=1,lx1
            zt(i) = 1.0 - zt(i)
         enddo
      endif

      if (ldim.eq.3) then
         do r=1,lx1
         do q=1,lx1
         do p=1,lx1
            b(p,q,r) = zr(p)*zs(q)*zt(r)
         enddo
         enddo
         enddo
      else
         do q=1,lx1
         do p=1,lx1
            b(p,q,1) = zr(p)*zs(q)
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
c     Interface subroutine starts here
c-----------------------------------------------------------------------
      subroutine semg_hmg_crs_solve(e,r)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! xxth
      include 'CTIMER'
      include 'TSTEP' ! ifield
      real e(1),r(1)

      if (icalld.eq.0) then ! timer info
         ncrsl=0
         tcrsl=0.0
      endif
      icalld = 1

      if (ifsync) call nekgsync()

      ncrsl  = ncrsl  + 1
      etime1=dnekclock()

      call crs_solve(xxth(ifield),e,r)

      tcrsl=tcrsl+dnekclock()-etime1

      return
      end
c-----------------------------------------------------------------------
      subroutine semg_hmg_rstr(r,ifdssum)
c     r = J^T r
      implicit none
      include 'SIZE'
      include 'REFINEMG'
      real r(1)
      real ro(lxyzc*lelt)
      logical ifdssum
      integer n

      n = lxyzc * hmg_nelv_o

      call hmg_interp_r2o(ro,r)
      call copy (r,ro,n)

      call col2(r,hmg_crs_cmlt,n)
      call fgslib_gs_op(hmg_crs_gsh,r,1,1,0)  !  "+"

      return
      end
c-----------------------------------------------------------------------
      subroutine semg_hmg_intp(w,e)
c     w = J e
      implicit none
      real w(1), e(1)

      call hmg_interp_o2r(w,e)

      return
      end
c-----------------------------------------------------------------------
      subroutine semg_hmg_mask(u)
      implicit none
      include 'SIZE'
      include 'REFINEMG'
      real u(1)
      integer n

      n = lxyzc * hmg_nelv_o
      call col2(u,hmg_crs_mask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine dbg_prt_crs_mat(a,ia,ja,se_to_gcrs,n,nel)
      implicit none
      include 'SIZE'
      include 'PARALLEL'

      real a(n,n,nel)
      integer ia(n,n,nel), ja(n,n,nel), ii,jj, i,j,n,e,nel
      integer*8 se_to_gcrs(1)

      integer er,egr,ego,ie_map_f2c,ie_map_c2f,ncut,nblk

      ncut = 2
      nblk = ncut**ldim

      do e=1,nel
         er = ie_map_c2f(e,nblk)
         egr = lglel(er)
         ego = ie_map_f2c(egr,nblk)
      do j=1,n
      do i=1,n
         ii = ia(i,j,e) + 1
         jj = ja(i,j,e) + 1
         write(*,*)'cm',nid,i,j,e,ego,'|',ii,jj,'|'
     $            ,se_to_gcrs(ii),se_to_gcrs(jj),a(i,j,e)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine test_mask(mask,n,nel,s6)
      implicit none
      include 'SIZE'
      character*6 s6
      real u(n),mask(1),avg
      integer n
      real m1,m2,m3,nn,glmin,glmax,glsum,nel

      call rone(u,n)
      m1 = glmin(u,n)
      m2 = glmax(u,n)
      m3 = glsum(u,n)
      nn = glsum(1.0*n,1)
      avg = m3 / nn
      if (nio.eq.0) write(*,*)'dbg msk1 ',s6,m1,m2,m3,avg,n,nn

c      call col2(u,mask,n)
      call h1mg_mask(u,mask,nel)

      m1 = glmin(u,n)
      m2 = glmax(u,n)
      m3 = glsum(u,n)
      nn = glsum(1.0*n,1)
      avg = m3 / nn
      if (nio.eq.0) write(*,*)'dbg msk2 ',s6,m1,m2,m3,avg,n,nn

      return
      end
c-----------------------------------------------------------------------
      subroutine test_mask2(s6)
      implicit none
      include 'SIZE'
      include 'REFINEMG'
      character*6 s6
      real u(lxyzc*lelv),avg
      integer n
      real m1,m2,m3,nn,glmin,glmax,glsum,nel

      n = lxyzc * hmg_nelv_o
      call rone(u,n)

      m1 = glmin(u,n)
      m2 = glmax(u,n)
      m3 = glsum(u,n)
      nn = glsum(1.0*n,1)
      avg = m3 / nn
      if (nio.eq.0) write(*,*)'dbg msk1 ',s6,m1,m2,m3,avg,n,nn

      call semg_hmg_mask(u)

      m1 = glmin(u,n)
      m2 = glmax(u,n)
      m3 = glsum(u,n)
      nn = glsum(1.0*n,1)
      avg = m3 / nn
      if (nio.eq.0) write(*,*)'dbg msk2 ',s6,m1,m2,m3,avg,n,nn

      return
      end
c-----------------------------------------------------------------------
