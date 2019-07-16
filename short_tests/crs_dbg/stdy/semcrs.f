c-----------------------------------------------------------------------
      subroutine get_local_crs_galerkin_t(a,ncl,nxc, h1,h2,h3,
     $                                          ifconv,w1,w2)
c     This routine generates Nelt submatrices of order ncl using
c     Galerkin projection

      include 'SIZE'
      include 'TOTAL'

      real    a(ncl,ncl,*), h1(*),h2(*),h3(*)
      logical ifconv
      real    w1(lx1*ly1*lz1,*),w2(lx1*ly1*lz1,*),w3(lx1*ly1*lz1,lelt)

      parameter (lcrd=lx1**ldim)
      common /ctmp1z/ b(lcrd,8)

      integer e, nel
      logical ifuf,ifcf

      nel=nelt
      if (ifield.eq.1) nel=nelv
      
      do j=1,ncl
         call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
      enddo

      isd  = 1
      imsh = 2
      if (ifield.eq.1) imsh=1

      ifuf = .false.
      ifcf = .true.

      nxyz = lx1*ly1*lz1
      do j = 1,ncl
         do e = 1,nel
            call copy(w1(1,e),b(1,j),nxyz)
         enddo
         call rzero(w2,nxyz*nelt) ! always nelt
         if (ifconv) call convect_new(w2,w1,ifuf,vxd,vyd,vzd,ifcf) ! C^e * bj
         call axhelm_new(w2,w1,h1,h2,h3) ! A^e * bj
         do e = 1,nel
         do i = 1,ncl
            a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * A^e * bj
         enddo ! i
         enddo ! e
      enddo ! j

      return
      end
c-----------------------------------------------------------------------
c
      subroutine set_up_adf_crs(a, h1,h2,h3, ifconv,w1,w2)

      include 'SIZE'
      include 'GEOM'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex

      integer gs_handle
      integer null_space,e

      character*3 cb
      common /scrxxt/ cmlt(lcr,lelt),mask(lcr,lelt)
      common /scrxxti/ ia(lcr,lcr,lelt), ja(lcr,lcr,lelt)
      real cmlt, mask
      integer ia,ja
      real z

      integer key(2),aa(2)
      common /scrch/ iwork(2,lx1*ly1*lz1*lelt)
      common /scrns/ w(7*lx1*ly1*lz1*lelt)
      integer w
      real wr(1)
      equivalence (wr,w)

      real a(*)
      real h1(*), h2(*), h3(*), w1(*), w2(*)
      logical ifconv

      common /ingv/ ngv
      integer*8 ngv

      integer n, nel

      t0 = dnekclock()

c     nxc is order of coarse grid space + 1, nxc=2, linear, 3=quad,etc.
c     nxc=param(82)
c     if (nxc.gt.lxc) then
c        nxc=lxc
c        write(6,*) 'WARNING :: coarse grid space too large',nxc,lxc 
c     endif
c     if (nxc.lt.2) nxc=2

      nxc     = 2
      nx_crs  = nxc

      nel = nelt
      if (ifield.eq.1) nel=nelv

      if(nio.eq.0) write(6,*) 'setup adf coarse grid, nx_crs=', nx_crs

      ncr     = nxc**ldim
      nxyz_c  = ncr
c
c     Set SEM_to_GLOB
c
      call get_vertex
      call set_vert(se_to_gcrs,ngv,nxc,nel,vertex,.true.)

c     Set mask
      z=0
      ntot=nel*nxyz_c
      nzc=1
      if (if3d) nzc=nxc
      call rone(mask,ntot)
      call rone(cmlt,ntot)
      nfaces=2*ldim
c     ifield=1			!c? avo: set in set_overlap through 'TSTEP'?

      if (ifield.eq.1) then
         do ie=1,nel
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'V  '  .or.  cb.eq.'v  ' .or.
     $          cb.eq.'W  '  .or.  cb.eq.'w  ' )
     $           call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
         enddo
         enddo
      elseif (ifield.gt.1) then
         do ie=1,nel
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'T  '  .or.  cb.eq.'t  ')
     $           call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
         enddo
         enddo
      elseif (ifield.eq.ifldmhd) then   ! no ifmhd ?avo?
         do ie=1,nelv
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'ndd'  .or.  cb.eq.'dnd'  .or.  cb.eq.'ddn')
     $          call facev(mask,ie,iface,z,nxc,nxc,nzc)
         enddo
         enddo
      endif

c     Set global index of dirichlet nodes to zero; xxt will ignore them

      call fgslib_gs_setup(gs_handle,se_to_gcrs,ntot,nekcomm,mp)
      call fgslib_gs_op   (gs_handle,mask,1,2,0)  !  "*"
      call fgslib_gs_op   (gs_handle,cmlt,1,1,0)  !  "+"
      call fgslib_gs_free (gs_handle)
      call set_jl_crs_mask(ntot,mask,se_to_gcrs)

      call invcol1(cmlt,ntot)

c     Setup local SEM-based Neumann operators (for now, just full...)

c      if (param(51).eq.1) then     ! old coarse grid
c         nxyz1=lx1*ly1*lz1
c         lda = 27*nxyz1*lelt
c         ldw =  7*nxyz1*lelt
c         call get_local_crs(a,lda,nxc,h1,h2,w,ldw)
c      else
c        NOTE: a(),h1,...,w2() must all be large enough
         n = lx1*ly1*lz1*nel
         call get_local_crs_galerkin_t(a,ncr,nxc,h1,h2,h3,ifconv,w1,w2)
c      endif

      call set_mat_ij(ia,ja,ncr,nel)
      null_space=0
      if (ifield.eq.1) then
         if (ifvcor)  null_space=1
      elseif (ifield.eq.ifldmhd) then
         if (ifbcor)  null_space=1
      endif

c     if(nio.eq.0) write(*,*) 'nullspace = ', null_space

      nz=ncr*ncr*nel
      isolver = param(40)

      call fgslib_crs_setup(xxth(ifield),isolver,nekcomm,mp,ntot,
     $     se_to_gcrs,nz,ia,ja,a, null_space, crs_param)

      t0 = dnekclock()-t0
      if (nio.eq.0) then
         write(6,*) 'done :: setup adf coarse grid ',t0, ' sec'
         write(6,*) ' '
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine adf_coarse_solve(x,r)
c
c     Given an input vector r, this generates the H1 coarse-grid solution
c
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'PARALLEL'
      include 'CTIMER'
      include 'TSTEP'

      real x(*), r(*)

      integer icalld1
      save    icalld1
      data    icalld1 /0/
      
      if (icalld1.eq.0) then ! timer info
         ncrsl=0
         tcrsl=0.0
         icalld1=1
      endif
      ncrsl = ncrsl  + 1

#ifdef TIMER
      etime1=dnekclock()
#endif
c      call crs_lu_solv(x,r)
      call fgslib_crs_solve(xxth(ifield),x,r)
#ifdef TIMER
      tcrsl=tcrsl+dnekclock()-etime1
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine yale_ij(ncl,nel,ia,ja,vert)
      implicit none
      integer ncl,nel

      integer ia(ncl,ncl,*)
      integer ja(ncl,ncl,*)
      integer vert(ncl,*)

      integer e,i,j,k

      do e=1,nel
         do i=1,ncl
         k = vert(i,e)
         do j=1,ncl
            ia(i,j,e) = k
            ja(j,i,e) = k
         enddo ! j
         enddo ! i
      enddo ! e
      return
      end
c-----------------------------------------------------------------------
      subroutine yale2full(ncl,nel,n,afull,aa,vert)
      implicit none
      integer ncl,nel,n
      real afull(n,*)
      real aa(ncl,ncl,*)
      integer vert(ncl,*)
      integer e,i,j, ii,jj

      call rzero(afull,n*n)

      do e=1,nel
         do j=1,ncl
            jj = vert(j,e)
            do i=1,ncl
               ii = vert(i,e)
               afull(ii,jj) = afull(ii,jj) + aa(i,j,e)
            enddo ! i
         enddo ! j
      enddo ! e
      return
      end
c-----------------------------------------------------------------------
      subroutine full2csr(n,a,nnz,ap,aj,acsr,tol)
      implicit none
      integer n,nnz
      real a(n,*)
      integer ap(*), aj(*)
      real acsr(*)
      real tol

      integer i,j,p
  
      p = 1
      do i=1,n
         ap(i) = p
         do j=1,n
         if(abs(a(i,j)).gt.tol) then
            aj(p) = j
            acsr(p) = a(i,j)
            p = p+1
         endif
      enddo ! j
      enddo ! i
      nnz = p-1
      return
      end
c-----------------------------------------------------------------------
      subroutine mask_blk(ncl,nel,a,mask,cmlt)
      implicit none
      integer ncl, nel
      real a(ncl,ncl,*), mask(ncl,*), cmlt(ncl,*)

      integer e,i

      do e=1,nel
      do i=1,ncl
         call dscal(ncl,mask(i,e),a(i,1,e),ncl)
         call dscal(ncl,mask(i,e),a(1,i,e),1)
         if(mask(i,e).eq.0.0E0) a(i,i,e) = cmlt(i,e)
      enddo ! k
      enddo ! e
      return
      end
c-----------------------------------------------------------------------
      subroutine mask_mat(ncl,nel,n,a,vertex,mask)
      implicit none
      integer ncl,nel,n
      real a(n,*), mask(ncl,*)
      integer vertex(ncl,*)
      integer e,i,k

      do e=1,nel
      do k=1,ncl
         i = vertex(k,e)
         call dscal(n,mask(k,e),a(i,1),n)
         call dscal(n,mask(k,e),a(1,i),1)
         if(mask(k,e).eq.0.0E0) a(i,i) = 1.0E0
      enddo ! k
      enddo ! e
      return
      end
c-----------------------------------------------------------------------
      subroutine crs_lu_fact( h1,h2,h3,ifconv,w1,w2)
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'TOTAL'


      real h1(*),h2(*),h3(*)
      logical ifconv
      real w1(*),w2(*)

      integer lg ! TODO common block
      parameter (lg = lcr*lelg)
      common /rcrslu/alu
      common /icrslu/n,ipiv
      real alu(lg*lg)
      integer n,ipiv(lg)
      

      real aaglob(lg*lg)

      integer ap(lcr*lg), aj(lcr*lg)

      integer ncl,nxc,nel,mel
      integer ncr,nnz
      real    aa(lcr,lcr,lelt)

      common /iglcrs/ vglob
      integer vglob(2**ldim,lelg)

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex

      common /scrxxt/ cmlt(lcr,lelt),mask(lcr,lelt)
      real cmlt, mask

      integer info
      integer iglmax, iglsum

      nxc = 2
      ncl = nxc**ndim
      nel = nelt
      if(ifield.eq.1) nel = nelv
      nnz = ncl*ncl*nel
      ncr = ncl*nel

      mel = nelgt

      call set_up_adf_crs(aa, h1,h2,h3, ifconv,w1,w2)
      call mask_blk(ncl,nel,aa,mask,cmlt)

      call icomm(vglob,vertex,ncr)
      call rcomm(aaglob,aa,nnz)

      n = iglmax(vertex,ncr)
      nnz = iglsum(nnz,1)

      if(nio.eq.0) write(*,*) 'crs_lu_fact :: n = ',n 
      call yale2full(ncl,mel,n,alu,aaglob,vglob)
      call dgetrf(n,n,alu,n,ipiv,info)
      if(info.ne.0) then
         write(*,*) 'dgetrf broken, info = ',info
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine crs_lu_solv(x,r)
c     Redundant LU solve        
      implicit none

      include 'SIZE'
      include 'DOMAIN'
      include 'TOTAL'

      real x(lcr,*), r(lcr,*)
      
      integer lg ! TODO common block
      parameter (lg = lcr*lelg)
      common /rcrslu/alu
      common /icrslu/n,ipiv
      real alu(lg*lg)
      integer n,ipiv(lg)


      common /ivrtx/ vertex(2**ldim,lelt)
      integer vertex

      common /iglcrs/ vglob
      integer vglob(2**ldim,lelg)

      common /scrxxt/ cmlt(lcr,lelt),mask(lcr,lelt)
      real cmlt, mask
      
      real xglob(lcr*lelg)
      real rglob(lcr,lelg)

      integer nxc,ncl,ncr,nel,mel
      integer i,k,e
      integer info

      nxc = 2
      ncl = nxc**ndim
      mel = nelg(ifield)
      nel = nelt
      if(ifield.eq.1) nel = nelv
      ncr = ncl*nel

      if(nio.eq.0) write(*,*) 'crs_lu_solv :: n = ',n

      call rcomm(rglob,r,ncr)
      call rzero(xglob,n)

c     Gather
      do e=1,mel
         do k=1,ncl
            i = vglob(k,e)
            xglob(i) = xglob(i) + rglob(k,e) 
         enddo ! k
      enddo ! e

c     Solve      
      call dgetrs('N',n,1,alu,n,ipiv,xglob,n,info)
      if(info.ne.0) then
         write(*,*) 'dgetrs broken, info = ',info
         call exitt
      endif

c     Scatter
      do e=1,nel
         do k=1,ncl
            i = vertex(k,e)
            x(k,e) = xglob(i)
         enddo ! k
      enddo ! e

      call col2(x,mask,ncr)
      return
      end
c-----------------------------------------------------------------------
      subroutine test_comm
c     Test collect_rarray and collect_iarray
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'

      integer e,eg,e0,i,n,nn
      real nels(lelg)

      real rar_in(lelt*3),rar_out(lelg*3)
      integer iar_in(lelt*3),iar_out(lelg*3)

      ! test case 1
      ifield = 2
      n=3 ! test case size

      do e=1,nelfld(ifield)
      do i=1,n
        eg = lglel(e)
        rar_in((e-1)*n+i) = real(eg)**2+i-1 + 0.5
        iar_in((e-1)*n+i) = real(eg)**2+i-1
        e0 = (e-1)*n+i
        write(*,60) 'aa',nid,'  arr=',i,e,rar_in(e0),iar_in(e0)
      enddo
      enddo

      call rcomm1(rar_out,rar_in,n)
      call icomm1(iar_out,iar_in,n)

      do e=1,nelg(ifield)
      do i=1,n
        e0 = (e-1)*n+i
        write(*,60) 'bb',nid,'  arr=',i,e,rar_out(e0),iar_out(e0)
      enddo
      enddo


      ! test case 2
      n=4
      if (nid.gt.0) n=5

      do i=1,n
        rar_in(i) = real(n*nid+i) + 0.5
        iar_in(i) = n*nid+i
        write(*,60) 'cc',nid,'  arr=',i,n,rar_in(i),iar_in(i)
      enddo

      call rcomm(rar_out,rar_in,n)
      call icomm(iar_out,iar_in,n)

      nn = n
      if (nid.gt.0) nn = n-1
      nn = nn+(np-1)*(nn+1)
      do i=1,nn
        write(*,60) 'dd',nid,'  arr=',i,nn,rar_out(i),iar_out(i)
      enddo

   60 format(A,i1,A,i8,i6,f14.2,i14)

      write(*,*)'before exitt',nid
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine rcomm(arr_out,arr_in,n)
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'

      real arr_out(*),arr_in(*)
      integer nn,n,mtype,nns(0:np)
      integer ii,id,iglsum

      mtype = 0 ! don't know what's this doing...

      id = 1
      ! Get nel 
      do ii=0,np-1
        call nekgsync()
        if (ii.eq.0) then
          nns(ii) = n
          call copy(arr_out(id),arr_in,n)
          id = id + nns(ii)
        elseif (np.gt.1) then
          if (nid.eq.0) then
            call crecv(mtype,nns(ii),4)
          endif
          if (ii.eq.nid) then
            call csend(mtype,n,4,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,arr_out(id),8*nns(ii))
            id = id + nns(ii)
          endif
          if (ii.eq.nid) then
            call csend(mtype,arr_in,8*n,0,0)
          endif
        endif
      enddo

      ! broadcast
      nn =iglsum(n,1)
      if (np.gt.1) call bcast(arr_out,nn*8)

      return
      end
c-----------------------------------------------------------------------
      subroutine icomm(arr_out,arr_in,n)
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'

      integer arr_out(*),arr_in(*)
      integer nn,n,mtype,nns(0:np)
      integer ii,id,iglsum

      mtype = 0 ! don't know what's this doing...

      id = 1
      ! Get nel 
      do ii=0,np-1
        call nekgsync()
        if (ii.eq.0) then
          nns(ii) = n
          call icopy(arr_out(id),arr_in,n)
          id = id + nns(ii)
        elseif (np.gt.1) then
          if (nid.eq.0) then
            call crecv(mtype,nns(ii),4)
          endif
          if (ii.eq.nid) then
            call csend(mtype,n,4,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,arr_out(id),4*nns(ii))
            id = id + nns(ii)
          endif
          if (ii.eq.nid) then
            call csend(mtype,arr_in,4*n,0,0)
          endif
        endif
      enddo

      ! broadcast
      nn =iglsum(n,1)
      if (np.gt.1) call bcast(arr_out,nn*4)

      return
      end
c-----------------------------------------------------------------------
 
      subroutine rcomm1(arr_out,arr_in,n)
c     arr_in is a local array of size n * nelt
c     arr_out = {arr_in}_{ii=1,np} + reordering w.r.t. global element index
c     ifcpall = T, broadcast arr_out to all processors
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'

      real arr_out(lelg*n),arr_in(lelt*n),arr_col(n*lelt,0:np)
      integer n,nel,nnel,mtype
      integer eg_nid,nels(0:np),egs(0:np),lglel_col(lelt,0:np)
      integer e0,e1,e,eg,e_now,ii,i,j

      nel = nelfld(ifield)
      nnel= nelg(ifield)

      mtype = 0 ! don't know what's this doing...

      ! Get nel 
      do ii=0,np-1
        call nekgsync()
        if (ii.eq.0) then
          nels(ii) = nel
          call icopy(lglel_col(1,ii),lglel,nel)
          call copy(arr_col(1,ii),arr_in,n*nel)
        elseif (np.gt.1) then
          if (nid.eq.0) then
            call crecv(mtype,nels(ii),4)
          endif
          if (ii.eq.nid) then
            call csend(mtype,nelt,4,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,lglel_col(1,ii),4*nels(ii)*n)
          endif
          if (ii.eq.nid) then
            call csend(mtype,lglel,4*nel*n,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,arr_col(1,ii),8*nels(ii)*n)
          endif
          if (ii.eq.nid) then
            call csend(mtype,arr_in,8*nel*n,0,0)
          endif
        endif
      enddo

      ! re-ordering
      if (nid.eq.0) then
        call izero(egs(0),np)
        do ii=0,np-1
        do e=1,nels(ii)
          eg = lglel_col(e,ii)
          eg_nid = gllnid(eg)

          e0 = 1+(eg-1)*n
          e1 = 1+(e-1)*n
          call copy(arr_out(e0),arr_col(e1,ii),n)
        enddo
        enddo
      endif

      ! broadcast
      if (np.gt.1) call bcast(arr_out,nnel*n*8)

      return
      end
c-----------------------------------------------------------------------
      subroutine icomm1(arr_out,arr_in,n)
c     arr_in is a local array of size n * nelt
c     arr_out = {arr_in}_{ii=1,np} + reordering w.r.t. global element index
c     ifcpall = T, broadcast arr_out to all processors
      implicit none
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'

      integer arr_out(lelg*n),arr_in(lelt*n),arr_col(n*lelt,0:np)
      integer n,nel,nnel,mtype
      integer eg_nid,nels(0:np),egs(0:np),lglel_col(lelt,0:np)
      integer e0,e1,e,eg,e_now,ii,i,j

      nel = nelfld(ifield)
      nnel= nelg(ifield)

      mtype = 0 ! don't know what's this doing...

      ! Get nel 
      do ii=0,np-1
        call nekgsync()
        if (ii.eq.0) then
          nels(ii) = nel
          call icopy(lglel_col(1,ii),lglel,nel)
          call icopy(arr_col(1,ii),arr_in,n*nel)
        elseif (np.gt.1) then
          if (nid.eq.0) then
            call crecv(mtype,nels(ii),4)
          endif
          if (ii.eq.nid) then
            call csend(mtype,nelt,4,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,lglel_col(1,ii),4*nels(ii)*n)
          endif
          if (ii.eq.nid) then
            call csend(mtype,lglel,4*nel*n,0,0)
          endif

          if (nid.eq.0) then
            call crecv(mtype,arr_col(1,ii),4*nels(ii)*n)
          endif
          if (ii.eq.nid) then
            call csend(mtype,arr_in,4*nel*n,0,0)
          endif
        endif
      enddo

      ! re-ordering
      if (nid.eq.0) then
        call izero(egs(0),np)
        do ii=0,np-1
        do e=1,nels(ii)
          eg = lglel_col(e,ii)
          eg_nid = gllnid(eg)

          e0 = 1+(eg-1)*n
          e1 = 1+(e-1)*n
          call copy(arr_out(e0),arr_col(e1,ii),n)
        enddo
        enddo
      endif

      ! broadcast
      if (np.gt.1) call bcast(arr_out,nnel*n*4)

      return
      end
c-----------------------------------------------------------------------
