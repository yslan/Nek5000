c-----------------------------------------------------------------------
c---- MG Prec for a scalar field
c-----------------------------------------------------------------------
      subroutine scprec(z,w,ifrefc) ! z = M w

!     assume one scalar field, with nelt

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL' ! TOTAL does not have HSMG

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)
      real     z(*), w(*)
      logical ifrefc

      mxmg = 1
      isw = 0
      ifield = 2

      call adfmn_mg_schwz_tens(z,w,h1,h2,h3,isw,mxmg,tmult,ifrefc)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_tens_genall() ! tensor setup
      include 'SIZE'
      include 'HSMGL'
      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)
      call tens_fact(mg_h1,mg_h2,h3) ! from nk_prec
      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_schwzmg_tens_solve(x,r,l,nx)
      real    x(*), r(*) 
      call tens_solv(x,r,l,nx) ! from nk_prec
      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_mg_schwz_tens(x,res,h1,h2,h3,iter,maxit,wt
     $                             ,ifrefc)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL' ! TOTAL does not have HSMGL

      real     x(*), res(*), wt(*)
      real     h1(*),h2(*),h3(*)
      integer  iter, maxit, lvmin
      logical  ifrefc

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      real     tol
      common /adfscr/ tmp(lt), r0(lt)
      real            tmp    , r0
      common /adfexc/ uex(lt)
      real            uex
      integer  icld
      data     icld / 0 /
      save     icld

      tol = 1.e-12
      ifield = 2 ! temp
      n = nx1*ny1*nz1*nelt
      if(icld.eq.0.or.ifrefc) then

        lv = max(nint(param(114)),2) ! Read number of levels from rea
        ! check if lv valid, Lan
        if (param(113).gt.0.5) then ! lx1/2
          lvmin = log2(nint(real(lx1)/(param(115)+1)))+nint(param(115))
          lv = min(lv,lvmin)
        else
          lvmin = nint((real(lx1)-(param(115)+1))/2.)+nint(param(115))
          lv = min(lv,lvmin)
        endif

        if(nio.eq.0)
     $    write(6,*) 'adfmn_mg_schwz_tens: into setup_tens, levels=',lv
        call adfmn_schwzmg_setup_tens(lv)
        icld = icld + 1
      endif
      resn = sqrt(glsc3(res,res,wt,n))
      call copy (r0, res, n)
      call rzero(x,n)

      iter = 0
      do it=1,maxit
        call adfmn_schwzmg_vcycle_tens(tmp,r0,icrs) ! v-cycle
        iter = iter + icrs
        call add2(x,tmp,n)
        call adfax (r0,x,h1,h2,h3)   !
        call sub2  (r0,res,n)
        call chsign(r0,n)      ! r0 = res - A x
        rn0 = sqrt(glsc3(r0,r0,wt,n))

        if(rn0.le.tol) then
          if(nio.eq.0) write(6,10) it, maxit, rn0, resn, tol
          goto 88
        endif
        if(nio.eq.0)   write(6, 9) it, maxit, rn0, resn, tol
      enddo
   9  format(2i5,1p3e12.4, ' mgshwz')
   10 format(2i5,1p3e12.4, ' mgshwzb')
   88 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmn_schwzmg_setup_tens(lv)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      include 'SEMHAT'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)
      real w1(lt),w2(lt),w3(lt)
      logical ifconv

      parameter (lxyz = lx1*ly1*lz1, 
     $           lxxe = lx1*lx1*lelt, lyye  = ly1*ly1*lelt,
     $           lzze = lz1*lz1*lelt, lxyze = lxyz*lelt    )

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical         ifdfrm, iffast, ifh2, ifsolv
      common /adfco/  co
      real            co(lmg_g*(ldim)*lelt) ! big enough?
      common /adfcoi/ pco
      integer         pco(lmgx,ldimt1) ! big enough?
      common /logmn/ ifrstr
      logical        ifrstr

      integer p_h1,p_h2,p_g,p_b,p_msk,e,p_c,p_el,p_ev
      integer icall
      save    icall
      data    icall /0/

      if (nio.eq.0) write(6,*) 'adf: famg_setup',icall,ifield
      icall = icall+1

      ifrstr = .true.  ! RAS
      ifrstr = .false. ! WAS
      if(ifrstr) then
        if(nio.eq.0) write(6,*) 'adfmn_schwzmg: Restricted AS'
      else
        if(nio.eq.0) write(6,*) 'adfmn_schwzmg: Weighted   AS'
      endif

      n = lx1*ly1*lz1*nelt

      mg_fld = ifield
      call adfmg_index_0 ! initialize index sets
      call adfmg_setup_nx(lv)! Sets the level schedules
      call adfmg_setup_semhat ! SEM hat matrices for each level
      call adfmg_setup_intp   ! Interpolation operators
      call adfmg_setup_dssum  ! set direct stiffness summation handles
      call adfmg_setup_wtmask ! set restriction weight matrices and bc masks
      call adfmn_setup_schwzmn_wt(.false.)
      call adfmg_setup_solve  ! set up the solver
      
      l=mg_h1_lmax ! == mg_lmax
      call adfmg_set_h1  (p_h1 ,l)
      call adfmg_set_h2  (p_h2 ,l)
      call adfmg_set_gb  (p_g  ,p_b,l)
      call adfmg_set_msk (p_msk,l) ! pointers to imask
      call adfmg_set_co  (p_c,l)   ! c coefficients
      call adfmg_set_cmlt
      call adfmg_tens_genall()

      call rzero(w3,n) ! corner mesh - Pablo
      if(mg_nh(1).eq.2) call crs_lu_fact(mg_h1,mg_h2,w3,.true.,w1,w2)
      return
      end
c-----------------------------------------------------------------------
      subroutine adfmn_schwzmg_vcycle_tens(z,rhs,icrs) ! V-cycle

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'HSMGL'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)

      real z(lt),rhs(lt)                        ! z = M \ rhs
      common /zcrcq/ e(4*lt), r(4*lt),w(lt)

      common /adfexc/ uex
      real            uex(lt)
      real            err(lt)

      logical ifdssum

      integer p_msk,p_b,icrs
      integer icall
      save    icall
      data    icall /0/

      ifield=2
      nel = nelfld(ifield)
      nel = nelt
      ntot=lx1*ly1*lz1*nel !

      if (icall.eq.0) then
         icall=1
      endif

      sigma =  1.0
      sigma =  0.5
      op    =  1.0                           ! Coefficients for residual
      om    = -1.0                           !    evaluation
      o0    =  0.0                           !    evaluation

      nsmth = 1 ! 

      do iv=1,1 ! only one V-cycle as preconditioner

       l  = mg_h1_lmax
       n  = mg_h1_n(l, mg_fld)
       is = 1                                    ! Solve index
       if (iv.eq.1) then
        call rzero(z,n)                          ! make zero
        call copy(r,rhs,n)                       ! r := rhs
        do ism=1,nsmth
        call adfmn_schwz_l_tens(e,r,sigma,l)     ! z := S r
        call add2(z,e,n)
        call adfmn_axmg(r,e,op,om,l,w)           ! r := r - Ae
        enddo
       else
        call copy(r,rhs,n)                       ! r := rhs
        call adfmg_axm(r,z,op,om,l,w)            ! r := r - Az
        do ism=1,nsmth
        call adfmn_schwz_l_tens(e,r,sigma,l)     ! e := S r
        call add2(z,e,n)                         ! z := z+e
        call adfmn_axmg(r,e,op,om,l,w)           ! r := r - Ae
        enddo
       endif

       ifdssum = .true.
       do l = mg_h1_lmax - 1, 2, -1               ! Restrict down to
         is = is + n                              ! coarse level
         n  = mg_h1_n(l, mg_fld)                  !       T
         call adfmg_rstr(r, l, ifdssum)           ! r := J  r
         call adfmn_schwz_l_tens(e(is),r,sigma,l) ! e := sigma W S r
         call adfmn_axmg(r,e(is),op,om,l,w)       ! r := r - A e
       enddo
       
       is = is + n                                
       l  = 1                                     ! lowest level
       ifdssum = .false.                          !        T
       call adfmg_rstr(r, l, ifdssum)             ! r  := J  r
       
       p_msk=p_mg_msk(l, mg_fld)
       call adfmg_mask(r, mg_imask(p_msk), nel)   !       -1    Mask and
       call adfmn_coarse_solve(e(is), r,icrs)     ! e := A  r   solve at 
       call adfmg_mask(e(is),mg_imask(p_msk),nel) !         coarse level

       do l = 2, mg_h1_lmax-1               ! Prolongate to finest level
         n  = mg_h1_n(l,mg_fld)
         im = is
         is = is - n
         call adfmg_intp (w,e(im),l-1)            ! w  :=  J e
         call add2      (e(is),w,n)               ! e  := e + w
       enddo

       l = mg_h1_lmax
       n = mg_h1_n(l,mg_fld)
       im = is
       call adfmg_intp(w,e(im),l-1)               ! w :=  J e
       call add2     (z,w,n)                      ! e  := e  + w
       call col2     (z,tmask,n)                  ! mask z !FIXME: check this?
       call dsavg    (z)                          ! ensure continuous z
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_schwz_l_tens(e,r,sigma,l) !h1mg_schwarz
      include 'SIZE'
      include 'HSMGL'
      include 'SOLN'

      real    e(*),r(*)
      real    sigma
      integer l
      common /logmn/ ifrstr
      logical        ifrstr

      n  = mg_h1_n(l,mg_fld)
      nh = mg_nh(l)

      call adfmn_schwz_l_tens_part1 (e,r,l) ! !
      if(ifrstr) then
        ! do nothing
        ! write(6,*) 'doing nothing in schwz l RAS path'
      else
        call adfmg_schwarz_wt  (e,l)          ! e  := W e
      endif
      call cmult               (e,sigma,n)    !  l       l
      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_schwz_l_tens_part1(e,r,l)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)

      integer enx,eny,enz,pmsk
      real    tmp(lx1,ly1,lz1,lelt)
      common /logmn/ ifrstr
      logical        ifrstr

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pmsk = p_mg_msk(l,mg_fld)
      nx = mg_nh(l)

      ifield = 2
      call adfmg_mask (r,mg_imask(pmsk),nelfld(ifield))
      call adfmn_schwzmg_tens_solve(e,r,l,nx) !

      if(ifrstr) then
        ! write(6,*) 'doing masking in schwz l part 1 RAS path'
        call adfmn_schwz_rs(e,l) !
      endif

      call adfmg_dssum(e,l)                           ! sum border nodes
      call adfmg_mask (e,mg_imask(pmsk),nelfld(ifield)) ! apply mask 

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_coarse_solve(e,r,icrs)
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'GEOM'
      include 'SOLN'
      include 'PARALLEL'
      include 'HSMGL'
      include 'CTIMER'
      include 'INPUT'
      include 'TSTEP'
      real    e(1),r(1)
      real    w(lx1*ly1*lz1*lelt), tmp(lx1*ly1*lz1*lelt)
      integer l
      real    eu(lx1*ly1*lz1*lelt)
c
      integer n_crs_tot
      save    n_crs_tot
      data    n_crs_tot /0/
      real    r0(lx1*ly1*lz1*lelt)
      integer p_msk

      real  is
c     NOTE coarse solver takes care of dssum - Pablo

      if (icalld.eq.0) then ! timer info
         ncrsl=0
         tcrsl=0.0
      endif
      icalld = 1

      if (ifsync) call nekgsync()

      ncrsl  = ncrsl  + 1
      etime1=dnekclock()

      l  = 1 ! lowest level
      nh = mg_nh  (l)
      nr = mg_h1_n(l,mg_fld)

c     nh = -1 ! uncomment to disable crs_solve

      if (nh.eq.2) then ! corner mesh is always solved wih crs_solve
         call rzero(e,nr)
         call copy (r0,r,nr)
         call adf_coarse_solve(e,r0)
         icrs = 1
      else ! repeated Schwarz until converges
         call adfmg_dssum(r,l)
         mxc = 200
         mxc = min(lgmres,mxc)
         call copy (r0,r,nr)
         call rzero(e,nr)
         p_msk=p_mg_msk(l, mg_fld)
         call adfmn_proj(e,r,icrs,mxc,l)
      endif

      tcrsl=tcrsl+dnekclock()-etime1

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmn_proj(x,res,iter,maxit,l)

c     Solve A t = r, via custome projection scheme

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real    res(1), x(1)
      integer iter, maxit, l

      common /adfwt/ cmlt(lx1*ly1*lz1*lelt)
      real cmlt ! multiplicity at coarse level only, very careful!

      real r(lt),w(lt),d(lt),v(lt),z(lt)
      common /dcrspj/ d       ! scratch
      common /ccrspj/ r,w     ! WORK ARRAYS
      common /zcrsp1/ v
      common /zcrsp2/ z
      common /zcrsps/ p,q
      real            p(lt,lgmres), q(lt,lgmres)
      real    x0(lt)
      real    proj_tol
      real    alpha, beta

      common /adfexc/ uex
      real            uex(lt)

      integer icalld
      data    icalld / 0 /
      save    icalld

      ngo = -99        ! TURN OFF VERBOSE OUTPUT
      n = mg_h1_n(l,mg_fld) ! length for multigrid(< ntot)
      nh = mg_nh  (l)
      proj_tol = 1.e-10
      if (istep.le.20.or.mod(istep,10).eq.0) ngo = nio

      if(icalld.eq.0) then
        icalld = icalld + 1
      endif

      m   = lgmres
      mxs = 1 ! maximum schwarz iter - once?
      iconv = 0 ! if converged
      sigma = 0.50
      sigma = 0.25
      zr = 0.
      op = 1.

      ! initial guess
c     call rzero(x,n)    ! comment out to allow input initial guess

      rn0 = sqrt(glsc3(res,res,cmlt,n))
      rn = rn0*1.e-5 ! rodcht
      rn = rn0*1.e-2 ! singlerod  - Pablo
      proj_tol = max(rn,proj_tol)

      iter = 1
      do while (iter.lt.maxit)    ! Main proj loop
         if(iter.eq.0) then                      ! first step :: r = res - A 0
            call copy   (r,res,n)  ! allow input initial guess
         else                                     ! second up  :: r = res - A x
            call adfmn_axmg(w,x,zr,op,l,d)        ! w = A x
            call sub3   (r,res,w,n)               ! r = r - w
         endif

         ! precoditioner solve
                                                  !       -1
c        call copy  (z,r,n)                       ! z  = M  v
                                                  !  j       j
         call adfmn_schwz_l_tens(z,r,sigma,l)
         call adfmn_axmg(w,z,zr,op,l,d)           ! w = A z
                                                  !        j
         do j=1,(iter-1)
           beta = glsc3(q(1,j),w,cmlt,n)
           call add2s2(z,p(1,j),-beta,n)
           call add2s2(w,q(1,j),-beta,n)
         enddo
         beta = sqrt(glsc3(w,w,cmlt,n))
         if(beta.lt.proj_tol) goto 900

         betai = 1./beta
         call copy (p(1,iter),z,n)
         call cmult(p(1,iter),betai,n)
         call copy (q(1,iter),w,n)
         call cmult(q(1,iter),betai,n)

         alpha = glsc3(q(1,iter),r,cmlt,n)
         call add2s2(x,p(1,iter),alpha,n)

         if (ngo.eq.0) write(6,9)
     $         n,iter,beta,alpha,proj_tol,dt
    9    format(i9,i5,1p4e12.4,' mgprj')
         iter = iter + 1
      enddo
  900 continue
      maxit = iter
      if (nio.eq.0) write(6,8)
     $   n,maxit,beta,proj_tol,dt
    8    format(i9,i5,1p3e12.4,' mgprjb')

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmn_axmg(w,p,aw,ap,l,wk)
c
c     w  := aw*w + ap*H*p, level l, with mask and dssum
c
c     Hu := div. h1 grad u + h2 u
c
c        ~= h1 A u + h2 B u
c
c     Here, we assume that pointers into g() and h1() and h2() have
c     been established
c
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'  ! nelfld

      real w(1),p(1),wk(1)

      integer p_h1,p_h2,p_g,p_b,p_msk, p_c
      logical ifh2
      common /adfco/  co
      real            co(lmg_g*(ldim)*lelt) ! big enough?
      common /adfcoi/ pco
      integer         pco(lmgx,ldimt1) ! big enough?

      real aw, ap

      p_h1  = p_mg_h1  (l,mg_fld)
      p_h2  = p_mg_h2  (l,mg_fld)
      p_g   = p_mg_g   (l,mg_fld)
      p_b   = p_mg_b   (l,mg_fld)
      p_msk = p_mg_msk (l,mg_fld)
c     write(6,*) 'inside adfmg axm:: p_msk',p_msk,l
      p_c   = pco      (l,mg_fld) ! right


c     Shouldn't these be set only at the beginning >:( - Pablo
c     if (p_h1 .eq.0) call adfmg_set_h1  (p_h1 ,l)
c     if (p_h2 .eq.0) call adfmg_set_h2  (p_h2 ,l)
c     if (p_g  .eq.0) call adfmg_set_gb  (p_g,p_b,l)
c     if (p_msk.eq.0) call adfmg_set_msk (p_msk,l)
c     if (p_c  .eq.0) call adfmg_set_co  (p_c,l)

      ifh2 = .true.

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)
      ng = 3*ldim-3
      nc = ldim


      call adfmg_axml (wk,p
     $    ,mg_h1(p_h1),mg_h2(p_h2),nx,ny,nz,nelfld(ifield)
     $    ,mg_g (p_g) , ng ,mg_b(p_b), mg_imask(p_msk),ifh2
     $    ,co(p_c), nc)

      ! wk = C p + A p

c     n = nx*ny*nz*nelfld(ifield)
      n = nx*ny*nz*nelfld(mg_fld)
      call adfmg_dssum  (wk,l) ! hsmg_dssum

      ! w = aw * w + ap * wk = 1 * w - 1 * wk = w - ( C p + A p)
      call add2sxy    (w,aw,wk,ap,n)

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_schwzmn_wt_1(wt,l,ifsqrt)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMGL'

      real wt(1),work(1)
      logical ifsqrt

      integer enx,eny,enz,pm

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pm = p_mg_msk(l,mg_fld)

      enx=mg_nh(l) ! no extrusion
      eny=mg_nh(l)
      enz=mg_nh(l)
      if(.not.if3d) enz=1
      ns = enx*eny*enz*nelfld(mg_fld)

      call rone(mg_work,ns)
      call adfmg_dssum(mg_work,l)                           ! sum border nodes

c     write(6,*) n,l,'n, l'
c     write(6,*) 'wt 1',(mg_work(j),j=1,n)

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nh(l)
      if (.not.if3d) nz=1
      nxyz = nx*ny*nz
      k    = 1
      do ie=1,nelfld(ifield)
c        call outmat(mg_work(k),nx,ny,'NEW WT',ie)
         call adfmg_setup_schwzmn_wt_2(wt,ie,nx,mg_work(k),ifsqrt)
         k = k+nxyz
      enddo

c     do ie=1,nelt
c       j0 = 1 + (ie-1)*4*2*nx
c       write(6,*) 'wt 2',(wt(j),j=j0,j0+4*nx*2-1)
c     enddo
c     stop

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_setup_schwzmn_wt_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      real wt(1),work(1)
      logical ifsqrt

      if(ldim.eq.2) call adfmg_setup_schwzmn_wt2d_2(wt,ie,n,work
     $                                             ,ifsqrt)
      if(ldim.eq.3) call adfmg_setup_schwzmn_wt3d_2(wt,ie,n,work
     $                                             ,ifsqrt)

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_schwzmn_wt2d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,4,2,nelt)
      real work(n,n)
      
      integer ie,i,j
      do j=1,n
         wt(j,1,1,ie)=1.0/work(1,j)
         wt(j,2,1,ie)=1.0/work(2,j)
         wt(j,3,1,ie)=1.0/work(n-1,j)
         wt(j,4,1,ie)=1.0/work(n,j)
      enddo
      do i=1,n
         wt(i,1,2,ie)=1.0/work(i,1)
         wt(i,2,2,ie)=1.0/work(i,2)
         wt(i,3,2,ie)=1.0/work(i,n-1)
         wt(i,4,2,ie)=1.0/work(i,n)
      enddo
      if(ifsqrt) then
         do ii=1,2
         do j=1,4
         do i=1,n
            wt(i,j,ii,ie)=sqrt(wt(i,j,ii,ie))
         enddo
         enddo
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_schwzmn_wt3d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,n,4,3,nelt)
      real work(n,n,n)
      
      integer ie,i,j,k
      integer lbr,rbr,lbs,rbs,lbt,rbt

      ierr = 0
      do k=1,n
      do j=1,n
         wt(j,k,1,1,ie)=1.0/work(1,j,k)
         wt(j,k,2,1,ie)=1.0/work(2,j,k)
         wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
         wt(j,k,4,1,ie)=1.0/work(n,j,k)
      enddo
      enddo
      do k=1,n
      do i=1,n
         wt(i,k,1,2,ie)=1.0/work(i,1,k)
         wt(i,k,2,2,ie)=1.0/work(i,2,k)
         wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
         wt(i,k,4,2,ie)=1.0/work(i,n,k)
      enddo
      enddo
      do j=1,n
      do i=1,n
         wt(i,j,1,3,ie)=1.0/work(i,j,1)
         wt(i,j,2,3,ie)=1.0/work(i,j,2)
         wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
         wt(i,j,4,3,ie)=1.0/work(i,j,n)
      enddo
      enddo
      if(ifsqrt) then
         do ii=1,3
         do k=1,4
         do j=1,4
         do i=1,n
            wt(i,j,k,ii,ie)=sqrt(wt(i,j,k,ii,ie))
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_schwz_rs(u,l) ! apply the weights
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      parameter(lt=lx1*ly1*lz1*lelt)

      common /adfrsm/ ras_mask, wrk
      real            ras_mask(lx1*ly1*lz1*lelt,lmgn)
     $              , wrk(lx1*ly1*lz1*lelt)

      real u(1)
      integer l

      nl = mg_nh(l)
      call col2(u,ras_mask(1,l),nl**ldim)
      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_setup_schwzmn_wt(ifsqrt)
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      
      integer l,i,nl,nlz,itmp
      common /adfrsm/ ras_mask, wrk
      real            ras_mask(lx1*ly1*lz1*lelt,lmgn)
     $              , wrk(lx1*ly1*lz1*lelt)

      common /logmn/ ifrstr
      logical        ifrstr
      logical ifsqrt

      i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)

      do l=1,mg_lmax

         mg_schwarz_wt_index(l,mg_fld)=i
         nl  = mg_nh(l)
         nlz = mg_nhz(l)
         i   = i+nl*nlz*4*ldim*nelt

         if (i .gt. lmg_swt*4*ldim*lelt) then
            itmp = i/(4*ldim*lelt)
            write(6,*) 'lmg_swt too small li',i,itmp,lmg_swt,l
            call exitt
         endif

         call adfmg_setup_schwzmn_wt_1(
     $      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)

         if(ifrstr) then
           write(6,*) 'creating ras mask in setup wt',l
           nx = mg_nh(l)
           ny = mg_nh(l)
           nz = mg_nhz(l)
           call adfmn_setup_mask_rt(ras_mask(1,l),nx,ny,nz,l,wrk)
         endif

c        if(nid.eq.0) write(6,*) 'write out ras mask',l
c        call adfmg_outpost(ras_mask(1,l),l)
      enddo

      mg_schwarz_wt_index(l,mg_fld)=i
c     call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_setup_mask_rt_3(wv,nx,ny,nz,l,w) ! No 3D
      include 'SIZE'
      include 'INPUT'
      include 'SOLN' ! for outpost debugging only, remove later
      include 'PARALLEL' ! for outpost debugging only, remove later
      integer nx,ny,nz,l
      real  w(nx,ny,nz,nelt)
      real wv(nx,ny,nz,nelt)

      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      common /logmn/ ifrstr
      logical        ifrstr

      n = nx*ny*nz*nelt
      call rone(w,n)
      call rone(wv,n)

c     set neumann nodes to zero
      ierr = 0
      two  = 2
      do ie=1,nelt
        ! general case
        call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
        if (ierr.ne.0) then
           ierr = -1
           call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
        endif

          if(lbr.eq.2) then !
            iy0 = 2
            if(lbs.eq.2) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.2) iy1 = ny
            iz0=2
            if(lbt.eq.2) iz0 = 1
            iz1=nz-1
            if(rbt.eq.2) iz1 = nz
            do k=iz0,iz1
            do j=iy0,iy1
              wv(1,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbr.eq.2) then
            iy0 = 2
            if(lbs.eq.2) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.2) iy1 = ny
            iz0=2
            if(lbt.eq.2) iz0 = 1
            iz1=nz-1
            if(rbt.eq.2) iz1 = nz
            do k=iz0,iz1
            do j=iy0,iy1
              wv(nx,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(lbs.eq.2) then
            ix0 = 2
            if(lbr.eq.2) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.2) ix1 = nx
            iz0=2
            if(lbt.eq.2) iz0 = 1
            iz1=nz-1
            if(rbt.eq.2) iz1 = nz
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,1,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbs.eq.2) then
            ix0 = 2
            if(lbr.eq.2) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.2) ix1 = nx
            iz0=2
            if(lbt.eq.2) iz0 = 1
            iz1=nz-1
            if(rbt.eq.2) iz1 = nz
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,ny,k,ie) = 0.
            enddo
            enddo
          endif
          if(lbt.eq.2) then
            ix0=2
            if(lbr.eq.2) ix0=1
            ix1=nx-1
            if(lbr.eq.2) ix1=nx
            iy0=2
            if(lbr.eq.2) iy0=1
            iy1=ny-1
            if(lbr.eq.2) iy1=ny
            do j=iy0,iy1
            do i=ix0,ix1
              wv(i,j,1,ie) = 0.
            enddo
            enddo
          endif
          if(rbt.eq.2) then
            ix0=2
            if(lbr.eq.2) ix0=1
            ix1=nx-1
            if(lbr.eq.2) ix1=nx
            iy0=2
            if(lbr.eq.2) iy0=1
            iy1=ny-1
            if(lbr.eq.2) iy1=ny
            do j=iy0,iy1
            do i=ix0,ix1
              wv(i,j,nz,ie) = 0.
            enddo
            enddo
          endif

          ieg = lglel(ie)
          if(ieg.gt.18) then
            write(6,*) 'ie,ieg  ',ie,ieg
            write(6,*) 'wv  ',(wv(i,1,1,ie),i=1,nx*ny*nz)
          endif
      enddo


      call copy(w, wv, n)

c     if(nid.eq.0) write(6,*) 'AA : write out ras mask',l
c     call adfmg_outpost(w,l)

      call adfmg_dssum(w,l) ! summed up

c     if(nid.eq.0) write(6,*) 'AB : write out ras mask',l
c     call adfmg_outpost(w,l)

      do i=1,nx*ny*nz*nelt
        if(w(i,1,1,1).lt.1e-7) then  ! one elem, right=outflow(t field)
          write(6,*) 'Error in rst mask',i,w(i,1,1,1),', set to 1.'
          w (i,1,1,1) = 1.
          wv(i,1,1,1) = 1.
        endif
        w(i,1,1,1) = 1./w(i,1,1,1)
      enddo

c     if(nid.eq.0) write(6,*) 'AC : write out ras mask',l
c     call adfmg_outpost(w,l)

      call col2(wv,w,n) ! zero remains zero, otherwize 1/multiplic

c     if(nid.eq.0) write(6,*) 'AD : write out ras mask',l
c     call adfmg_outpost(w,l)
c     if(l.eq.2) call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_setup_mask_rt_2(wv,nx,ny,nz,l,w) ! No 3D
      include 'SIZE'
      include 'INPUT'
      include 'SOLN' ! for outpost debugging only, remove later!
      integer nx,ny,nz,l
      real  w(nx,ny,nz,nelt)
      real wv(nx,ny,nz,nelt)

      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      common /logmn/ ifrstr
      logical        ifrstr

      n = nx*ny*nz*nelt
      call rone(w,n)
      call rone(wv,n)

c     set neumann nodes to zero
      ierr = 0
      two  = 2
      do ie=1,nelt
        ! general case
        call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
        if (ierr.ne.0) then
           ierr = -1
           call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
        endif
        sumbc = lbr+rbr+lbs+rbs
        if(sumbc.gt.4)
     $  write(6,*) ie,lbr,rbr,lbs,rbs,' fast bcs'

          iz0 = 1
          iz1 = nz
          if(lbr.eq.2) then !
            iy0 = 2
            if(lbs.eq.2) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.2) iy1 = ny
            do k=iz0,iz1
            do j=iy0,iy1
              wv(1,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbr.eq.2) then
            iy0 = 2
            if(lbs.eq.2) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.2) iy1 = ny
            do k=iz0,iz1
            do j=iy0,iy1
              wv(nx,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(lbs.eq.2) then
            ix0 = 2
            if(lbr.eq.2) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.2) ix1 = nx
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,1,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbs.eq.2) then
            ix0 = 2
            if(lbr.eq.2) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.2) ix1 = nx
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,ny,k,ie) = 0.
            enddo
            enddo
          endif
      enddo
      call copy(w, wv, n)
      call adfmg_dssum(w,l) ! summed up

      do i=1,nx*ny*nz*nelt
        if(w(i,1,1,1).lt.1e-7) then  ! one elem, right=outflow(t field)
          write(6,*) 'Error in rst mask',i,w(i,1,1,1),', set to 1.'
          w (i,1,1,1) = 1.
          wv(i,1,1,1) = 1.
        endif
        w(i,1,1,1) = 1./w(i,1,1,1)
      enddo

      call col2(wv,w,n) ! zero remains zero, otherwize 1/multiplic

      return
      end
c----------------------------------------------------------------------
      subroutine adfmn_setup_mask_rt(wv,nx,ny,nz,l,w) ! No 3D
      include 'SIZE'
      include 'INPUT'
      include 'SOLN' ! for outpost debugging only, remove later!
      integer nx,ny,nz,l
      real    w (nx,ny,nz,nelt)
      real    wv(nx,ny,nz,nelt)
      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      if(ldim.eq.2) call adfmn_setup_mask_rt_2(wv,nx,ny,nz,l,w)
      if(ldim.eq.3) call adfmn_setup_mask_rt_3(wv,nx,ny,nz,l,w)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfax(ax,xx,h1,h2,h3)
c
c     apply advection diffusion operator, ax = A xx
c     This routine takes care of setting the convection field to be
c     just the velocity from SOLN, vxd is then used by tens.f
c
      include 'SIZE'
      include 'TOTAL'
      include 'ANISO' ! cx, kx, kxd

      real     ax(*), xx(*), h1(*), h2(*), h3(*)


      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV, ifuf, ifcf


      integer icald
      data    icald / 0 /
      save    icald

      nt = nx1*ny1*nz1*nelt

      ! not dealiased
c     call convop(ax,xx) ! adv, physical space
c     call col2  (ax,bm1,nt)

      if (icald.eq.0) then 
         call copy(cx,vx,nt)
         call copy(cy,vy,nt)
         call copy(cz,vz,nt)
         call set_convect_new(vxd,vyd,vzd,cx,cy,cz)
         call set_anisotr_new(kxd,kyd,kzd,kx,ky,kz)
         icald = 1
      endif

      ifuf = .false.
      ifcf = .true.
      call rzero(ax,nt) !Lan
      call convect_new(ax,xx,ifuf,vxd,vyd,vzd,ifcf)

c     call rzero (ax,nt) ! turn on so that only do Laplacian

c     do ie=1,nelt
c       write(6,*) ie,ifdfrm(ie), iffast(ie),'logic ie,'
c     enddo

      imesh = 2
      isd = 1

      call axhelm_new(ax,xx,h1,h2,h3)
      call dssum (ax,lx1,ly1,lz1)
      call col2  (ax,tmask(1,1,1,1,1),nt)
      return
      end
c-----------------------------------------------------------------------
c---- Util functions
c-----------------------------------------------------------------------
      subroutine adf_proj(x,res,h1,h2,h3,iter,maxit,wt)

c     Solve A t = r, via custome projection scheme

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      parameter (lt=lx1*ly1*lz1*lelt)
      real x(*), res(*), h1(*) , h2(*), h3(*)
      integer iter, maxit, itall

      real r(lt),w(lt),v(lt),z(lt), wt(lt) ! tmult
      common /cadfpj/ r,w     ! WORK ARRAYS
      common /zadfp1/ v
      common /zadfp2/ z
      common /zadfps/ p,q
      real            p(lt,lgmres), q(lt,lgmres)

      real   x0(lt), er(lt)

      real    proj_tol
      real    alpha, beta

      common /adfexc/ uex
      real            uex(lt)


      ngo = -99        ! TURN OFF VERBOSE OUTPUT
      n = nx1*ny1*nz1*nelt
      proj_tol = 1.e-10 ! staganates at 5e-9
c     proj_tol = 1.e-12 ! staganates at 5e-9

      if (istep.le.20.or.mod(istep,10).eq.0) ngo = nio


      m     = lgmres
      mxs = 1 ! maximum schwarz iter - once?
      mxmg = 500 ! max MG
      iconv  = 0 ! if converged
      sigma = 0.5

      ! initial guess
      call rzero(x,n)                     ! Solution
      call rzero(x0,n)           ! all 0 initial guess

      iter = 1
      itall = 0

      do while (iter.lt.maxit)             ! Main proj loop
         if(iter.eq.1) then                ! first step :: r = res - A 0
            call copy   (r,res,n)
         else                              ! second up  :: r = res - A x
            call adfax  (w,x,h1,h2,h3)     ! w = A x
            call sub3   (r,res,w,n)        ! r = r - w
         endif

         ! preconditioner solve
                                           !       -1
c        call copy  (z,r,n)                ! z  = M  r

         img = 0
         mxmg = 1
         if (param(51).eq.0) then
           call col3(z,bintm1,r,n)
c          call copy(z,r,n)
         else
           call adfmn_mg_schwz_tens(z,r,h1,h2,h3,img,mxmg,wt,.false.)
c          call copy(z,r,n)  ! turn off preconditioner
         endif
         itall = itall + img

         call adfax    (w,z,h1,h2,h3)      !w = A z
         do j=1,(iter-1)
           beta = glsc3(q(1,j),w,wt,n)
           call add2s2(z,p(1,j),-beta,n)
           call add2s2(w,q(1,j),-beta,n)
         enddo
         beta = sqrt(glsc3(w,w,wt,n))
         if(beta.lt.proj_tol) goto 900
         betai = 1./beta
         call copy (p(1,iter),z,n)
         call cmult(p(1,iter),betai,n)
         call copy (q(1,iter),w,n)
         call cmult(q(1,iter),betai,n)

         alpha = glsc3(q(1,iter),r,wt,n)
         call add2s2(x,p(1,iter),alpha,n)

         if (ngo.eq.0) write(6,9)
     $         n,iter,beta,alpha,proj_tol,dt
    9    format(i9,i5,1p4e12.4,' tproj')
         istep = iter
         time  = real(iter)
         iter = iter + 1
         tmin = glmin(x,n)
         tmax = glmax(x,n)
         if (nio.eq.0) write (6,*) 'iter,tmin,tmax:',iter,tmin,tmax
         if(mod(iter,max(iostep,1)).eq.0) then
c          call sub3(er,x,uex,n)
c          call col2(er,tmask,n)
           call outpost(vx,vy,vz,pr,x,'   ')
c          call exitt
         endif

      enddo
  900 continue
      maxit = iter
      if (nio.eq.0) write(6,8)
     $   n,maxit,itall,beta,proj_tol,dt
    8    format(i9,i5,i8,1p3e12.4,' tprojb')

      istep = iter
      time  = real(iter)
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_gmres(x,res,h1,h2,h3,iter,maxit,wt)

c     Solve A t = r, via gmres

      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      real x(*), res(*), h1(*), h2(*), h3(*)
      integer iter, maxit

      real x0(lt)

      real r(lt),w(lt),d(lt),v(lt,lgmres),z(lt,lgmres)
     $   , wt(lt) ! multiplicity
      common /dgmres/ d       ! Preconditioner
      common /cgmres/ r,w     ! WORK ARRAYS
      common /egmre1/ v
      common /egmre2/ z

      real    gmtol_in
      real    alpha, l, temp
      common /gctmp/  h(lgmres,lgmres), gamma(lgmres+1)
     $              , c(lgmres), s(lgmres), wk1(lgmres)

      integer icalld
      save    icalld
      data    icalld / 0 / ! only do this once

      ngo = -99        ! TURN OFF VERBOSE OUTPUT
      nt = nx1*ny1*nz1*nelt
      gmtol_in = 1.e-12 !

      if (istep.le.20.or.mod(istep,10).eq.0) ngo = nio

      if (icalld.lt.0) call rone(d,nt)
      if (icalld.eq.0) then ! formulate preconditioner once
         call rone (d, nt)
      endif
      icalld=icalld+1

      m     = lgmres
      mxs = 1 ! maximum schwarz iter - once?
      iconv  = 0 ! if converged
      call rzero(x,nt)                     ! Solution
      call rzero(x0,nt)                     ! Solution

      do while (iconv.eq.0.and.iter.lt.maxit)    ! Main GMRES loop
         if(iter.eq.0) then                      ! first step :: r = res - A 0
            call copy   (r,res,nt)
         else                                     ! second up  :: r = res - A x
            call adfax  (w,x,h1,h2,h3)            ! w = A x
            call sub3   (r,res,w,nt)              ! r = r - w
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc3(r,r,wt,nt))        ! gamma  = \/ (r,r)

         if (gamma(1).eq.0.) goto 9000            ! check for lucky convergence

         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,nt)            ! v  = r / gamma
                                                  !  1            1
         do j=1,m ! Krylov space
            iter = iter+1
                                                  !       -1
            call col3  (z(1,j),d,v(1,j),nt)       ! z  = M  v
                                                  !  j       j
            isw = 0
c           call adf_schwz(z(1,j),x0,v(1,j),h1,isw,mxs,wt) ! works

            call adfax(w,z(1,j),h1,h2,h3)         !w = A z
                                                  !       j
            do i=1,j                              ! Gram-Schmidt
               h(i,j)=vlsc3(w,v(1,i),wt,nt)       ! h    = (w,v )
            enddo                                 !  i,j       i
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs
            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),nt)   ! w = w - h    v
            enddo                                 !          i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                  !            ______
            alpha = sqrt(glsc3(w,w,wt,nt))        ! alpha =  \/ (w,w)

            if (alpha.eq.0.) goto 900             ! converged

            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l

            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)
            rnorm = abs(gamma(j+1))

            if (rnorm .lt. gmtol_in) goto 900     !converged
            ! if (alpha .lt. gmtol_in) goto 900     !converged

            if (j.eq.m) goto 1000                 ! not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,nt)       ! v    = w / alpha

            if (ngo.eq.0) write(6,9)
     $         istep,iter,j,rnorm,gmtol_in,dt
    9       format(i9,2i5,1p3e12.4,' tgmrs')
         enddo
  900    iconv = 1
 1000    continue

         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,z(1,i),c(i),nt)  ! x = x + c  z
                                           !          i  i
         enddo
      enddo
 9000 continue
      maxit = iter
      if (nio.eq.0) write(6,8)
     $   istep,maxit,j,rnorm,gmtol_in,dt
    8    format(i9,2i5,1p3e12.4,' tgmrsb')

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_setup_wtmask ! h1mg_setup_wtmask
      include 'SIZE'
      include 'HSMGL'
      integer i,l
      i = mg_mask_index(mg_lmax,mg_fld-1)
      do l=1,mg_lmax
         mg_rstr_wt_index(l,mg_fld)=i
         mg_mask_index   (l,mg_fld)=i
         i=i+mg_nh(l)*mg_nhz(l)*2*ldim*nelt
         if(i .gt. lmgs*lmg_rwt*2*ldim*lelt) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
            call exitt
         endif
         call adfmg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
         call adfmg_setup_mask_hs(
     $           mg_mask(mg_mask_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      mg_mask_index(l,mg_fld)=i
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_rstr_wt(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelt)
      real wt(nx,nz,2,ldim,nelt)
      
      integer ie
      !init border nodes to 1
      call rzero(w,nx*ny*nz*nelt)
c     print *, 'Setup rstr wt: ',nx,ny,nz,nelt
      if (.not.if3d) then
         do ie=1,nelt
            do i=1,nx
               w(i,1,1,ie)=1.0
               w(i,ny,1,ie)=1.0
            enddo
            do j=1,ny
               w(1,j,1,ie)=1.0
               w(nx,j,1,ie)=1.0
            enddo
         enddo
      else
         do ie=1,nelt
            do j=1,ny
            do i=1,nx
               w(i,j,1,ie)=1.0
               w(i,j,nz,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=1.0
               w(i,ny,k,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=1.0
               w(nx,j,k,ie)=1.0
            enddo
            enddo
         enddo
      endif
      call adfmg_dssum(w,l)
      !invert the count w to get the weight wt
      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               wt(j,1,1,1,ie)=1.0/w(1,j,1,ie)
               wt(j,1,2,1,ie)=1.0/w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=1.0/w(i,1,1,ie)
               wt(i,1,2,2,ie)=1.0/w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=1.0/w(1,j,k,ie)
               wt(j,k,2,1,ie)=1.0/w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=1.0/w(i,1,k,ie)
               wt(i,k,2,2,ie)=1.0/w(i,ny,k,ie)
            enddo
            enddo
            do j=1,ny
            do i=1,nx
               wt(i,j,1,3,ie)=1.0/w(i,j,1,ie)
               wt(i,j,2,3,ie)=1.0/w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_mask_hs(wt,nx,ny,nz,l,w) ! hs version
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelt)
      real wt(nx,nz,2,ldim,nelt)
      
      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
c     init everything to 1

      n = nx*ny*nz*nelt
      call rone(w,n)

c     set dirichlet nodes to zero
      ierr = 0
      two  = 2
      do ie=1,nelt
         call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
         if (ierr.ne.0) then
            ierr = -1
            call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
         endif

         if(lbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(nx,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(lbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,ny,k,ie)=0.0
            enddo
            enddo
         endif
         if(if3d) then
            if(lbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,1,ie)=0.0
               enddo
               enddo
            endif
            if(rbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,nz,ie)=0.0
               enddo
               enddo
            endif
         endif
      enddo
c     do direct stiffness multiply
      call adfmg_dsprod(w,l)

c     store weight
      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               wt(j,1,1,1,ie)=w(1,j,1,ie)
               wt(j,1,2,1,ie)=w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=w(i,1,1,ie)
               wt(i,1,2,2,ie)=w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=w(1,j,k,ie)
               wt(j,k,2,1,ie)=w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=w(i,1,k,ie)
               wt(i,k,2,2,ie)=w(i,ny,k,ie)
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               wt(j,k,1,3,ie)=w(i,j,1,ie)
               wt(j,k,2,3,ie)=w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('B INVALID BC FOUND in genfast$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_index_0 ! initialize index sets
      include 'SIZE'
      include 'HSMGL'

      n = lmgn*(lmgs+1)

      call izero( mg_rstr_wt_index      , n )
      call izero( mg_mask_index         , n )
      call izero( mg_solve_index        , n )
      call izero( mg_fast_s_index       , n )
      call izero( mg_fast_d_index       , n )
      call izero( mg_schwarz_wt_index   , n )
      
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_solve
      include 'SIZE'
      include 'HSMGL'
      
      integer l,i,nl,nlz
      i = mg_solve_index(mg_lmax+1,mg_fld-1)
      do l=1,mg_lmax
         mg_solve_index(l,mg_fld)=i
         i=i+mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
         if(i .gt. lmg_solve*lelt) then
            itmp = i/lelv
            write(6,*) 'lmg_solve too small',i,itmp,lmg_solve,l
            call exitt
         endif
      enddo
      mg_solve_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_dssum
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMGL'
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      common /c_is1/ glo_num(lxyz*lelt)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      integer l
      
      ncrnr = 2**ldim
      call get_vert

      do l=1,mg_lmax  ! set up direct stiffness summation for each level
         nx=mg_nh(l)
         ny=mg_nh(l)
         nz=mg_nhz(l)
         call setupds(mg_gsh_handle(l,mg_fld),nx,ny,nz
     $                ,nelt,nelgt,vertex,glo_num)
c        nx=nx+2
c        ny=ny+2
c        nz=nz+2
c        if(.not.if3d) nz=1
c        call setupds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz
c    $                ,nelt,nelgt,vertex,glo_num)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_intp
      include 'SIZE'
      include 'HSMGL'
      include 'SEMHAT'
      integer l,nf,nc

      do l=1,mg_lmax-1

         nf=mg_nh(l+1)
         nc=mg_nh(l)

!        Standard multigrid coarse-to-fine interpolation
         call adfmg_setup_intpm(
     $           mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
         call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)

!        Fine-to-coarse interpolation for variable-coefficient operators
         call adfmg_setup_intpm(
     $           mg_jhfc(1,l),mg_zh(1,l),mg_zh(1,l+1),nc,nf)
         call transpose(mg_jhfct(1,l),nf,mg_jhfc(1,l),nc)
c        call outmat(mg_jhfc(1,l),nc,nf,'MG_JHFC',l)

      enddo
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_intpm(jh,zf,zc,nf,nc)
      integer nf,nc
      real jh(nf,nc),zf(1),zc(1)
      include 'SIZE'
      real w(2*lx1+2)
      do i=1,nf
         call fd_weights_full(zf(i),zc,nc-1,1,w)
         do j=1,nc
            jh(i,j)=w(j)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_semhat ! SEM hat matrices for each level
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      include 'SEMHAT'

      do l=1,mg_h1_lmax
         n = mg_nx(l)     ! polynomial order
         call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zglhat,dgl,jgl,n,wh)
         call copy(mg_ah(1,l),ah,(n+1)*(n+1))
         call copy(mg_bh(1,l),bh,n+1)
c        write(6,*) 'bh  ',(bh(i),i=1,n+1)
c        write(6,*) 'mgbh',(mg_bh(i,l),i=1,n+1)
         call copy(mg_dh(1,l),dh,(n+1)*(n+1))
         call transpose(mg_dht(1,l),n+1,dh,n+1)
         call copy(mg_zh(1,l),zh,n+1)
         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1
      enddo
      end
c-----------------------------------------------------------------------
      subroutine adfmg_set_cmlt
      implicit none
      include 'SIZE'
      include 'HSMGL'

      common /adfwt/ cmlt(lx1*ly1*lz1*lelt)
      real cmlt

      integer l,n

      l = 1 ! coarse level
      n = mg_h1_n(l,mg_fld)
      call rone(cmlt,n)
      call adfmg_dssum(cmlt,l)
      call invcol1(cmlt,n)
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_set_vec(co, cx,cy,cz) ! advection coefficients
      implicit none
      include 'SIZE'
      include 'HSMGL'
      include 'MASS'     ! bm1
      include 'TSTEP'    ! nelfld
      include 'PARALLEL' ! nid

      common /ctmp1/ w
      real w(lx1*ly1*lz1*lelt*2)

      common /adfcoi/ pco
      integer pco(lmgx,ldimt1) ! big enough?

      real co(*), cx(*), cy(*), cz(*)

      integer e,l
      integer n,nc
      integer nx,ny,nz,nxyz
      integer nxl,nyl,nzl
      integer p_c, l_c

c     write(6,*) 'adfmg set co: ',mg_lmax
      l                 = mg_lmax
      pco    (l,mg_fld) = 1 ! NOT 0 ! ! ! !
      n                 = mg_h1_n(l,mg_fld)

      nc = ldim  ! 2 or 3 elements to adv cr,cs,ct tensor

      do l=mg_lmax-1,1,-1
         pco(l,mg_fld) = pco    (l+1,mg_fld) + n*nc
         n             = mg_h1_n(l  ,mg_fld)
      enddo

      do e=1,nelfld(mg_fld)
       do l=mg_lmax,1,-1

         nx = mg_nh(l)
         ny = mg_nh(l)
         nz = mg_nhz(l)
         nxyz = nx*ny*nz
         p_c = pco(l,mg_fld) + nc*nx*ny*nz*(e-1)

         if (l.eq.mg_lmax) then
           ! fine grid, evaluate cr=cx rx + cy ry + cz rz
           call pack_rst(e,nc, co(p_c), cx,cy,cz)

           call adfmg_scale_mass_co     ! Divide out Weights for intp
     $          (co(p_c),mg_bh(1,l),nc,nx,ny,nz,w,.true.)
         else
c           Generate G and B by interpolating their continous counterparts onto
c           the coarse grid and collocating with coarse-grid quadrature weights
            call adfmg_intp_cfc_e
     $           (co(p_c),co(l_c),nc,nx,ny,nz,nxl,nyl,nzl,e,l,w)
            call adfmg_scale_mass_co       ! Multiply by weights
     $           (co(l_c),mg_bh(1,l+1),nc,nxl,nyl,nzl,w,.false.)
         endif
         l_c = p_c ! lag
         nxl = nx
         nyl = ny
         nzl = nz
       enddo
       call adfmg_scale_mass_co     ! Multiply by weights
     $      (co(l_c),mg_bh(1,1),nc,nxl,nyl,nzl,w,.false.)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_set_co(p_c,l0) ! advection coefficients
c     Assumes that cx,cy,cz have been set manually or from a previous 
c     call to adfax, in which that case were set to vx,vy,vz
      implicit none
      include 'SIZE'
      include 'HSMGL'
      include 'ANISO'

      common /adfco/  co
      real            co(lmg_g*ldim*lelt) ! TODO mg_co

      common /adfcoi/ pco
      integer pco(lmgx,ldimt1) ! big enough?

      integer p_c, l0

      call adfmg_set_vec(co,cx,cy,cz)
      p_c  = pco(l0,mg_fld)
      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_set_gb  (p_g,p_b,l0)
      include 'SIZE'
      include 'HSMGL'
      include 'MASS'   ! bm1
      include 'TSTEP'  ! nelfld

      include 'PARALLEL'  ! nid

      integer p_g,p_b,e
      common /ctmp1/ w(lx1*ly1*lz1*lelt*2)

      l                 = mg_h1_lmax
      p_mg_b (l,mg_fld) = 0
      p_mg_g (l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)


      ng = 3*(ldim-1)  ! 3 or 6 elements to symm dxd tensor

      do l=mg_h1_lmax-1,1,-1
         p_mg_b (l,mg_fld) = p_mg_b (l+1,mg_fld) + n
         p_mg_g (l,mg_fld) = p_mg_g (l+1,mg_fld) + n*ng
         n                 = mg_h1_n(l  ,mg_fld)
      enddo

c     write(6,*) 'adfmg set gb :: ifld, nel',ifield,nelfld(ifield)
      do e=1,nelfld(ifield)
       do l=mg_h1_lmax,1,-1

         nx = mg_nh(l)
         ny = mg_nh(l)
         nz = mg_nhz(l)
         nxyz = nx*ny*nz
         p_g = p_mg_g (l,mg_fld) + ng*nx*ny*nz*(e-1)
         p_b = p_mg_b (l,mg_fld) +    nx*ny*nz*(e-1)
         if (l.eq.mg_h1_lmax) then
            call adf_gxfer_e(mg_g(p_g),ng,e)            ! Fine grid=original G
            call copy(mg_b(p_b) ,bm1(1,1,1,e),nxyz)     ! Fine grid=original B
            call adfmg_scale_mass                       ! Divide out Wghts
     $         (mg_b(p_b),mg_g(p_g),mg_bh(1,l),ng,nx,ny,nz,w,.true.)
         else
c        Generate G and B by interpolating their continous counterparts onto
c        the coarse grid and collocating with coarse-grid quadrature weights

            call adfmg_intp_gfc_e
     $            (mg_g(p_g),mg_g(l_g),ng,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call adfmg_intp_fc_e
     $            (mg_b(p_b),mg_b(l_b)   ,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call adfmg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,l+1),ng,nxl,nyl,nzl,w,.false.)

         endif

         l_b = p_b
         l_g = p_g

         nxl = nx
         nyl = ny
         nzl = nz

       enddo

       call adfmg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,1),ng,nxl,nyl,nzl,w,.false.)

      enddo

      p_b  = p_mg_b (l0,mg_fld)
      p_g  = p_mg_g (l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_intp_fc_e(uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMGL'

      real uf(nxf,nyf,nzf),uc(nxc,nyc,nzc),w(1)

      if (if3d) then

         n1=nxf*nyf
         n2=nzf
         n3=nzc
         call mxm(uf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()
         lc=1 + n1*n3
         lc0=lc

         n1=nxf
         n2=nyf
         n3=nyc

         do k=1,nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,w(lc),n3)
            lf = lf + n1*n2
            lc = lc + n1*n3
         enddo

         lf=lc0  ! Rewind fine pointer to start of coarse data
         n1=nxc
         n2=nxf
         n3=nyc*nzc
         call mxm(mg_jhfc(1,l),n1,w(lf),n2,uc,n3)

      else ! 2D

         n1=nxf
         n2=nyf
         n3=nyc
         call mxm(uf,n1,mg_jhfct(1,l),n2,w,n3)

         n1=nxc
         n2=nxf
         n3=nyc
         call mxm(mg_jhfc(1,l),n1,w,n2,uc,n3)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_intp_gfc_e(gc,gf,ng,nxc,nyc,nzc,nxf,nyf,nzf
     $                           ,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMGL'

      real gf(ng,nxf,nyf,nzf),gc(ng,nxc,nyc,nzc),w(1)


      if (if3d) then

         n1=ng*nxf*nyf
         n2=nzf
         n3=nzc
         call mxm(gf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()
         lc=1 + n1*n3
         lc0=lc

         n1=ng*nxf
         n2=nyf
         n3=nyc

         do k=1,nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,w(lc),n3)
            lf = lf + n1*n2
            lc = lc + n1*n3
         enddo

         lf=lc0  ! Rewind fine pointer to start of coarse data
         n1=ng
         n2=nxf
         n3=nxc

         do k=1,nyc*nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,gc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      else ! 2D

         n1=ng*nxf
         n2=nyf
         n3=nyc
         call mxm(gf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()

         n1=ng
         n2=nxf
         n3=nxc

         do k=1,nyc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,gc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_scale_mass_co (co,wt,nc,nx,ny,nz,wk,ifinv)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'HSMGL'

      real    co(nc,1),wt(1),wk(1)
      logical ifinv

      common /ctmp0/ wi(2*lx1+4)

      n = nx*ny*nz

      if (nx.le.2*lx1) then
        if (ifinv) then ! divide weights from b, g
            call invers2(wi,wt,nx)
         else           ! multiply weights into b, g
            call copy(wi,wt,nx)
         endif
      else
         call exitti('adfmg_scale_mass_co: wi too small$',nx)
      endif

      if (ldim.eq.3) then
         l=0
         do k=1,nz
         do j=1,ny
            wjk=wi(j)*wi(k)
            do i=1,nx
               l=l+1
               wk(l) = wjk*wi(i)
            enddo
         enddo
         enddo

         do k=1,n
            co(1,k) = wk(k)*co(1,k)
            co(2,k) = wk(k)*co(2,k)
            co(3,k) = wk(k)*co(3,k)
         enddo

      else      ! 2D
         l=0
         do j=1,ny
         do i=1,nx
            l=l+1
            wk(l) = wi(i)*wi(j)
         enddo
         enddo

         do k=1,n
            co(1,k) = wk(k)*co(1,k)
            co(2,k) = wk(k)*co(2,k)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_scale_mass (b,g,wt,ng,nx,ny,nz,wk,ifinv)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'HSMGL'

      real b(1),g(ng,1),wt(1),wk(1)
      logical ifinv

      common /ctmp0/ wi(2*lx1+4)

      n = nx*ny*nz

      if (nx.le.2*lx1) then
        if (ifinv) then ! divide weights from b, g
            call invers2(wi,wt,nx)
         else           ! multiply weights into b, g
            call copy(wi,wt,nx)
         endif
      else
         call exitti('mg_scale_mass: wi too small$',nx)
      endif

      if (if3d) then
         l=0
         do k=1,nz
         do j=1,ny
            wjk=wi(j)*wi(k)
            do i=1,nx
               l=l+1
               wk(l) = wjk*wi(i)
            enddo
         enddo
         enddo

         do k=1,n
            b(k)   = wk(k)*b(k)
            g(1,k) = wk(k)*g(1,k)
            g(2,k) = wk(k)*g(2,k)
            g(3,k) = wk(k)*g(3,k)
            g(4,k) = wk(k)*g(4,k)
            g(5,k) = wk(k)*g(5,k)
            g(6,k) = wk(k)*g(6,k)
         enddo

      else      ! 2D
         l=0
         do j=1,ny
         do i=1,nx
            l=l+1
            wk(l) = wi(i)*wi(j)
         enddo
         enddo

         do k=1,n
            b(k)   = wk(k)*b(k)
            g(1,k) = wk(k)*g(1,k)
            g(2,k) = wk(k)*g(2,k)
            g(3,k) = wk(k)*g(3,k)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_intp_cfc_e
     $           (cc,cf,nc,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMGL'

      real cf(nc,nxf,nyf,nzf),cc(nc,nxc,nyc,nzc),w(1)

      if (if3d) then

         n1=nc*nxf*nyf
         n2=nzf
         n3=nzc
         call mxm(cf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()
         lc=1 + n1*n3
         lc0=lc

         n1=nc*nxf
         n2=nyf
         n3=nyc

         do k=1,nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,w(lc),n3)
            lf = lf + n1*n2
            lc = lc + n1*n3
         enddo

         lf=lc0  ! Rewind fine pointer to start of coarse data
         n1=nc
         n2=nxf
         n3=nxc

         do k=1,nyc*nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,cc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      else ! 2D

         n1=nc*nxf
         n2=nyf
         n3=nyc
         call mxm(cf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()

         n1=nc
         n2=nxf
         n3=nxc

         do k=1,nyc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,cc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pack_rst(e,nc,cf,cx,cy,cz)  ! mesh 1
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'WZ'

      integer e, nc
      real cf(nc,*)
      real cx(lx1,ly1,lz1,*)
      real cy(lx1,ly1,lz1,*)
      real cz(lx1,ly1,lz1,*)

      integer i,nxyz
      real ux,uy,uz

      nxyz = lx1*ly1*lz1
      if(ldim.eq.3) then
       do i=1,nxyz
        ux = cx(i,1,1,e)
        uy = cy(i,1,1,e)
        uz = cz(i,1,1,e)
        cf(1,i) = w3m1(i,1,1) *
     $          (ux*rxm1(i,1,1,e) + uy*rym1(i,1,1,e) + uz*rzm1(i,1,1,e))
        cf(2,i) = w3m1(i,1,1) *
     $          (ux*sxm1(i,1,1,e) + uy*sym1(i,1,1,e) + uz*szm1(i,1,1,e))
        cf(3,i) = w3m1(i,1,1) *
     $          (ux*txm1(i,1,1,e) + uy*tym1(i,1,1,e) + uz*tzm1(i,1,1,e))
       enddo
      else
       do i=1,nxyz
        ux = cx(i,1,1,e)
        uy = cy(i,1,1,e)
        cf(1,i) = w3m1(i,1,1) *(ux*rxm1(i,1,1,e) + uy*rym1(i,1,1,e))
        cf(2,i) = w3m1(i,1,1) *(ux*sxm1(i,1,1,e) + uy*sym1(i,1,1,e))
       enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_gxfer_e (g,ng,e) 
      include 'SIZE'
      include 'TOTAL'
      include 'ANISO'

      common /adfcoef/ h1,h2,h3
      real h1(lx1,ly1,lz1,lelt)
      real h2(lx1,ly1,lz1,lelt)
      real h3(lx1,ly1,lz1,lelt) ! FIXME should get h3 as input
      
      real g(ng,*)
      integer e

      nxyz = lx1*ly1*lz1

c     ifdfrm(e) = .true.  ! TOO LATE

      if (if3d) then
         do i=1,nxyz
            w1 = kxd(i,1,1,e)
            w2 = kyd(i,1,1,e)
            w3 = kzd(i,1,1,e)
            w0 = bm1(i,1,1,e)*h3(i,1,1,e)/h1(i,1,1,e)
            g(1,i) = g1m1(i,1,1,e) + w0*w1*w1
            g(2,i) = g2m1(i,1,1,e) + w0*w2*w2
            g(3,i) = g3m1(i,1,1,e) + w0*w3*w3
            g(4,i) = g4m1(i,1,1,e) + w0*w1*w2
            g(5,i) = g5m1(i,1,1,e) + w0*w1*w3
            g(6,i) = g6m1(i,1,1,e) + w0*w2*w3
         enddo
      else
         do i=1,nxyz
            w1 = kxd(i,1,1,e)
            w2 = kyd(i,1,1,e)
            w0 = bm1(i,1,1,e)*h3(i,1,1,e)/h1(i,1,1,e)
            g(1,i) = g1m1(i,1,1,e) + w0*w1*w1
            g(2,i) = g2m1(i,1,1,e) + w0*w2*w2
            g(3,i) = g4m1(i,1,1,e) + w0*w1*w2
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_set_msk(p_msk ,l0) ! boundary mask
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'
      integer p_msk

      l                  = mg_h1_lmax
      p_mg_msk(l,mg_fld) = 0
      n                  = mg_h1_n(l,mg_fld)


      do l=mg_h1_lmax,1,-1
         nx = mg_nh  (l)
         ny = mg_nh  (l)
         nz = mg_nhz (l)

         p_msk = p_mg_msk(l,mg_fld)

c        write(6,*) 'inside adfmg_set_msk:::'
         call adfmg_setup_mask ! h1mg_setup_mask
     $     (mg_imask(p_msk),nm,nx,ny,nz,nelfld(mg_fld),l,mg_work)

         if (l.gt.1) p_mg_msk(l-1,mg_fld)=p_mg_msk(l,mg_fld)+nm

c        write(6,*) p_msk,p_mg_msk(l-1,mg_fld),p_mg_msk(l,mg_fld),'M'

      enddo

      p_msk = p_mg_msk(l0,mg_fld)

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_setup_mask(mask,nm,nx,ny,nz,nel,l,w) ! h1 version
      include 'SIZE'
      include 'INPUT'        ! if3d
      include 'TSTEP'        ! ifield

      integer mask(1)        ! Pointer to Dirichlet BCs
      integer nx,ny,nz,l
      real w(nx,ny,nz,nel)
      
      integer e,count,ptr
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      zero = 0
      nxyz = nx*ny*nz
      n    = nx*ny*nz*nel

      call rone(w,n)   ! Init everything to 1

      ierr = 0
      ierrmx = 0       ! BC verification
      two    = 2
      ifield = 2
      do e=1,nel       ! Set dirichlet nodes to zero

         call adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,two,ierr)

c        write(6,6) e,lbr,rbr,lbs,rbs,ierr,nx
c   6    format(i5,2x,4i3,2x,i2,3x,i5,'  lbr,rbr,lbs')

         if (lbr.eq.1) call facev(w,e,4,zero,nx,ny,nz)
         if (rbr.eq.1) call facev(w,e,2,zero,nx,ny,nz)
         if (lbs.eq.1) call facev(w,e,1,zero,nx,ny,nz)
         if (rbs.eq.1) call facev(w,e,3,zero,nx,ny,nz)
         if (if3d) then
            if (lbt.eq.1) call facev(w,e,5,zero,nx,ny,nz)
            if (rbt.eq.1) call facev(w,e,6,zero,nx,ny,nz)
         endif
         ierrmx = max(ierrmx,ierr)
      enddo

c     call hsmg_dsprod(w,l)    ! direct stiffness multiply
      call adfmg_dsprod(w,l)    ! direct stiffness multiply

c     write(6,*) 'w l ',(w(i,1,1,1),i=1,nx*ny*nel)
c
c     Prototypical mask layout, nel=5:
c
c    e=1 ...             10
c      1  2  3  4  5 ... 10 | 11 12 13 14 | 15 | 16 |
c     11 15 16 ...          |  3 p1 p2 p3 |  0 |  0 | ...
c                              ^
c                              |
c                              |_count for e=1
c

      nm  = 1                  ! store mask
      do e=1,nel

         mask(e) = nel+nm
         count   = 0          ! # Dirchlet points on element e
         ptr     = mask(e)

         do i=1,nxyz
            if (w(i,1,1,e).eq.0) then
               nm    = nm   +1
               count = count+1
               ptr   = ptr  +1
               mask(ptr) = i + nxyz*(e-1)   ! where I mask on element e 
            endif
         enddo


         ptr       = mask(e)
         mask(ptr) = count

         nm        = nm+1     ! bump pointer to hold next count

      enddo

      nm = nel + nm-1 ! Return total number of mask pointers/counters

      ierrmx = iglmax(ierrmx,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nio,ierr,' BC FAIL h1'
         call exitti('D INVALID BC FOUND in adfmg_setup_mask$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_set_h1  (p_h1 ,l0)
      include 'SIZE'
      include 'HSMGL'
      include 'INPUT'

      integer pf,pc
      integer p_h1

      l                 = mg_h1_lmax
      p_mg_h1(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

c     call copy(mg_h1,h1,n)   ! Fine grid is just original h1
      call svisc(mg_h1,2)     ! Override with prec h1 Lan


      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h1(l,mg_fld) = p_mg_h1(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h1(l+1,mg_fld)
         pc                = p_mg_h1(l  ,mg_fld)

         call adfmg_intp_fc (mg_h1(pc),mg_h1(pf),l)

      enddo

      p_h1 = p_mg_h1(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_set_h2  (p_h2 ,l0)
      include 'SIZE'
      include 'HSMGL'

c     As a first pass, rely on the cheesy common-block interface to get h2
      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)

      integer p_h2,pf,pc

      l                 = mg_h1_lmax
      p_mg_h2(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

      call copy (mg_h2,h2,n)   ! Fine grid is just original h2

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h2(l,mg_fld) = p_mg_h2(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h2(l+1,mg_fld)
         pc                = p_mg_h2(l  ,mg_fld)

         call adfmg_intp_fc (mg_h2(pc),mg_h2(pf),l)

      enddo

      p_h2 = p_mg_h2(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_setup_nx(lv)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      mg_fld = ifield                    ! FIELD pointer

      mg_h1_lmax        = lv
      mg_lmax           = mg_h1_lmax

      ! test two level
c     mg_h1_lmax        = 2
c     mg_lmax           = mg_h1_lmax

      if  ( param(113).gt.0.5 ) then
         ! Schedule: divide by 2
         mg_nx(mg_h1_lmax) = lx1-1        ! Polynomial degree 
         do k = mg_h1_lmax-1, 1,-1       ! of multigrid levels
            mg_nx(k) = (mg_nx(k+1)+1)/2-1
         enddo
      else
         ! Schedule: minus 2
         mg_nx(mg_h1_lmax) = lx1-1        ! Polynomial degree 
         do k = mg_h1_lmax-1, 1,-1       ! of multigrid levels
            mg_nx(k) = mg_nx(k+1)-2
         enddo
      endif
      if (nint(param(115)).gt.0.5) mg_nx(1)=1
      if (nint(param(115)).gt.1.5.AND.lv.gt.2) mg_nx(2)=3
c      mg_nx(2) = 2 ! test
c      mg_nx(1) = 1 ! Pablo: Always corner grid as coarsest level, Lan
      do k = mg_h1_lmax, 1,-1
        if(nio.eq.0) write(*,*)'prec-mg schedule',k,mg_nx(k)
      enddo

      call icopy(mg_ny, mg_nx, mg_h1_lmax)
      call icopy(mg_nz, mg_nx, mg_h1_lmax)
      if (ldim.eq.2) call izero(mg_nz, mg_h1_lmax)

      if (nio.eq.0) then
         write(6, *) 'adf_mg_nx:', (mg_nx(i), i = 1, mg_h1_lmax)
         write(6, *) 'adf_mg_ny:', (mg_ny(i), i = 1, mg_h1_lmax)
         write(6, *) 'adf_mg_nz:', (mg_nz(i), i = 1, mg_h1_lmax)
      endif

      do ifld=1,ldimt1
      do l=1,mg_lmax
         mg_h1_n(l, ifld) = (mg_nx(l) + 1) * 
     $                      (mg_ny(l) + 1) *
     $                      (mg_nz(l) + 1) * 
     $                      nelfld(ifld) ! vector length on each level
      enddo
      enddo
c     if(nio.eq.0) then
c     do l=1,mg_lmax
c     write(6,*) 'adf h1 n:',mg_h1_n(l,2),'=(',mg_nx(l),'+1)^2 *'
c    $           ,nelfld(2),' ?'
c     enddo
c     endif

c     S-matrix size
c      nxxs  = 0
c      nxyzs = 0
c      do l=1,mg_lmax
c         nxx  = (mg_nx(l)+1)**2
c         nxyz = (mg_nx(l)+1)*(mg_ny(l)+1)*(mg_nz(l)+1)
c         nxxs = nxxs + nxx
c         nxyzs= nxyzs + nxyz
c      enddo
c
c      if (nio.eq.0) then
c       write(6,*) nxxs ,lmg_fasts,'  mg: S-matrix sizes (used, alloc)'
c       write(6,*) nxyzs,lmg_g    ,'  mg: G-matrix sizes (used, alloc)'
c      endif
c      if(nxyzs.gt.lmg_g)call exitti('ERROR: Increase lmg_g to:$',nxyzs)
c      if(nxxs.gt.lmg_fasts)call exitti('Increase lmg_fasts to:$',nxxs)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_axm(w,p,aw,ap,l,wk)
c
c     w  := aw*w + ap*H*p, level l, with mask and dssum
c
c     Hu := div. h1 grad u + h2 u
c
c        ~= h1 A u + h2 B u
c
c     Here, we assume that pointers into g() and h1() and h2() have
c     been established
c
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'  ! nelfld

      real w(1),p(1),wk(1)

      integer p_h1,p_h2,p_g,p_b,p_msk, p_c
      logical ifh2, ifh3

      common /adfco/  co
      real            co(lmg_g*(ldim)*lelt) ! big enough?
      common /adfcoi/ pco
      integer         pco(lmgx,ldimt1) ! big enough?

      real aw, ap

      p_h1  = p_mg_h1  (l,mg_fld)
      p_h2  = p_mg_h2  (l,mg_fld)
      p_g   = p_mg_g   (l,mg_fld)
      p_b   = p_mg_b   (l,mg_fld)
      p_msk = p_mg_msk (l,mg_fld)
c     write(6,*) 'inside adfmg axm:: p_msk',p_msk,l
      p_c   = pco      (l,mg_fld) ! right

      if (p_h1 .eq.0) call adfmg_set_h1  (p_h1 ,l)
      if (p_h2 .eq.0) call adfmg_set_h2  (p_h2 ,l)
      if (p_g  .eq.0) call adfmg_set_gb  (p_g,p_b,l)
      if (p_msk.eq.0) call adfmg_set_msk (p_msk,l)
      if (p_c  .eq.0) call adfmg_set_co  (p_c,l)

      ifh2 = .true.
      ifh3 = .false.

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)
      ng = 3*ldim-3
      nc = ldim


      call adfmg_axml (wk,p
     $    ,mg_h1(p_h1),mg_h2(p_h2),nx,ny,nz,nelfld(ifield)
     $    ,mg_g (p_g) , ng ,mg_b(p_b), mg_imask(p_msk),ifh2
     $    ,co(p_c), nc)

      ! wk = C p + A p

c     n = nx*ny*nz*nelfld(ifield)
      n = nx*ny*nz*nelfld(mg_fld)
      call adfmg_dssum  (wk,l) ! hsmg_dssum
      ! w = aw * w + ap * wk = 1 * w - 1 * wk = w - ( C p + A p)
      call add2sxy    (w,aw,wk,ap,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_axml
     $  (w,p,h1,h2,nx,ny,nz,nel,g,ng,b,mask,ifh2,c,nc)
c
c     w  := aw*w + ap*H*p, level l, with mask and dssum
c
c     Hu := div. h1 grad u + h2 u
c
c        ~= h1 A u + h2 B u
c

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      real w (nx*ny*nz,nel), p (nx*ny*nz,nel)
     $   , h1(nx*ny*nz,nel), h2(nx*ny*nz,nel)
     $   , b (nx*ny*nz,nel), g (ng*nx*ny*nz,nel)
      real c (nc*nx*ny*nz,nel)
      integer nc
      integer mask(1)

      logical ifh2

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp0/ ur(lxyz),us(lxyz),ut(lxyz)

      integer e

      do e=1,nel
         call adf_axe(w(1,e),p(1,e),h1(1,e),h2(1,e),g(1,e),ng,b(1,e)
     $            ,nx,ny,nz,ur,us,ut,ifh2,e,c(1,e),nc)
   
         im = mask(e)
         call adfmg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_axe
     $     (w,p,h1,h2,g,ng,b,nx,ny,nz,ur,us,ut,ifh2,e,c,nc) !

      include 'SIZE'
      include 'INPUT'   ! if3d
      logical ifh2

      real w (*), p (*)
     $   , h1(*), h2(*) !,  h3(*)
     $   , b (*), g (ng,*)
     $   , ur(*), us(*), ut(*)
      real c (nc,*)
      integer e

      real tmp(lx1*ly1*lz1) ! largest size for one e

      nxyz = nx*ny*nz

      call gradl_rst(ur,us,ut,p,nx,if3d)

      if(ldim.eq.3) then
        do i=1,nxyz
          tmp(i) = c(1,i)*ur(i) + c(2,i)*us(i) + c(3,i)*ut(i)
        enddo
      else ! 2D
        do i=1,nxyz
          tmp(i) = c(1,i)*ur(i) + c(2,i)*us(i)
        enddo
      endif

      if (if3d) then
         do i=1,nxyz
            wr = g(1,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
            ws = g(4,i)*ur(i) + g(2,i)*us(i) + g(6,i)*ut(i)
            wt = g(5,i)*ur(i) + g(6,i)*us(i) + g(3,i)*ut(i)
            ur(i) = wr*h1(i)
            us(i) = ws*h1(i)
            ut(i) = wt*h1(i)
         enddo
      else ! if2d
         do i=1,nxyz
            wr = g(1,i)*ur(i) + g(3,i)*us(i)
            ws = g(3,i)*ur(i) + g(2,i)*us(i)
            ur(i) = wr*h1(i)
            us(i) = ws*h1(i)
         enddo
      endif

      call gradl_rst_t(w,ur,us,ut,nx,if3d)

      ! add in advection
      call add2(w,tmp,nxyz)

      ! add in reaction
      if(ifh2) call addcol4(w,h2,b,p,nxyz)

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_schwarz_wt(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'

      real e(*)
      integer l
      
      if(.not.if3d) call adfmg_schwarz_wt2d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      if(if3d) call adfmg_schwarz_wt3d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_schwarz_wt2d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,nelt)
      real wt(n,4,2,nelt)
      
      integer ie,i,j
      do ie=1,nelt
         do j=1,n
            e(1  ,j,ie)=e(1  ,j,ie)*wt(j,1,1,ie)
            e(2  ,j,ie)=e(2  ,j,ie)*wt(j,2,1,ie)
            e(n-1,j,ie)=e(n-1,j,ie)*wt(j,3,1,ie)
            e(n  ,j,ie)=e(n  ,j,ie)*wt(j,4,1,ie)
         enddo
         do i=3,n-2
            e(i,1  ,ie)=e(i,1  ,ie)*wt(i,1,2,ie)
            e(i,2  ,ie)=e(i,2  ,ie)*wt(i,2,2,ie)
            e(i,n-1,ie)=e(i,n-1,ie)*wt(i,3,2,ie)
            e(i,n  ,ie)=e(i,n  ,ie)*wt(i,4,2,ie)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_schwarz_wt3d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelt)
      real wt(n,n,4,3,nelt)
      
      integer ie,i,j,k
      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
         enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_dssum(u,l)
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'

      include 'INPUT'

      real u(1)

      if (ifsync) call nekgsync()
      etime1=dnekclock()

      call fgslib_gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)
      tdadd =tdadd + dnekclock()-etime1
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_dsprod(u,l)
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'
      real u(1)

      if (ifsync) call nekgsync()

      call fgslib_gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)
      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_mask(w,mask,nel)
      include 'SIZE'

      real    w   (1)
      integer mask(1)        ! Pointer to Dirichlet BCs
      integer e
      
      do e=1,nel
         im = mask(e)
c        write(6,*) e,im,mask(im),'im in adfmg mask'
         call adfmg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine adfmg_intp_fc(uc,uf,l) ! l is coarse level

      include 'SIZE'
      include 'HSMGL'

      real uc(1),uf(1)


      nc = mg_nh(l)
      nf = mg_nh(l+1)
      call adfmg_tnsr(uc,nc,uf,nf,mg_jhfc(1,l),mg_jhfct(1,l))

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_intp(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMGL'
      call adfmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      return
      end
c------------------------------------------   T  -----------------------
      subroutine adfmg_rstr(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMGL'
      logical ifdssum

      real r(1)
      integer l

      call adfmg_do_wt(r,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

      call adfmg_tnsr1(r,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))

      if (ifdssum) call adfmg_dssum(r,l)

      return
      end
c----------------------------------------------------------------------
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
      subroutine adfmg_tnsr(v,nv,u,nu,A,At)
      integer nv,nu
      real v(1),u(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call adfmg_tnsr2d(v,nv,u,nu,A,At)
      else
         call adfmg_tnsr3d(v,nv,u,nu,A,At,At)
      endif
      return
      end
c----------------------------------------------------------------------
c     computes
c              T
c     v = A u B
      subroutine adfmg_tnsr2d(v,nv,u,nu,A,Bt)
      integer nv,nu
      real v(nv*nv,nelt),u(nu*nu,nelt),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work((lx1+2)*(lx1+2))
      integer ie
      do ie=1,nelt
         call mxm(A,nv,u(1,ie),nu,work,nu)
         call mxm(work,nv,Bt,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine adfmg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv,nelt),u(nu*nu*nu,nelt),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer ie, i
      do ie=1,nelt
         call mxm(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_tnsr1(v,nv,nu,A,At)
c
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
c
      integer nv,nu
      real v(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call adfmg_tnsr1_2d(v,nv,nu,A,At)
      else
         call adfmg_tnsr1_3d(v,nv,nu,A,At,At)
      endif
      return
      end
c-------------------------------------------------------T--------------
      subroutine adfmg_tnsr1_2d(v,nv,nu,A,Bt) ! u = A u B
      integer nv,nu
      real v(1),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work(lx1*lx1)
      integer e

      nv2 = nv*nv
      nu2 = nu*nu

      if (nv.le.nu) then
         iv=1
         iu=1
         do e=1,nelt
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
            iv = iv + nv2
            iu = iu + nu2
         enddo
      else
         do e=nelt,1,-1
            iu=1+nu2*(e-1)
            iv=1+nv2*(e-1)
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_tnsr1_3d(v,nv,nu,A,Bt,Ct) ! v = [C (x) B (x) A] u
      integer nv,nu
      real v(1),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer e,e0,ee,es

      e0=1
      es=1
      ee=nelt

      if (nv.gt.nu) then
         e0=nelt
         es=-1
         ee=1
      endif

      nu3 = nu**3
      nv3 = nv**3

      do e=e0,ee,es
         iu = 1 + (e-1)*nu3
         iv = 1 + (e-1)*nv3
         call mxm(A,nv,v(iu),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
      enddo

      return
      end
c----------------------------------------------------------------------
c     u = wt .* u
      subroutine adfmg_do_wt(u,wt,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz
      real u(nx,ny,nz,nelt)
      real wt(nx,nz,2,ldim,nelt)
      
      integer e

c     if (nx.eq.2) then
c        do e=1,nelt
c           call outmat(wt(1,1,1,1,e),nx,nz,'wt 1-1',e)
c           call outmat(wt(1,1,2,1,e),nx,nz,'wt 2-1',e)
c           call outmat(wt(1,1,1,2,e),nx,nz,'wt 1-2',e)
c           call outmat(wt(1,1,2,2,e),nx,nz,'wt 2-2',e)
c        enddo
c        call exitti('hsmg_do_wt quit$',nelt)
c     endif

      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               u( 1,j,1,ie)=u( 1,j,1,ie)*wt(j,1,1,1,ie)
               u(nx,j,1,ie)=u(nx,j,1,ie)*wt(j,1,2,1,ie)
            enddo
            do i=2,nx-1
               u(i, 1,1,ie)=u(i, 1,1,ie)*wt(i,1,1,2,ie)
               u(i,ny,1,ie)=u(i,ny,1,ie)*wt(i,1,2,2,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
               u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
            enddo
            enddo
            do k=1,nz
            do i=2,nx-1
               u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
               u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
            enddo
            enddo
            do j=2,ny-1
            do i=2,nx-1
               u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
               u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_mask_e(w,mask) ! Zero out Dirichlet conditions
      include 'SIZE'
      real w(1)
      integer mask(0:1)

      n=mask(0)
      do i=1,n
c        write(6,*) i,mask(i),n,' MG_MASK'
         w(mask(i)) = 0.
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine adfmg_outpost(b,l)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      integer n
      real    b(1)
      
      common  /adfoutp/ tmp
      real              tmp(lx1*ly1*lz1*lelt*lmgn)
      integer i,j,ie

      if(l.eq.mg_h1_lmax) then
        call outpost(vx,vy,vz,pr,b,'   ')
        return
      endif

      is = 1
      do il=l,mg_h1_lmax-1
        if(il.eq.l) then ! first pass
          call adfmg_intp(tmp(is),b,il)
        else
          call adfmg_intp(tmp(is),tmp(im),il)
        endif
        im = is
        n = mg_h1_n(il+1,mg_fld)
        is = is + n
      enddo
      if(nio.eq.0) write(6,*) 'adfmg outpost: orig field level',l
      call outpost(vx,vy,vz,pr,tmp(im),'   ')

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ier)
      integer                    lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ier
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'TSTEP'
      include 'SOLN' ! vx, vy
c
      integer fbc(6)
      logical ifr

      real    w(lx1,ly1,lz1)
      common /scradf/ w

      common /logmn/ ifrstr
      logical        ifrstr
c
c     ibc = 0  <==>  Dirichlet, with extension
c     ibc = 1  <==>  Dirichlet, no extension, for boundaries
c     ibc = 2  <==>  Neumann,  domain boundary
c     ibc = 3  <==>  Neumann,  element boundary

      ! make them all Dirichlet as a start?

      do iface=1,2*ldim
         ied = eface(iface)
         ibc = -1

         if (ifmhd) call mhd_bc_dn(ibc,iface,e) ! can be overwritten by 'mvn'

         if (cbc(ied,e,ifield).eq.'   ' .or.
     $       cbc(ied,e,ifield).eq.'E  ' .or.
     $       cbc(ied,e,ifield).eq.'P  ' .or.
     $       cbc(ied,e,ifield).eq.'p  ')then
             ibc = 0
             if(ifrstr) then ! ifr = .true. : turn on Neumann
             fl = 0.
             call surface_flux(fl,vx,vy,vz,e,ied,w) ! face # ied or iface?
               if(fl.gt. 1.e-7) then
c                write(6,*) e,ied,fl,iface,'outflow face'
                 ibc = 2 ! fake neumann
               endif
             endif
         endif

         if (cbc(ied,e,ifield).eq.'T  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'t  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'O  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'o  ') ibc = 2

         if (cbc(ied,e,ifield).eq.'f  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'I  ') ibc = 2

         fbc(iface) = ibc

         if (ier.eq.-1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
  1      format(2i3,i8,i3,2x,a3,'  adf_get_fast_bc_error')

      enddo

      if (ier.eq.-1) call exitti('Error A adf_get_fast_bc$',e)

      lbr = fbc(1)
      rbr = fbc(2)
      lbs = fbc(3)
      rbs = fbc(4)
      lbt = fbc(5)
      rbt = fbc(6)

      ier = 0 
      if (ibc.lt.0) ier = lglel(e)

c     write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
c   6 format(i5,2x,4i3,3x,4(1x,a3),'  adf_get_fast_bc')

      return
      end
c----------------------------------------------------------------------
      subroutine printmax_nt(u,n,txt8)  ! 2+3
      include 'SIZE'
      include 'TOTAL'
      real umx, umn, ume
      real u(1)
      integer n
      character*8 txt8

      umx = glmax(u,n)
      umn = glmin(u,n)
      ! track mean
      ume = glsum(u,n)
      ume = ume / real(n)
      if(nid .eq. 0) then
          write(6,19) n,time,umx,umn,ume,txt8
   19     format(i9,1p4e14.6,a10)
      endif

      return
      end
c-----------------------------------------------------------------------
