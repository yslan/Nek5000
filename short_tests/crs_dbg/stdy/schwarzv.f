c-----------------------------------------------------------------------
c---- Varying level minimal overlap Schwarz for velo
c-----------------------------------------------------------------------
      subroutine flprec(z,w,ifrefc) ! z = M w

!     assume ndim velocity fields, with nelv

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real     z(lt,1), w(lt,1)
      logical  ifrefc

      common /flcoef/ h1(lx1*ly1*lz1*lelt)

      mxmg = 1
      isw = 0
      ifield = 1

      if (ifrefc) then
         call svisc(h1,1) ! set viscosity
         if (nio.eq.0) write (6,*)
     $      'flprec: changing viscosity in',1,h1(1)
      endif

c     write(6,*) nid,'flprec: into mg schwz'
      if(ldim.eq.3.or.param(104).gt.0.5) then ! use tensor solve
        call flmn_mg_schwz_tens(z,w,isw,mxmg,vmult,ifrefc)
      else
        call flmn_mg_schwz(z,w,isw,mxmg,vmult,ifrefc)
      endif

c     write(6,*) nid,'flprec: returned from mg schwz'

      return
      end
c-----------------------------------------------------------------------
c---- tensor path
c-----------------------------------------------------------------------
      subroutine flmg_tens_genall()
      include 'SIZE'
      include 'HSMGL'
      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)
      integer n
      n = nx1*ny1*nz1*nelv
      call rzero(h3,n)
      call tens_fact(mg_h1,mg_h2,h3) ! from nk_prec
      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwzmg_tens_solve(x,r,l,nx,id) !
      real    x(*), r(*) 
      call tens_solv(x,r,l,nx) ! from nk_prec
      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_mg_schwz_tens(x,res,iter,maxit,wt,ifr) ! multi level
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real     x(lt,1), res(lt,1), wt(1)
      integer  iter, maxit
      real     r0(lt,ldim)

      real      tmp(lt,ldim)
      real      tol
      integer   icld
      data      icld / 0 /
      save      icld
      integer   lv
      character name*4
      logical   ifr

      logical  ifcpld  ! if coupled solve

      real w(2*lt), d(2*lt)

      tol = 1.e-12
      ifield = 1
      lv = 2
      lv = max(nint(param(81)),2) ! Tensor had problem with lx=4? - Li
      if(icld.eq.0) then
        call flmn_schwzmg_setup(lv) !!
        icld = icld + 1
        if(nid.eq.0) write(6,*) 'MG prec: tensor uncoupled solve'
      endif

      if(ifr) then ! refactor
        call flmg_tens_genall()
      endif

      ifield = 1 ! vel
      nel = nelfld(ifield) ! nel = nelv
      n = nx1*ny1*nz1*nel

      call opcopy(r0,r0(1,2),r0(1,ldim),res,res(1,2),res(1,ldim))
      call rzero (x(1,1),n)
      call rzero (x(1,2),n)
      if(ldim.eq.3) call rzero (x(1,3),n)

      call opcolv (r0,r0(1,2),r0(1,ldim),bm1) ! mass

      call dssum(r0(1,1),lx1,ly1,lz1)
      call col2 (r0(1,1),v1mask,n)
      call dssum(r0(1,2),lx1,ly1,lz1)
      call col2 (r0(1,2),v2mask,n)
      if(ldim.eq.3) then
        call dssum(r0(1,3),lx1,ly1,lz1)
        call col2 (r0(1,3),v3mask,n)
      endif

      do id=1,ldim
        rn0 = sqrt(glsc2(r0(1,id),r0(1,id),n))
        if(nid.eq.0) write(6,*) 'fl mg tens: init res',rn0,id,nel,nx1
      enddo

      iter = 0

      do id=1,ldim
        if(id.eq.1) name='velx'
        if(id.eq.2) name='vely'
        if(id.eq.3) name='velz'
        call flmn_schwzmg_solve_v_tens(tmp(1,id),r0(1,id),name,icrs) ! v-cycle
        call add2  (x(1,id),tmp(1,id),n)
        iter = iter + icrs
      enddo

      if(ldim.eq.2) then ! apply pressure projection
        call incomprn(x(1,1),x(1,2),tmp,tmp(1,2))
      else
        call incomprn(x(1,1),x(1,2),x(1,3),tmp(1,2))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_schwzmg_solve_v_tens(z,rhs,name,icrs) ! V-cycle
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'HSMGL'
      parameter(lt=lx1*ly1*lz1*lelt)
      real    z(1),rhs(1)                  ! z = M \ rhs
      common /zcrhi/ h2inv(lt)                  ! Work arrays
      common /zcrvh/ h2(lt)
      common /zcrcq/  e(4*lt), r(4*lt),w(lt)
      real            err(lt)

      integer p_msk,p_b,icrs
      integer icall
      save    icall
      data    icall /0/
      character name*4

      ifield=1
      nel = nelfld(ifield)
      nel = nelv
      ntot=lx1*ly1*lz1*nel !

      id=0
      if(name.eq.'velx') id=1
      if(name.eq.'vely') id=2
      if(name.eq.'velz') id=3

      sigma = 0.25
      op    =  1.0                           ! Coefficients for residual
      om    = -1.0                           !    evaluation
      o0    =  0.0                           !    evaluation
      nsmth = 2 !

      do iv=1,1

       l  = mg_h1_lmax
       n  = mg_h1_n(l, mg_fld)
       is = 1                                   ! Solve index

       if (iv.eq.1) then
        call rzero        (z,n) ! make zero
        call copy         (r,rhs,n)             ! r := rhs

        do ism=1,nsmth
        call flmn_schwz_l_tens(e,r,sigma,l,id)         ! z := S r
        call add2          (z,e,n)
        call flmn_axmg     (r,e,op,om,l,w,id)       ! r := r - Ae
        enddo
       else
        call copy         (r,rhs,n)             ! r := rhs
        call flmn_axmg    (r,z,op,om,l,w,id)       ! r := r - Az
        do ism=1,nsmth
        call flmn_schwz_l_tens (e,r,sigma,l,id)         ! e := S r
        call add2         (z,e,n)               ! z := z+e
        call flmn_axmg    (r,e,op,om,l,w,id)       ! r := r - Ae
        enddo
       endif

       do l = mg_h1_lmax - 1, 2, -1             ! Restrict down to
         is = is + n                            ! coarse level
         n  = mg_h1_n(l, mg_fld)                !       T
         call flmg_rstr(r, l, .true.)           ! r := J  r
         call flmn_schwz_l_tens(e(is),r,sigma,l,id)     ! e := sigma W S r
         call flmn_axmg(r,e(is),op,om,l,w,id)       ! r := r - A e
        write(6,*) 'should not be in here for two level'
       enddo

       ! previous lowest level
       is = is + n
       nlow = mg_h1_n(1,mg_fld)
       l  = 1                                   !        T
       call flmg_rstr(r, l, .true.)            ! r  := J  r

       ip = id + (mg_fld-1)*ldim
       p_msk=p_mg_msk(l, ip)

       call flmg_mask(r, mg_imask(p_msk), nel)  !       -1    Mask and
       call flmn_coarse_solve(e(is), r,name,icrs)    ! e := A  r   solve at
       call flmg_mask(e(is),mg_imask(p_msk),nel)!  1    1  1  coarse level

       do l = 2, mg_h1_lmax-1                   ! Prolongate to finest level
         n  = mg_h1_n(l,mg_fld)
         im = is
         is = is - n
         call flmg_intp (w,e(im),l-1)           ! w  :=  J e
         call add2      (e(is),w,n)             ! e  := e  + w

         write(6,*) 'should not be in here for two level: loc b'
       enddo

       l = mg_h1_lmax
       n = mg_h1_n(l,mg_fld)
       im = is
       call flmg_intp(w,e(im),l-1)              ! w :=  J e
       call add2     (z,w,n)                    ! e  := e  + w
       if(id.eq.1) then
         call col2     (z,v1mask,n)                ! mask z
       elseif(id.eq.2) then
         call col2     (z,v2mask,n)                ! mask z
       elseif(id.eq.3) then
         call col2     (z,v3mask,n)                ! mask z
       endif
       call dsavg    (z)                        ! ensure continuous z
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l_tens(e,r,sigma,l,id)
      include 'SIZE'
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)
      real    sigma
      integer l
      common /flogmn/ ifrstr
      logical         ifrstr

      n  = mg_h1_n(l,mg_fld)
      nh = mg_nh(l)

      call flmn_schwz_l_tens_part1 (e,r,l,id) ! !
      if(ifrstr) then
        ! do nothing
      else
        call flmg_schwarz_wt  (e,l)          ! e  := W e
      endif
      call cmult               (e,sigma,n)    !  l       l

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l_tens_part1(e,r,l,id)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)
      integer pmsk
      common /flogmn/ ifrstr
      logical         ifrstr

      n    = mg_h1_n (l,mg_fld)

      ip = id + (mg_fld-1)*ldim
      pmsk = p_mg_msk(l,ip)

      nx   = mg_nh(l)
      ifield = 1
      call flmg_mask (r,mg_imask(pmsk),nelfld(ifield))  ! Zero Dirichlet nodes

      call flmn_schwzmg_tens_solve(e,r,l,nx,id) ! Do the local solves

      if(ifrstr) then
        call flmn_schwz_rs(e,l) !
      endif
      call flmg_dssum(e,l)                           ! sum border nodes
      call flmg_mask (e,mg_imask(pmsk),nelfld(ifield)) ! apply mask 

      return
      end
c-----------------------------------------------------------------------
c---- LU path
c-----------------------------------------------------------------------
      subroutine flmn_mg_schwz(x,res,iter,maxit,wt,ifr) ! multi level
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real     x(lt,1), res(lt,1), wt(1)
      integer  iter, maxit
      real     r0(lt,ldim)
      common /fldcoef/ h1(lx1*ly1*lz1*lelt)
      real             h1
      common /flco/  co
      real           co(lmg_g*ldim*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?
      common /flcp/  cp
      real           cp(lmg_g*(ldim*ldim)*lelt) ! big enough?
      common /flcpi/ pcp
      integer        pcp(lmgx,ldimt1) ! big enough?
      integer        p_c, p_p

      real     tmp(lt,ldim)
      real     tol
      integer  icld
      data     icld / 0 /
      save     icld
      integer  lv
      character name*4
      logical  ifr

      logical  ifcpld  ! if coupled solve

      real w(2*lt), d(2*lt)
c     real ldbg(4*lt*lt)

      integer p_msk

      ifcpld = .true.  ! coupled solve, param(51) = 1
      if(param(51).ge.1.5) ifcpld = .false.

      tol = 1.e-12
      ifield = 1 ! temp
      lv = 4
      lv = 3 ! best choice
      lv = 2
      lv = max(nint(param(81)),2)
      if(icld.eq.0) then
        call flmn_schwzmg_setup(lv) !!
        icld = icld + 1
        if(ifcpld) then
          if(nid.eq.0) write(6,*) 'MG prec: Coupled solve'
        else
          if(nid.eq.0) write(6,*) 'MG prec: Un-coupled solve'
        endif
      endif

      if(ifr) then
c       write(6,*) nid,lv,'Refactoring'
        call flmg_set_co  (p_c,lv) ! c coefficients
        call flmg_set_cp  (p_p,lv) ! for C prime operator
        if(ifcpld) then
        ! coupled, one big system for ndim velocities
c       write(6,*) nid,lv,'Refactoring: into build Coupled'
          call flmg_schwzmn_buildlb_cpled() ! build element local matrices
        else
        ! not coupled, same system for ndim velocities
c       write(6,*) nid,lv,'Refactoring: into build Uncoupled'
          call flmg_schwzmn_buildlb() ! build element local matrices
        endif
      endif

      ifield = 1 ! vel
      nel = nelfld(ifield) ! nel = nelv
      n = nx1*ny1*nz1*nel

      call opcopy(r0,r0(1,2),r0(1,ldim),res,res(1,2),res(1,ldim))
      call rzero (x(1,1),n)
      call rzero (x(1,2),n)
      if(ldim.eq.3) call rzero (x(1,3),n)

      call opcolv (r0,r0(1,2),r0(1,ldim),bm1) ! mass

      call dssum(r0(1,1),lx1,ly1,lz1)
      call col2 (r0(1,1),v1mask,n)
      call dssum(r0(1,2),lx1,ly1,lz1)
      call col2 (r0(1,2),v2mask,n)
      if(ldim.eq.3) then
        call dssum(r0(1,3),lx1,ly1,lz1)
        call col2 (r0(1,3),v3mask,n)
      endif

c     do id=1,ldim
c       call dssum(r0(1,id),lx1,ly1,lz1)
c       call col2(r0(1,id),v1mask,n)
c     enddo

      do id=1,ldim
c       resn = sqrt(glsc2(res(1,id),res(1,id),n))
        rn0 = sqrt(glsc2(r0(1,id),r0(1,id),n))
        if(nid.eq.0) write(6,*) 'fl schwz ml: init res',rn0,id,nel,nx1
      enddo

      umx = glmax(vx,n)
      umn = glmin(vx,n)
      vmx = glmax(vy,n)
      vmn = glmin(vy,n)
      if(nid.eq.0) write(6,*) 'fl schwz ml: init vel',umx,umn,vmx,vmn

c     write(6,*) 'r0 1',(r0(i,1),i=1,n)
c     write(6,*) 'r0 2',(r0(i,2),i=1,n)

      iter = 0 

      if(ifcpld) then
        ! coupled solve
c       if(nid.eq.0) write(6,*) 'Into: solve v coupled'
        ifield = 1
        call flmn_schwzmg_solve_v_cpled(tmp(1,1),r0(1,1)) ! v-cycle
        call opadd2(x,x(1,2),x(1,ldim),tmp,tmp(1,2),tmp(1,ldim))
c       call outpost(x(1,1),x(1,2),x(1,2),pr,t,'   ')
      else
        ! uncoupled solve
c       if(nid.eq.0) write(6,*) 'Into: solve v un coupled'
        do id=1,ldim
          if(id.eq.1) name='velx'
          if(id.eq.2) name='vely'
          if(id.eq.3) name='velz'
          call flmn_schwzmg_solve_v(tmp(1,id),r0(1,id),name,icrs) ! v-cycle
          iter = iter + icrs
          call add2  (x(1,id),tmp(1,id),n)
        enddo
      endif

      if(ldim.eq.2) then ! apply pressure projection
c       if(nid.eq.0) write(6,*) 'Into: incomprn pres proj'
c       call rzero(tmp(1,2),n)
        call incomprn(x(1,1),x(1,2),tmp,tmp(1,2))
      else
c       call rzero(tmp(1,2),n)
        call incomprn(x(1,1),x(1,2),x(1,3),tmp(1,2))
      endif
c     call outpost(x(1,1),x(1,2),x(1,2),pr,t,'   ')
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_schwzmg_solve_v_cpled(z,rhs) ! V-cycle
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'HSMGL'
      parameter(lt=lx1*ly1*lz1*lelt)
      real    z(lt,1),rhs(lt,1)                    ! z = M \ rhs
      common /flcoef/ h1(lx1*ly1*lz1*lelt)
      real            h1
      common /zcrcq/ e(ldim*4*lt), r(ldim*4*lt),w(ldim*lt), zw(ldim*lt)
      real           e, r, w, zw
      real           err(lt)
      integer p_msk,p_b

      ifield=1
      nel = nelfld(ifield)
      nel = nelv
      ntot=ldim*lx1*ly1*lz1*nel !

      sigma =  0.5
      sigma =  0.25
      op    =  1.0                           ! Coefficients for residual
      om    = -1.0                           !    evaluation
      o0    =  0.0                           !    evaluation
      nsmth = 5 !
      nsmth = 1 !

c     n  = mg_h1_n(l, mg_fld)
c     call rzero        (z(1,1),n) ! make zero, not contiguous
c     call rzero        (z(1,2),n) ! make zero, not contiguous
c     if(ldim.eq.3) call rzero(z(1,3),n) ! make zero, not contiguous

      do iv=1,1

       l  = mg_h1_lmax
       n  = mg_h1_n(l, mg_fld)
       is = 1                                   ! Solve index

       if (iv.eq.1) then
        call rzero(zw,ntot)            ! contiguous
        call fieldgap(r, rhs, .false.) ! remove gap

c      write(6,*) 'r  a',(r(i),i=1,ldim*n)
        do ism=1,nsmth
        call flmn_schwz_l_cpl(e,r,sigma,l)         ! z := S r
        call add2         (zw,e,ntot)
        call flmn_axmg_cpl(r,e,op,om,l,w)       ! r := r - Ae
c      write(6,*) 'zw1b',(zw(i),i=1,ldim*n)
c      write(6,*) 'r  b',(r(i),i=1,ldim*n)
c      write(6,*) 'r2 b',(r(i+n),i=1,n)
c      call exitt
        enddo
       else
c       call copy         (r,rhs,ntot)             ! r := rhs
        call fieldgap(r, rhs, .false.)
        call flmn_axmg_cpl(r,zw,op,om,l,w)       ! r := r - Az
        do ism=1,nsmth
        call flmn_schwz_l_cpl (e,r,sigma,l)         ! e := S r
        call add2         (zw,e,ntot)               ! z := z+e
        call flmn_axmg_cpl(r,e,op,om,l,w)       ! r := r - Ae
        enddo
       endif

       do l = mg_h1_lmax - 1, 2, -1             ! Restrict down to
         is = is + n*ldim                       ! coarse level
         n  = mg_h1_n(l, mg_fld)                !       T
         call flmg_rstr_cpl(r, l, .true.)           ! r := J  r
         call flmn_schwz_l_cpl(e(is),r,sigma,l)     ! e := sigma W S r
         call flmn_axmg_cpl(r,e(is),op,om,l,w)       ! r := r - A e
        write(6,*) 'should not be in here for two level'
       enddo

       ! previous lowest level
       is = is + n*ldim
       nlow = mg_h1_n(1,mg_fld)
       l  = 1                                   !        T
c      write(6,*) 'r br',(r(i),i=1,2*n)
       call flmg_rstr_cpl(r, l, .true.)            ! r  := J  r
c      write(6,*) 'r ar',(r(i),i=1,2*nlow)
c      p_msk=p_mg_msk(l, mg_fld)
c      rmn = glmin(r,2*nlow)
c      rmx = glmax(r,2*nlow)
c      write(6,*) rmn,rmx,'r crs min max a'
c      write(6,*) 'rcmp',(r(i),i=1,2*nlow)
       call flmg_mask_cpl(r,nel,l)!       -1    Mask and
c      rmn = glmin(r,2*nlow)
c      rmx = glmax(r,2*nlow)
c      write(6,*) rmn,rmx,'r crs min max b'
c      write(6,*) 'rcmp',(r(i),i=1,2*nlow)
       call flmn_coarse_solve_cpl(e(is), r)    ! e := A  r   solve at 
c      emn = glmin(e(is),2*nlow)
c      emx = glmax(e(is),2*nlow)
c      write(6,*) emn,emx,'e crs min max c'
c      write(6,*) 'ecmp',(e(is+i-1),i=1,2*nlow)
       call flmg_mask_cpl(e(is),nel,l)! 1    1  1  coarse level
c      emn = glmin(e(is),2*nlow)
c      emx = glmax(e(is),2*nlow)
c      write(6,*) emn,emx,'e crs min max d'
c      write(6,*) 'ecmp',(e(is+i-1),i=1,2*nlow)
c      call exitt

       do l = 2, mg_h1_lmax-1                   ! Prolongate to finest level
         n  = mg_h1_n(l,mg_fld)
         im = is
         is = is - n*ldim
         call flmg_intp_cpl (w,e(im),l-1)           ! w  :=  J e
         call add2      (e(is),w,n*ldim)             ! e  := e  + w

         write(6,*) 'should not be in here for two level: loc b'
       enddo

       l = mg_h1_lmax
       n = mg_h1_n(l,mg_fld)
       im = is
       call flmg_intp_cpl(w,e(im),l-1)              ! w :=  J e
c      write(6,*) ' w  ',(w(i),i=1,2*n)

c      write(6,*) 'zw a',(zw(i),i=1,2*n)
       call add2(zw,w,ntot)
c      write(6,*) 'zw b',(zw(i),i=1,2*n)

       call col2     (zw(1  ),v1mask,n)             ! mask z
       call col2     (zw(1+n),v2mask,n)             ! mask z
       if(ldim.eq.3) call col2 (zw(1+2*n),v3mask,n) ! mask z
       call dsavg    (zw(1  ))                      ! ensure continuous z
       call dsavg    (zw(1+n))                      ! ensure continuous z
       if(ldim.eq.3) call dsavg(zw(1+2*n))          ! ensure continuous z
c      zmn = glmin(zw,2*n)
c      zmx = glmax(zw,2*n)
c      write(6,*) zmn,zmx,'zw soln'
c      write(6,*) 'zw c',(zw(i),i=1,2*n)
      enddo

      call fieldgap (z, zw, .true.) ! add in gap
c     call outpost(z(1,1),z(1,2),z(1,2),pr,t,'   ')
c      write(6,*) 'z  1',(z(i,1),i=1,n)
c      write(6,*) 'z  2',(z(i,2),i=1,n)
c     call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_coarse_solve_cpl(e,r)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      real    e(1),r(1)
      integer l
      integer p_msk

      l   = 1 ! lowest level
      nh  = mg_nh  (l)
      nr  = mg_h1_n(l,mg_fld)
      mxc = 300
      mxc = 20
      mxc = min(mxc,lgmres) ! cannot be greater than that
      nel = nelfld(mg_fld)
      rn0 = sqrt(glsc2(r,r,ldim*nr))
      if(nid.eq.0) write(6,*) 'coarse cpl: initial norm of r',rn0,nr
      call flmn_crs_proj_cpl(e,r,mxc,l)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_crs_proj_cpl(x,res,maxit,l)

c     Solve A t = r, via custome projection scheme

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lnt=ldim*lx1*ly1*lz1*lelt)
      real    res(1), x(1)
      integer iter, maxit, l

      real r(lnt),w(lnt),d(lnt),v(lnt),z(lnt), wt(lnt) ! multiplicity
      common /dcrspj/ d       ! scratch
      common /ccrspj/ r,w     ! WORK ARRAYS
      common /zcrsp1/ v
      common /zcrsp2/ z
      common /zcrsps/ p,q
      real            p(lnt,lgmres), q(lnt,lgmres)
      real    proj_tol
      real    alpha, beta

      ngo = -99        ! TURN OFF VERBOSE OUTPUT
      n = mg_h1_n(l,mg_fld) ! length for multigrid(< ntot)
      nh = mg_nh  (l)
      nt = n*ldim

      proj_tol = 1.e-10
      if (istep.le.20.or.mod(istep,10).eq.0) ngo = nio

      m   = lgmres
      mxs = 1 ! maximum schwarz iter - once?
      iconv = 0 ! if converged
      sigma = 0.5
      sigma = 0.25
      zr = 0.
      op = 1.

      ! initial guess
      call rzero(x,nt)                     ! Solution
      rn0 = sqrt(glsc2(res,res,nt))
      rn = rn0*1.e-2
c     proj_tol = rn ! relative tolerance

      iter = 1
      do while (iter.lt.maxit)    ! Main proj loop
         if(iter.eq.1) then                      ! first step :: r = res - A 0
            call copy   (r,res,nt)
         else                                     ! second up  :: r = res - A x
            call flmn_axmg_cpl(w,x,zr,op,l,d)        ! w = A x
            call sub3   (r,res,w,nt)               ! r = r - w
         endif

         ! precoditioner solve
                                                  !       -1
c        call copy  (z,r,nt)                       ! z  = M  v
                                                  !  j       j
         call flmn_schwz_l_cpl(z,r,sigma,l)          ! z  = M  v
                                                  !  j       j
         call flmn_axmg_cpl(w,z,zr,op,l,d)           ! w = A z
                                                  !       j
         do j=1,(iter-1)
           beta = glsc2(q(1,j),w,nt)
           call add2s2(z,p(1,j),-beta,nt)
           call add2s2(w,q(1,j),-beta,nt)
         enddo
         beta = sqrt(glsc2(w,w,nt))

         if(beta.lt.proj_tol) goto 900

         betai = 1./beta
         call copy (p(1,iter),z,nt)
         call cmult(p(1,iter),betai,nt)
         call copy (q(1,iter),w,nt)
         call cmult(q(1,iter),betai,nt)

         alpha = glsc2(q(1,iter),r,nt)
         call add2s2(x,p(1,iter),alpha,nt)

         if (ngo.eq.0) write(6,9)
     $         nt,iter,beta,alpha,proj_tol,dt,'Vcpl'
    9    format(i9,i5,1p4e12.4,' mgprj',1x,a4)
         iter = iter + 1

      enddo
  900 continue
      if (nio.eq.0) write(6,8)
     $   nt,iter,maxit,beta,proj_tol,dt,'Vcpl'
    8    format(i9,i5,i5,1p3e12.4,' mgprjb',1x,a4)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_mask_cpl(w,nel,l)
      include 'SIZE'
      include 'HSMGL'

      real    w(1)
c     integer m(1)        ! Pointer to Dirichlet BCs
      integer nel, l
      integer pmsk

      n = mg_h1_n(l,mg_fld)

c     call flmg_mask(w(1), m, nel)
c     call flmg_mask(w(1+n), m, nel)
c     if(ldim.eq.3) call flmg_mask(w(1+2*n), m, nel)

      do id=1,ldim
        ip = id + (mg_fld-1)*ldim
        pmsk = p_mg_msk(l,ip)
        call flmg_mask (w(1+(id-1)*n),mg_imask(pmsk),nel)  ! Zero Dir nodes
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_intp_cpl(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMGL'

      nc = mg_h1_n(l,mg_fld)
      nf = mg_h1_n(l+1,mg_fld)

      call flmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      call flmg_tnsr(uf(1+nf),mg_nh(l+1),uc(1+nc),mg_nh(l)
     $                       ,mg_jh(1,l),mg_jht(1,l))
      if(ldim.eq.3) then
        call flmg_tnsr(uf(1+2*nf),mg_nh(l+1),uc(1+2*nc),mg_nh(l)
     $                           ,mg_jh(1,l),mg_jht(1,l))
      endif

      return
      end
c------------------------------------------   T  -----------------------
      subroutine flmg_rstr_cpl(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMGL'
      logical ifdssum

      real r(1)
      integer l
      real           w(lx1*ly1*lz1*lelt)
      common /flscr/ w

      nf = mg_h1_n(l+1,mg_fld) !
      nc = mg_h1_n(l  ,mg_fld) ! nc < nf, must

      call flmg_do_wt(r(1),mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call flmg_tnsr1(r(1),mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      if (ifdssum) call flmg_dssum(r(1),l)
      ! fill up 1 to nc

      call copy(w,r(1+nf),nf) ! entire field

      call flmg_do_wt(w,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call flmg_tnsr1(w,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      if (ifdssum) call flmg_dssum(w,l)

      call copy(r(1+nc),w,nc) ! contiguous
      ! fill up (1+nc) to 2*nc

      if(ldim.eq.3) then

        call copy(w,r(1+2*nf),nf) ! entire field

        call flmg_do_wt(w
     $                 ,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                 ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
        call flmg_tnsr1(w,mg_nh(l),mg_nh(l+1)
     $                 ,mg_jht(1,l),mg_jh(1,l))
        if (ifdssum) call flmg_dssum(w,l)

        call copy(r(1+2*nc),w,nc) ! contiguous

      endif

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l_cpl(e,r,sigma,l) !h1mg_schwarz
      include 'SIZE'
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)
      real    sigma
      integer l
      common /flogmn/ ifrstr
      logical         ifrstr

      n  = mg_h1_n(l,mg_fld)
      nh = mg_nh(l)

      call flmn_schwz_l_part1_cpl (e,r,l) ! !
      if(ifrstr) then
        ! do nothing
      else
        call flmg_schwarz_wt  (e(1  ),l)          ! e  := W e
        call flmg_schwarz_wt  (e(1+n),l)          ! e  := W e
        if(ldim.eq.3) call flmg_schwarz_wt(e(1+2*n),l)! e  := W e
      endif
      call cmult               (e(1  ),sigma,n)    !  l       l
      call cmult               (e(1+n),sigma,n)    !  l       l
      if(ldim.eq.3) call cmult (e(1+2*n),sigma,n)! e  := W e

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l_part1_cpl(e,r,l)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield, nelfld
      include 'HSMGL'
      include 'SOLN'

      parameter(lt=lx1*ly1*lz1*lelt)
      real    e(1),r(1)
      integer pmsk
      common /flogmn/ ifrstr
      logical         ifrstr
      common /flscr/ w, tmp
      real           w(ldim*lt), tmp(ldim*lt)

      n    = mg_h1_n (l,mg_fld)
      nx   = mg_nh(l)
      ifield = 1
      nel = nelfld(ifield)

      do id=1,ldim
        ip = id + (mg_fld-1)*ldim
        pmsk = p_mg_msk(l,ip)
        call flmg_mask (r(1+(id-1)*n),mg_imask(pmsk),nel)  ! Zero Dir nodes
      enddo

      ! re-ordering r : put into a new array, to keep the INPUT r
      ! unchanged
      call reorder_field(tmp,r,w,nx,nel,.false.)

      call flmn_schwzmg_lu_cpled(e,tmp,l,nx) ! Do the local solves

      ! re-ordering e
      call reorder_field(e,e,w,nx,nel,.true.)

      if(ifrstr) then
        call flmn_schwz_rs(e(1    ),l) !
        call flmn_schwz_rs(e(1+n  ),l) !
        if(ldim.eq.3) call flmn_schwz_rs(e(1+2*n),l) !
      endif

      call flmg_dssum(e(1  ),l)                           ! sum border nodes
      call flmg_dssum(e(1+n),l)                           ! sum border nodes
      if(ldim.eq.3) call flmg_dssum(e(1+2*n),l) !

      do id=1,ldim
        ip = id + (mg_fld-1)*ldim
        pmsk = p_mg_msk(l,ip)
        call flmg_mask (e(1+(id-1)*n),mg_imask(pmsk),nel) ! apply mask 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_axmg_cpl(w,p,aw,ap,l,wk)
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
      real    w(1),p(1),wk(1)
      integer p_h1,p_h2,p_g,p_b,p_msk, p_c, p_p
      logical ifh2
      common /flco/  co
      real           co(lmg_g*(ldim)*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?
      common /flcp/  cp
      real           cp(lmg_g*(ldim*ldim)*lelt) ! big enough?
      common /flcpi/ pcp
      integer        pcp(lmgx,ldimt1) ! big enough?
      common /flogmn/ ifrstr
      logical         ifrstr
      real aw, ap

      p_h1  = p_mg_h1  (l,mg_fld)
      p_h2  = p_mg_h2  (l,mg_fld)
      p_g   = p_mg_g   (l,mg_fld)
      p_b   = p_mg_b   (l,mg_fld)
      p_c   = pco      (l,mg_fld) ! right
      p_p   = pcp      (l,mg_fld) ! right
      if (p_c  .eq.0) call flmg_set_co  (p_c,l)
      if (p_p  .eq.0) call flmg_set_cp  (p_p,l)

      ifh2 = .false.

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)
      ng = 3*ldim-3     ! 3, 6
      nc = ldim         ! 2, 3
      ncp = ldim*ldim   ! 4, 9

      call flmg_axml_cpl(wk,p
     $               ,mg_h1(p_h1),mg_h2(p_h2),nx,ny,nz,nelfld(mg_fld)
     $               ,mg_g (p_g) , ng ,mg_b(p_b), ifh2
     $               , co(p_c), nc,cp(p_p),ncp,l)

c     n = nx*ny*nz*nelfld(ifield) !
      n = nx*ny*nz*nelfld(mg_fld)

      call flmg_dssum  (wk(1),l) ! hsmg_dssum
      call flmg_dssum  (wk(1+n),l) ! hsmg_dssum
      if(ldim.eq.3) call flmg_dssum  (wk(1+2*n),l) ! hsmg_dssum

      ! w = aw * w + ap * wk = 1 * w - 1 * wk = w - ( C p + A p)
      call add2sxy    (w(1),aw,wk(1),ap,n)
      call add2sxy    (w(1+n),aw,wk(1+n),ap,n)
      if(ldim.eq.3) call add2sxy(w(1+2*n),aw,wk(1+2*n),ap,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_axml_cpl
     $  (w,p,h1,h2,nx,ny,nz,nel,g,ng,b,ifh2,c,nc,cp,ncp,l)
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

      real w (nx*ny*nz,nel,ldim), p(nx*ny*nz,nel,ldim)
     $   , h1(nx*ny*nz,nel), h2(nx*ny*nz,nel)
     $   , b (nx*ny*nz,nel), g (ng*nx*ny*nz,nel)
      real c (nc*nx*ny*nz,nel)
      real cp(ncp*nx*ny*nz,nel)
      integer nc,ncp
c     integer mask(1)

      logical ifh2
      parameter(lxyz=lx1*ly1*lz1)
      common /ctmp0/ ur(lxyz),us(lxyz),ut(lxyz)
      real           ur, us, ut
      common /flscr/ tmp, tmp2
      real           tmp(lxyz*ldim*lelt),tmp2(lxyz*ldim*lelt)
      integer e,p_msk

      n = nx*ny*nz
      nt = n*ndim

c     write(6,*) 'p a ',(p(i,1,1),i=1,n*ndim*nel)
      call reorder_field(tmp,p,tmp2,nx,nel,.false.) ! p->tmp2->tmp
c     write(6,*) 'p b ',(tmp(i),i=1,n*ndim*nel)
      do e=1,nel ! fills up <nt> numbers for each elem
         ile = 1+(e-1)*(nt)
         call flmg_axe_cpl(tmp2(ile),tmp(ile),h1(1,e),h2(1,e)
     $ ,g(1,e),ng,b(1,e),nx,ny,nz,ur,us,ut,ifh2,e,c(1,e),nc,cp(1,e),ncp)
      enddo

c     write(6,*) 'w a ',(tmp2(i),i=1,n*ndim*nel)
      call reorder_field(w,tmp2,tmp,nx,nel,.true.) ! tmp2->tmp->w
c     write(6,*) 'w b ',(w(i,1,1),i=1,n*ndim*nel)

      do id=1,ldim
        ip = id + (mg_fld-1)*ldim
        p_msk = p_mg_msk(l,ip)
        do e=1,nel
           ! im = mask(e)
           im = mg_imask(p_msk+e-1)
           call flmg_mask_e(w(1,1,id),mg_imask(p_msk+im-1)) ! Zero out Dir.
c          call flmg_mask_e(w(1,1,2),mask(im)) ! Zero out Dir. conditions
        enddo
      enddo

c     write(6,*) 'w e ',(w(i,1,1),i=1,n*ndim*nel)
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_axe_cpl
     $     (w,p,h1,h2,g,ng,b,nx,ny,nz,ur,us,ut,ifh2,e,c,nc,cp,ncp) !

      include 'SIZE'
      include 'INPUT'   ! if3d

      real w (nx*ny*nz*ldim), p (nx*ny*nz*ldim)
     $   , h1(nx*ny*nz), h2(nx*ny*nz)
     $   , b (nx*ny*nz), g (ng,nx*ny*nz)
     $   , ur(nx*ny*nz), us(nx*ny*nz), ut(nx*ny*nz)
      real c (nc,nx*ny*nz)
      real cp(ncp,nx*ny*nz)
      integer e

      real tmp(lx1*ly1*lz1) ! largest size for one e

      nxyz = nx*ny*nz

      if(if3d) then
        do id=1,ldim
          call gradl_rst(ur,us,ut,p(1+(id-1)*nxyz),nx,if3d)
          do i=1,nxyz
            tmp(i) = c(1,i)*ur(i) + c(2,i)*us(i) + c(3,i)*ut(i)
            wr = g(1,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
            ws = g(4,i)*ur(i) + g(2,i)*us(i) + g(6,i)*ut(i)
            wt = g(5,i)*ur(i) + g(6,i)*us(i) + g(3,i)*ut(i)
            ur(i) = wr*h1(i)
            us(i) = ws*h1(i)
            ut(i) = wt*h1(i)
          enddo
          call gradl_rst_t(w(1+(id-1)*nxyz),ur,us,ut,nx,if3d)
          call add2(w(1+(id-1)*nxyz),tmp,nxyz)
        enddo
        do i=1,nxyz        ! do coupled terms
          id1 = i
          id2 = i + (1)*nxyz
          id3 = i + (2)*nxyz
          w(id1) = w(id1) + cp(1,i)*p(id1)
     $           + cp(2,i)*p(id2) + cp(3,i)*p(id3)
          w(id2) = w(id2) + cp(4,i)*p(id1)
     $           + cp(5,i)*p(id2) + cp(6,i)*p(id3)
          w(id3) = w(id3) + cp(7,i)*p(id1)
     $           + cp(8,i)*p(id2) + cp(9,i)*p(id3)
        enddo
      else ! 2D
        do id=1,ldim
          call gradl_rst(ur,us,ut,p(1+(id-1)*nxyz),nx,if3d)
          do i=1,nxyz
            tmp(i) = c(1,i)*ur(i) + c(2,i)*us(i)
            wr = g(1,i)*ur(i) + g(3,i)*us(i)
            ws = g(3,i)*ur(i) + g(2,i)*us(i)
            ur(i) = wr*h1(i)
            us(i) = ws*h1(i)
          enddo
          call gradl_rst_t(w(1+(id-1)*nxyz),ur,us,ut,nx,if3d)
          call add2(w(1+(id-1)*nxyz),tmp,nxyz)
        enddo
        do i=1,nxyz        ! do coupled terms
          id1 = i
          id2 = i + (1)*nxyz
          w(id1) = w(id1) + cp(1,i)*p(id1)
     $                    + cp(2,i)*p(id2)
          w(id2) = w(id2) + cp(3,i)*p(id1)
     $                    + cp(4,i)*p(id2)
          ! write(6,*) i,id1,id2,cp(1,i),cp(2,i),cp(3,i),cp(4,i),'test'
          ! write(6,*) i,id1,id2,w(id1),w(id2),p(id1),p(id2),'test'
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine reorder_field(uo,ui,w,nx,nel,ifdim) !
      include 'SIZE'
      include 'TOTAL'

      logical ifdim ! if last dimension in output is dim

c     Both input and output are contiguous
      ! ifdim .eq. true
c     real    uo(nx**ldim,1:nel,1:ldim), ui(nx**ldim,1:ldim,1:nel)
      ! ifdim .eq. false
c     real    uo(nx**ldim,1:ldim,1:nel), ui(nx**ldim,1:nel,1:ldim)

      real    uo(1), ui(1), w(1) ! use w in case ui uo are the same
      integer nx, nel

      ne = nx**ldim  ! for each elem
      if(ifdim) then !
        do id=1,ldim
        do ie=1,nel
          ii = 1 + (id-1)*ne + (ie-1)*ne*ldim
          io = 1 + (ie-1)*ne + (id-1)*ne*nel
          call copy(w(io),ui(ii),ne)
        enddo
        enddo
      else            !
        do ie=1,nel
        do id=1,ldim
          ii = 1 + (ie-1)*ne + (id-1)*ne*nel
          io = 1 + (id-1)*ne + (ie-1)*ne*ldim
          call copy(w(io),ui(ii),ne)
        enddo
        enddo
      endif
      call copy(uo, w, ne*nel*ldim)

      return
      end
c-----------------------------------------------------------------------
      subroutine fieldgap(uo,ui,ifexpd) !
      include 'SIZE'
      include 'TOTAL'
      parameter(lt=lx1*ly1*lz1*lelt)
      logical ifexpd
      real    ui(1), uo(1)

      na = nx1*ny1*nz1*nelv ! actual size
      ii = 1
      io = 1
      if(ifexpd) then       ! expand, assume ui is contiguous
        ni = na
        no = lt
      else                  ! condense, uo is contiguous
        ni = lt
        no = na
      endif

      do id=1,ldim
        call copy (uo(io),ui(ii),na)
        ii = ii + ni
        io = io + no
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_schwzmg_solve_v(z,rhs,name,icrs) ! V-cycle
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'HSMGL'
      parameter(lt=lx1*ly1*lz1*lelt)
      real    z(1),rhs(1)                  ! z = M \ rhs
      common /zcrhi/ h2inv(lt)                  ! Work arrays
      common /zcrvh/ h2(lt)
      common /flcoef/ h1(lx1*ly1*lz1*lelt)
      real             h1
      common /zcrcq/ e(4*lt), r(4*lt),w(lt)
      real            err(lt)

      integer p_msk,p_b
      integer icall
      save    icall
      data    icall /0/
      character name*4

      ifield=1
      nel = nelfld(ifield)
      nel = nelv
      ntot=lx1*ly1*lz1*nel !

c     if (icall.eq.0) then
c        icall=1
c        call rzero (h2   , ntot)
c        call rzero (h2inv, ntot)
c     endif
      id=0
      if(name.eq.'velx') id=1
      if(name.eq.'vely') id=2
      if(name.eq.'velz') id=3

      sigma =  0.5
      sigma = 0.25
      op    =  1.0                           ! Coefficients for residual
      om    = -1.0                           !    evaluation
      o0    =  0.0                           !    evaluation
      nsmth = 5 !
      nsmth = 1 !

      do iv=1,1

       l  = mg_h1_lmax
       n  = mg_h1_n(l, mg_fld)
       is = 1                                   ! Solve index

       if (iv.eq.1) then
        call rzero        (z,n) ! make zero
        call copy         (r,rhs,n)             ! r := rhs

c     write(6,*) 'r a ',(r(i),i=1,n)

        do ism=1,nsmth
        call flmn_schwz_l (e,r,sigma,l,id)         ! z := S r
        call add2          (z,e,n)
        call flmn_axmg    (r,e,op,om,l,w,id)       ! r := r - Ae

c     write(6,*) 'z b ',(z(i),i=1,n)
c     write(6,*) 'r b ',(r(i),i=1,n)
c     if(name.eq.'vely') call exitt

        enddo
       else
        call copy         (r,rhs,n)             ! r := rhs
        call flmn_axmg    (r,z,op,om,l,w,id)       ! r := r - Az
        do ism=1,nsmth
        call flmn_schwz_l (e,r,sigma,l,id)         ! e := S r
        call add2         (z,e,n)               ! z := z+e
        call flmn_axmg    (r,e,op,om,l,w,id)       ! r := r - Ae
        enddo
       endif

       do l = mg_h1_lmax - 1, 2, -1             ! Restrict down to
         is = is + n                            ! coarse level
         n  = mg_h1_n(l, mg_fld)                !       T
         call flmg_rstr(r, l, .true.)           ! r := J  r
         call flmn_schwz_l(e(is),r,sigma,l,id)     ! e := sigma W S r
         call flmn_axmg(r,e(is),op,om,l,w,id)       ! r := r - A e
        write(6,*) 'should not be in here for two level'
       enddo

       ! previous lowest level
       is = is + n
       nlow = mg_h1_n(1,mg_fld)
       l  = 1                                   !        T
c      write(6,*) 'r br',(r(i),i=1,n)
       call flmg_rstr(r, l, .true.)            ! r  := J  r
c      write(6,*) 'r ar',(r(i),i=1,nlow)

       ip = id + (mg_fld-1)*ldim
       p_msk=p_mg_msk(l, ip)

       call flmg_mask(r, mg_imask(p_msk), nel)  !       -1    Mask and
c      write(6,*) 'rcmp',(r(i),i=1,nlow)
       call flmn_coarse_solve(e(is), r,name,icrs)    ! e := A  r   solve at
c      write(6,*) 'ecmp',(e(is+i-1),i=1,nlow)
       call flmg_mask(e(is),mg_imask(p_msk),nel)!  1    1  1  coarse level
c      write(6,*) 'ecmp',(e(is+i-1),i=1,nlow)
c      if(name.eq.'vely') call exitt

       do l = 2, mg_h1_lmax-1                   ! Prolongate to finest level
         n  = mg_h1_n(l,mg_fld)
         im = is
         is = is - n
         call flmg_intp (w,e(im),l-1)           ! w  :=  J e
         call add2      (e(is),w,n)             ! e  := e  + w

         write(6,*) 'should not be in here for two level: loc b'
       enddo

       l = mg_h1_lmax
       n = mg_h1_n(l,mg_fld)
       im = is
       call flmg_intp(w,e(im),l-1)              ! w :=  J e
       call add2     (z,w,n)                    ! e  := e  + w
       if(id.eq.1) then
         call col2     (z,v1mask,n)                ! mask z
       elseif(id.eq.2) then
         call col2     (z,v2mask,n)                ! mask z
       elseif(id.eq.3) then
         call col2     (z,v3mask,n)                ! mask z
       endif
       call dsavg    (z)                        ! ensure continuous z
      enddo

c     write(6,*) 'z   ',(z(i),i=1,n)
c     if(name.eq.'vely') call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l(e,r,sigma,l,id)
      include 'SIZE'
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)
      real    sigma
      integer l
      common /flogmn/ ifrstr
      logical         ifrstr

      n  = mg_h1_n(l,mg_fld)
      nh = mg_nh(l)

      call flmn_schwz_l_part1 (e,r,l,id) ! !
      if(ifrstr) then
        ! do nothing
      else
        call flmg_schwarz_wt  (e,l)          ! e  := W e
      endif
      call cmult               (e,sigma,n)    !  l       l

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwz_l_part1(e,r,l,id)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMGL'
      include 'SOLN'

      real    e(1),r(1)
      integer pmsk
      common /flogmn/ ifrstr
      logical         ifrstr

      n    = mg_h1_n (l,mg_fld)

      ip = id + (mg_fld-1)*ldim
      pmsk = p_mg_msk(l,ip)

      nx   = mg_nh(l)
      ifield = 1
      call flmg_mask (r,mg_imask(pmsk),nelfld(ifield))  ! Zero Dirichlet nodes

      call flmn_schwzmg_lu(e,r,l,nx,id) ! Do the local solves

      if(ifrstr) then
        call flmn_schwz_rs(e,l) !
      endif
      call flmg_dssum(e,l)                           ! sum border nodes
      call flmg_mask (e,mg_imask(pmsk),nelfld(ifield)) ! apply mask 

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_axmg(w,p,aw,ap,l,wk,id)
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
      common /flco/  co
      real           co(lmg_g*(ldim)*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?

      common /flogmn/ ifrstr
      logical         ifrstr
      real aw, ap

      p_h1  = p_mg_h1  (l,mg_fld)
      p_h2  = p_mg_h2  (l,mg_fld)
      p_g   = p_mg_g   (l,mg_fld)
      p_b   = p_mg_b   (l,mg_fld)
      p_c   = pco      (l,mg_fld) ! right
      if (p_c  .eq.0) call flmg_set_co  (p_c,l)

      ip = id + (mg_fld-1)*ldim
      p_msk = p_mg_msk (l,ip)

      ifh2 = .false.

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)
      ng = 3*ldim-3
      nc = ldim

      call flmg_axml (wk,p
     $               ,mg_h1(p_h1),mg_h2(p_h2),nx,ny,nz,nelfld(mg_fld)
     $               ,mg_g (p_g) , ng ,mg_b(p_b), mg_imask(p_msk),ifh2
     $               , co(p_c), nc)

c     n = nx*ny*nz*nelfld(ifield) !
      n = nx*ny*nz*nelfld(mg_fld)

c     write(6,*) l,nx,ny,nz,'what'
c     write(6,*) 'wk a',(wk(i),i=1,n)
      call flmg_dssum  (wk,l) ! hsmg_dssum
      ! w = aw * w + ap * wk = 1 * w - 1 * wk = w - ( C p + A p)
      call add2sxy    (w,aw,wk,ap,n)
c     write(6,*) 'w  b',(w(i),i=1,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_axml
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

         call flmg_axe(w(1,e),p(1,e),h1(1,e),h2(1,e),g(1,e),ng,b(1,e)
     $            ,nx,ny,nz,ur,us,ut,ifh2,e,c(1,e),nc)
   
         im = mask(e)
         call flmg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_axe
     $     (w,p,h1,h2,g,ng,b,nx,ny,nz,ur,us,ut,ifh2,e,c,nc) !

      include 'SIZE'
      include 'INPUT'   ! if3d
      logical ifh2

      real w (nx*ny*nz), p (nx*ny*nz)
     $   , h1(nx*ny*nz), h2(nx*ny*nz)
     $   , b (nx*ny*nz), g (ng,nx*ny*nz)
     $   , ur(nx*ny*nz), us(nx*ny*nz), ut(nx*ny*nz)
      real c (nc,nx*ny*nz)
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
      call add2(w,tmp,nxyz)

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_coarse_solve(e,r,name,icrs)
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
      integer l
      integer p_msk
      character name*4

      l   = 1 ! lowest level
      nh  = mg_nh  (l)
      nr  = mg_h1_n(l,mg_fld)
      mxc = 300
      mxc = 20
      nel = nelfld(mg_fld)
      rn0 = sqrt(glsc2(r,r,nr))
      if(nid.eq.0) write(6,*) 'coarse: initial norm of r',rn0,nr
c     if(nid.eq.0) write(6,*) r(1),r(2),r(1000)
      call rzero(e,nr)
      call flmn_crs_proj(e,r,icrs,mxc,l,name)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_crs_proj(x,res,iter,maxit,l,name)

c     Solve A t = r, via custome projection scheme

      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real    res(1), x(1)
      integer iter, maxit, l

      real r(lt),w(lt),d(lt),v(lt),z(lt), wt(lt) ! multiplicity
      common /dcrspj/ d       ! scratch
      common /ccrspj/ r,w     ! WORK ARRAYS
      common /zcrsp1/ v
      common /zcrsp2/ z
      common /zcrsps/ p,q
      real            p(lt,lgmres), q(lt,lgmres)
      real    proj_tol
      real    alpha, beta
      character name*4

      logical iflu
      iflu = .true.
      if(param(104) .gt. 0.5) iflu = .false.
      if(ldim.eq.3) iflu = .false.
      if(nid.eq.0) write(6,*) 'flmn crs proj: if using LU:',iflu

      ngo = -99        ! TURN OFF VERBOSE OUTPUT
      n = mg_h1_n(l,mg_fld) ! length for multigrid(< ntot)
      nh = mg_nh  (l)
      proj_tol = 1.e-10
      if (istep.le.20.or.mod(istep,10).eq.0) ngo = nio

      m   = lgmres
      mxs = 1 ! maximum schwarz iter - once?
      iconv = 0 ! if converged
      sigma = 0.25
      zr = 0.
      op = 1.

      id = 0
      if(name.eq.'velx') id=1
      if(name.eq.'vely') id=2
      if(name.eq.'velz') id=3

      ! initial guess
      call rzero(x,n)                     ! Solution
      rn0 = sqrt(glsc2(res,res,n))
      rn = rn0*1.e-2
c     proj_tol = rn ! relative tolerance

      iter = 1
      do while (iter.lt.maxit)    ! Main proj loop
         if(iter.eq.1) then                      ! first step :: r = res - A 0
            call copy   (r,res,n)
         else                                     ! second up  :: r = res - A x
            call flmn_axmg(w,x,zr,op,l,d,id)        ! w = A x
            call sub3   (r,res,w,n)               ! r = r - w
         endif

         ! precoditioner solve
                                                  !       -1
         call copy  (z,r,n)                       ! z  = M  v
                                                  !  j       j
c        if(iflu) then                            !
c          call flmn_schwz_l(z,r,sigma,l,id)        ! z  = M  v
c        else
c          call flmn_schwz_l_tens(z,r,sigma,l,id)   
c        endif
                                                  !  j       j
         call flmn_axmg(w,z,zr,op,l,d,id)           ! w = A z
                                                  !       j
         do j=1,(iter-1)
           beta = glsc2(q(1,j),w,n)
           call add2s2(z,p(1,j),-beta,n)
           call add2s2(w,q(1,j),-beta,n)
         enddo
         beta = sqrt(glsc2(w,w,n))

         if(beta.lt.proj_tol) goto 900

         betai = 1./beta
         call copy (p(1,iter),z,n)
         call cmult(p(1,iter),betai,n)
         call copy (q(1,iter),w,n)
         call cmult(q(1,iter),betai,n)

         alpha = glsc2(q(1,iter),r,n)
         call add2s2(x,p(1,iter),alpha,n)

         if (ngo.eq.0) write(6,9)
     $         n,iter,beta,alpha,proj_tol,dt,name
    9    format(i9,i5,1p4e12.4,' mgprj',1x,a4)
         iter = iter + 1

c        if(mod(iter,5).eq.0) then
c          call adfmg_debug_coarse(x,nh)
c          call adfmg_debug_coarse(r,nh)
c        endif
      enddo
  900 continue
      if (nio.eq.0) write(6,8)
     $   n,iter,maxit,beta,proj_tol,dt,name
    8    format(i9,i5,i5,1p3e12.4,' mgprjb',1x,a4)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_schwzmg_setup(lv)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      include 'SEMHAT'

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical  ifdfrm, iffast, ifh2, ifsolv

      common /flogmn/ ifrstr
      logical         ifrstr

      integer p_h1,p_h2,p_g,p_b,p_u,e, p_c
      integer p_msk(3),p_tmp

      integer icall
      save    icall
      data    icall /0/

      integer lv

      ifield = 1
      write(6,*) 'fluid: famg_setup',icall,ifield,lv
      icall = icall+1
      if (icall.gt.10) stop

      ifrstr = .true.  ! RAS
      ifrstr = .false. ! WAS

      if(ifrstr) then
        if(nio.eq.0) write(6,*) 'flmn_schwzmg: Restricted AS'
      else
        if(nio.eq.0) write(6,*) 'flmn_schwzmg: Weighted   AS'
      endif

      n = lx1*ly1*lz1*nelv
      param(59) = 1
      do e=1,nelv
         ifdfrm(e)=.true.
      enddo
      call geom_reset(1)  ! Recompute g1m1 etc. with deformed only

      ifield = 1
      mg_fld = ifield    ! rely on mg_fld for pointers to fields
      call flmg_index_0 ! initialize index sets
      call flamg_setup_nx(lv)! Sets the level schedules
      call flmg_setup_semhat ! SEM hat matrices for each level
      call flmg_setup_intp   ! Interpolation operators
      call flmg_setup_dssum  ! set direct stiffness summation handles
      call flmg_setup_wtmask ! set restriction weight matrices
                             ! mg_rstr_wt_index, mg_rstr_wt
      call flmn_setup_schwzmn_wt(.false.) ! ras_mask???
      call flmg_setup_solve  ! set up the solver

      l=mg_h1_lmax ! == mg_lmax
      call flmg_set_h1  (p_h1 ,l)
      call flmg_set_h2  (p_h2 ,l)
      call flmg_set_gb  (p_g  ,p_b,l)
      call flmn_set_msk (p_msk,p_tmp,l,ifrstr) ! mg_imask, p_mg_msk

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_index_0 ! initialize index sets
      include 'SIZE'
      include 'HSMGL'

      n = lmgn*(lmgs+1) ! lmgn: levels, lmgs: fields
c     ns = lmgn

      call izero( mg_rstr_wt_index     , n )
      call izero( mg_mask_index        , n )
      call izero( mg_solve_index       , n )
      call izero( mg_fast_s_index      , n )
      call izero( mg_fast_d_index      , n )
      call izero( mg_schwarz_wt_index  , n )
      
      return
      end
c-----------------------------------------------------------------------
      subroutine flamg_setup_nx(lv)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      integer lv

      if(lv.eq.3) then
        mg_h1_lmax        = 3  ! test 3 levels
        mg_lmax           = mg_h1_lmax
      elseif(lv.eq.4) then
        mg_h1_lmax        = 4  ! test 4 levels
        mg_lmax           = mg_h1_lmax
      else ! lv == 2 or anything else
        mg_h1_lmax        = 2  ! test 2 level
        mg_lmax           = mg_h1_lmax
      endif

      mg_nx(mg_h1_lmax) = lx1-1        ! Polynomial degree 
      do k = mg_h1_lmax-1, 1,-1       ! of multigrid levels
         mg_nx(k) = mg_nx(k+1)-2
      enddo

      if(mg_nx(1).eq.1) then
        mg_nx(1) = 2 ! hard fix for N=1 ?
      endif

      call icopy(mg_ny, mg_nx, mg_h1_lmax)
      call icopy(mg_nz, mg_nx, mg_h1_lmax)
      if (ldim.eq.2) call izero(mg_nz, mg_h1_lmax)

      if (nio.eq.0) then
         write(6, *) 'fl_mg_nx:', (mg_nx(i), i = 1, mg_h1_lmax)
         write(6, *) 'fl_mg_ny:', (mg_ny(i), i = 1, mg_h1_lmax)
         write(6, *) 'fl_mg_nz:', (mg_nz(i), i = 1, mg_h1_lmax)
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

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_semhat ! SEM hat matrices for each level
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      include 'SEMHAT'

      do l=1,mg_h1_lmax
         n = mg_nx(l)     ! polynomial order
         call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zglhat,dgl,jgl,n,wh)
         call copy(mg_ah(1,l),ah,(n+1)*(n+1))
         call copy(mg_bh(1,l),bh,n+1)
         call copy(mg_dh(1,l),dh,(n+1)*(n+1))
         call transpose(mg_dht(1,l),n+1,dh,n+1)
         call copy(mg_zh(1,l),zh,n+1)
         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1
      enddo
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_intp
      include 'SIZE'
      include 'HSMGL'
      include 'SEMHAT'
      integer l,nf,nc

      do l=1,mg_lmax-1

         nf=mg_nh(l+1)
         nc=mg_nh(l)

!        Standard multigrid coarse-to-fine interpolation
         call flmg_setup_intpm(
     $           mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
         call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)

!        Fine-to-coarse interpolation for variable-coefficient operators
         call flmg_setup_intpm(
     $           mg_jhfc(1,l),mg_zh(1,l),mg_zh(1,l+1),nc,nf)
         call transpose(mg_jhfct(1,l),nf,mg_jhfc(1,l),nc)
c        call outmat(mg_jhfc(1,l),nc,nf,'MG_JHFC',l)

      enddo
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_intpm(jh,zf,zc,nf,nc)
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
      subroutine flmg_setup_dssum ! rely on mg_fld
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMGL'
      include 'TSTEP' ! for nelfld
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      common /c_is1/ glo_num(lxyz*lelt)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      integer l
      
      ncrnr = 2**ldim
      call get_vert

      nel = nelfld(mg_fld) ! !
      nelgfld = nelg(mg_fld) ! Trust this to be correct? otherwise
                             ! just use nelgv
      do l=1,mg_lmax       ! direct stiffness summation for each level
         nx=mg_nh(l)
         ny=mg_nh(l)
         nz=mg_nhz(l)
         call setupds(mg_gsh_handle(l,mg_fld),nx,ny,nz
     $                ,nel,nelgfld,vertex,glo_num)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_wtmask ! mg_fld
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP' ! nelfld
      integer i,l
      ! i = mg_mask_index(mg_lmax,mg_fld-1)
      i = 0
      do l=1,mg_lmax
         mg_rstr_wt_index(l,mg_fld)=i
         nel = nelfld(mg_fld)
         i = i+mg_nh(l)*mg_nhz(l)*2*ldim*nel
         call flmg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_schwzmn_wt_1(wt,l,ifsqrt)
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

      enx=mg_nh(l) ! no extrusion
      eny=mg_nh(l)
      enz=mg_nh(l)
      if(.not.if3d) enz=1
      ns = enx*eny*enz*nelfld(mg_fld)

      call rone(mg_work,ns)
      call flmg_dssum(mg_work,l)                           ! sum border nodes

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nh(l)
      if (.not.if3d) nz=1
      nxyz = nx*ny*nz
      k    = 1
      do ie=1,nelfld(mg_fld)
c        call outmat(mg_work(k),nx,ny,'NEW WT',ie)
         call flmg_setup_schwzmn_wt_2(wt,ie,nx,mg_work(k),ifsqrt)
         k = k+nxyz
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_setup_schwzmn_wt_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      real wt(1),work(1)
      logical ifsqrt

      if(ldim.eq.2) call flmg_setup_schwzmn_wt2d_2(wt,ie,n,work
     $                                             ,ifsqrt)
      if(ldim.eq.3) call flmg_setup_schwzmn_wt3d_2(wt,ie,n,work
     $                                             ,ifsqrt)

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_schwzmn_wt2d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,4,2,1)
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
      subroutine flmg_setup_schwzmn_wt3d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,n,4,3,1)
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
      subroutine flmg_setup_solve
      include 'SIZE'
      include 'HSMGL'
      
      integer l,i,nl,nlz
      i = mg_solve_index(mg_lmax+1,mg_fld-1)
      do l=1,mg_lmax
         mg_solve_index(l,mg_fld)=i
         i=i+mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
         if(i .gt. lmg_solve*lelv) then
            itmp = i/lelv
            write(6,*) 'lmg_solve too small',i,itmp,lmg_solve,l
            call exitt
         endif
      enddo
      mg_solve_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_set_h1  (p_h1 ,l0)
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'
      integer pf,pc

c     Rely on h1 array in flcoef already filled in

      common /scrvh/ h2    (lx1,ly1,lz1,lelv)
     $             , h2inv (lx1,ly1,lz1,lelv)

      common /flcoef/ h1(lx1*ly1*lz1*lelt)
      real            h1

      integer p_h1

      l                 = mg_h1_lmax
      p_mg_h1(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

      ifield = mg_fld
      write (6,*) 'mg_fld=',mg_fld

      call copy (mg_h1,h1,n)   ! Fine grid is just original h1
      if (nio.eq.0) write (6,*) 'viscosity in flmg_set_h1: ', h1(1)

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h1(l,mg_fld) = p_mg_h1(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h1(l+1,mg_fld)
         pc                = p_mg_h1(l  ,mg_fld)

         call flmg_intp_fc (mg_h1(pc),mg_h1(pf),l)

      enddo

      p_h1 = p_mg_h1(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_set_h2  (p_h2 ,l0)
      include 'SIZE'
      include 'HSMGL'

c     As a first pass, rely on the cheesy common-block interface to get h2

      common /scrvh/ h2    (lx1,ly1,lz1,lelv)
     $             , h2inv (lx1,ly1,lz1,lelv)

      common /adfcoef/ h1(lx1*ly1*lz1*lelt)
      real             h1

      integer p_h2,pf,pc

      l                 = mg_h1_lmax
      p_mg_h2(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

      call rzero(h2   ,n)
      call rzero(h2inv,n)

      call copy (mg_h2,h2,n)   ! Fine grid is just original h2

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h2(l,mg_fld) = p_mg_h2(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h2(l+1,mg_fld)
         pc                = p_mg_h2(l  ,mg_fld)

         call flmg_intp_fc (mg_h2(pc),mg_h2(pf),l)

      enddo

      p_h2 = p_mg_h2(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_set_gb  (p_g,p_b,l0)
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

      do e=1,nelfld(mg_fld)
       do l=mg_h1_lmax,1,-1

         nx = mg_nh(l)
         ny = mg_nh(l)
         nz = mg_nhz(l)
         nxyz = nx*ny*nz

         p_g = p_mg_g (l,mg_fld) + ng*nx*ny*nz*(e-1)
         p_b = p_mg_b (l,mg_fld) +    nx*ny*nz*(e-1)

         if (l.eq.mg_h1_lmax) then
            call fl_gxfer_e (mg_g(p_g) ,ng,e  ) ! Fine grid=original G
            call copy    (mg_b(p_b) ,bm1(1,1,1,e),nxyz) ! Fine grid=original B
            call flmg_scale_mass                          ! Divide out Wghts
     $         (mg_b(p_b),mg_g(p_g),mg_bh(1,l),ng,nx,ny,nz,w,.true.)
         else
c           Generate G and B by interpolating their continous counterparts onto
c           the coarse grid and collocating with coarse-grid quadrature weights

            call flmg_intp_gfc_e
     $            (mg_g(p_g),mg_g(l_g),ng,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call flmg_intp_fc_e
     $            (mg_b(p_b),mg_b(l_b)   ,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call flmg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,l+1),ng,nxl,nyl,nzl,w,.false.)

         endif

         l_b = p_b
         l_g = p_g

         nxl = nx
         nyl = ny
         nzl = nz

       enddo

       call flmg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,1),ng,nxl,nyl,nzl,w,.false.)

      enddo

      p_b  = p_mg_b (l0,mg_fld)
      p_g  = p_mg_g (l0,mg_fld)

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_set_msk(p_msk, p_tmp ,l0,ifrs) ! mg_set_msk
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP' ! nelfld
      integer p_tmp, p_msk(ldim)
      logical ifrs

      integer pi_last

      ! make it dimension-dependent
      pi_last = 0

      do id=1,ldim
      l              = mg_h1_lmax
      ip = id + (mg_fld-1)*ldim ! 1,2,3 for 3D fluid
                                ! 1,2   for 2D fluid
      n              = mg_h1_n(l,mg_fld) ! same n?

      p_mg_msk(l,ip) = pi_last

      do l=mg_h1_lmax,1,-1
         nx = mg_nh  (l)
         ny = mg_nh  (l)
         nz = mg_nhz (l)
         p_tmp = p_mg_msk(l,ip)

         call flmn_setup_mask ! nm is an output
     $  (mg_imask(p_tmp),nm,nx,ny,nz,nelfld(mg_fld),l,mg_work,ifrs,id)

         if (l.gt.1) p_mg_msk(l-1,ip)=p_mg_msk(l,ip)+nm

      enddo

      p_msk(id) = p_mg_msk(l0,ip)
      pi_last = p_mg_msk(1,ip) + nm ! save the end point from prev dir

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_setup_mask(mask,nm,nx,ny,nz,nel,l,w,ifr,id) ! h1 version
      include 'SIZE'
      include 'INPUT'        ! if3d
      include 'TSTEP'        ! ifield
      include 'HSMGL'         ! mg_fld

      integer mask(1)        ! Pointer to Dirichlet BCs
      integer nx,ny,nz,l
      real w(nx,ny,nz,nel)
      logical ifr
      
      integer e,count,ptr
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
      logical ifalgn,ifnorx,ifnory,ifnorz

      zero = 0
      nxyz = nx*ny*nz
      n    = nx*ny*nz*nel

      call rone(w,n)   ! Init everything to 1

      ierr = 0
      ierrmx = 0       ! BC verification
      two    = 2
      ifield = mg_fld
      do e=1,nel       ! Set dirichlet nodes to zero

         call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,two,ierr,ifr)

         if (lbr.eq.1) call facev(w,e,4,zero,nx,ny,nz)
         if (rbr.eq.1) call facev(w,e,2,zero,nx,ny,nz)
         if (lbs.eq.1) call facev(w,e,1,zero,nx,ny,nz)
         if (rbs.eq.1) call facev(w,e,3,zero,nx,ny,nz)
         if (if3d) then
            if (lbt.eq.1) call facev(w,e,5,zero,nx,ny,nz)
            if (rbt.eq.1) call facev(w,e,6,zero,nx,ny,nz)
         endif
         ! SYM
         if((id.eq.1).and.((lbr.eq.4).or.(rbr.eq.4))) then
           if (lbr.eq.4) then
             ifc=4
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
           if (rbr.eq.4) then
             ifc=2
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
         endif
         if((id.eq.2).and.((lbs.eq.4).or.(rbs.eq.4))) then
           if(lbs.eq.4) then
             ifc = 1
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
           if (rbs.eq.4) then
             ifc = 3
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
         endif
         if((id.eq.3).and.((lbt.eq.4).or.(rbt.eq.4))) then
           if (lbt.eq.4) then
             ifc = 5
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
           if (rbt.eq.4) then
             ifc = 6
             call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,e) ! bdry
             call facev(w,e,ifc,zero,nx,ny,nz)
           endif
         endif
         ierrmx = max(ierrmx,ierr)
      enddo

      call flmg_dsprod(w,l)    ! need mg_fld

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
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL h1'
         call exitti('D INVALID BC FOUND in flmg_setup_mask$',ierrmx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_set_cp  (p_p,l0) ! advection prime coefficients
      include 'SIZE'
      include 'HSMGL'
      include 'MASS'   ! bm1
      include 'TSTEP'  ! nelfld
      include 'PARALLEL' ! nid

      common /flcp/  cp
      real           cp(lmg_g*(ldim*ldim)*lelt) ! big enough?
      common /flcpi/ pcp
      integer        pcp(lmgx,ldimt1) ! big enough?

      integer p_p,e, l0
      common /ctmp1/ w(lx1*ly1*lz1*lelt*2)

      ncp = ldim*ldim  ! 4 or 9 elements to adv cr,cs,ct tensor

c     write(6,*) 'adfmg set cp: ',mg_lmax
      l                 = mg_lmax
      pcp    (l,mg_fld) = 1 ! NOT 0 ! ! ! !
      n                 = mg_h1_n(l,mg_fld)

      do l=mg_lmax-1,1,-1
         pcp(l,mg_fld) = pcp    (l+1,mg_fld) + n*ncp
         n             = mg_h1_n(l  ,mg_fld)
      enddo

      do e=1,nelfld(mg_fld)
       do l=mg_lmax,1,-1

         nx = mg_nh(l)
         ny = mg_nh(l)
         nz = mg_nhz(l)
         nxyz = nx*ny*nz
         p_p = pcp(l,mg_fld) + ncp*nx*ny*nz*(e-1)

         if (l.eq.mg_lmax) then
           ! fine grid, evaluate cr=cx rx + cy ry + cz rz
           call cp_e(cp(p_p),ncp,e)

c          write(6,*) 'cp a',(cp(p_p+i-1),i=1,4*nxyz)

           call flmg_scale_mass_cp     ! Divide out Wghts for intp
     $         (cp(p_p),mg_bh(1,l),ncp,nx,ny,nz,w,.true.)

c          write(6,*) 'cp b',(cp(p_p+i-1),i=1,4*nxyz)
         else
c           Generate c by interpolating continous parts onto
c           the coarse grid and collocating with coarse-grid quadrature weights
            call flmg_intp_cfc_e !
     $            (cp(p_p),cp(l_p),ncp,nx,ny,nz,nxl,nyl,nzl,e,l,w)

c          write(6,*) 'cp c',(cp(p_p+i-1),i=1,4*nxyz)

            call flmg_scale_mass_cp     ! Multiply by wghts
     $         (cp(l_p),mg_bh(1,l+1),ncp,nxl,nyl,nzl,w,.false.)

c          write(6,*) 'cp d',(cp(l_p+i-1),i=1,4*nxl*nyl*nzl)
         endif
         l_p = p_p ! lag
         nxl = nx
         nyl = ny
         nzl = nz
       enddo
       call flmg_scale_mass_cp     ! Multiply by wghts
     $          (cp(l_p),mg_bh(1,1),ncp,nxl,nyl,nzl,w,.false.)
c          write(6,*) 'cp e',(cp(l_p+i-1),i=1,4*nxl*nyl*nzl)
      enddo
      p_p  = pcp(l0,mg_fld)
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine cp_e(cp,ncp,e)  ! mesh 1
      include 'SIZE'
      include 'TOTAL'

      real    cp(ncp,1)
      integer e
      real    ur,us,ut,vr,vs,vt,wr,ws,wt
      common  /flscr/ ur(lx1*ly1*lz1), us(lx1*ly1*lz1)
     $              , ut(lx1*ly1*lz1)
     $              , vr(lx1*ly1*lz1), vs(lx1*ly1*lz1)
     $              , vt(lx1*ly1*lz1)
     $              , wr(lx1*ly1*lz1), ws(lx1*ly1*lz1)
     $              , wt(lx1*ly1*lz1)

      nxyz = lx1*ly1*lz1
      N = lx1 - 1

      if(ldim.eq.3) then
        call local_grad3(ur,us,ut,vx,N,e,dxm1,dxtm1)
        call local_grad3(vr,vs,vt,vy,N,e,dym1,dytm1)
        call local_grad3(wr,ws,wt,vz,N,e,dzm1,dztm1)
        do i=1,nxyz
           cp(1,i) = w3m1(i,1,1)*(ur(i)*rxm1(i,1,1,e)
     $                          + us(i)*sxm1(i,1,1,e)
     $                          + ut(i)*txm1(i,1,1,e))
           cp(2,i) = w3m1(i,1,1)*(ur(i)*rym1(i,1,1,e)
     $                          + us(i)*sym1(i,1,1,e)
     $                          + ut(i)*tym1(i,1,1,e))
           cp(3,i) = w3m1(i,1,1)*(ur(i)*rzm1(i,1,1,e)
     $                          + us(i)*szm1(i,1,1,e)
     $                          + ut(i)*tzm1(i,1,1,e))
           cp(4,i) = w3m1(i,1,1)*(vr(i)*rxm1(i,1,1,e)
     $                          + vs(i)*sxm1(i,1,1,e)
     $                          + vt(i)*txm1(i,1,1,e))
           cp(5,i) = w3m1(i,1,1)*(vr(i)*rym1(i,1,1,e)
     $                          + vs(i)*sym1(i,1,1,e)
     $                          + vt(i)*tym1(i,1,1,e))
           cp(6,i) = w3m1(i,1,1)*(vr(i)*rzm1(i,1,1,e)
     $                          + vs(i)*szm1(i,1,1,e)
     $                          + vt(i)*tzm1(i,1,1,e))
           cp(7,i) = w3m1(i,1,1)*(wr(i)*rxm1(i,1,1,e)
     $                          + ws(i)*sxm1(i,1,1,e)
     $                          + wt(i)*txm1(i,1,1,e))
           cp(8,i) = w3m1(i,1,1)*(wr(i)*rym1(i,1,1,e)
     $                          + ws(i)*sym1(i,1,1,e)
     $                          + wt(i)*tym1(i,1,1,e))
           cp(9,i) = w3m1(i,1,1)*(wr(i)*rzm1(i,1,1,e)
     $                          + ws(i)*szm1(i,1,1,e)
     $                          + wt(i)*tzm1(i,1,1,e))
        enddo
      else
        call local_grad2(ur,us,vx,N,e,dxm1,dxtm1)
        call local_grad2(vr,vs,vy,N,e,dym1,dytm1)

        do i=1,nxyz
           cp(1,i) = w3m1(i,1,1)*(ur(i)*rxm1(i,1,1,e)
     $                          + us(i)*sxm1(i,1,1,e))
           cp(2,i) = w3m1(i,1,1)*(ur(i)*rym1(i,1,1,e)
     $                          + us(i)*sym1(i,1,1,e))
           cp(3,i) = w3m1(i,1,1)*(vr(i)*rxm1(i,1,1,e)
     $                          + vs(i)*sxm1(i,1,1,e))
           cp(4,i) = w3m1(i,1,1)*(vr(i)*rym1(i,1,1,e)
     $                          + vs(i)*sym1(i,1,1,e))

c       write(6,*) i,vx(i,1,1,1),vy(i,1,1,1),(cp(j,i),j=1,4), 'cp'
        enddo
c       call rzero(cp,4*nxyz) ! debug -> agree with uncoupled
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_intp_cfc_e(cc,cf,nc,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMGL'

      real cf(nc,nxf,nyf,nzf),cc(nc,nxc,nyc,nzc),w(1)
      integer e,l

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
      subroutine flmg_scale_mass_cp (cp,wt,ncp,nx,ny,nz,wk,ifinv)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'HSMGL'

      real    cp(ncp,1),wt(1),wk(1)
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
            cp(1,k) = wk(k)*cp(1,k)
            cp(2,k) = wk(k)*cp(2,k)
            cp(3,k) = wk(k)*cp(3,k)
            cp(4,k) = wk(k)*cp(4,k)
            cp(5,k) = wk(k)*cp(5,k)
            cp(6,k) = wk(k)*cp(6,k)
            cp(7,k) = wk(k)*cp(7,k)
            cp(8,k) = wk(k)*cp(8,k)
            cp(9,k) = wk(k)*cp(9,k)
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
            cp(1,k) = wk(k)*cp(1,k)
            cp(2,k) = wk(k)*cp(2,k)
            cp(3,k) = wk(k)*cp(3,k)
            cp(4,k) = wk(k)*cp(4,k)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_set_co  (p_c,l0) ! advection coefficients
      include 'SIZE'
      include 'HSMGL'
      include 'MASS'   ! bm1
      include 'SOLN'   ! vx,vy,vz
      include 'TSTEP'  ! nelfld
      include 'PARALLEL' ! nid

      integer p_c, e, l0
      common /ctmp1/ w(lx1*ly1*lz1*lelt*2)

      common /flco/  co
      real           co(lmg_g*ldim*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?

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
           call pack_rst(e,nc, co(p_c), vx,vy,vz)

           call flmg_scale_mass_co     ! Divide out Wghts for intp
     $         (co(p_c),mg_bh(1,l),nc,nx,ny,nz,w,.true.)
         else
c           Generate c by interpolating continous parts onto
c           the coarse grid and collocating with coarse-grid quadrature weights
            call flmg_intp_cfc_e
     $            (co(p_c),co(l_c),nc,nx,ny,nz,nxl,nyl,nzl,e,l,w)
            call flmg_scale_mass_co     ! Multiply by wghts
     $         (co(l_c),mg_bh(1,l+1),nc,nxl,nyl,nzl,w,.false.)
         endif
         l_c = p_c ! lag
         nxl = nx
         nyl = ny
         nzl = nz
       enddo
       call flmg_scale_mass_co     ! Multiply by wghts
     $          (co(l_c),mg_bh(1,1),nc,nxl,nyl,nzl,w,.false.)
      enddo
      p_c  = pco(l0,mg_fld)

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwzmn_buildlb_cpled() ! h1mg_setup_fdm
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      
      integer l,i,j,nl

      common /fllngc/ lb ! system to solve
      real            lb(2*(lx1**2),2*(lx1**2),lelt,lmgn) ! large, (Pablo made this small)

      ifield=1
      two = 2
      ierr = 0
      do l=1,mg_lmax ! does not start from lowest one??
         nl = mg_nh(l)
         write(6,*) 'buildlb cpled:',l,nl,mg_nh(l),mg_lmax,'levels'
         call flmg_schwzmn_build_l_cpled(lb(1,1,1,l),nl,l) ! level l
      enddo

c     call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwzmn_build_l_cpled(lb,nl,l) ! minimal overlap
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      include 'SOLN' ! to vx, vy, vz
      include 'GEOM'
      integer l,i,j,nl
      real    lb(2*(nl**2),2*(nl**2),1) ! Pablo made this small
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
      integer ierr
      integer ph1, pg, pb, pc, pp
      common /flco/  co
      real           co(lmg_g*(ldim)*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?

      common /flcp/  cp
      real           cp(lmg_g*(ldim*ldim)*lelt) ! big enough?
      common /flcpi/ pcp
      integer        pcp(lmgx,ldimt1) ! big enough?

      common /flogmn/ ifrstr
      logical         ifrstr

      common /flilngc/ ipiv ! integer array
      integer          ipiv(2*(lx1**2),lelt,lmgn)   ! pivot index

      nx   = mg_nh(l)
      ny   = mg_nh(l)
      nz   = 1          ! 2D
      if(ldim.eq.3) nz = nx
      nxyz = nx*ny*nz

      ph1  = p_mg_h1(l,mg_fld)
      pb   = p_mg_b (l,mg_fld)
      pg   = p_mg_g (l,mg_fld)
      pc   = pco    (l,mg_fld)
      pp   = pcp    (l,mg_fld)

      nc = ldim
      ng = 3*(ldim-1)
      ncp = ldim*ldim
      ns = ldim*nxyz ! actual shape for the coupled system

      two  = 2
      ierr = 0
      i = 1
      do ie=1,nelv
         ifield = mg_fld
         call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
c        write(6,*) lbr,rbr,lbs,rbs,ie,nelv,'what bcs',ifrstr
         call flmn_build_lbe_l_cpled(lb(1,1,ie),lbr,rbr,lbs,rbs,lbt,rbt
     $   ,mg_h1(ph1),mg_g(pg),mg_b(pb),co(pc),cp(pp),ng,nc,ncp,l,nl,ie)

c     if(l.eq.2) then
c       write(6,*) 'ie=',ie
c       write(6,*) 'lb 2',(lb(k,1,ie),k=1,ns*ns)
c     endif

c     if(l.eq.1) then
c       write(6,*) 'ie=',ie
c       write(6,*) 'lb 1',(lb(k,1,ie),k=1,ns*ns)
c     endif

         call dgetrf(ns,ns,lb(1,1,ie),ns,ipiv(1,ie,l),info) ! LU
         if(info.ne.0) write(6,*) 'issue in dgetrf',ie,info
         ph1  = ph1 + nxyz    ! element count
         pb   = pb  + nxyz    ! element count
         pc   = pc  + nc*nxyz
         pg   = pg  + ng*nxyz ! element count
         pp   = pp  + ncp*nxyz
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_l_cpled(lbe,lbx,rbx,lby,rby,lbz,rbz
     $                            ,h1,gf,bm,co,cp,ng,nc,ncp,l,nl,ie)
      include 'SIZE'

      integer lbx,rbx,lby,rby,lbz,rbz
      real    lbe(1) ! build array
      real    h1(1), gf(1), bm(1), co(1), cp(1)
      real    llx,lmx,lrx,lly,lmy,lry,llz,lmz,lrz

      if(ldim.eq.3) then
        call flmn_build_lbe_3d_l_cpled(lbe,lbx,rbx,lby,rby,lbz,rbz
     $                            ,h1,gf,bm,co,cp,ng,nc,ncp,l,nl,ie)
      else
        call flmn_build_lbe_2d_l_cpled(lbe,lbx,rbx,lby,rby
     $                            ,h1,gf,bm,co,cp,ng,nc,ncp,l,nl,ie)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_3d_l_cpled(Lbe,lbx,rbx,lby,rby,lbz,rbz
     $                               ,h1,gf,bm,co,cp,ng,nc,ncp,l,nl,ie)
      include 'SIZE'

      integer lbx,rbx,lby,rby,lbz,rbz
      real    Lbe(3*nl*nl*nl,3*nl*nl*nl) ! build array, 3D
      real    h1(1), gf(ng,1), bm(1), co(nc,1), cp(ncp,1)
      integer l,nl,ie

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_2d_l_cpled(lbe,lbx,rbx,lby,rby
     $                              ,h1,gf,bm,co,cp,ng,nc,ncp,l,nl,ie)
      include 'SIZE'
      include 'GEOM' ! rxm1, jacm1
      include 'WZ'   ! w3m1
      include 'DXYZ' ! dxm1
      include 'HSMGL' ! mg_dh

      integer l,nl,nc,ng,ie,i,j,p,q
      integer lbx,rbx,lby,rby
      real    lbe(nl,nl,nl,nl,2,2) ! build array

      real    h1(nl,nl), bm(nl,nl)
      real    gf(ng,nl,nl), co(nc,nl,nl), cp(ncp,nl,nl)

      common /flscr/ dr, drt, ds, dst, tmpr, tmps, wrk, wrk2
      real    dr  (lx1*lx1*lx1*lx1) , drt (lx1*lx1*lx1*lx1)
     $      , ds  (lx1*lx1*lx1*lx1) , dst (lx1*lx1*lx1*lx1)
     $      , tmpr(lx1*lx1*lx1*lx1) , tmps(lx1*lx1*lx1*lx1)    
     $      , wrk (lx1*lx1*lx1*lx1)
     $      , wrk2(lx1,lx1,lx1,lx1), wrk4(lx1*lx1,lx1*lx1,4)

      ! coefficients

      n    = nl    ! 1D points
      nz   = 1
      nf   = n*n   ! number of points

c     call gen_iden(wrk,n)
c     call kron2(dr,wrk,n,n,mg_dh(1,l),n,n) ! Dr
c     call kron2(ds,mg_dh(1,l),n,n,wrk,n,n) ! Ds
      call rzero(lbe,n*n*n*n*4)
      do q=1,n
      do p=1,n
      do j=1,n
      do i=1,n
        hv = h1(i,j)
        ij = i + (j-1)*n
        iv = i + (j-1)*n + (p-1)*n*n + (q-1)*n*n*n
        g1 = gf(1,i,j)
        g2 = gf(2,i,j)
        g3 = gf(3,i,j)
        tmpr(iv) = hv* (g1*dr(iv) + g3*ds(iv))
        tmps(iv) = hv* (g3*dr(iv) + g2*ds(iv))
      enddo
      enddo
      enddo
      enddo

      call transpose(drt,nf,dr,nf)
      call transpose(dst,nf,ds,nf)
      call mxm(drt,nf,tmpr,nf,wrk2,nf) ! nf x nf
      call mxm(dst,nf,tmps,nf,wrk,nf) !
      call add2(wrk2,wrk,nf*nf)        ! add to lbe

      do q=1,n
      do p=1,n
      do j=1,n
      do i=1,n
        cr = co(1,i,j)
        cs = co(2,i,j)
        iv = i + (j-1)*n + (p-1)*n*n + (q-1)*n*n*n
        wrk2(iv,1,1,1) = wrk2(iv,1,1,1) + cr * dr(iv) + cs * ds(iv)
        tmp = cr*dr(iv)+cs*ds(iv)
      enddo
      enddo
      enddo
      enddo

      call copy (lbe(1,1,1,1,1,1),wrk2,nf*nf)
      call copy (lbe(1,1,1,1,2,2),wrk2,nf*nf) ! A + C

      call rzero(lbe(1,1,1,1,1,2),nf*nf) ! off diag blocks
      call rzero(lbe(1,1,1,1,2,1),nf*nf) ! off diag blocks
      ! only diagonal
      do j=1,n
      do i=1,n
        c11 = cp(1,i,j)
        c12 = cp(2,i,j)
        c21 = cp(3,i,j)
        c22 = cp(4,i,j)
        lbe(i,j,i,j,1,1) = lbe(i,j,i,j,1,1) + c11
        lbe(i,j,i,j,1,2) = lbe(i,j,i,j,1,2) + c12
        lbe(i,j,i,j,2,1) = lbe(i,j,i,j,2,1) + c21
        lbe(i,j,i,j,2,2) = lbe(i,j,i,j,2,2) + c22
        ! write(6,*) i,j,c11,c12,c21,c22,'i,j,c'
      enddo
      enddo

      ! figure out boundary conditions
      call flbe_bdry_2d_l_cpled(lbe,lbx,rbx,lby,rby
     $                  ,h1,gf,bm,co,cp,dr,ds,drt,dst,ng,nc,ncp,nl,l,ie)

c     write(6,*) 'lbea',(lbe(i,1,1,1,1,1),i=1,nf*nf*4)
      call blocktomem(lbe,lbe,wrk4,nf,ldim)
c     write(6,*) 'lbeb',(lbe(i,1,1,1,1,1),i=1,nf*nf*4)
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine blocktomem(uo,ui,w,n,nd)
      include 'SIZE'
c     set row ia of matrix a to c
      integer n, nd
      integer i,j
      real    uo(1), ui(1), w(1)

      if(nd.eq.3) then
        call blocktomem_3d(uo,ui,w,n)
      else
        call blocktomem_2d(uo,ui,w,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine blocktomem_2d(uo,ui,w,n)
      include 'SIZE'
c     re-order matrix, from block to memory-first
      integer n, nd
      integer i,j
      real    uo(2*n,2*n) ! memory access form
      real    ui(n,n,2,2) ! block form
      real    w (n,n,2,2)
      real    c

      call copy(w,ui,n*n*4) ! in case ui and uo are the same

      do jt=1,2*n
        jb = int((jt-1)/n) ! 0: jt<=n, 1: jt>=n+1
        js = jt - n*jb     ! jt: jt<=n, jt-n: jt>=n+1
        jb = jb + 1        ! 1: jt<=n, 2: jt>=n+1
        call copy(uo(1  ,jt),w(1,js,1,jb),n)
        call copy(uo(1+n,jt),w(1,js,2,jb),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine blocktomem_3d(uo,ui,w,n)
      include 'SIZE'
c     re-order matrix
      integer n, nd
      integer i,j
      real    uo(3*n,3*n)
      real    ui(n,n,3,3)
      real    w (n,n,3,3)

      return
      end
c-----------------------------------------------------------------------
      subroutine flbe_bdry_2d_l_cpled(lbe,lbr,rbr,lbs,rbs
     $                  ,h1,gf,bm,co,cp,dr,ds,drt,dst,ng,nc,ncp,n,l,ie)
      ! assemble
      include 'SIZE'
      include 'GEOM'
      include 'WZ'
      include 'HSMGL'
      integer ng, nc, n, l, ie
      real    lbe(n,n,n,n,2,2)
      real    tmp(n*n,n*n)
      real    co(nc,n,n), h1(n,n), gf(ng,n,n), bm(n*n)
      real    cp(ncp,n,n)
      real    dr(n,n,n,n), ds(n,n,n,n)
     $      , drt(n,n,n,n), dst(n,n,n,n)
      integer lbr, rbr, lbs, rbs
      integer icase(8) ! save this into one array
      real    tmpr(n*n), tmps(n*n)
      integer pm

      logical ifalgn,ifnorx,ifnory,ifnorz

      ! in 3D this will be a nightmare
      call izero(icase,8) ! 8 cases, 4 sides, 4 points
      if(lbr.eq.0) icase(1) = 1
      if(rbr.eq.0) icase(2) = 1
      if(lbs.eq.0) icase(3) = 1
      if(rbs.eq.0) icase(4) = 1
      icase(5) = icase(1)*icase(3) ! 1 when both sides are 1
      icase(6) = icase(1)*icase(4) ! 1 when both sides are 1
      icase(7) = icase(2)*icase(3) ! 1 when both sides are 1
      icase(8) = icase(2)*icase(4) ! 1 when both sides are 1

      nf = n*n
      call rzero(tmp,nf*nf)

      do ni=1,4
        if(icase(ni).eq.1) then ! need assemble
          if(ni.eq.1) then
            istart = 1
            jstart = 1
            iskip  = n
            jskip  = n
            idsrc  = n - 1
            jdsrc  = n - 1
          elseif(ni.eq.2) then
            istart = n
            jstart = n
            iskip  = n
            jskip  = n
            idsrc  = - n + 1
            jdsrc  = - n + 1
          elseif(ni.eq.3) then
            istart = 1
            jstart = 1
            iskip  = 1
            jskip  = 1
            idsrc  =   nf - n
            jdsrc  =   nf - n
          elseif(ni.eq.4) then
            istart = nf-n+1
            jstart = nf-n+1
            iskip  = 1
            jskip  = 1
            idsrc  = - nf + n
            jdsrc  = - nf + n
          endif
          do jj=1,n
          do ii=1,n
            i = istart + (ii-1)*iskip
            j = jstart + (jj-1)*jskip
            isrc = i + idsrc
            jsrc = j + jdsrc
            cr = co(1,i,1)
            cs = co(2,i,1)
            tmp(i,j) = tmp(i,j) + cr*dr(isrc,1,jsrc,1)
     $                          + cs*ds(isrc,1,jsrc,1)
            do k=1,nf
              tmpr(k) = ( gf(1,k,1)*dr(k,1,jsrc,1)
     $                +   gf(3,k,1)*ds(k,1,jsrc,1) ) * h1(i,1)
              tmps(k) = ( gf(3,k,1)*dr(k,1,jsrc,1)
     $                +   gf(2,k,1)*ds(k,1,jsrc,1) ) * h1(i,1)
            enddo
            vdr = vlsc2(dr(1,1,isrc,1),tmpr,nf)
            vds = vlsc2(ds(1,1,isrc,1),tmps,nf)
            tmp(i,j) = tmp(i,j) + (vdr + vds)
            if(i.eq.j) then
              lbe(i,1,i,1,1,1) = lbe(i,1,i,1,1,1) + cp(1,i,1)
              lbe(i,1,i,1,1,2) = lbe(i,1,i,1,1,2) + cp(2,i,1)
              lbe(i,1,i,1,2,1) = lbe(i,1,i,1,2,1) + cp(3,i,1)
              lbe(i,1,i,1,2,2) = lbe(i,1,i,1,2,2) + cp(4,i,1)
            endif
          enddo
          enddo
        endif
      enddo

      do ni=5,8
      if(icase(ni).eq.1) then !
        if(ni.eq.5) then
          i = 1
        elseif(ni.eq.6) then
          i = nf - n + 1
        elseif(ni.eq.7) then
          i = n
        elseif(ni.eq.8) then
          i = nf
        endif
        isrc = nf - i + 1
        cr = co(1,i,1)
        cs = co(2,i,1)
        tmp(i,i) = tmp(i,i) + cr*dr(isrc,1,isrc,1)
     $                      + cs*ds(isrc,1,isrc,1)
        do k=1,nf
          tmpr(k) = ( gf(1,k,1)*dr(k,1,isrc,1)
     $            +   gf(3,k,1)*ds(k,1,isrc,1) ) * h1(i,1)
          tmps(k) = ( gf(3,k,1)*dr(k,1,isrc,1)
     $            +   gf(2,k,1)*ds(k,1,isrc,1) ) * h1(i,1)
        enddo
        vdr = vlsc2(dr(1,1,isrc,1),tmpr,nf)
        vds = vlsc2(ds(1,1,isrc,1),tmps,nf)
        tmp(i,i) = tmp(i,i) + (vdr + vds)
        lbe(i,1,i,1,1,1) = lbe(i,1,i,1,1,1) + cp(1,i,1) ! *iden(isrc,jsrc)
        lbe(i,1,i,1,1,2) = lbe(i,1,i,1,1,2) + cp(2,i,1) ! *iden(isrc,jsrc)
        lbe(i,1,i,1,2,1) = lbe(i,1,i,1,2,1) + cp(3,i,1) ! *iden(isrc,jsrc)
        lbe(i,1,i,1,2,2) = lbe(i,1,i,1,2,2) + cp(4,i,1) ! *iden(isrc,jsrc)
      endif
      enddo

      call add2(lbe(1,1,1,1,1,1),tmp,nf*nf) ! big
      call add2(lbe(1,1,1,1,2,2),tmp,nf*nf) ! big

      ! make use of the existing mask
      zr = 0.
      id = 1
      ip = id + (mg_fld-1)*ldim
      pm = p_mg_msk(l,ip) ! starting loc for each elem
      im = mg_imask(pm+ie-1)  ! move by elem number
      nm = mg_imask(pm+im-1)  ! number of Dirichlet noes
      do i=1,nm
        nr = mg_imask(pm+im+i-1) - (ie-1)*nf
c       call mat_setrow(lbe(1,1,1,1,1,1),nf,nr,zr)
c       call mat_setrow(lbe(1,1,1,1,1,2),nf,nr,zr)
c       call mat_setcol(lbe(1,1,1,1,1,1),nf,nr,zr)
c       call mat_setcol(lbe(1,1,1,1,2,1),nf,nr,zr)
        lbe(nr,1,nr,1,1,1) = 1.
      enddo

      id = 2
      ip = id + (mg_fld-1)*ldim

      pm = p_mg_msk(l,ip) ! starting loc for each elem
      im = mg_imask(pm+ie-1)  ! move by elem number
      nm = mg_imask(pm+im-1)  ! number of Dirichlet noes
      do i=1,nm
        nr = mg_imask(pm+im+i-1) - (ie-1)*nf
c       call mat_setrow(lbe(1,1,1,1,2,1),nf,nr,zr)
c       call mat_setrow(lbe(1,1,1,1,2,2),nf,nr,zr)
c       call mat_setcol(lbe(1,1,1,1,1,2),nf,nr,zr)
c       call mat_setcol(lbe(1,1,1,1,2,2),nf,nr,zr)
        lbe(nr,1,nr,1,2,2) = 1.
      enddo

      ! treat SYM
c     if((lbr.eq.4).or.(rbr.eq.4).or.(lbs.eq.4).or.(rbs.eq.4)) then
c       ! hack: only one face can be SYM
c       if(lbr.eq.4) ifc = 4
c       if(rbr.eq.4) ifc = 2
c       if(lbs.eq.4) ifc = 1
c       if(rbs.eq.4) ifc = 3
c       call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,ie) ! bdry

c       call facind (i0,i1,j0,j1,k0,k1,n,n,1,ifc)
c       if(.not.ifalgn .or. ifnorx) then ! zero out x
c         ! write(6,*) nid,ifalgn,ifnorx,ifnory,n,ie,'x sym face'
c         do k=k0,k1
c         do j=j0,j1
c         do i=i0,i1
c           call oper_ind(ir,i,j,k,n,n,1) ! find ind in the operator
c           call mat_setrow(lbe(1,1,1,1,1,1),nf,ir,zr)
c           call mat_setcol(lbe(1,1,1,1,1,1),nf,ir,zr)
c           call mat_setrow(lbe(1,1,1,1,1,2),nf,ir,zr)
c           call mat_setcol(lbe(1,1,1,1,2,1),nf,ir,zr)
c           lbe(ir,1,ir,1,1,1) =  1.
c         enddo
c         enddo
c         enddo
c       endif
c       if(ifnory) then ! zero out y
c         ! write(6,*) nid,ifalgn,ifnorx,ifnory,n,ie,'y sym face'
c         do k=k0,k1
c         do j=j0,j1
c         do i=i0,i1
c           call oper_ind(ir,i,j,k,n,n,1) ! find ind in the operator
c           call mat_setrow(lbe(1,1,1,1,2,2),nf,ir,zr)
c           call mat_setcol(lbe(1,1,1,1,2,2),nf,ir,zr)
c           call mat_setrow(lbe(1,1,1,1,2,1),nf,ir,zr)
c           call mat_setcol(lbe(1,1,1,1,1,2),nf,ir,zr)
c           lbe(ir,1,ir,1,2,2) =  1.
c         enddo
c         enddo
c         enddo
c       endif
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine oper_ind(ir,i,j,k,nx,ny,nz)
c      ifcase in preprocessor notation
      integer ir, i, j, k, nx, ny, nz

      ir = i + (j-1)*nx + (k-1)*nx*ny ! 2D: k==1

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwzmn_buildlb() ! h1mg_setup_fdm
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      
      integer l,i,j,nl

      common /fllng/ lb ! system to solve
      real           lb(lx1**2,lx1**2,ldim,lelt,lmgn) ! so large

      ifield=1
      two = 2
      ierr = 0
      do l=1,mg_lmax ! does not start from lowest one??
         nl = mg_nh(l)
         if(nid.eq.0) write(6,*) 'buildlb:',l,nl,mg_lmax,'levels'
         call flmg_schwzmn_build_l(lb(1,1,1,1,l),nl,l) ! level l
      enddo
c     call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwzmn_build_l(lb,nl,l) ! minimal overlap
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      include 'SOLN' ! to vx, vy, vz
      include 'GEOM'
      integer l,i,j,nl
      real    lb(nl**ldim,nl**ldim,ldim,1)
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
      integer ierr
      integer ph1, pg, pb, pc
      common /flco/  co
      real           co(lmg_g*(ldim)*lelt) ! big enough?
      common /flcoi/ pco
      integer        pco(lmgx,ldimt1) ! big enough?
      common /flogmn/ ifrstr
      logical         ifrstr

      common /flilng/ ipiv ! integer array
      integer         ipiv(lx1**ldim,ldim,lelt,lmgn)   ! pivot index

      nx   = mg_nh(l)
      ny   = mg_nh(l)
      nz   = 1          ! 2D
      if(ldim.eq.3) nz = nx
      nxyz = nx*ny*nz

      ph1  = p_mg_h1(l,mg_fld)
      pb   = p_mg_b (l,mg_fld)
      pg   = p_mg_g (l,mg_fld)
      pc   = pco    (l,mg_fld)

      nc = ldim
      ng = 3*(ldim-1)
      ns = nxyz

      two  = 2
      ierr = 0
      i = 1
      do ie=1,nelv
         ifield = mg_fld
         call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
c        write(6,*) lbr,rbr,lbs,rbs,ie,nelv,'what bcs',ifrstr
         do id=1,ldim
          call flmn_build_lbe_l(lb(1,1,id,ie),lbr,rbr,lbs,rbs,lbt,rbt
     $        ,mg_h1(ph1),mg_g(pg),mg_b(pb),co(pc),ng,nc,l,nl,id,ie)

          call dgetrf(ns,ns,lb(1,1,id,ie),ns,ipiv(1,id,ie,l),info) ! LU
          if(info.ne.0) write(6,*) 'issue in dgetrf',ie,info
         enddo

         ph1  = ph1 + nxyz    ! element count
         pb   = pb  + nxyz    ! element count
         pc   = pc  + nc*nxyz
         pg   = pg  + ng*nxyz ! element count
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_l(lbe,lbx,rbx,lby,rby,lbz,rbz
     $                            ,h1,gf,bm,co,ng,nc,l,nl,id,ie)
      include 'SIZE'

      integer lbx,rbx,lby,rby,lbz,rbz
      real    lbe(1) ! build array
      real    h1(1), gf(1), bm(1), co(1)
      real    llx,lmx,lrx,lly,lmy,lry,llz,lmz,lrz

      if(ldim.eq.3) then
        call flmn_build_lbe_3d_l(lbe,lbx,rbx,lby,rby,lbz,rbz
     $                            ,h1,gf,bm,co,ng,nc,l,nl,id,ie)
      else
        call flmn_build_lbe_2d_l(lbe,lbx,rbx,lby,rby
     $                            ,h1,gf,bm,co,ng,nc,l,nl,id,ie)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_3d_l(Lbe,lbx,rbx,lby,rby,lbz,rbz
     $                               ,h1,gf,bm,co,ng,nc,l,nl,id,ie)
      include 'SIZE'

      integer lbx,rbx,lby,rby,lbz,rbz
      real    Lbe(nl,nl,nl,nl,nl,nl) ! build array, 3D
      real    h1(1), gf(ng,1), bm(1), co(nc,1)
      integer l,nl,ie

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_build_lbe_2d_l(lbe,lbx,rbx,lby,rby
     $                              ,h1,gf,bm,co,ng,nc,l,nl,id,ie)
      include 'SIZE'
      include 'GEOM' ! rxm1, jacm1
      include 'WZ'   ! w3m1
      include 'DXYZ' ! dxm1
      include 'HSMGL' ! mg_dh

      integer l,nl,nc,ng,ie,i,j,p,q
      integer lbx,rbx,lby,rby
      real    lbe(nl,nl,nl,nl) ! build array

      real    h1(nl,nl), bm(nl,nl)
      real    gf(ng,nl,nl), co(nc,nl,nl)

      common /flscr/ dr, drt, ds, dst, tmpr, tmps, wrk
      real    dr  (lx1*lx1*lx1*lx1) , drt (lx1*lx1*lx1*lx1)
     $      , ds  (lx1*lx1*lx1*lx1) , dst (lx1*lx1*lx1*lx1)
     $      , tmpr(lx1*lx1*lx1*lx1) , tmps(lx1*lx1*lx1*lx1)    
     $      , wrk (lx1*lx1*lx1*lx1)

      ! coefficients
      real    mas

      n    = nl    ! 1D points
      nz   = 1
      np   = n*n ! number of points

c     call gen_iden(wrk,n)
c     call kron2(dr,wrk,n,n,mg_dh(1,l),n,n) ! Dr
c     call kron2(ds,mg_dh(1,l),n,n,wrk,n,n) ! Ds

      call rzero(lbe,n*n*n*n)
      do q=1,n
      do p=1,n
      do j=1,n
      do i=1,n
        hv = h1(i,j)
        ij = i + (j-1)*n
        iv = i + (j-1)*n + (p-1)*n*n + (q-1)*n*n*n
        g1 = gf(1,i,j)
        g2 = gf(2,i,j)
        g3 = gf(3,i,j)
        tmpr(iv) = hv* (g1*dr(iv) + g3*ds(iv))
        tmps(iv) = hv* (g3*dr(iv) + g2*ds(iv))
      enddo
      enddo
      enddo
      enddo

      call transpose(drt,np,dr,np)
      call transpose(dst,np,ds,np)
      call mxm(drt,np,tmpr,np,lbe,np) ! np x np
      call mxm(dst,np,tmps,np,wrk,np) !
      call add2(lbe,wrk,np*np)        ! add to lbe

      ! keep the advection as is for now
      do q=1,n
      do p=1,n
      do j=1,n
      do i=1,n
        cr = co(1,i,j)
        cs = co(2,i,j)
        iv = i + (j-1)*n + (p-1)*n*n + (q-1)*n*n*n
        lbe(i,j,p,q) = lbe(i,j,p,q)
     $   + cr * dr(iv) + cs * ds(iv)
      enddo
      enddo
      enddo
      enddo

      ! figure out boundary conditions
      call flbe_bdry_2d_l(lbe,lbx,rbx,lby,rby
     $                  ,h1,gf,bm,co,dr,ds,drt,dst,ng,nc,nl,l,id,ie)

      return
      end
c-----------------------------------------------------------------------
      subroutine flbe_bdry_2d_l(lbe,lbr,rbr,lbs,rbs
     $                  ,h1,gf,bm,co,dr,ds,drt,dst,ng,nc,n,l,id,ie)
      ! assemble
      include 'SIZE'
      include 'GEOM'
      include 'WZ'
      include 'HSMGL'
      integer ng, nc, n, l, ie
      real    lbe(n,n,n,n)
      real    tmp(n*n,n*n)
      real    co(nc,n,n), h1(n,n), gf(ng,n,n), bm(n*n)
      real    dr(n,n,n,n), ds(n,n,n,n)
     $      , drt(n,n,n,n), dst(n,n,n,n)
      integer lbr, rbr, lbs, rbs
      integer icase(8) ! save this into one array
      real    mas
      real    tmpr(n*n), tmps(n*n)
      integer pm

      logical ifalgn,ifnorx,ifnory,ifnorz

      ! in 3D this will be a nightmare
      call izero(icase,8) ! 8 cases, 4 sides, 4 points
      if(lbr.eq.0) icase(1) = 1
      if(rbr.eq.0) icase(2) = 1
      if(lbs.eq.0) icase(3) = 1
      if(rbs.eq.0) icase(4) = 1
      icase(5) = icase(1)*icase(3) ! 1 when both sides are 1
      icase(6) = icase(1)*icase(4) ! 1 when both sides are 1
      icase(7) = icase(2)*icase(3) ! 1 when both sides are 1
      icase(8) = icase(2)*icase(4) ! 1 when both sides are 1

      nl = n*n
      zr = 0.
      call rzero(tmp,nl*nl)

      do ni=1,4
        if(icase(ni).eq.1) then ! need assemble
          if(ni.eq.1) then
            istart = 1
            jstart = 1
            iskip  = n
            jskip  = n
            idsrc  = n - 1
            jdsrc  = n - 1
          elseif(ni.eq.2) then
            istart = n
            jstart = n
            iskip  = n
            jskip  = n
            idsrc  = - n + 1
            jdsrc  = - n + 1
          elseif(ni.eq.3) then
            istart = 1
            jstart = 1
            iskip  = 1
            jskip  = 1
            idsrc  =   nl - n
            jdsrc  =   nl - n
          elseif(ni.eq.4) then
            istart = nl-n+1
            jstart = nl-n+1
            iskip  = 1
            jskip  = 1
            idsrc  = - nl + n
            jdsrc  = - nl + n
          endif
          do jj=1,n
          do ii=1,n
            i = istart + (ii-1)*iskip
            j = jstart + (jj-1)*jskip
            isrc = i + idsrc
            jsrc = j + jdsrc
            cr = co(1,i,1)
            cs = co(2,i,1)
            tmp(i,j) = tmp(i,j) + cr*dr(isrc,1,jsrc,1)
     $                          + cs*ds(isrc,1,jsrc,1)
            do k=1,nl
              tmpr(k) = ( gf(1,k,1)*dr(k,1,jsrc,1)
     $                +   gf(3,k,1)*ds(k,1,jsrc,1) ) * h1(i,1)
              tmps(k) = ( gf(3,k,1)*dr(k,1,jsrc,1)
     $                +   gf(2,k,1)*ds(k,1,jsrc,1) ) * h1(i,1)
            enddo
            vdr = vlsc2(dr(1,1,isrc,1),tmpr,nl)
            vds = vlsc2(ds(1,1,isrc,1),tmps,nl)
            tmp(i,j) = tmp(i,j) + (vdr + vds)
          enddo
          enddo
        endif
      enddo

      do ni=5,8
      if(icase(ni).eq.1) then !
        if(ni.eq.5) then
          i = 1
        elseif(ni.eq.6) then
          i = nl - n + 1
        elseif(ni.eq.7) then
          i = n
        elseif(ni.eq.8) then
          i = nl
        endif
        isrc = nl - i + 1
        cr = co(1,i,1)
        cs = co(2,i,1)
        tmp(i,i) = tmp(i,i) + cr*dr(isrc,1,isrc,1)
     $                      + cs*ds(isrc,1,isrc,1)
        do k=1,nl
          tmpr(k) = ( gf(1,k,1)*dr(k,1,isrc,1)
     $            +   gf(3,k,1)*ds(k,1,isrc,1) ) * h1(i,1)
          tmps(k) = ( gf(3,k,1)*dr(k,1,isrc,1)
     $            +   gf(2,k,1)*ds(k,1,isrc,1) ) * h1(i,1)
        enddo
        vdr = vlsc2(dr(1,1,isrc,1),tmpr,nl)
        vds = vlsc2(ds(1,1,isrc,1),tmps,nl)
        tmp(i,i) = tmp(i,i) + (vdr + vds)
      endif
      enddo

      call add2(lbe,tmp,nl*nl) ! big

      ! make use of the existing mask
      ip = id + (mg_fld-1)*ldim
      pm = p_mg_msk(l,ip)     ! starting loc for each elem
      im = mg_imask(pm+ie-1)  ! move by elem number
      nm = mg_imask(pm+im-1)  ! number of Dirichlet noes
      do i=1,nm
        nr = mg_imask(pm+im+i-1) - (ie-1)*nl
c       call mat_setrow(lbe,nl,nr,zr)
c       call mat_setcol(lbe,nl,nr,zr)
        lbe(nr,1,nr,1) = 1.
      enddo

      ! sym - taken care of in mask setup
c     if((lbr.eq.4).or.(rbr.eq.4).or.(lbs.eq.4).or.(rbs.eq.4)) then
c       if(lbr.eq.4) ifc = 4
c       if(rbr.eq.4) ifc = 2
c       if(lbs.eq.4) ifc = 1
c       if(rbs.eq.4) ifc = 3

c       call chknord (ifalgn,ifnorx,ifnory,ifnorz,ifc,ie) ! bdry

c       call facind (i0,i1,j0,j1,k0,k1,n,n,1,ifc)
c       if(id.eq.1) then
c         if(.not.ifalgn .or. ifnorx) then ! zero out x
c           write(6,*) nid,ifalgn,ifnorx,ifnory,n,ie,'x sym face, uncpl'
c           do k=k0,k1
c           do j=j0,j1
c           do i=i0,i1
c             call oper_ind(ir,i,j,k,n,n,1) ! find ind in the operator
c             call mat_setrow(lbe,nl,ir,zr)
c             call mat_setcol(lbe,nl,ir,zr)
c             lbe(ir,1,ir,1) =  1.
c           enddo
c           enddo
c           enddo
c         endif
c       elseif(id.eq.2) then
c         if(ifnory) then ! zero out y
c           write(6,*) nid,ifalgn,ifnorx,ifnory,n,ie,'y sym face, uncpl'
c           do k=k0,k1
c           do j=j0,j1
c           do i=i0,i1
c             call oper_ind(ir,i,j,k,n,n,1) ! find ind in the operator
c             call mat_setrow(lbe,nl,ir,zr)
c             call mat_setcol(lbe,nl,ir,zr)
c             lbe(ir,1,ir,1) =  1.
c           enddo
c           enddo
c           enddo
c         endif
c       endif
c     endif

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_setup_rstr_wt(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'HSMGL'
      integer nx,ny,nz,l
      real w(nx,ny,nz,1)
      real wt(nx,nz,2,ldim,1)
      
      integer ie
      !init border nodes to 1
      nel = nelfld(mg_fld)
      call rzero(w,nx*ny*nz*nel)
      if (.not.if3d) then
         do ie=1,nel
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
         do ie=1,nel
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
      call flmg_dssum(w,l)
      !invert the count w to get the weight wt
      if (.not. if3d) then
         do ie=1,nel
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
         do ie=1,nel
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
      subroutine flmg_setup_mask_hs(wt,nx,ny,nz,l,w) ! hs version
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'HSMGL'
      integer nx,ny,nz,l
      real w(nx,ny,nz,1)
      real wt(nx,nz,2,ldim,1)
      
      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
      common /flogmn/ ifrstr
      logical         ifrstr

      nel = nelfld(mg_fld)
      n = nx*ny*nz*nel
      call rone(w,n)

c     set dirichlet nodes to zero
      ierr = 0
      two  = 2
      ifield = mg_fld
      do ie=1,nel
         call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
         if (ierr.ne.0) then
            ierr = -1
            call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
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
         do ie=1,nel
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
         do ie=1,nel
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
c-----------------------------------------------------------------------
      subroutine flmg_intp_fc_e(uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
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
      subroutine flmg_intp_gfc_e(gc,gf,ng,nxc,nyc,nzc,nxf,nyf,nzf
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
      subroutine flmg_scale_mass_co (co,wt,nc,nx,ny,nz,wk,ifinv)
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
      subroutine flmg_scale_mass (b,g,wt,ng,nx,ny,nz,wk,ifinv)
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
      subroutine fl_gxfer_e (g,ng,e)
      include 'SIZE'
      include 'TOTAL'

      real g(ng,1)
      integer e

      nxyz = lx1*ly1*lz1

c     ifdfrm(e) = .true.  ! TOO LATE

      if (if3d) then
         do i=1,nxyz
            g(1,i) = g1m1(i,1,1,e)
            g(2,i) = g2m1(i,1,1,e)
            g(3,i) = g3m1(i,1,1,e)
            g(4,i) = g4m1(i,1,1,e)
            g(5,i) = g5m1(i,1,1,e)
            g(6,i) = g6m1(i,1,1,e)
         enddo
      else
         do i=1,nxyz
            g(1,i) = g1m1(i,1,1,e)
            g(2,i) = g2m1(i,1,1,e)
            g(3,i) = g4m1(i,1,1,e)
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwarz_wt(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      
      if(.not.if3d) call flmg_schwarz_wt2d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      if(if3d) call flmg_schwarz_wt3d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end
c----------------------------------------------------------------------
      subroutine flmg_schwarz_wt2d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,nelv)
      real wt(n,4,2,nelv)
      
      integer ie,i,j
      do ie=1,nelv
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
      subroutine flmg_schwarz_wt3d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelv)
      real wt(n,n,4,3,nelv)
      
      integer ie,i,j,k
      do ie=1,nelv
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
      subroutine flmg_dssum(u,l)
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'
      real u(1)

      if (ifsync) call nekgsync()
      etime1=dnekclock()

c     write(6,*) l,mg_fld,mg_gsh_handle(l,mgfld),'in adfmg dssum'
      call fgslib_gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)
      tdadd =tdadd + dnekclock()-etime1

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_dsprod(u,l)
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'
      real u(1)

      if (ifsync) call nekgsync()

      call fgslib_gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)
      return
      end
c-----------------------------------------------------------------------
      subroutine flmg_mask(w,mask,nel)
      include 'SIZE'

      real    w   (1)
      integer mask(1)        ! Pointer to Dirichlet BCs
      integer e
      
      do e=1,nel
         im = mask(e)
c        write(6,*) e,im,mask(im),'im in adfmg mask'
         call flmg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_mask_e(w,mask) ! Zero out Dirichlet conditions
!
!      this assumes that the field is (nx*ny*nz*nel)
!
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
c-----------------------------------------------------------------------
      subroutine flmg_intp_fc(uc,uf,l) ! l is coarse level

      include 'SIZE'
      include 'HSMGL'

      real uc(1),uf(1)

      nc = mg_nh(l)
      nf = mg_nh(l+1)
      call adfmg_tnsr(uc,nc,uf,nf,mg_jhfc(1,l),mg_jhfct(1,l))

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_intp(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMGL'
      call flmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      return
      end
c------------------------------------------   T  -----------------------
      subroutine flmg_rstr(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMGL'
      logical ifdssum

      real r(1)
      integer l

      call flmg_do_wt(r,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

      call flmg_tnsr1(r,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))

      if (ifdssum) call flmg_dssum(r,l)

      return
      end
c----------------------------------------------------------------------
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
      subroutine flmg_tnsr(v,nv,u,nu,A,At)
      integer nv,nu
      real v(1),u(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call flmg_tnsr2d(v,nv,u,nu,A,At)
      else
         call flmg_tnsr3d(v,nv,u,nu,A,At,At)
      endif
      return
      end
c----------------------------------------------------------------------
c     computes
c              T
c     v = A u B
      subroutine flmg_tnsr2d(v,nv,u,nu,A,Bt)
      integer nv,nu
      real v(nv*nv,nelv),u(nu*nu,nelv),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work((lx1+2)*(lx1+2))
      integer ie
      do ie=1,nelv
         call mxm(A,nv,u(1,ie),nu,work,nu)
         call mxm(work,nv,Bt,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine flmg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv,nelv),u(nu*nu*nu,nelv),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer ie, i
      do ie=1,nelv
         call mxm(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine flmg_tnsr1(v,nv,nu,A,At)
c
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
c
      integer nv,nu
      real v(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call flmg_tnsr1_2d(v,nv,nu,A,At)
      else
         call flmg_tnsr1_3d(v,nv,nu,A,At,At)
      endif
      return
      end
c-------------------------------------------------------T--------------
      subroutine flmg_tnsr1_2d(v,nv,nu,A,Bt) ! u = A u B
      integer nv,nu
      real v(1),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work(lx1*lx1)
      integer e

      nv2 = nv*nv
      nu2 = nu*nu

      nel = nelv
      if (nv.le.nu) then
         iv=1
         iu=1
         do e=1,nel
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
            iv = iv + nv2
            iu = iu + nu2
         enddo
      else
         do e=nel,1,-1
            iu=1+nu2*(e-1)
            iv=1+nv2*(e-1)
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine flmg_tnsr1_3d(v,nv,nu,A,Bt,Ct) ! v = [C (x) B (x) A] u
      integer nv,nu
      real v(1),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer e,e0,ee,es

      e0=1
      es=1
      ee=nelv

      if (nv.gt.nu) then
         e0=nelv
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
      subroutine flmg_do_wt(u,wt,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz
      real u(nx,ny,nz,nelv)
      real wt(nx,nz,2,ldim,nelv)
      
      integer e

c     if (nx.eq.2) then
c        do e=1,nelv
c           call outmat(wt(1,1,1,1,e),nx,nz,'wt 1-1',e)
c           call outmat(wt(1,1,2,1,e),nx,nz,'wt 2-1',e)
c           call outmat(wt(1,1,1,2,e),nx,nz,'wt 1-2',e)
c           call outmat(wt(1,1,2,2,e),nx,nz,'wt 2-2',e)
c        enddo
c        call exitti('hsmg_do_wt quit$',nelv)
c     endif

      if (.not. if3d) then
         do ie=1,nelv
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
         do ie=1,nelv
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
      subroutine flmn_schwz_rs(u,l) ! apply the weights
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      parameter(lt=lx1*ly1*lz1*lelt)

      common /flrsm/ ras_mask, wrk
      real           ras_mask(lx1*ly1*lz1*lelt,lmgn,lmgs)

      real u(1)
      integer l

      nl = mg_nh(l)

      if(ldim.eq.3) then
        call fl_schwarz_rs3d(u,ras_mask(1,l,mg_fld),nl)
      else
        call fl_schwarz_rs2d(u,ras_mask(1,l,mg_fld),nl)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_setup_schwzmn_wt(ifsqrt)
      logical ifsqrt
      include 'SIZE'
      include 'INPUT'
      include 'HSMGL'
      
      integer l,i,nl,nlz
      common /flrsm/ ras_mask, wrk
      real           ras_mask(lx1*ly1*lz1*lelt,lmgn,lmgs)
     $              , wrk(lx1*ly1*lz1*lelt)

      common /flogmn/ ifrstr
      logical         ifrstr

      i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)

c     do l=2,mg_lmax
      do l=1,mg_lmax

         mg_schwarz_wt_index(l,mg_fld)=i
         nl  = mg_nh(l)
         nlz = mg_nhz(l)
         i   = i+nl*nlz*4*ldim*nelv

         call flmg_setup_schwzmn_wt_1(
     $      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)

         if(ifrstr) then
           write(6,*) 'creating ras mask in setup wt',l
           nx = mg_nh(l)
           ny = mg_nh(l)
           nz = mg_nhz(l)
           call flmn_setup_mask_rt(ras_mask(1,l,mg_fld),nx,ny,nz,l,wrk)
         endif

      enddo

      mg_schwarz_wt_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_setup_mask_rt(wv,nx,ny,nz,l,w) ! No 3D
      include 'SIZE'
      include 'INPUT'
      include 'SOLN' ! for outpost debugging only, remove later!
      include 'TSTEP'
      include 'HSMGL'
      integer nx,ny,nz,l
      real  w(nx,ny,nz,1)
      real wv(nx,ny,nz,1)

      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      common /flogmn/ ifrstr
      logical         ifrstr

      nel = nelfld(mg_fld)
      ifield = mg_fld ! bc routine requires this
      n = nx*ny*nz*nel
      call rone(w,n)
      call rone(wv,n)

c     set neumann nodes to zero
      ierr = 0
      two  = 2
      do ie=1,nel
        ! general case
        call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
        if (ierr.ne.0) then
           ierr = -1
           call flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr,ifrstr)
        endif
        sumbc = lbr+rbr+lbs+rbs
        if(ldim.eq.3) then ! 3D NOT TESTED !
          iz0=2
          if(lbt.eq.2) iz0 = 1
          iz1=nz-1
          if(rbt.eq.2) iz1 = nz
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
          ix0=2
          if(lbr.eq.2) ix0=1
          ix1=nx-1
          if(lbr.eq.2) ix1=nx
          iy0=2
          if(lbr.eq.2) iy0=1
          iy1=ny-1
          if(lbr.eq.2) iy1=ny
          if(lbt.eq.2) then
            do j=iy0,iy1
            do i=ix0,ix1
              wv(i,j,1,ie) = 0.
            enddo
            enddo
          endif
          if(rbt.eq.2) then
            do j=iy0,iy1
            do i=ix0,ix1
              wv(i,j,nz,ie) = 0.
            enddo
            enddo
          endif
        else ! 2D
          iz0 = 1
          iz1 = nz
          if(lbr.eq.3) then !
            iy0 = 2
            if(lbs.eq.3) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.3) iy1 = ny
            ! ! dir-neum corner becomes neum
            ! iy0 = 1
            ! iy1 = ny
            do k=iz0,iz1
            do j=iy0,iy1
              wv(1,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbr.eq.3) then
            iy0 = 2
            if(lbs.eq.3) iy0 = 1
            iy1 = ny-1
            if(rbs.eq.3) iy1 = ny
            do k=iz0,iz1
            do j=iy0,iy1
              wv(nx,j,k,ie) = 0.
            enddo
            enddo
          endif
          if(lbs.eq.3) then
            ix0 = 2
            if(lbr.eq.3) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.3) ix1 = nx
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,1,k,ie) = 0.
            enddo
            enddo
          endif
          if(rbs.eq.3) then
            ix0 = 2
            if(lbr.eq.3) ix0 = 1
            ix1 = nx-1
            if(rbr.eq.3) ix1 = nx
            do k=iz0,iz1
            do i=ix0,ix1
              wv(i,ny,k,ie) = 0.
            enddo
            enddo
          endif
        endif
      enddo
      call copy(w, wv, n)
      call flmg_dssum(w,l) ! summed up

      do i=1,nx*ny*nz*nel
        if(w(i,1,1,1).lt.1e-7) then  ! one elem, right=outflow(t field)
          write(6,*) 'Error in rst mask',i,w(i,1,1,1),', set to 1.'
          w (i,1,1,1) = 1.
          wv(i,1,1,1) = 1.
        endif
        w(i,1,1,1) = 1./w(i,1,1,1)
      enddo

      call col2(wv,w,n) ! zero remains zero, otherwize 1/multiplic
c     if(l.eq.1) then
c       call adfmg_debug_coarse(wv,nx)
c     elseif(l.eq.2) then
c       call outpost(vx,vy,vz,pr,wv,'   ')
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flmn_get_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ier,ifr)
      integer                lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ier
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
      common /scrfl/ w
c
c     ibc = 0  <==>  Dirichlet, with extension
c     ibc = 1  <==>  Dirichlet, no extension, for boundaries
c     ibc = 2  <==>  Neumann,  domain boundary
c     ibc = 3  <==>  Neumann,  element boundary
c     ibc = 4  <==>  SYM

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
c            ibc = 3 ! more Neumann faces
             if(ifr) then ! ifr = .true. : turn on Neumann
             fl = 0.
             call surface_flux(fl,vx,vy,vz,e,ied,w) ! face # ied or iface?
c              if(fl.le.-1.e-7) then
c                write(6,*) e,ied,iface,'inflow  face'
c                ibc = 0 ! only turn to Dir when inflow
c              endif
               if(fl.gt. 1.e-7) then
c                write(6,*) e,ied,fl,iface,'outflow face'
                 ibc = 3 ! fake neumann
               endif
             endif
         endif

         if (cbc(ied,e,ifield).eq.'V  ') ibc = 1 !
         if (cbc(ied,e,ifield).eq.'v  ') ibc = 1 !
         if (cbc(ied,e,ifield).eq.'W  ') ibc = 1 !
         if (cbc(ied,e,ifield).eq.'w  ') ibc = 1 !
         if (cbc(ied,e,ifield).eq.'O  ') ibc = 2 !
         if (cbc(ied,e,ifield).eq.'o  ') ibc = 2 !
         if (cbc(ied,e,ifield).eq.'SYM') ibc = 4 !

         fbc(iface) = ibc

         if (ier.eq.-1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
  1      format(2i3,i8,i3,2x,a3,'  flmn_get_bc_error')

      enddo

      if (ier.eq.-1) call exitti('Error A flmn_get_bc$',e)

      lbr = fbc(1)
      rbr = fbc(2)
      lbs = fbc(3)
      rbs = fbc(4)
      lbt = fbc(5)
      rbt = fbc(6)

      ier = 0 
      if (ibc.lt.0) ier = lglel(e)

c     write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
c   6 format(i5,2x,4i3,3x,4(1x,a3),'  flmn_get_bc')

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwzmg_lu_cpled(x,r,l,nx) !
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'

      common /fllngc/ lb ! large system to solve
      real            lb(2*lx1**2,2*lx1**2,lelt,lmgn) ! Pablo made this small
      common /flilngc/ ipiv ! integer array
      integer          ipiv(2*(lx1**2),lelt,lmgn)   ! pivot index

      integer info       ! return information
      real    x(ldim*(nx**ldim),1), r(ldim*(nx**ldim),1) ! 3D declaration
      character*1 tr

      ns = ldim*(nx**ldim) ! coupled system

      ! only factor once
      tr = 'N' ! no transpose
      do ie=1,nelv
        ile = 1 + (ie-1)*ns*ns
        call copy(x(1,ie),r(1,ie),ns)
        call dgetrs(tr,ns,1,lb(ile,1,1,l),ns,ipiv(1,ie,l),x(1,ie)
     $                                                   ,ns,info)
        if(info.ne.0) write(6,*) 'issue in dgetrs',ie,info,l
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flmn_schwzmg_lu(x,r,l,nx,id) !
      include 'SIZE'
      include 'HSMGL'
      include 'CTIMER'

      common /fllng/  lb
      real            lb(lx1**2,lx1**2,ldim,lelt,lmgn)  ! system to solve
      common /flilng/ ipiv ! integer array
      integer         ipiv(lx1**ldim,ldim,lelt,lmgn)   ! pivot index

      integer info       ! return information
      real    x(nx**ldim,1), r(nx**ldim,1) ! 3D declaration
      character*1 tr

      ns = nx**ldim

      ! only factor once
      tr = 'N' ! no transpose
      do ie=1,nelv

        ile = 1 + (ie-1)*ns*ns*ldim + (id-1)*ns*ns

        call copy(x(1,ie),r(1,ie),ns)

        call dgetrs(tr,ns,1,lb(ile,1,1,1,l),ns,ipiv(1,id,ie,l),x(1,ie)
     $                                                   ,ns,info)

        if(info.ne.0) write(6,*) 'issue in dgetrs',ie,info,l
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine fl_schwarz_rs2d(u,mk,n) !
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'
      integer n
      real u(n,n,1)
      real mk(n,n,1)
      integer ie,i,j
      integer ntot

      ntot = n*n*nelfld(mg_fld)
      call col2(u,mk,ntot)

      return
      end
c----------------------------------------------------------------------
      subroutine fl_schwarz_rs3d(u,mk,n) !
      include 'SIZE'
      include 'HSMGL'
      include 'TSTEP'
      integer n
      real u (n,n,n,1)
      real mk(n,n,n,1)
      integer ie,i,j,k
      integer ntot

      ntot = n*n*n*nelfld(mg_fld)
      call col2(u,mk,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine flax(ax,xx) ! only on mesh 1
c
c     apply advection diffusion operator, ax = A xx
c
      include 'SIZE'
      include 'TOTAL'
      parameter(lt=lx1*ly1*lz1*lelt)
      real     ax(lt), xx(lt), h1(lt)
      common   /flmgtmp/ h2, tmp
      real     h2(lt), tmp(lt)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer icald
      data    icald / 0 /
      save    icald

      nt = nx1*ny1*nz1*nelv
      if(icald.eq.0) then
        if(nid.eq.0) write(6,*) 'flax: setup dealiased vel. fields'
        call setup_convect(2)
        icald = 1
      endif

      ifield = 1
      imesh = 2      ! 2 for temperature mesh
      call convop(ax,xx) ! adv, physical space
      call col2  (ax,bm1,nt)

c     call rzero (ax,nt) ! turn on so that only do Laplacian

      isd = 1
      imesh = 2      ! 2 for temperature mesh
      ifield = 1
      call rzero (h2,nt)
      call axhelm(tmp,xx,h1,h2,imesh,isd)
      call add2  (ax,tmp,nt)
      call dssum (ax,lx1,ly1,lz1)
      call col2  (ax,v1mask,nt)

      return
      end
c----------------------------------------------------------------------
      subroutine svisc(h1,jfield)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'

      real h1(lx1,ly1,lz1,lelt)

      n = lx1*ly1*lz1*nelv
      if (jfield.gt.1) n = lx1*ly1*lz1*nelt

      if (param(3).eq.0.) then
c        call vprops
         call copy(h1,vdiff(1,1,1,1,jfield),n) ! ==1, vel
c        call copy(h1,vdiff(1,1,1,1,1),n) ! ==1, vel FIXME
      else
         tmp = param(3)
         if (tmp.lt.0) tmp = 1./abs(tmp)
         call cfill(h1,tmp,n)
      endif

      return
      end
