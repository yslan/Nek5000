c     This is for setting up prec
c     Compute F: compute_fk_nt(SS.f)
c        comp_Fs, done
c        comp_Ff
c        comp_visc_t, done
c        comp_visc_v, done
c        comp_Frhs_q, done
c        comp_Frhs_f, done
c     Compute J: JacobiMatVec_loc(SS.f)
c        comp_J_rhs_q, done
c        comp_J_rhs_f, done
c        comp_J_visc_t, done
c        comp_J_visc_v, done
c        comp_JFf
c        comp_JFs, done
c     Compute prec h: ss_set_prec_h(SS.f)
c        set_prec_h
c
c-----------------------------------------------------------------------
      subroutine comp_Fs(rk,ifld) ! compute the kth residual (scalar)
c     input:  t (global
c     output: rk = f(rk) (there is a outer routine doing binv and tmask
c     F(t) = q(t) - v dot grad t + div[ r(t) grad t ]
c          = bm1*q - C*t - A_r*t - Robin_BC_conv + flux_BC
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lt,ltd,ifld,n
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      real rk(lt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt)

      real c1v,c2v,c3v,u1v,u2v,u3v
      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)
      logical ifpseudo, ifbouss
      common /snsflags/ ifpseudo, ifbouss
      real tau,taui,alpha,taulag
      common /snsvars/ alpha, tau, taui, taulag ! ToDo: Temp

      n = lx1*ly1*lz1*nelt

      call rzero(rk,lt)
      call rzero(t1,lt)

      ifield = ifld
      imesh = 2
      if (ifield.eq.1) imesh = 1

      call setprop
      call sethlm(helm1,helm2,0)
      call bcneusc(helm2,-1)  ! Robin conv

      call copy(t2,t(1,1,1,1,ifield-1),n)
      call convect_new(rk,t2,.false.,c1v,c2v,c3v,.true.) ! Conv. term
      call axhelm(t1,t2,helm1,helm2,2,1)                 ! Diff. term
      call add2(rk,t1,n)

      call chsign(rk,n)
      call bcneusc(t2,1)      ! Neumann/Robin flux
      call add2(rk,t2,n)

      call setqvol(bq)
      call col2(bq,bm1,n)
      call add2(rk,bq,n)

c      if (ifpseudo) then
c         call admcol3(rk,t2,bm1,-taui,n) ! rk = rk - t*bm1/tau
c      endif

      call dssum(rk,lx1,ly1,lz1)
      call col2(rk,tmask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_JFs(jsk,sk) ! compute J(u_k) * s_k (coupled)
c     input:    sk
c     output: J*sk
c     
c     For F(u) = q(u) - v dot grad u + div [r(u) grad(u)]
c         J*p  = J_q*p- v dot grad p + div [r(u) grad(p)] + div [J_f * p grad u]
c              = J_q*p- C*p - A_r*p - A_{(J_f*p)}*u
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lt, ltd
      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)


      real jsk(lt), sk(lt)
      real helm1(lt),helm2(lt)
      real t1(lt), t2(lt), t_ref(lt)
      real ceps,cinveps

      real alpha,tau,taui,taulag
      logical ifpseudo, ifbouss
      real c1v,c2v,c3v,u1v,u2v,u3v
      common /snsvars/ alpha, tau, taui, taulag
      common /snsflags/ ifpseudo, ifbouss
      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      integer n, nel, ifld
      real Jeps

      Jeps = 1.E-6
  
      nel = nelt
      if (ifield.eq.1) nel = nelv
      n = lx1*ly1*lz1*nel
 
      call outpost(vx,vy,vz,pr,sk,'js1')
 
      call rzero(jsk,n)
      call rzero(t1,n)
 
      call copy(t_ref,t(1,1,1,1,ifield-1),n)
 
      ifield = 2 ! ToDo: merge outside
      ifld = ifield
 
      ! Set helm1 helm2
      call setprop
      call sethlm(helm1,helm2,0)
      call bcneusc(helm2,-1)     ! Robin conv
 
      ceps = 1.E-6
      cinveps = 1./ceps
 
      ! convection + (linear) diffusion
      call convect_new(jsk,sk,.false.,c1v,c2v,c3v,.true.) ! do dealais before
      call axhelm(t2,sk,helm1,helm2,2,1)
      call add2(jsk,t2,n)
 
      ! nonlinear diffusion
      call comp_J_visc_t(t1,sk,Jeps,vdiff(1,1,1,1,ifld),ifld,.false.) ! t1 = J_r*p
      call axhelm(t2,t_ref,t1,helm2,2,1)
      call add2(jsk,t2,n)
 
      ! nonlinear diffusive induced convection
      call comp_J_rhs_q(t1,sk,Jeps,t_ref,ifld,.true.)
      call add2(jsk,t1,n)
 
      call dssum(jsk,lx1,ly1,lz1)
      call chsign(jsk,n)
 
      if (ifpseudo) then
         call admcol3(jsk,sk,bm1,-taui,n) ! jsk = sk - sk*bm1/tau
      endif

      call tt_cpws(jsk)
 
      ifield = ifld
      return
      end
c-----------------------------------------------------------------------
c     input:  v_nt
c     output  z_nt
      subroutine tmp_tt_scprec(z,w,ifrefc) ! z = M w
      implicit none
      include 'SIZE'
      include 'TOTAL'
!     assume one scalar field, with nelt

      logical ifrefc
      integer ntot
      real     z(1), w(1)

      ifield = 2
      ntot=nx1*ny1*nz1*nelt

c      call col2(w,bm1,ntot)
      call scprec(z,w,ifrefc)

c      call col2(w,binvm1,ntot)
c      call col2(z,binvm1,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine tmp_tt_set_jac_h
c     Set up Jac coef for prec
c        J = - div r(t) grad t - v dot frad t - Robin_conv
c     nonlinear deendence on t of viscosity 
c        h1, h2, v  ->  h1, h2+div C, v+C
c        C = J_r grad u

      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'TOTAL'
      include 'ANISO'
      include 'STDY' ! ifpstime
      include 'NEWTON_loc' ! dtNT

      integer lt,n,ifld
      parameter(lt=lx1*ly1*lz1*lelt)
      common /adfcoef/ h1,h2,h3
      real h1(lt),h2(lt),h3(lt)

      real tau,taui,alpha,taulag
      common /snsvars/ alpha, tau, taui, taulag ! ToDo: Temp
      logical ifpseudo, ifbouss
      common /snsflags/ ifpseudo, ifbouss
      real acrs(lcr,lcr,lelt)

      real ctemp(lt),h2temp(lt),r1(lt),visc_ref(lt)
      real JJeps

      logical ifconv

      n = nx1*ny1*nz1*nelt
      JJeps = 1.E-6

      ifield = 2  ! FIXME, merge this later, Lan
      ifld = ifield

      if (ifpstime) then
        ifpseudo = .true.
        tau = dtNT
        taui= 1./tau
      endif

      call setprop
      call sethlm(h1,h2,0)
      call bcneusc(h2,-1) ! Robin conv
      call rzero(h3,n)
      call rone(kx,n)
      call rone(ky,n)
      call rone(kz,n)
      call rone(r1,n)

      if (ifpseudo) call cadd(h2,-taui,n)

      ! h2=h2+div(C), vel=vel+C, C = r'(u) grad u
      call gradm1(cx,cy,cz,t)
      call comp_visc(visc_ref,t)

      ! d/dx
      call comp_J_visc_t(ctemp,cx,JJeps,visc_ref,ifld,.false.)
      call add3(cx,vx,ctemp,n)

      call dudxyz (h2temp,ctemp,rxm1,sxm1,txm1,jacm1,1,1)
      call add2(h2,h2temp,n)

      ! d/dy
      call comp_J_visc_t(ctemp,cy,JJeps,visc_ref,ifld,.false.)
      call add3(cy,vy,ctemp,n)

      call dudxyz (h2temp,ctemp,rym1,sym1,tym1,jacm1,1,1)
      call add2(h2,h2temp,n)

      ! d/dz
      if (ndim.eq.3) then
      call comp_J_visc_t(ctemp,cz,JJeps,visc_ref,ifld,.false.)
      call add3(cz,vz,ctemp,n)

      call dudxyz (h2temp,ctemp,rzm1,szm1,tzm1,jacm1,1,1)
      call add2(h2,h2temp,n)
      endif


      ! rhs q
      call comp_J_rhs_q(h2temp,r1,JJeps,t(1,1,1,1,ifld-1)
     $ ,ifield,.false.)
      call add2(h2,h2temp,n)

      call set_convect_new(vxd,vyd,vzd,cx,cy,cz)
      call set_anisotr_new(kxd,kyd,kzd,kx,ky,kz) ! aniso


      ! coarse grid solver
      ifconv = .true.
      call rzero(h2temp,n)
c     call set_up_adf_crs(acrs,h1,h2,h2temp,ifconv,ctemp,visc_ref)
      
      ! Set LU solver
      call crs_lu_fact(h1,h2,h2temp,ifconv,ctemp,visc_ref)
      return
      end
c-----------------------------------------------------------------------
      subroutine comp_visc_t(visc,t_in,ifld)
c     input: t
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      real t_in(lt)
     $   , visc(lt)
     $   , t_bak(lt)
      real vdiff_bak(lt,ldimt1)

      integer n,ie,mfield,ifld,ifld_bak

      n = nx1*ny1*nz1*nelt

      ! baskup
      ifld_bak = ifield
      ifield = ifld
      call copy(t_bak,t(1,1,1,1,ifield-1),n)
      call copy(vdiff_bak,vdiff,lt*ldimt1)

      call copy(t(1,1,1,1,ifield-1),t_in,n)
      do ie=1,nelt
        call nekuvp(ie)
      enddo
      call copy(visc,vdiff(1,1,1,1,ifield),n)

      ! restore
      call copy(vdiff,vdiff_bak,lt*ldimt1)
      call copy(t(1,1,1,1,ifield-1),t_bak,n)
      ifield = ifld_bak

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_visc_v(visc,vx_in,vy_in,vz_in) ! I don't think we need this - Pablo
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real visc(lt)
      real vx_in(lt), vx_bak(lt)
     $   , vy_in(lt), vy_bak(lt)
     $   , vz_in(lt), vz_bak(lt)
     $   , rhsfx(lt), bfx_bak(lt)
     $   , rhsfy(lt), bfy_bak(lt)
     $   , rhsfz(lt), bfz_bak(lt)
      real vdiff_bak(lt,ldimt1)
      integer ie

      ! backup
      call opcopy(vx_bak,vy_bak,vz_bak,vx,vy,vz)
      call copy(vdiff_bak,vdiff,lt*ldimt1)

      call opcopy(vx,vy,vz,vx_in,vy_in,vz_in)

      ifield = 1
      do ie=1,nelt
        call nekuvp(ie)
      enddo
      call copy(visc,vdiff,nx1*ny1*nz1*nelv)

      ! restore
      call copy(vdiff,vdiff_bak,lt*ldimt1)
      call opcopy(vx,vy,vz,vx_bak,vy_bak,vz_bak)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_Frhs_q(rhsq,t_in,ifld,ifwt) ! userq
c     nekuq
c        in:  t = t_in (ifield)
c        out: rhsq
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)
      real t_in(lt)
     $   , rhsq(lt)
     $   , t_bak(lt)

      logical ifwt
      integer ifld,ifld_bak,n,ie

      n = nx1*ny1*nz1*nelt

      ! backup
      ifld_bak = ifield ! ToDo: remove ifield dependency
      call copy(t_bak,t,n)
      if (.not.ifcvfld(ifield)) time = time-dt ! Set time to t^n-1 for user function

      ifield = ifld
      call copy(t,t_in,n) ! input

      call setqvol(rhsq)  ! output
      if (ifwt) call col2(rhsq,bm1,n)

      ! restore
      if (.not.ifcvfld(ifield)) time = time+dt
      call copy(t,t_bak,n)
      ifield = ifld_bak

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_Frhs_f(rhsf,v_in,ifwt) ! userf
c     makeuf
c        in:  vx = vx_in (ifield)
c        out: rhsfx
c     ToDo: pr? stress? ref: makef
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real v_in(lt,ldim), v_bak(lt,ldim)
     $   , rhsf(lt,ldim)
      logical ifwt

      ! backup
      call opcopy(v_bak,v_bak(1,2),v_bak(1,ndim),vx,vy,vz)
      time = time - dt

      call opcopy(vx,vy,vz,v_in,v_in(1,2),v_in(1,ndim))
      call nekuf(rhsf,rhsf(1,2),rhsf(1,ndim))
      if (ifwt) call opcolv(rhsf,rhsf(1,2),rhsf(1,ndim),bm1)

      ! restore
      time = time + dt
      call opcopy(vx,vy,vz,v_bak,v_bak(1,2),v_bak(1,ndim))

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_J_rhs_q(Jacp,p,Jeps,t_ref,ifld,ifwt)
c     Jq*p = [q(t_ref+eps*p)-q(t_ref)]/Jeps
      implicit none
      include 'SIZE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real Jacp(lt)
     $   , p(lt)
     $   , t_ref(lt)
     $   , tep(lt)
     $   , q_ref(lt)
     $   , q_tep(lt)
      real Jeps
      integer ifld,n
      logical ifwt

      n = nx1*ny1*nz1*nelt

      call add3s2(tep,t_ref,p,1.0,Jeps,n)

      call comp_Frhs_q(q_ref,t_ref,ifld,ifwt)
      call comp_Frhs_q(q_tep,tep,ifld,ifwt)

      call sub3(Jacp,q_tep,q_ref,n)
      call cmult(Jacp,1./Jeps,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_J_rhs_f(Jacp,p,Jeps,v_ref,ifwt)
c     Jf*p = [f(v_ref+eps*p)-f(v_ref)]/Jeps
c     ToDo: coordinates the size of working arrays
      implicit none
      include 'SIZE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real Jacp(lt,ldim)
     $   , p(lt,ldim)
     $   , v_ref(lt,ldim)
     $   , vep(lt,ldim)
     $   , f_ref(lt,ldim)
     $   , f_vep(lt,ldim)
      real Jeps
      integer ifld,n
      logical ifwt

      n = nx1*ny1*nz1*nelv

      call add3s2(vep(1,1),v_ref(1,1),p(1,1),1.0,Jeps,n)
      call add3s2(vep(1,2),v_ref(1,2),p(1,2),1.0,Jeps,n)
      if (ndim.eq.3)
     $   call add3s2(vep(1,ndim),v_ref(1,ndim),p(1,ndim),1.0,Jeps,n)

      call comp_Frhs_f(f_ref,v_ref,ifwt)
      call comp_Frhs_f(f_vep,vep,ifwt)

      call sub3(Jacp(1,1),f_vep(1,1),f_ref(1,1),n)
      call sub3(Jacp(1,2),f_vep(1,2),f_ref(1,2),n)
      call cmult(Jacp(1,1),1./Jeps,n)
      call cmult(Jacp(1,2),1./Jeps,n)

      if (ndim.eq.3) then
        call sub3(Jacp(1,ndim),f_vep(1,ndim),f_ref(1,ndim),n)
        call cmult(Jacp(1,ndim),1./Jeps,n)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine comp_J_visc_t(Jacp,p,Jeps,visc_ref,ifld,ifwt)
      implicit none
      include 'SIZE'
      include 'SOLN'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real Jacp(lt)
     $   , p(lt)
     $   , tep(lt)
     $   , visc_ref(lt)
     $   , visc_uep(lt)
      real Jeps
      integer ifld,n
      logical ifwt

      n = nx1*ny1*nz1*nelt

      call add3s2(tep,t(1,1,1,1,ifld-1),p,1.0,Jeps,n)
      call comp_visc_t(visc_uep,tep,ifld)
      call sub3(Jacp,visc_uep,visc_ref,n)
      call cmult(Jacp,1./Jeps,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_J_visc_v(Jacp,p,Jeps,visc_ref,ifld,ifwt)
      implicit none
      include 'SIZE'
      include 'SOLN'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelv)
      real Jacp(lt)
     $   , p(lt,ldim)
     $   , vxep(lt)
     $   , vyep(lt)
     $   , vzep(lt)
     $   , visc_ref(lt)
     $   , visc_uep(lt)
      real Jeps
      integer ifld,n
      logical ifwt

      n = nx1*ny1*nz1*nelv

      call add3s2(vxep,vx,p(1,1),1.0,Jeps,n)
      call add3s2(vyep,vy,p(1,2),1.0,Jeps,n)
      if(ndim.eq.3) call add3s2(vzep,vz,p(1,ndim),1.0,Jeps,n)
      call comp_visc_v(visc_uep,vx,vy,vz)
      call sub3(Jacp,visc_uep,visc_ref,n)
      call cmult(Jacp,1./Jeps,n)

      return
      end
c-----------------------------------------------------------------------
c     ToDo: move ten to somewhere
c-----------------------------------------------------------------------
      subroutine tens_fact(h1,h2,h3)
c     Computes Schur decomposition, assumes vxd,vyd,vzd
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'HSMGL'
      include 'TOTAL'

      real h1(*), h2(*), h3(*)

      common /flogmn/ ifrstr
      logical         ifrstr

      common /coarse/ acrs
      real acrs(lcr,lcr,lelt)

      common /cschur/ arank
      character*3 arank(lelt,lmgn,ldimt1)

      integer lwork,lsch
      parameter (
     $   lsch = ldimt1*lmgn*lelt*(lzm**2+(2*lzm+1)*(lxm**2+lym**2)),
     $   lwork = 2*(lzm+4)*(lxm+4)*(lym+4) )

      common /rschur/ S, QR, QL, WORK
      real S(lsch),QR(lsch),QL(lsch),WORK(lwork)

      common /ischur/ iptrs,ps
      integer iptrs(lelt,lmgn,ldimt1),ps

      integer IWORK(2*(lxm+lym+lzm))
      logical BWORK(2*(lxm+lym+lzm))

      real A(6*lxb*lxb), B(6*lyb*lyb), C(6*lzb*lzb)

      integer bcs(6,lelt)

      integer nx,ny,nz
      integer mx,my,mz
      integer lda,ldb,ldc

      integer n,nel
      integer e,l,fld
      integer p1,p2,slen,lmin
      integer two,ierr

      nel = nelt
      if(ifield.eq.1) nel=nelv
      n = nx1*ny1*nz1*nel

      if(ifield.eq.1) then
       do e=1,nel
         two  = 2
         ierr = 0
         call flmn_get_bc(    bcs(1,e) ,bcs(2,e)
     $                       ,bcs(3,e) ,bcs(4,e)
     $                       ,bcs(5,e) ,bcs(6,e)
     $                       ,e,two,ierr,ifrstr)
       enddo ! e
      else
       do e=1,nel
         two  = 2
         ierr = 0
         call adf_get_fast_bc(bcs(1,e) ,bcs(2,e)
     $                       ,bcs(3,e) ,bcs(4,e)
     $                       ,bcs(5,e) ,bcs(6,e)
     $                       ,e,two,ierr)
       enddo ! e
      endif

      call tens_init(h1,h2,bcs)

c     Assuming that setup starts at ifield = 1 or 2 
      if( (ifield.eq.1.and.ifflow) .or. 
     $    (ifield.eq.2.and.(.not.ifflow)) .or.
     $    (ifield.eq.2.and.iptrs(1,1,1).eq.0) ) ps=1
   
      lmin = 1
c      if( mg_nh(1).eq.2 ) lmin = 2 ! Corner mesh never done with tensor
      do l=mg_lmax,lmin,-1
      
c     Interior GLL nodes at current level
      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)
c     Basis functions in omega_bar (before Schwarz assembly)
      lda = nx + lxb-lx1
      ldb = ny + lyb-ly1
      ldc = nz + lzb-lz1
c     Extended DOFs in omega_bar (after Schwarz assembly)
      mx = max(nx, lda-2)
      my = max(ny, ldb-2)
      mz = max(nz, ldc-2)

      call tens_gen_basis(nx,ny,nz)
      
      if(nio.eq.0) write(6,*) 'Factoring Schur, level ',l,mx,my,mz
      if(ndim.eq.2) then
       slen = mx*mx + my*my
       do e=1,nel
        iptrs(e,l,ifield) = ps
        p1 = ps
        p2 = p1 + mx*mx
        call adf_svd(e,nx,ny,A,lda,B,ldb,h1,h2,h3,bcs(1,e),WORK)
        call schfact22(mx,my,S(p1),S(p2),QR(p1),QR(p2),QL(p1),QL(p2),
     $                 A,B,WORK,IWORK,BWORK)
        ps = ps + slen
       enddo ! e
      elseif(ndim.eq.3) then
       do e=1,nel
        iptrs(e,l,ifield) = ps
        p1 = ps
        call adf_als(e,arank(e,l,ifield),nx,ny,nz,A,lda,B,ldb,C,ldc,
     $                                   h1,h2,h3,bcs(1,e),WORK)
        call schfact3d(arank(e,l,ifield),mx,my,mz,S(p1),QR(p1),QL(p1),
     $                                   slen,A,B,C,WORK,IWORK,BWORK)
        ps = ps + slen
       enddo ! e
      endif ! ndim 
      enddo ! l
      return
      end
c-----------------------------------------------------------------------
      subroutine tens_solv(x,r,l,nx)
      implicit none
      include 'SIZE'
      include 'HSMGL'
      include 'TOTAL'
      include 'CTIMER'

      real x(nx**ndim,*), r(nx**ndim,*) ! 3D declaration
      integer l,nx

      common /cschur/ arank
      character*3 arank(lelt,lmgn,ldimt1)

      integer lwork,lsch
      parameter (
     $   lsch = ldimt1*lmgn*lelt*(lzm**2+(2*lzm+1)*(lxm**2+lym**2)),
     $   lwork = 2*(lzm+4)*(lxm+4)*(lym+4) )

      common /rschur/ S, QR, QL, WORK
      real S(lsch),QR(lsch),QL(lsch),WORK(lwork)

      common /ischur/ iptrs,ps
      integer iptrs(lelt,lmgn,ldimt1),ps

      integer p1,p2
      integer ny,nz
      integer mx,my,mz
      integer lda,ldb,ldc
      integer e,nel

      nel = nelt
      if(ifield.eq.1) nel = nelv

c     Interior GLL nodes at current level
      ny = mg_nh(l)
      nz = mg_nhz(l)
         
c     Basis functions in omega_bar (before Schwarz assembly)
      lda = nx + lxb-lx1
      ldb = ny + lyb-ly1
      ldc = nz + lzb-lz1

c     Extended DOFs in omega_bar (after Schwarz assembly)
      mx = max(nx, lda-2)
      my = max(ny, ldb-2)
      mz = max(nz, ldc-2)

      if (ndim.eq.2) then
       do e=1,nel
        p1 = iptrs(e,l,ifield)
        p2 = p1 + mx*mx
        call copy(x(1,e),r(1,e),mx*my)
        call schsolv22(mx,my,x(1,e),S(p1),S(p2),
     $       QR(p1),QR(p2),QL(p1),QL(p2),WORK)
       enddo ! e

      elseif(ndim.eq.3) then
       do e=1,nel
        p1 = iptrs(e,l,ifield)
        call schsolv3d(arank(e,l,ifield),mx,my,mz,x(1,e),
     $                        S(p1),QR(p1),QL(p1),r(1,e),WORK)
       enddo ! e
      endif ! ndim
      return
      end
c-----------------------------------------------------------------------
 
