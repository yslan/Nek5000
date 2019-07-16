c-----------------------------------------------------------------------
c     list of helper subroutines:
c     cask - compute A s_k (dssummed)
c     ccsk - compute C(u_k) s_k (+ C(s_k) u_k) (dssummed)
c     cjskf_pr - compute J(u_k) s_k with incomprn
c-----------------------------------------------------------------------
      subroutine sns_common ! all common blocks used

      logical ifeo, ifpsio

      common /snsvars/ alpha, tau, taui, taulag
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsio/ ifeo, ifpsio

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_init ! initialize variables to sensible values

      include 'SIZE'
      include 'TOTAL'

      logical ifeo, ifpsio

      common /snsvars/ alpha
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsio/ ifeo, ifpsio
      common /snsclock/ snstime, snstime_stp, ptime

      time = 0.
      istep = 0.
      ptime = 0.

      ! Li: why do we not make these param from rea?
      atoln = 1.e-20
      rtoln = 1.e-3

      atolp = 1.e-9

      atolg = 1.e-8
      rtolg = 1.e-6
      rtolg = 1.e-4

      atols = 1.e-9
      rtols = 1.

      itmaxg = 500
      itmaxn = 10
      itmaxp = 10000

      alpha = 1.

      ifield = 1

      ifeo = .false.
      ifpsio = .false.

      if (ifflow) then
         call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)
         call incomprn(vx,vy,vz,pr)
      endif

      if (ifheat) then
         ifield = 2
         call bcdirsc(t)
         ifield = 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_init_max ! initialize variables to sensible values

      include 'SIZE'
      include 'TOTAL'

      logical ifeo, ifpsio

      common /snsvars/ alpha
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsio/ ifeo, ifpsio

      time = 0.
      istep = 0.

      atoln = 1.e-12
      rtoln = 1.e-6

      atolp = 1.e-18

      rtolg = 1.e-8
      atolg = 1.e-12

      itmaxg = 500
      itmaxn = 20
      itmaxp = 100000

      alpha = 1.

      ifield = 1

      ifeo = .false.
      ifpsio = .false.

      if (ifflow) then
         call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)
         call incomprn(vx,vy,vz,pr)
      endif

      if (ifheat) then
         ifield = 2
         call bcdirsc(t)
         ifield = 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_fluid

      ! new JFNK method specifically for incompressible, steady N-S

      ! First compute r_k based on current solution, u_k
      ! Then solve J(u_k) s_k = -F_k = -P B^-1 r_k
      ! Use s_k to update: u_k+1 = u_k + alpha * s_k

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)
      parameter(l1=lx1*ly1*lz1)

      real sk(lt,ldima), rk(lt,ldima), fk(lt,ldima), ud(lt,ldima)

      logical ifpseudo, ifprint

      common /cprint/ ifprint
      common /snsvars/ alpha, tau, taui
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsl2/ el2,fl2,sl2
      common /snsflags/ ifpseudo
      common /snspvars/ up(lt,ldima)

      ifprint = .false.

      itn = min(max(itmaxn,10),10000)

      nel = nelv
      if (ifield.eq.2) nel = nelt

      n = lx1*ly1*lz1*nel
      isns=0
      tau = dt
      taui = 1 / tau

      fl2 = tol_newt*1.1

      ! zeroth step

      call opzero(sk,sk(1,2),sk(1,3))


      call cfk(rk)

      call checke
      call checkf(rk)
      call checks(sk)

      isns = 0
      istep = 0

      call sns_out

      ! end zeroth step

      do isns=1,itn
         istep = isns
         iter_gmres = 0

         call opchsgn(rk,rk(1,2),rk(1,3))
         call sns_gmres(sk,rk,vmult,tmult,it)
         iter_gmres = iter_gmres + it
         call opmask(sk,sk(1,2),sk(1,3))

         call opadds(vx,vy,vz,sk,sk(1,2),sk(1,3),alpha,n,2)
         call incomprn(vx,vy,vz,pr)

         call cfk(rk)

         call checke
         call checkf(rk)
         call checks(sk)

         call sns_out

         if (nio.eq.0) then
            if (fl2.le.atoln) write (6,*) 'reached atoln', fl2
            if (isns.eq.itn)  write (6,*) 'reached max itn', isns
            if (sl2.le.atols) write (6,*) 'reached atols', sl2
         endif

         if (fl2.le.atoln.or.isns.eq.itn.or.sl2.le.atols) then
            if (nio.eq.0) write (6,*) 'done with fluid - newton'
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cask(ask,sk) ! compute A s_k

      ! assumes scnv and setprop has been called

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ask(lt,ldim),sk(lt,ldim)
      real helm1(lt),helm2(lt)

      call sethlm(helm1,helm2,0)

      call axhelm(ask,sk,helm1,helm2,1,1)
      call axhelm(ask(1,2),sk(1,2),helm1,helm2,1,2)
      if (ldim.eq.3) call axhelm(ask(1,3),sk(1,3),helm1,helm2,1,3)

      return
      end
c-----------------------------------------------------------------------
      subroutine crkf(rk) ! compute the kth residual (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsflags/ ifpseudo, ifbouss

      real rk(lt,ldim), g(lt,ldim)
      real helm1(lt),helm2(lt)
      real t1(lt,ldim)

      n = lx1*ly1*lz1*nelv

      ifield = 1

      call rzero(rk,lt*ldim)
      call rzero(t1,lt)

      call setprop
      call sethlm(helm1,helm2,0)

      call scnv

      call convect_new(rk,u1v,.true.,c1v,c2v,c3v,.true.)
      call convect_new(rk(1,2),u2v,.true.,c1v,c2v,c3v,.true.)
      if (ldim.eq.3)
     $ call convect_new(rk(1,3),u3v,.true.,c1v,c2v,c3v,.true.)

      call axhelm(t1,vx,helm1,helm2,1,1)
      call axhelm(t1(1,2),vy,helm1,helm2,1,2)
      if (ldim.eq.3) call axhelm(t1(1,3),vz,helm1,helm2,1,3)

      do i=1,ldim
         call add2(rk(1,i),t1(1,i),n)
         call chsign(rk(1,i),n)
      enddo

      call makeuf

      call add2(rk,bfx,n)
      call add2(rk(1,2),bfy,n)
      if (ldim.eq.3) call add2(rk(1,3),bfz,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine scnv

      include 'SIZE'
      include 'TOTAL'

      parameter (ltd=lxd*lyd*lzd*lelt)

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      call set_convect_new(c1v,c2v,c3v,vx,vy,vz)

      call intp_rstd_all(u1v,vx,nelv)
      call intp_rstd_all(u2v,vy,nelv)
      if (ldim.eq.3) call intp_rstd_all(u3v,vz,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine cfkf(rk)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)
      parameter (lt1=lx1*ly1*lz1*lelt)
      parameter (lt2=lx2*ly2*lz2*lelt)
      parameter (ldima=ldim+ldimt)

      real rk(lt1,ldima), pt(lt2)

      n2 = lx2*ly2*lz2*nelv

      call scnv

      call copy(pt,pr,n2)

      call crkf(rk)
      call cpwf(rk)

c     call dssum(rk,lx1,ly1,lz1)
c     call dssum(rk(1,2),lx1,ly1,lz1)
c     if (ldim.eq.3) call dssum(rk(1,3),lx1,ly1,lz1)

      call copy(pr,pt,n2)


      return
      end
c-----------------------------------------------------------------------
      subroutine cfks(rk)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real rk(lt,ldimt)

      call crks(rk)
      call cpws(rk)

      return
      end
c-----------------------------------------------------------------------
      subroutine cfk(rk)

      include 'SIZE'
      include 'INPUT'

      real rk(lx1*ly1*lz1*lelt,ldim+ldimt)

      itemp = 1
      if (ifflow) itemp = ldim+1

      if (ifflow) call cfkf(rk)
      if (.not.ifflow) call scnv
      if (ifheat) call cfks(rk(1,itemp))

      return
      end
c-----------------------------------------------------------------------
      subroutine cjsk(jsk,sk)

      include 'SIZE'
c     include 'INPUT'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima = ldim + ldimt)

      common /snsfields/ fk(lt,ldima)
      common /snsvars/ alpha, tau, taui, taulag
      common /snstemp/ tk(lt,ldima)

      real jsk(lt,ldima), sk(lt,ldima)

      call checkvec(sk,sl2)

      eps = sqrt(1.e-16) / sl2

      n = lx1*ly1*lz1*nelv

      call snscopyto(tk)

      call snscopy(jsk,tk)
      call snsadd2s2(jsk,sk,eps)
      call snscopyfrom(jsk)

      call cfk(jsk)
      call snsadd2s2(jsk,fk,-1.)
      call snscmult2(jsk,jsk,1/eps)

      call snsadd2s2(jsk,sk,-1/tau)

      call snscopyfrom(tk)

      return
      end
c-----------------------------------------------------------------------
      subroutine cgk(rk)

      include 'SIZE'
      include 'INPUT'

      real rk(lx1*ly1*lz1*lelt,ldim+ldimt)

      itemp = 1
      if (ifflow) itemp = ldim+1

      if (ifflow) call cgkf(rk)
      if (ifheat) call cgks(rk(1,itemp))

      return
      end
c-----------------------------------------------------------------------
      subroutine cbuoy(buoy,therm)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /rayleigh_r/ rapr,ta2pr

      real buoy(lt,ldim), therm(lt), g(lt,ldim)

      n = lx1*ly1*lz1*nelv

      do i=1,ldim
         call rzero(buoy(1,i),n)
      enddo

      call add2s2(buoy(1,ldim),therm,rapr,n)
      call col2(buoy(1,ldim),bm1,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine cpwf(rk) ! compute the P_w W^-1 on r_k

      include 'SIZE'
      include 'SOLN'

      parameter(lt=lx1*ly1*lz1*lelt)

      real rk(lt,3), zero(lt)

      ifield = 1
      imesh = 1

      nit = nio
      nio = -1

      ! dssum is inside opbinv1
      call opbinv1(rk,rk(1,2),rk(1,3),rk,rk(1,2),rk(1,3),1.)

      if (ldim.eq.2) then
         call rzero(zero,lt)
         call incomprn(rk,rk(1,2),zero,pr)
      else
         call incomprn(rk,rk(1,2),rk(1,3),pr)
      endif

      call opmask(rk,rk(1,2),rk(1,ldim))

      nio = nit

      return
      end
c------------------------------------------------------------------------
      subroutine cwf(rk) ! compute the W^-1 on r_k

      include 'SIZE'
      include 'SOLN'

      parameter(lt=lx1*ly1*lz1*lelt)

      real rk(lt,ldim)

      ifield = 1
      imesh = 1

      nit = nio
      nio = -1

      ! dssum is inside opbinv1
      call opbinv1(rk,rk(1,2),rk(1,ldim),rk,rk(1,2),rk(1,ldim),1.)

      call opcolv(rk,rk(1,2),rk(1,ldim),vmult)

      call opmask(rk,rk(1,2),rk(1,ldim))

      nio = nit

      return
      end
c------------------------------------------------------------------------
      subroutine cjskf_pr(jsk,sk) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop

      call snscopy(t2,sk)
      call snsmask(t2)

      call update_conv(t2)

      call ccsk(jsk,t2,param(53).eq.0.)
      call cask(t3,t2)

      do i=1,ldim
         call add2(jsk(1,i),t3(1,i),nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      if (ifbouss) then
         ! Jacobian application result for the temperature field
         call convect_new(t1,t2(1,itemp),.false.,
     $                    c1v,c2v,c3v,.true.)
         call copy(jsk(1,itemp),t1,nv)
         call convect_new(t1,t,.false.,c1s,c2s,c3s,.true.)
         call add2(jsk(1,itemp),t1,nv)

         ifield = 2
         call sethlm(helm1,helm2,0)

         call axhelm(t1,t2(1,itemp),helm1,helm2,2,1)
         call add2(jsk(1,itemp),t1,nt)

         call dssum(jsk(1,itemp),lx1,ly1,lz1)
         call col2(jsk(1,itemp),tmask,nt)
         call chsign(jsk(1,itemp),nt)

         ifield = 1

         ! add buoyancy contribution
         call cbuoy(t3,t2(1,itemp))
         do i=1,ldim
            call add2(jsk(1,i),t3(1,i),nv)
         enddo

      endif

      call cpwf(jsk)

      if (ifbouss) call cpws(jsk(1,itemp))

      if (ifpseudo) then
         call opadds(jsk,jsk(1,2),jsk(1,3),
     $               t2,t2(1,2),t2(1,3),-taui,nv,2)
         if (ifbouss) call add2s2(jsk(1,itemp),t2(1,itemp),-taui,nt)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cjskf(jsk,sk) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call opcopy(t2,t2(1,2),t2(1,3),sk,sk(1,2),sk(1,3))
      call opmask(t2,t2(1,2),t2(1,3))

      if (ifbouss) then
         call copy(t2(1,itemp),sk(1,itemp),nt)
         call col2(t2(1,itemp),tmask,nt)
      endif

      do i=1,ldim
         call rzero(jsk(1,i),nv)
      enddo

      if (ifbouss) call rzero(jsk(1,itemp),nt)

      call update_conv(t2)
      call ccpsk(jsk,t2)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      if (ifbouss) then
         ! Jacobian application result for the temperature field
         call convect_new(t1,t2(1,itemp),.false.,
     $                    c1v,c2v,c3v,.true.)
         call copy(jsk(1,itemp),t1,nv)
         call convect_new(t1,t,.false.,c1s,c2s,c3s,.true.)
         call add2(jsk(1,itemp),t1,nv)

         ifield = 2
         call sethlm(helm1,helm2,0)

         call axhelm(t1,t2(1,itemp),helm1,helm2,2,1)
         call add2(jsk(1,itemp),t1,nt)

         call dssum(jsk(1,itemp),lx1,ly1,lz1)
         call col2(jsk(1,itemp),tmask,nt)
         call chsign(jsk(1,itemp),nt)

         ifield = 1

         ! add buoyancy contribution
         call cbuoy(t3,t2(1,itemp))
         do i=1,ldim
            call add2(jsk(1,i),t3(1,i),nv)
         enddo

      endif

c     call cwf(jsk)
      call cpwf(jsk)

      if (ifbouss) call cpws(jsk(1,itemp))

      if (ifpseudo) then
         call opadds(jsk,jsk(1,2),jsk(1,3),
     $               t2,t2(1,2),t2(1,3),-taui,nv,2)
         if (ifbouss) call add2s2(jsk(1,itemp),t2(1,itemp),-taui,nv)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_scalar

      ! new JFNK method specifically for incompressible, steady N-S

      ! First compute r_k based on current solution, u_k
      ! Then solve J(u_k) s_k = -F_k = -P B^-1 r_k
      ! Use s_k to update: u_k+1 = u_k + alpha * s_k

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real sk(lt,ldimt), rk(lt,ldimt)

      common /snsvars/ alpha
      common /snsl2/ el2, fl2
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg

      if (nio.eq.0) write (6,*) 'inside sns_scalar'

      itmax = min(max(10,itmaxn),10000)
      call setup_adf_visc

      n = lx1*ly1*lz1*nelt

      call setfq

      do isns=1,itmax
         istep = isns
         do i=1,1
            ifield=i+1
            if (nio.eq.0) write (6,*) 'inside loop'
            call cfk(rk(1,i))
            call sns_check(rk(1,i))
            call snschsgn(rk(1,i))
            call sns_gmres(
     $           sk(1,i),rk(1,i),vmult,tmult(1,1,1,1,i),iter_gmres)

            ! formulation: li
c           call adf_proj(sk(1,i),rk(1,i),h1,it,mxit,tmult(1,1,1,1,i))

            call col2(sk(1,i),tmask(1,1,1,1,i),n)
            call add2s2(t(1,1,1,1,i),sk(1,i),alpha,n)
         enddo

         call sns_out
         if (fl2.lt.atoln) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine crks(rk) ! compute the kth residual (scalar)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      real rk(lt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt)

      logical ifpseudo, ifbouss

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)
      common /snsflags/ ifpseudo, ifbouss

      n = lx1*ly1*lz1*nelt

      call rzero(rk,lt)
      call rzero(t1,lt)

      ifield = 2
      imesh = 2

      call setprop
      call sethlm(helm1,helm2,0)

      call copy(t2,t,n)

      call convect_new(rk,t2,.false.,c1v,c2v,c3v,.true.)

      call axhelm(t1,t2,helm1,helm2,2,1)

      call add2(rk,t1,n)

      call chsign(rk,n)

      call setqvol(bq)
      call col2(bq,bm1,n)

      call add2(rk,bq,n) ! Todo investigate non-zero bq
      call dssum(rk,lx1,ly1,lz1)
      call col2(rk,tmask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine cpws(rk) ! compute the P_w W^-1 on r_k

      include 'SIZE'
      include 'TOTAL'

      real rk(lx1*ly1*lz1*lelt)

      n = lx1*ly1*lz1*nelt

      call col2(rk,bintm1,n)
      call col2(rk,tmask,n)

      return
      end
c------------------------------------------------------------------------
      subroutine cjsks(jsk,sk) ! compute J(u_k) * s_k (coupled)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt), sk(lt)
      real helm1(lt),helm2(lt)
      real t1(lt), t2(lt)

      common /snsvars/ alpha, tau, taui
      common /snsflags/ ifpseudo, ifbouss
      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      n = lx1*ly1*lz1*nelt

c     call outpost(vx,vy,vz,pr,sk,'js1')

      call rzero(jsk,n)
      call rzero(t1,n)

      ifld = ifield
      ifield = 2

      call setprop
      call sethlm(helm1,helm2,0)

      call col3(t1,sk,tmask(1,1,1,1,ifield-1),n)

      nel = nelv
      if (ifield.eq.2) nel = nelt

      call convect_new(jsk,t1,.false.,c1v,c2v,c3v,.true.)

      call axhelm(t2,t1,helm1,helm2,2,1)
      call add2(jsk,t2,n)

      call dssum(jsk,lx1,ly1,lz1)
      call col2(jsk,tmask(1,1,1,1,ifield-1),n)
      call chsign(jsk,n)

      call cpws(jsk)

      if (ifpseudo) then
         call add2s2(jsk,t1,-taui,n)
      endif

      ifield = ifld

      return
      end
c------------------------------------------------------------------------
      subroutine setfq

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e

      real fx(lx1,ly1,lz1,lelt)
      real fy(lx1,ly1,lz1,lelt)
      real fz(lx1,ly1,lz1,lelt)
      real q(lx1,ly1,lz1,lelt)

      lxyz=lx1*ly1*lz1

      nv = lxyz*nelv
      nt = lxyz*nelt

      ! Todo replace loops with makeuf and setqvol

      do e=1,nelv
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         call userf(i,j,k,lglel(e))
         fx(i,j,k,e) = ffx
         fy(i,j,k,e) = ffy
         fz(i,j,k,e) = ffz
         bfx(i,j,k,e) = ffx * bm1(i,j,k,e)
         bfy(i,j,k,e) = ffy * bm1(i,j,k,e)
         bfz(i,j,k,e) = ffz * bm1(i,j,k,e)
      enddo
      enddo
      enddo
      enddo

      do e=1,nelt
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         call userq(i,j,k,lglel(e))
         q(i,j,k,e) = qvol
         bq(i,j,k,e,1) = qvol * bm1(i,j,k,e)
      enddo
      enddo
      enddo
      enddo

c     call outpost(fx,fy,fz,pr,q,'ffq')

      return
      end
c------------------------------------------------------------------------
      subroutine cgkf(fk)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      common /snspvars/ up(lt,ldim)
      common /snsvars/ alpha, tau, taui

      real fk(lt,ldim)
      real t1(lt,ldim)

      n = lx1*ly1*lz1*nelv

      call opsub3(t1,t1(1,2),t1(1,3),vxlag,vylag,vzlag,vx,vy,vz)

      do i=1,ldim
         call add2s2(fk(1,i),t1(1,i),taui,n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_pseudo

      ! new JFNK method specifically for incompressible, steady N-S

      ! First compute r_k based on current solution, u_k
      ! Then solve J(u_k) s_k = -F_k = -P B^-1 r_k
      ! Use s_k to update: u_k+1 = u_k + alpha * s_k

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      real sk(lt,ldima), gk(lt,ldima), ud(lt,ldima)
      real tk(lt,ldima), pk(lt,ldima)

      common /snsfields/ fk(lt,ldima)

      logical ifpseudo, ifprint, ifbouss, ifeo, ifpsio, ifprec

      common /snsvars/ alpha, tau, taui, taulag
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsl2/ el2,fl2,sl2,gl2
      common /snsflags/ ifpseudo, ifbouss
      common /snsio/ ifeo, ifpsio

      common /cprint/ ifprint

      n = lx1*ly1*lz1*nelv
      ifprint = .false.
      isns=0
      tau = dt
      taui = 1 / tau

      itn = min(max(itmaxn,10),1000000)
      itp = min(max(itmaxp,100),1000000)
      nel = nelv

      ifield = 1
      if (ifheat) ifield = 2
      if (ifheat.and.ifflow) ifbouss = .true.

      ifpseudo = .true.
      istep = 0

      call snslag

      call snszero(sk)

      ! zeroth step
      nt = lx1*ly1*lz1*nelt

      call cfk(fk)
      call snscopy(gk,fk)
      call cgk(gk)

      fl2lag = fl2
      call sns_check(fk)

      call sns_sett
      call sns_out

      istep = 1

      do ipseudo=1,itp
         call sns_sett
      do isns=1,itn
         call snschsgn(gk)

         call sns_gmres(sk,gk,vmult,tmult,it)
         iter_gmres = iter_gmres + it
         call snsmask(sk)

         call snscopyto(tk)
         call snsadd2s2(tk,sk,alpha)
         call snscopyfrom(tk)

         if (ifflow) call incomprn(vx,vy,vz,pr)

         ! create gk
         call cfk(fk)
         call snscopy(gk,fk)
         call cgk(gk)

         call checkg(gk)

         if (gl2/fl2.lt.rtoln) goto 100
c        if (gl2/fl2.lt.fl2/fl2lag) goto 100
c        if (gl2.lt.atoln) goto 100
      enddo
  100    continue
         time = time + tau
         fl2lag = fl2
         call sns_check(fk)

         call sns_out

         istep = ipseudo+1
         iter_gmres = 0

         taulag = tau
         tau = tau * fl2lag / fl2
         taui = 1 / tau

         call snslag

         if (fl2.lt.atolp.and.sl2.lt.atols) then
            call outpost(vx,vy,vz,pr,t,'   ')
            if (ifpsio) call sns_psi_omega
            return
         endif
      enddo

      if (nio.eq.0) write (6,*) 'done with fluid - newton'

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_exact

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical ifeo

      common /snsexact/ uex(lt,ldim), tex(lt,ldimt)
      common /snsio/ ifeo

      ifeo = .true.

      call opcopy(uex,uex(1,2),uex(1,3),vx,vy,vz)
      call copy(tex,t,lx1*ly1*lz1*nelt)
      call outpost(uex,uex(1,2),uex(1,3),pr,tex,'   ')

      return
      end
c-----------------------------------------------------------------------
      subroutine ccpsk(cpsk,sk) ! compute C'(u_k) * s_k

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real cpsk(lt,ldim), sk(lt,ldim)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real t2(lt,ldim)

      nel = nelv
      if (ifield.eq.2) nel = nelt

      call convect_new(cpsk(1,1),u1s,.true.,c1v,c2v,c3v,.true.)
      call convect_new(cpsk(1,2),u2s,.true.,c1v,c2v,c3v,.true.)
      if (ldim.eq.3)
     $ call convect_new(cpsk(1,3),u3s,.true.,c1v,c2v,c3v,.true.)

      call convect_new(t2,u1v,.true.,c1s,c2s,c3s,.true.)
      call convect_new(t2(1,2),u2v,.true.,c1s,c2s,c3s,.true.)
      if (ldim.eq.3)
     $ call convect_new(t2(1,3),u3v,.true.,c1s,c2s,c3s,.true.)

      call opadd2(cpsk,cpsk(1,2),cpsk(1,3),t2,t2(1,2),t2(1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine update_conv(sk)

      include 'SIZE'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real cpsk(lt,ldim), sk(lt,ldim)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      call set_convect_new(c1s,c2s,c3s,sk,sk(1,2),sk(1,3))

      call intp_rstd_all(u1s,sk,nelv)
      call intp_rstd_all(u2s,sk(1,2),nelv)
      if (ldim.eq.3) call intp_rstd_all(u3s,sk(1,3),nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd_all(uf,u,nel)

      include 'SIZE'
      include 'INPUT'

      parameter (lxyz1=lx1*ly1*lz1)
      parameter (lxyzd=lxd*lyd*lzd)

      real uf(lxyzd,lelt), u(lxyz1,lelt)

      do i=1,nel
         call intp_rstd(uf(1,i),u(1,i),lx1,lxd,if3d,0) ! 0 --> forward
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cpwb_pr(rk) ! compute the P_w W^-1 on r_k

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      real rk(lt,ldim)
      real zero(lt)

      nit = nio
      nio = -1

      call rzero(zero,lt)

      ifield = 1
      imesh = 1

      ! dssum is inside opbinv1
      if (ldim.eq.2) then
         call opbinv1(rk,rk(1,2),zero,rk,rk(1,2),zero,1.)
         call incomprn(rk,rk(1,2),zero,pr)
      else
         call opbinv1(rk,rk(1,2),rk(1,3),rk,rk(1,2),rk(1,3),1.)
         call incomprn(rk,rk(1,2),rk(1,3),pr)
      endif

      call opmask(rk,rk(1,2),rk(1,ldim))

      nio = nit

      return
      end
c------------------------------------------------------------------------
      subroutine cgks(fk)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      common /snsvars/ alpha, tau, taui

      real fk(lt,ldimt)
      real t1(lt,ldimt)

      nt = lx1*ly1*lz1*nelt

      call sub3(t1,tlag,t,nt)
      call add2s2(fk,t1,taui,nt)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_gmres(phi,res,wt,wtt,iter)

c     Solve the Helmholtz equation by non-preconditioned GMRES iteration.

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      logical ifpseudo, ifbouss, ifprec

      common /snsgmresc0/ x_NT(lt,ldima)
     $              , r_NT(lt,ldima), w_NT(lt,ldima)
     $              , h_NT(lgmres,lgmres), gamma_NT(lgmres+1)
     $              , c_NT(lgmres), s_NT(lgmres)
      common /snsgmresc1/ v_NT(lt,ldima,lgmres)
      common /snsgmresc2/ z_NT(lt,ldima,lgmres)
      common /snsflags/ ifpseudo, ifbouss, ifprec

      real x_NT,r_NT,w_NT,h_NT,gamma_NT,c_NT,s_NT,v_NT,z_NT

      common /snsspltprec/ ml_NT(lt,ldima)
     $              ,     mu_NT(lt,ldima)
      real ml_NT,mu_NT

      common /snsscr/ tmp(lt,ldima)
      real            tmp

c     w is a work vector
c     c and s store the Givens rotations
c     V stores the orthogonal Krylov subspace basis
c          -1
c     Z = M   V

      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsclock/ snstime, snstime_stp, ptime, ptimelag


      integer  n,outer,iflag,isep
      real     phi(lt,ldim+ldimt),res(lt,ldim+ldimt),wt(lt),wtt(lt)
      real     res_out,tol,alpha,l,temp
      real tmpk

      integer  i,j,k,iconv,iter,m
      real     rnorm,rinorm,tolpss


      n = lx1*ly1*lz1*nelv

      ifprec = .false.
c     if (param(51).gt.0.and.isns.gt.2) ifprec = .true.
      !  ! reason for keeping this? - Li
      if (param(51).gt.0) ifprec = .true.

      iter  = 0
      itemp = 1
      m     = lgmres

      tol = atolg
      rtol = rtolg

      tolps = tol
      tolpss= tolps
      iconv = 0

      itmax = itmaxg
      eps = 1.e-7

      if (nio.eq.0) write (6,*) 'starting Newton GMRES'
      if (nio.eq.0) write (6,*) 'tol,rtol',tol,rtol

      call rzero(x_NT,lt*ldima)
      call rzero(w_NT,lt*ldima)
      call rzero(h_NT,m*m)
      ptimelag=ptime

      if (ifflow) itemp = ldim + 1
      outer = 0
      do while (iconv.eq.0.and.iter.lt.itmax)
         outer = outer+1
         if (iter.eq.0) then
            call snscopy(r_NT,res) ! r = res
         else
            !update residual
            call snscopy(r_NT,res) ! r = res
            if (ifflow) then
               call cjskf_pr(w_nt,x_nt)
               if (ifheat) call cjsks(w_nt(1,itemp),x_nt(1,itemp))
            else
               call cjsks(w_nt,x_nt)
            endif

            call snsadd2s2(r_nt,w_nt,-1.)
         endif

         gamma_NT(1) = snsglsc3(r_NT,r_NT,wt,wtt)                ! gamma  = (r,r)
         gamma_NT(1) = sqrt(gamma_NT(1))                 ! gamma  = sqrt{ (r,r) }
         if(nid.eq.0) write (6,*) 'gamma_nt', gamma_nt(1)
         if(nid.eq.0) write (6,*) 'testing...'
         !check for lucky convergence
         rnorm = 0.
         if (gamma_NT(1) .eq. 0.) goto 9000
         temp = 1./gamma_NT(1)
         call snscmult2(v_NT,r_NT,temp)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1
            if (.not.ifprec) call snscopy(z_nt(1,1,j),v_nt(1,1,j))
            if (ifflow) then
               ! yes fluid preconditioner
               ifield = 1
               if (j.eq.1) then
                  if (ifprec) then
                     if (iter.gt.1) then
                        call snscopyto(tmp)
                        call opadd2(vx,vy,vz,x_nt,x_nt(1,2),x_nt(1,3))
                        call incomprn(vx,vy,vz,pr)
                     endif
                     tmpk = dnekclock()
                     call flprec(z_nt(1,1,j),v_nt(1,1,j),.true.)
                     ptime = ptime + dnekclock() - tmpk
                     if (iter.gt.1) call snscopyfrom(tmp)
                  endif
               else
                  if (ifprec) then
                     tmpk = dnekclock()
                     call flprec(z_nt(1,1,j),v_nt(1,1,j),.false.)
                     ptime = ptime + dnekclock() - tmpk
                  endif
               endif
c              if (iter.eq.10) call exitt0
               call cjskf_pr(w_nt,z_nt(1,1,j))
               if (ifheat) call cjsks(w_nt(1,itemp),z_nt(1,itemp,j))
            else
               ifield = 2
               if (ifprec) then
                  tmpk = dnekclock()
                  call scprec(z_nt(1,1,j),v_nt(1,1,j)) ! z = M v
                  ptime = ptime + dnekclock() - tmpk
               endif
               call cjsks(w_nt,z_nt(1,1,j))        ! w = A z
            endif

c           write (6,*) snsglsc3(w_nt,w_nt,wt,wtt)

            !modified Gram-Schmidt
            do i=1,j
               h_NT(i,j)=snsglsc3(w_NT,v_NT(1,1,i),wt,wtt)        ! h    = (w,v )
               call snsadd2s2(w_nt,v_nt(1,1,i),-h_nt(i,j))
            enddo                                 !         i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_NT(i,j)
               h_NT(i  ,j)=  c_NT(i)*temp + s_NT(i)*h_NT(i+1,j)
               h_NT(i+1,j)= -s_NT(i)*temp + c_NT(i)*h_NT(i+1,j)
            enddo
                                                    !            ______
            alpha = sqrt(snsglsc3(w_NT,w_NT,wt,wtt))     ! alpha =  \/ (w,w)
c           write (6,*) 'alpha',iter,alpha
            if (alpha.eq.0.) goto 900 !converged
            l = sqrt(h_NT(j,j)*h_NT(j,j)+alpha*alpha)
            temp = 1./l
            c_NT(j) = h_NT(j,j) * temp
            s_NT(j) = alpha  * temp
            h_NT(j,j) = l
            gamma_NT(j+1) = -s_NT(j) * gamma_NT(j)
            gamma_NT(j)   =  c_NT(j) * gamma_NT(j)

            rlag = rnorm
            rnorm = abs(gamma_NT(j+1))
            if (iter.eq.1) rinorm=abs(gamma_nt(1))

c            if ((nid.eq.0).and.(ipstep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,ipstep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')
            ratio = rnorm/rinorm
            if (nio.eq.0) then
c              write (6,1234) iflag, isns, iter, tol,
               write (6,1234) isns, iter, rnorm, rinorm,
     $         ratio, alpha,tol,rtol
            endif
 1234 format('nt_gmresf',i5,i6,1p6e12.4)

            if (rnorm .lt. tol) goto 900 !converged
            if (ratio .lt. rtol) goto 900 !converged
c           if (iter.gt.1.and.rnorm/rlag.gt.(1-eps)) goto 900 !converged

            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call snscmult2(v_NT(1,1,j+1),w_NT,temp)   ! v    = w / alpha
                                             !  j+1
         enddo ! end j=1,m loop
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
         res_out = rnorm   ! output |Ax-b|
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma_NT(k)
            do i=j,k+1,-1
               temp = temp - h_NT(k,i)*c_NT(i)
            enddo
            c_NT(k) = temp/h_NT(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
c           call snsadds(x_nt,x_nt(1,2),x_nt(1,3),
c    $                  v_nt(1,1,i),v_nt(1,2,i),v_nt(1,3,i),c_nt(i),n,2)
            call snsadd2s2(x_NT,z_NT(1,1,i),c_NT(i))     ! x = x + c  z
         enddo                                        !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo ! end while loop
 9000 continue

      call snscopy(phi,x_nt)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

c     if ((nid.eq.0).and. (mod(ipstep,iocomm).eq.0) ) then
c         write(6,9999) ipstep,iterNT,iter,tolpss
c     endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end
c------------------------------------------------------------------------
      subroutine snscopy(a,b)

      include 'SIZE'
      include 'INPUT'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      real a(lt,ldima),b(lt,ldima)

      itemp = 1

      nt=lx1*ly1*lz1*nelt

      if (ifflow) then
         itemp = ldim + 1
         call opcopy(a,a(1,2),a(1,ldim),b,b(1,2),b(1,ldim))
      endif

      if (ifheat) call copy(a(1,itemp),b(1,itemp),nt)

      return
      end
c------------------------------------------------------------------------
      subroutine snsmask(a,b)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      real a(lt,ldima),b(lt,ldima)

      nt=lx1*ly1*lz1*nelt

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opmask(a,a(1,2),a(1,ldim))
      endif

      if (ifheat) call col2(a(1,itemp),tmask,nt)

      return
      end
c------------------------------------------------------------------------
      subroutine snsadd2s2(a,b,c)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      real a(lt,ldima), b(lt,ldima), c

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opadds(a,a(1,2),a(1,ldim),b,b(1,2),b(1,ldim),c,nv,2)
      endif

      if (ifheat) call add2s2(a(1,itemp),b(1,itemp),c,nt)

      return
      end
c------------------------------------------------------------------------
      subroutine snscmult2(a,b,const)             !  v  = r / gamma

      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real a(lt,ldim+ldimt),b(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opcmult2(a,b,const)
      endif

      if (ifheat) call cmult2(a(1,itemp),b(1,itemp),const,lxyz*nelt)

      return
      end
c------------------------------------------------------------------------
      subroutine snschsgn(a)

      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real a(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opchsgn(a,a(1,2),a(1,ldim))
      endif

      if (ifheat) call chsign(a(1,itemp),lxyz*nelt)

      return
      end
c------------------------------------------------------------------------
      subroutine snslag

      include 'SIZE'
      include 'TOTAL'

      if (ifflow) call lagvel
      if (ifheat) then
         ifldt = ifield
         ifield = 2
         call lagscal
         ifield = ifldt
      endif

      return
      end
c------------------------------------------------------------------------
      subroutine snszero(a)

      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real a(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opzero(a,a(1,2),a(1,ldim))
      endif

      if (ifheat) call rzero(a(1,itemp),lxyz*nelt)

      return
      end
c------------------------------------------------------------------------
      subroutine snscopyto(a)

      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real a(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opcopy(a,a(1,2),a(1,ldim),vx,vy,vz)
      endif

      if (ifheat) call copy(a(1,itemp),t,lxyz*nelt)

      return
      end
c------------------------------------------------------------------------
      subroutine snscopyfrom(a)

      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real a(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         call opcopy(vx,vy,vz,a,a(1,2),a(1,ldim))
      endif

      if (ifheat) call copy(t,a(1,itemp),lxyz*nelt)

      return
      end
c------------------------------------------------------------------------
      subroutine checkf(sk)

      include 'SIZE'

      real sk(lx1*ly1*lz1*lelt,ldim+ldimt)

      common /snsl2/ el2,fl2,sl2

      call checkvec(sk,fl2)

      return
      end
c-----------------------------------------------------------------------
      subroutine checkg(sk)

      include 'SIZE'

      common /snsl2/ el2,fl2,sl2,gl2
      real sk(lx1*ly1*lz1*lelt,ldim+ldimt)

      call checkvec(sk,gl2)

      return
      end
c-----------------------------------------------------------------------
      subroutine checke

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      logical ifbouss, ifeo

      common /snsio/ ifeo
      common /snsl2/ el2,fl2
      common /snsexact/ uex(lt,ldim), tex(lt,ldimt)
      common /snsflags/ ifpseudo, ifbouss

      real ud(lt),vd(lt),wd(lt),td(lt)

      nt = lx1*ly1*lz1*nelt
      itemp = ldim + 1

      call opsub3(ud,vd,wd,uex,uex(1,2),uex(1,3),vx,vy,vz)
      call sub3(td,tex,t,nt)

      if (ifbouss) then
         call opnorm(el2,ud,vd,wd,'L2 ')
         el2 = sqrt(glsc3(td,td,bm1,nt)/voltm1 + el2*el2)
      else
         if (ifield.eq.1) call opnorm(el2,ud,vd,wd,'L2 ')
         if (ifield.eq.2) el2 = sqrt(glsc3(td,td,bm1,nt)/voltm1)
      endif

      if (ifeo) call outpost(ud,vd,wd,pr,td,'err')

      if (nio.eq.0) write (6,*) 'l2 error', istep, el2

      return
      end
c-----------------------------------------------------------------------
      subroutine checkd(dl2)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      logical ifbouss

      real ud(lt),vd(lt),wd(lt),td(lt)

      nt = lx1*ly1*lz1*nelt
      itemp = ldim + 1

      call opsub3(ud,vd,wd,vx,vy,vz,vxlag,vylag,vzlag)
      call sub3(td,t,tlag,nt)

      if (ifbouss) then
         call opnorm(dl2,ud,vd,wd,'L2 ')
         el2 = sqrt(glsc3(td,td,bm1,nt)/voltm1 + dl2*dl2)
      else
         if (ifield.eq.1) call opnorm(dl2,ud,vd,wd,'L2 ')
         if (ifield.eq.2) dl2 = sqrt(glsc3(td,td,bm1,nt)/voltm1)
      endif

      if (nio.eq.0) write (6,*) 'dl2', istep, time, dl2

      return
      end
c-----------------------------------------------------------------------
      subroutine checks(sk)

      include 'SIZE'

      common /snsl2/ el2,fl2,sl2
      real sk(lx1*ly1*lz1*lelt,ldim+ldimt)

      call checkvec(sk,sl2)

      return
      end
c-----------------------------------------------------------------------
      subroutine checkf_r(sk)

      include 'SIZE'
      include 'TOTAL'

      logical ifpseudo, ifbouss

      common /snsl2/ el2,fl2,sl2,fl2_r
      common /snsflags/ ifpseudo, ifbouss

      real sk(lx1*ly1*lz1*lelt,ldim+ldimt)

      itemp = ldim + 1

      if (ifbouss) then
         call opnorm(ul2,vx,vy,vz,'L2 ')
         call opnorm(fl2_r,sk,sk(1,2),sk(1,3),'L2 ')
         fl2_r = fl2_r / ul2
         fl2_r = fl2_r + sqrt(glsc3(sk(1,itemp),sk(1,itemp),bm1,nt)/
     $                        glsc3(t,t,bm1,nt))
         fl2_r = t1 + t2
      else
         if (ifield.eq.1) then
            call opnorm(fl2_r,sk,sk(1,2),sk(1,3),'L2 ')
            call opnorm(ul2,vx,vy,vz,'L2 ')
            fl2_r = fl2_r / ul2
         else if (ifield.eq.2) then
            fl2_r = sqrt(glsc3(sk,sk,bm1,nt)/glsc3(t,t,bm1,nt))
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine checkvec(sk,dl2)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      logical ifbouss

      common /snsflags/ ifpseudo, ifbouss

      real sk(lt,ldim+ldimt)

      nt = lx1*ly1*lz1*nelt
      itemp = 1

      dl2 = 0.

      if (ifflow) then
         itemp = ldim + 1
         call opnorm(dl2,sk,sk(1,2),sk(1,ldim),'L2 ')
      endif

      if (ifheat) dl2 = sqrt(glsc3(sk(1,itemp),sk(1,itemp),bm1,nt)/
     $            voltm1) + dl2

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_out

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      logical ifpseudo, ifeo, ifpsio

      common /snsl2/ el2,fl2,sl2,fl2_r
      common /snsivars/ isns, iter_gmres
      common /snsvars/ alpha, tau, taui
      common /snsflags/ ifpseudo
      common /snsio/ ifeo, ifpsio
      common /snsclock/ snstime, snstime_stp, ptime, ptimelag

      if (mod(istep,max(iostep,1)).eq.0) then
         call outpost(vx,vy,vz,pr,t,'   ')
         if (ifpsio) call sns_psi_omega
      endif

      if (nio.eq.0.and.ifpseudo) then
         write (6,4000) istep, time, tau, fl2, snstime, snstime_stp
         write (6,3000) istep, isns, iter_gmres, tau, el2, fl2, sl2
      else
         dpt=ptime-ptimelag
         if (nio.eq.0.and.ifield.eq.1) then
            write (6,1000) isns, iter_gmres,el2,fl2,sl2
            write (6,5000) isns, fl2, snstime, snstime_stp, ptime, dpt
         endif
         if (nio.eq.0.and.ifield.eq.2) then
            write (6,2000) isns, iter_gmres,el2,fl2,sl2
            write (6,5000) isns, fl2, snstime, snstime_stp, ptime, dpt
         endif
      endif

c     el2 L2 norm of the error compared to exact solution: uex and tex.
c     fl2 L2 norm of the residual vector.
c     sl2 L2 norm of the update vector.
c     
c     el2 scaled by the inverse of the L2 norm of the solution
c     fl2 and sl2 are scaled by the inverse of domain size.
      
 1000 format('sns_outf',i5,i8,1p3e12.4)
 2000 format('sns_outt',i5,i8,1p3e12.4)
 3000 format('sns_outp',i7,i5,i7,1p4e12.4)
 4000 format('Step',i7,', t=',1pe11.4,', Tau=',1pe11.4,
     $       ', F=',1p3e11.4)
 5000 format('Step',i7,', F=',1p5e11.4) 

      call check_ioinfo
      if (ioinfodmp.ne.0) then
         if (nio.eq.0) write (6,*) 'exit due to ioinfo...'
         call exitt0
      endif

      return
      end
c-----------------------------------------------------------------------
      function snsglsc3(a,b,mult,multt)
C
C     Perform inner-product in double precision
C
      include 'SIZE'
      include 'INPUT'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      real a(lt,ldima),b(lt,ldima),mult(lt)
      real tmp,work(1)

      tmp = 0.
      tmp2 = 0.

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         tmp = opglsc3(a,b,mult)
      endif

      if (ifheat) then
         tmp2 = glsc3(a(1,itemp),b(1,itemp),multt,lx1*ly1*lz1*nelt)
         tmp = tmp + tmp2
      endif

c     call gop(tmp,work,'+  ',1)
      snsglsc3 = tmp

      return
      end
c-----------------------------------------------------------------------
      function opglsc3(a,b,mult)
C
C     Perform inner-product in double precision
C
      include 'SIZE'

      parameter(lt=lx1*ly1*lz1*lelt)

      real a(lt,ldim),b(lt,ldim),mult(lt)
      real tmp,work(1)

      n = lx1*ly1*lz1*nelv

      tmp = 0.
      do i=1,n
         tmp = tmp + (a(i,1)*b(i,1)+a(i,2)*b(i,2))*mult(i)
         if (ldim.eq.3) tmp = tmp + a(i,3)*b(i,3)*mult(i)
      enddo

      call gop(tmp,work,'+  ',1)
      opglsc3 = tmp

      return
      end
c-----------------------------------------------------------------------
      subroutine opcmult2(a,b,const)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      real a(lt,ldim),b(lt,ldim)

      do i=1,lx1*ly1*lz1*nelv
         a(i,1)=b(i,1)*const
         a(i,2)=b(i,2)*const
         if (ldim.eq.3) a(i,3)=b(i,3)*const
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine snsdssum(sk)

      include 'SIZE'
      include 'INPUT'

      parameter (lt=lx1*ly1*lz1*lelt)

      real sk(lt,ldim+ldimt)

      real zero(lt)

      itemp = 1
      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      if (ifflow) then
         itemp = ldim + 1
         if (ldim.eq.2) then
            call opdssum(sk,sk(1,2),zero)
         else
            call opdssum(sk,sk(1,2),sk(1,3))
         endif
      endif

      if (ifheat) call dssum(sk(1,itemp),lx1,ly1,lz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_psi_omega

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iftemp

      real psi(lt), omega(lt,3), rhs(lt)
      real w1(lt),w2(lt),h1(lt),h2(lt)

      n = lx1*ly1*lz1*nelv

      call comp_vort3(omega,w1,w2,vx,vy,vz)
      call col3(rhs,bm1,omega,n)
      call rone(h1,n)
      call rzero(h2,n)
      tol = param(22)
c     call chsign(rhs,n)
      call hmholtz('psi ',psi,rhs,h1,h2,v1mask,vmult,1,tol,1000,1)
      call dsavg(psi)

      iftemp = ifxyo
      ifxyo = .true.
      call outpost(psi,omega,vz,pr,t,'psi')
      ifxyo = iftemp

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_check(fk)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)

      common /snsl2/ el2,fl2,sl2,gl2

      real ud(lt,ldima), fk(lt,ldima)

      call checke
      call checkf(fk)

      itemp = 1
      if (ifflow) then
         itemp = ldim + 1
         call opsub3(ud,ud(1,2),ud(1,ldim),vxlag,vylag,vzlag,vx,vy,vz)
      endif
      if (ifheat) call sub3(ud(1,itemp),tlag,t,lx1*ly1*lz1*nelt)

      call checks(ud)

      return
      end
c-----------------------------------------------------------------------
      subroutine cbad(jsk,sk) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call snscopy(t2,sk)
      call snsmask(t2)

      call rzero(jsk,nv)
      call rzero(jsk(1,2),nv)

      call update_conv(t2)
      call ccpsk_bad(jsk,t2)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      call cwf(jsk)

      return
      end
c-----------------------------------------------------------------------
      subroutine cjac(jsk,sk) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)
      real zero(lt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call snscopy(t2,sk)
      call snsmask(t2)

      call rzero(jsk,nv)
      call rzero(jsk(1,2),nv)

      call update_conv(t2)
      call ccpsk(jsk,t2)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      ! cpwf
      call rzero(vz,nv)
      call opbinv1(jsk,jsk(1,2),nv,jsk,jsk(1,2),nv,1.)

      if (ldim.eq.2) then
         call rzero(zero,lt)
         call incomprn(jsk,jsk(1,2),zero,pr)
      else
         call incomprn(jsk,jsk(1,2),jsk(1,3),pr)
      endif

      call opmask(jsk,jsk(1,2),jsk(1,ldim))

c     call dssum(jsk,lx1,ly1,lz1)
c     call dssum(jsk(1,2),lx1,ly1,lz1)

c     call cpwf(jsk)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_sett
C---------------------------------------------------------------------
C
C     No need to comment !!
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
      include 'CTIMER'

      real*8 eetime0,eetime1,eetime2
      save   eetime0,eetime1,eetime2
      data   eetime0,eetime1,eetime2 /0.0, 0.0, 0.0/

      common /snsclock/ snstime, snstime_stp

C     Only node zero makes comments.
      if (nio.ne.0) return

      if (istep.eq.0) eetime0=dnekclock()

      eetime1=eetime2
      eetime2=dnekclock()

      if (istep.gt.0) then
         snstime_stp = eetime2-eetime1   ! time per timestep
         snstime     = eetime2-eetime0   ! sum of all timesteps
         if (istep.eq.0) then
           snstime_stp = 0
           snstime     = 0
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine sns_checka

      include 'SIZE'

      real fk(lx1*ly1*lz1*lelt,ldim+ldimt)

      call cfk(fk)
      call sns_check(fk)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_checku

      include 'SIZE'
      include 'TSTEP'

      logical ifeo, ifpsio

      common /snsio/ ifeo, ifpsio
      common /snsl2/ el2, fl2, sl2

      if (mod(istep,max(iostep,1)).eq.0.or.istep.eq.nsteps) then
         call sns_checka
         if (ifpsio) call sns_psi_omega
         if (nio.eq.0) write (6,*) 'fdl2',istep,time,fl2,sl2
      endif

      return
      end
c-----------------------------------------------------------------------
      function l2dof(idof)

      ! returns an integer that corresponds to a real dof

      include 'SIZE'

      l2dof = idof + (idof-1) / (lx1-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setdof(idof,ndof)

      include 'SIZE'
      include 'TOTAL'

      integer idof(lx1*ly1*lz1,lelt)
      integer jdof(lx1*ly1*lz1,lelt)
      integer kdof(lx1*ly1*lz1,lelt)
c     integer*8 glo_num(lx1*ly1*lz1,lelt), ngv
      integer e

      ndof = 0

      n = lx1*ly1*lz1*nelt

      call izero(idof,n)
      call izero(jdof,n)
      call izero(kdof,n)

c     getgln(glo_num,ngv)

      do e=1,nelt
      do i=1,lx1*ly1*lz1
         jdof(i,e) = e
         write (6,*) e
      enddo
      enddo

      do i=1,lx1*ly1*lz1*nelt
         if (v1mask(i,1,1,1).eq.0.) jdof(i,1) = 0
      enddo

      write (6,*) 'end of initialization'

      call idsop(jdof,'m  ',lx1,ly1,lz1)

      icount = 0

      do e=1,nelt
      do i=1,lx1*ly1*lz1
         write (6,*) jdof(i,e)
         if (jdof(i,e).eq.e) then
            icount = icount + 1
            kdof(i,e) = icount
         endif
      enddo
      enddo

      write (6,*) 'end of setting'

      do i=1,lx1*ly1*lz1*nelt
         write (6,*) kdof(i,1)
      enddo

      write (6,*) 'before dssum'

      icount = 0

      do i=1,n
         icount = icount + 1
         if (kdof(i,1).ne.0) then
            idof(kdof(i,1),1) = icount
         endif
      enddo

      write (6,*) 'result'
      ndof = iglmax(kdof,n)

      do i=1,ndof
         write (6,*) idof(i,1)
      enddo

      write (6,*) 'ndof = ',ndof

      return
      end
c------------------------------------------------------------------------
      subroutine setdofp(idof,ndof)

      include 'SIZE'
      include 'TOTAL'

      integer idof(lx1*ly1*lz1*lelt)
      integer*8 jdof(lx1*ly1*lz1*lelt)
      integer kdof(lx1*ly1*lz1*lelt)

      ndof = 0

      n = lx1*ly1*lz1*nelt

      call izero(idof,n)
      call izero(kdof,n)

      call usrsetvert(jdof,1,lx1,ly1,lz1)

      do i=1,n
         if (v1mask(i,1,1,1).ne.0.) kdof(i) = jdof(i)
      enddo

      write (6,*) 'result'

      icount = 1
      jmark = 0

      do i=1,n
         if (kdof(i).gt.jmark) then
            idof(icount) = i
            jmark = kdof(i)
            icount = icount + 1
         endif
      enddo

      do i=1,lx1*ly1*lz1*nelt
         write (6,*) 'idof',i,idof(i)
      enddo

      ndof = icount - 1

      write (6,*) 'ndof = ',ndof

      return
      end
c------------------------------------------------------------------------
      subroutine idsop(u,op,nx,ny,nz)
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TSTEP'
      include  'CTIMER'

      integer u(1)
      character*3 op
      character*10 s1,s2
c
c     o gs recognized operations:
c
c             o "+" ==> addition.
c             o "*" ==> multiplication.
c             o "M" ==> maximum.
c             o "m" ==> minimum.
c             o "A" ==> (fabs(x)>fabs(y)) ? (x) : (y), ident=0.0.
c             o "a" ==> (fabs(x)<fabs(y)) ? (x) : (y), ident=MAX_DBL
c             o "e" ==> ((x)==0.0) ? (y) : (x),        ident=0.0.
c
c             o note: a binary function pointer flavor exists.
c
c
c     o gs level:
c
c             o level=0 ==> pure tree
c             o level>=num_nodes-1 ==> pure pairwise
c             o level = 1,...num_nodes-2 ==> mix tree/pairwise.
c
c
      ifldt = ifield
c     if (ifldt.eq.0)       ifldt = 1
      if (ifldt.eq.ifldmhd) ifldt = 1

c     if (nio.eq.0)
c    $   write(6,*) istep,' dsop: ',op,ifield,ifldt,gsh_fld(ifldt)

      if(ifsync) call nekgsync()

      if (op.eq.'+  ') call fgslib_gs_op(gsh_fld(ifldt),u,2,1,0)
      if (op.eq.'sum') call fgslib_gs_op(gsh_fld(ifldt),u,2,1,0)
      if (op.eq.'SUM') call fgslib_gs_op(gsh_fld(ifldt),u,2,1,0)

      if (op.eq.'*  ') call fgslib_gs_op(gsh_fld(ifldt),u,2,2,0)
      if (op.eq.'mul') call fgslib_gs_op(gsh_fld(ifldt),u,2,2,0)
      if (op.eq.'MUL') call fgslib_gs_op(gsh_fld(ifldt),u,2,2,0)

      if (op.eq.'m  ') call fgslib_gs_op(gsh_fld(ifldt),u,2,3,0)
      if (op.eq.'min') call fgslib_gs_op(gsh_fld(ifldt),u,2,3,0)
      if (op.eq.'mna') call fgslib_gs_op(gsh_fld(ifldt),u,2,3,0)
      if (op.eq.'MIN') call fgslib_gs_op(gsh_fld(ifldt),u,2,3,0)
      if (op.eq.'MNA') call fgslib_gs_op(gsh_fld(ifldt),u,2,3,0)

      if (op.eq.'M  ') call fgslib_gs_op(gsh_fld(ifldt),u,2,4,0)
      if (op.eq.'max') call fgslib_gs_op(gsh_fld(ifldt),u,2,4,0)
      if (op.eq.'mxa') call fgslib_gs_op(gsh_fld(ifldt),u,2,4,0)
      if (op.eq.'MAX') call fgslib_gs_op(gsh_fld(ifldt),u,2,4,0)
      if (op.eq.'MXA') call fgslib_gs_op(gsh_fld(ifldt),u,2,4,0)
c
      return
      end
c------------------------------------------------------------------------
      subroutine sns_extrap(ek)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /snsvars/ alpha, tau, taui, taulag

      real ek(lt,ldim+ldimt)

      itemp = 1

      if (ifflow) then
         itemp = ldim + 1
         nv = lx1*ly1*lz1*nelv
         call opcopy(ek,ek(1,2),ek(1,3),vxlag,vylag,vzlag)
         call opcmult(ek,ek(1,2),ek(1,3),-tau/taulag)
         call opadds(ek,ek(1,2),ek(1,3),vx,vy,vz,1+tau/taulag,nv,2)
      endif

      if (ifheat) then
         nt = lx1*ly1*lz1*nelt
         call add3s2(ek(1,itemp),t,tlag,1+tau/taulag,-tau/taulag,nt)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ccpsk_bad(cpsk,sk) ! compute C'(u_k) * s_k

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real cpsk(lt,ldim), sk(lt,ldim)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real t2(lt,ldim)

      nel = nelv
      if (ifield.eq.2) nel = nelt

      call convect_new(cpsk(1,1),u1s,.true.,c1v,c2v,c3v,.true.)
      call convect_new(cpsk(1,2),u2s,.true.,c1v,c2v,c3v,.true.)
      if (ldim.eq.3)
     $ call convect_new(cpsk(1,3),u3s,.true.,c1v,c2v,c3v,.true.)

      return
      end
c-----------------------------------------------------------------------
      subroutine cddd(jsk,sk) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call snscopy(t2,sk)
      call snsmask(t2)

      call rzero(jsk,nv)
      call rzero(jsk(1,2),nv)

      call update_conv(t2)
      call ccpsk(jsk,t2)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      call cwf(jsk)

      return
      end
c-----------------------------------------------------------------------
      subroutine ceee(jsk,sk,ifcorrect) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss, ifcorrect

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call snscopy(t2,sk)
      call snsmask(t2)

      call rzero(jsk,nv)
      call rzero(jsk(1,2),nv)

      call update_conv(t2)
      call ccpsk_bad(jsk,t2)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      call dssum(jsk,lx1,ly1,lz1)
      call dssum(jsk(1,2),lx1,ly1,lz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine cadv(jsk,sk,ifcorrect) ! compute J(u_k) * s_k (fluid)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(lt2=lx2*ly2*lz2*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      logical ifpseudo, ifbouss, ifcorrect

      real jsk(lt,ldim+ldimt), sk(lt,ldim+ldimt)
      real helm1(lt),helm2(lt)
      real t1(lt)
      real t2(lt,ldim+ldimt)
      real t3(lt,ldim+ldimt)
      real t4(lxd*lyd*lzd*lelt)

      common /snsflags/ ifpseudo, ifbouss
      common /snsvars/ alpha, tau, taui

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      ifield = 1
      itemp = ldim + 1

      call setprop
      call sethlm(helm1,helm2,0)

      call snscopy(t2,sk)
      call snsmask(t2)

      call rzero(jsk,nv)
      call rzero(jsk(1,2),nv)

      call update_conv(t2)
      call ccsk(jsk,t2,ifcorrect)

      do i=1,ldim
         call axhelm(t1,t2(1,i),helm1,helm2,1,i)
         call add2(jsk(1,i),t1,nv)
      enddo

      call opmask(jsk,jsk(1,2),jsk(1,3))
      call opchsgn(jsk,jsk(1,2),jsk(1,3))

      call dssum(jsk,lx1,ly1,lz1)
      call dssum(jsk(1,2),lx1,ly1,lz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine assemble(a,aloc,idof,n,nd,ndof)

      include 'SIZE'

      real a(n,nd,n,nd)
      real aloc(n,nd,n,nd)
      integer idof(ndof)

      do j2=1,nd
      do j=1,ndof
      do i2=1,nd
      do i=1,ndof
         a(i,i2,j,j2) = aloc(idof(i),i2,idof(j),j2)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine genmat2(aloc,imat,ifexact)

      include 'SIZE'
      include 'TOTAL'
      include 'PREC'

      parameter (l1=lx1*ly1*lz1)
      parameter (lt=l1*lelt)

      real aloc(lk*nelt,2,lk*nelt,2)
      real p3(lt,2)
      real p4(lt,2)
      real zero(lt)

      logical ifexact

      n = lx1*ly1*lz1*nelt

      do i2=1,2
      do i=1,n
         call rzero(p3,n)
         call rzero(p3(1,2),n)
         call rzero(zero,n)
         p3(i,i2) = 1.
         call opdssum(p3,p3(1,2),zero)
         call update_conv(p3)

         if (imat.eq.0) then
            call cask(p4,p3)
         else if(imat.eq.1) then
            call ccsk(p4,p3,ifexact)
         else if (imat.eq.2) then
            call cadv(p4,p3,ifexact)
         else if (imat.eq.3) then
            call cadv(p4,p3,ifexact)
            call col2(p4,binvm1,n)
            call col2(p4(1,2),binvm1,n)
         else if (imat.eq.4) then
            call rzero(vz,n)
            call cadv(p4,p3,ifexact)
            call opbinv1(p4,p4(1,2),zero,p4,p4(1,2),zero,1.)
            call incomprn(p4,p4(1,2),vz,pr)
         else if (imat.eq.5) then
            call rzero(vz,n)
            call incomprn(p3,p3(1,2),vz,pr)
            call opcopy(p4,p4(1,2),p4(1,2),p3,p3(1,2),p3(1,2))
         else if (imat.eq.6) then
            call col3(p4,p3,bm1,n)
            call col3(p4(1,2),p3(1,2),bm1,n)
         endif

         if (imat.le.3) call opdssum(p4,p4(1,2),zero)

         do j2=1,2
         do j=1,n
            aloc(j,j2,i,i2) = p4(j,j2)
         enddo
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compact2(a,b,ndof)

      include 'SIZE'
      include 'TOTAL'

      real b(lx1*ly1*lz1*nelt,2,lx1*ly1*lz1*nelt,2)
      real a(ndof,2,ndof,2)

      do i2=1,2
      do i=1,ndof
      do j2=1,2
      do j=1,ndof
         a(j,j2,i,i2) = b(j,j2,i,i2)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compactv(a,b,ndof)

      include 'SIZE'
      include 'TOTAL'

      real b(lx1*ly1*lz1*nelt,ldim)
      real a(ndof,ldim)

      do i2=1,ldim
      do i=1,ndof
         a(i,i2) = b(i,i2)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine matdump0(a,str,nd)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)

      real a(nd,nd)

      character*3 str

      open (unit=25,file=str)

      do j=1,nd
      do i=1,nd
         write (25,*) a(i,j)
      enddo
      enddo

      close (unit=25)

      return
      end
c-----------------------------------------------------------------------
      subroutine vecdump(a,str,nd)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)

      real a(nd)

      character*3 str

      open (unit=25,file=str)

      do i=1,nd
         write (25,*) a(i)
      enddo

      close (unit=25)

      return
      end
c-----------------------------------------------------------------------
      subroutine ivecdump(a,str,nd)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)

      integer a(nd)

      character*3 str

      open (unit=25,file=str)

      do i=1,nd
         write (25,*) a(i)
      enddo

      close (unit=25)

      return
      end
c-----------------------------------------------------------------------
      subroutine matdump(a,str,nd)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)

      real a(l1*nelt,l1*nelt)

      character*3 str

      open (unit=25,file=str)

      do j=1,nd
      do i=1,nd
         write (25,*) a(i,j)
      enddo
      enddo

      close (unit=25)

      return
      end
c-----------------------------------------------------------------------
      subroutine matdump2(a,str,nd)

      include 'SIZE'
      include 'TOTAL'

      parameter (l1=lx1*ly1*lz1)

      real a(l1*nelt,2,l1*nelt,2)

      character*3 str

      open (unit=25,file=str)

      do j2=1,2
      do j=1,nd
      do i2=1,2
      do i=1,nd
         write (25,*) a(i,i2,j,j2)
      enddo
      enddo
      enddo
      enddo

      close (unit=25)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_prec_apply(yp,mmat,pvec,idof,ndof)

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)

      real mmat(ndof,ldim,ndof,ldim)
      real yp(ndof,ldim)

      integer idof(ndof*ldim)
      integer pvec(ndof*ldim)
      character*1 tr

      tr = 'N'

      nla = ndof * ldim

      call dgetrs(tr,nla,1,mmat,nla,pvec,yp,nla,info)

      if (info.ne.0) then
         if (nio.eq.0) write (6,*) 'error in dgetrs in sns_pres_apply',
     $   info
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_prec(a,b) ! a = M b

      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (l1=lx1*ly1*lz1)


c     real mmat((l1*nelt*ldim)**2)
      real mmat((1*ldim)**2)
      integer pvec(l1*lelt*ldim)
      integer idof(l1*lelt)

      real zero(l1*lelt)

      common /snsprec/ mmat
      common /snsiprec/ pvec, idof, ndof

      integer icalld
      save    icalld
      data    icalld /0/

      real a(lt,ldim)
      real b(lt,ldim)

      real ap(lt,ldim)
      real app(lt,ldim)

      icalld = icalld + 1

      n = lx1*ly1*lz1*nelv

      call opzero(a,a(1,2),a(1,ldim))

      if (icalld.eq.1) call sns_prec_init(mmat,pvec,ndof,idof)
c     call sns_prec_init(mmat,pvec,ndof,idof)

      do j=1,ldim
      do i=1,ndof
         ap(i,j) = b(idof(i),j)
      enddo
      enddo

      call compactv(app,ap,ndof)
      call sns_prec_apply(app,mmat,pvec,idof,ndof)

      do j=1,ldim
      do i=1,ndof
         a(idof(i),j) = app(i,j)
      enddo
      enddo

      call dssum(a,lx1,ly1,lz1)
      call dssum(a(1,2),lx1,ly1,lz1)
      if (ldim.eq.3) call dssum(a(1,3),lx1,ly1,lz1)

c     call dssum(a,lx1,ly1,lz1)
c     call dssum(a(1,2),lx1,ly1,lz1)
c     if (ldim.eq.3) call dssum(a(1,3),lx1,ly1,lz1)

c     call rzero(zero,n)
c     call incomprn(a,a(1,2),zero,pr)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_prec_init(mmat,pvec,ndof,idof)

      include 'SIZE'
      include 'PREC' !lk

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (l1=lx1*ly1*lz1)

      real aloc(lk*nelt,2,lk*nelt,2)
      real amat(lk*nelt,2,lk*nelt,2)

      real mmat((lk*nelt*ldim)**2)
      integer pvec(lk*nelt*2)
      integer idof(lk*nelt)

      n = lk*nelt !lx1*ly1*lz1*nelt

      write (6,*) 'inside sns_prec_init'

      call setdof(idof,ndof)
      nla = ndof * 2

      call scnv
c     write (6,*) 'generating test matrix'
      call genmat2(aloc,2,.false.)
c     call matdump2(aloc,'uuu',l1*nelt)

      call assemble(amat,aloc,idof,n,2,ndof)
c     call matdump2(amat,'ttt',ndof)

      call compact2(mmat,amat,ndof)
c     call matdump0(mmat,'sss',ndof*2)

      call dgetrf(nla,nla,mmat,nla,pvec,info)

      if (info.ne.0) then
         if (nio.eq.0) write (6,*) 'error in dgetrf in sns_pres_init',
     $   info
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ccsk(csk,sk,ifexact) ! compute C(u_k) * s_k (+ C(s_k)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(ltd=lxd*lyd*lzd*lelt)

      real csk(lt,ldim), sk(lt,ldim)

      logical ifexact

      common /snsconvs/ c1s(ltd), c2s(ltd), c3s(ltd),
     $                  u1s(ltd), u2s(ltd), u3s(ltd)

      common /snsconv/ c1v(ltd), c2v(ltd), c3v(ltd),
     $                 u1v(ltd), u2v(ltd), u3v(ltd)

      real t2(lt,ldim)

      nel = nelv
      if (ifield.eq.2) nel = nelt

      call convect_new(csk(1,1),u1s,.true.,c1v,c2v,c3v,.true.)
      call convect_new(csk(1,2),u2s,.true.,c1v,c2v,c3v,.true.)
      if (ldim.eq.3)
     $ call convect_new(csk(1,3),u3s,.true.,c1v,c2v,c3v,.true.)

      if (ifexact) then
         call convect_new(t2,u1v,.true.,c1s,c2s,c3s,.true.)
         call convect_new(t2(1,2),u2v,.true.,c1s,c2s,c3s,.true.)
         if (ldim.eq.3)
     $   call convect_new(t2(1,3),u3v,.true.,c1s,c2s,c3s,.true.)

         call opadd2(csk,csk(1,2),csk(1,3),t2,t2(1,2),t2(1,3))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_cht

      ! new JFNK method specifically for incompressible, steady N-S

      ! First compute r_k based on current solution, u_k
      ! Then solve J(u_k) s_k = -F_k = -P B^-1 r_k
      ! Use s_k to update: u_k+1 = u_k + alpha * s_k

      include 'SIZE'
      include 'TOTAL'
      include 'PREC'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)
      parameter(l1=lx1*ly1*lz1)

      real aloc(lk*nelt,2,lk*nelt,2)

      real sk(lt,ldima), rk(lt,ldima), fk(lt,ldima), ud(lt,ldima)

      logical ifpseudo, ifprint, ifbouss

      common /cprint/ ifprint
      common /snsvars/ alpha, tau, taui
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsl2/ el2,fl2,sl2
      common /snsflags/ ifpseudo, ifbouss
      common /snspvars/ up(lt,ldima)

      ifprint = .false.
      ifbouss = .false.

      itn = min(max(itmaxn,10),10000)

c     nel = nelv
c     if (ifield.eq.2) nel = nelt

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      isns=0
      tau = dt
      taui = 1 / tau

      fl2 = tol_newt*1.1

      ! zeroth step

      call opzero(sk,sk(1,2),sk(1,3))
      call rzero(sk(1,3),nt)

      call cfk(rk)
      ifxyo = .true.
c     call outpost(rk,rk(1,2),vz,pr,rk(1,3),'rs1')

      call checke
      call checkf(rk)
      call checks(sk)

      isns = 0
      istep = 0

      call sns_out

      ! end zeroth step

      do isns=1,itn
         istep = isns
         iter_gmres = 0

c        call opchsgn(rk,rk(1,2),rk(1,3))
         call snschsgn(rk)
         call sns_gmres(sk,rk,vmult,tmult,it)
         iter_gmres = iter_gmres + it
         call opmask(sk,sk(1,2),sk(1,3))

         call opadds(vx,vy,vz,sk,sk(1,2),sk(1,3),alpha,nv,2)
         call add2s2(t,sk(1,3),alpha,nt)
         call incomprn(vx,vy,vz,pr)
         call outpost(vx,vy,vz,pr,t,'ttt')

         call cfk(rk)

         call checke
         call checkf(rk)
c        call exitt0
         call checks(sk)

         call sns_out

         if (nio.eq.0) then
            if (fl2.le.atoln) write (6,*) 'reached atoln', fl2
            if (isns.eq.itn)  write (6,*) 'reached max itn', isns
            if (sl2.le.atols) write (6,*) 'reached atols', sl2
         endif

         if (fl2.le.atoln.or.isns.eq.itn.or.sl2.le.atols) then
            if (nio.eq.0) write (6,*) 'done with fluid - newton'
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine checkvec_(sk,dl2)

      include 'SIZE'
      include 'TOTAL'

      parameter(lt=lx1*ly1*lz1*lelt)

      logical ifbouss

      common /snsflags/ ifpseudo, ifbouss

      real sk(lt,ldim+ldimt)

      nt = lx1*ly1*lz1*nelt
      itemp = 1

      dl2 = 0.

      if (ifflow) then
         itemp = ldim + 1
         call opnorm(dl2,sk,sk(1,2),sk(1,ldim),'L2 ')
      endif

      if (ifheat) dl2 = sqrt(glsc3(sk(1,itemp),sk(1,itemp),bm1,nt)/
     $            voltm1+dl2*dl2)

      return
      end
c-----------------------------------------------------------------------
      subroutine sns_newton

      ! new Newton-Krylov method specifically for incompressible, steady N-S

      ! First compute r_k based on current solution, u_k
      ! Then solve J(u_k) s_k = -F_k = -P B^-1 r_k
      ! Use s_k to update: u_k+1 = u_k + alpha * s_k

      include 'SIZE'
      include 'TOTAL'
      include 'PREC'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ldima=ldim+ldimt)
      parameter(l1=lx1*ly1*lz1)

      real aloc(lk*nelt,2,lk*nelt,2)

      real sk(lt,ldima), rk(lt,ldima), fk(lt,ldima), ud(lt,ldima)

      logical ifpseudo, ifprint

      common /cprint/ ifprint
      common /snsvars/ alpha, tau, taui
      common /snstol/ rtoln, atoln, rtolp, atolp, rtolg, atolg, etolg
     $                rtols, atols
      common /snsivars/ isns, iter_gmres, itmaxn, itmaxp, itmaxg
      common /snsl2/ el2,fl2,sl2
      common /snsflags/ ifpseudo
      common /snspvars/ up(lt,ldima)

      ifprint = .false.

      itn = min(max(itmaxn,1),10000)

      itemp = 1
      if (ifflow) itemp = ldim + 1

      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

      isns=0
      tau = dt
      taui = 1 / tau

      fl2 = tol_newt*1.1

      call setup_adf_visc

      call setfq

      ! zeroth step

c     call opzero(sk,sk(1,2),sk(1,3))
      call snszero(sk)

      call cfk(rk)

      call checke
      call checkf(rk)
      call checks(sk)

      isns = 0
      istep = 0

      call sns_sett
      call sns_out

      ! end zeroth step

      do isns=1,itn
         istep = isns
         iter_gmres = 0

c        call opchsgn(rk,rk(1,2),rk(1,3))
         call snschsgn(rk)
         call sns_gmres(sk,rk,vmult,tmult,it)
         iter_gmres = iter_gmres + it
c        call opmask(sk,sk(1,2),sk(1,3))
         call snsmask(sk)
         ifxyo = .true.
         call outpost(sk,sk(1,2),vz,pr,t,'sss')

         if (ifflow) then
            call opadds(vx,vy,vz,sk,sk(1,2),sk(1,3),alpha,nv,2)
            call incomprn(vx,vy,vz,pr)
         endif

c        call q_filter(param(103))

         if (ifheat) call add2s2(t,sk(1,itemp),alpha,nt)

         call cfk(rk)

         call checke
         call checkf(rk)
         call checks(sk)

         call sns_sett
         call sns_out

         if (nio.eq.0) then
            if (fl2.le.atoln) write (6,*) 'reached atoln', fl2
            if (isns.eq.itn)  write (6,*) 'reached max itn', isns
            if (sl2.le.atols) write (6,*) 'reached atols', sl2
         endif

         if (fl2.le.atoln.or.isns.eq.itn.or.sl2.le.atols) then
            if (nio.eq.0) write (6,*) 'done with fluid - newton'
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_adf_visc

      include 'SIZE'
      include 'TOTAL'

      common /adfcoef/ h1(lx1*ly1*lz1*lelt)
      real             h1

      jmesh=imesh
      jfield=ifield

      imesh=2
      ifield=2
      it = 0
      mxit = itmaxg
      mxit = min(itmaxg,lgmres)
      call vprops

      n = lx1*ly1*lz1*nelt
      call copy(h1,vdiff(1,1,1,1,2),n)

      hmn = glmin(h1,n)
      hmx = glmax(h1,n)
      write(6,*) hmn,hmx,'h min max'

      imesh=jmesh
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
