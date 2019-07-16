c     Steady state solvers for Nek5000
c      Methods to steady in Nek5000
c        BDFxEXTy (Nek5000 default time-dependent solver)
c        SS methods (main)
c          Newton
c            NK: Newton-Krylov, currently only support PNPN-2 for Flow
c            JFNK-af: Jacobian-free NK with analytic f
c            JFNK-bf: JFNK w/ BDF f
c          Pseudo-transient
c            psNK
c            psJFNKaf
c            psJFNKbf
c     Math symbol: 
c      target   0 = du/dt = f(u)
c      NT solve 0 = g = (u-uc)/dtau - f(u) at ref uc, ifpstime or 1/dtau= 0
c      GMRES solve J \cdot del u = -g(u)
c     Currently support:
c        - dt fixed
c        - dtNT: SER
c        - fully-couple f and Jacobian
c     User interface:
c        - only have to control params at userchk, see steady_ui
c
c     ref:
c        Yu-Hsiang Lan (seanlan1994@gmail.com)
c        Ping-Hsuan Tsai, Kento Kaneko, Li Lu, Pablo Brubeck
c
c     Dev log:
c        dependencies to other files
c           compute_f
c              scnv (K)
c              cfkf (K)
c              cfks (K)
c           JacobiMatVec_loc
c              cjskf_pr (K)
c              cjsks (K)
c           gmres_newton
c              tt_set_jac_h (P)
c           apply_prec
c              flprec (L)
c              scprec (L)
c-----------------------------------------------------------------------
      subroutine steady_ui
c     In userchk:
c       call ss_default ! todo: move to nek_init
c       ! set steady param at here
c       call steady_ui 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON_loc'
      include 'STDY'

      ! user will do something at here = userchk(0)

      if (istep.eq.0) then

        call ss_set_order
        call ss_set_length

        call copy_nek2newton(bdf_saved(1,1,1))   ! save initial state

        call userchk_nt ! userchk(0)

      else

        ! Steady Solvers
        if (.true.) then

          call copy_newton2nek(bdf_saved(1,1,1)) ! recover initial state
          call ss_drive

          if (nio.eq.0)write(*,*)'Error: nested jfnk, exitt!'
          stop

        endif

        ! BDF
        call userchk_nt

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_default
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'SOLN'
      include 'STDY'

      ifpstime = .false.  ! T: pseudo-time, F: Newton
      ifjfnk   = .true.   ! F: NK, T: JFNK
      ifbdfF   = .true.   ! T: bdf f, F: analytic f

      ifusrodr = .false.
      ifprec   = .false.
      !if (param(51).gt.0) ifprec=.true. ! todo, put later

c     pseudo time
      dtau = 5.*dt
      tol_ps = 1.E-10
c     SER
      SER_op = 0
      SER_max= 1.2
      SER_min= 0.8
      tau_max= 1.E12

c     newton
      alpha_nt = 1.0 ! inexact NT: u = u + alpha * del_u  
      Jeps = 1.E-5   ! Jeps based line, 1.E-6
      Jeps_op = 0    ! 0: const, 1: sqrt( (||u||+1)*eps )/||p||, 2: (average_i(||u_i||)/||v||+1)*sqrt(eps)
      atol_gmres = 1.E-12!1.E-14 !1.E-8
      rtol_gmres = 1.E-10!1.E-16 !1.E-6
      atol_nt = 1.E-12
      rtol_nt = 1.E-10
      maxit_ps= 1
      maxit_nt= 200
      maxit_g = 500
      ifntcf = .false.
      ifntcg = .false.

c     gmres

c     f
      dt_op = 0 ! bdff


c     Debug options
      ifdbg   = .false.
      ifdbgPS = .false.
      ifdbgNT = .false.
      ifdbgF  = .false.
      if (ifdbgPS.or.ifdbgPS.or.ifdbgPS) ifdbg = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_drive
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'STDY'

      call ss_init ! update param, zero flds

      if (ifpstime) then
        call ps_drive
      else 
        call nt_drive
      endif
cc      call sns_drive ! Kento's NK, set in user


      if (nio.eq.0) write(*,*)'End of SS drive'
      call exitt0

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_init
c     reset logic and param after user specifying
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'
      integer ic

c     set up logic and params

      if (.not.ifjfnk) ifbdfF=.false.

      if (.not.ifpstime) dtau = 1.E12 ! not used, but will be shown in print

      if (nfldc.gt.toteq) then ! check toteq is large enough for storage of copying
        write(*,*)'toteq is not large enough!'
        write(*,*)'nfldc (default)  = ',ndim+1+ldimt
        write(*,*)'nfldc (user def) = ',nfldc
        write(*,*)'toteq            = ',toteq
        call exitt
      endif

      if (.not.ifjfnk.AND..not.ifbdff.AND.lx1.eq.lxd) then
        write(*,*)'PNPN is not supported in NK(af) for now'
        call exitt
      endif

c     zeros
      call newton_init ! todo

c     print flags, params... ! todo clean up with uniform format
      if (nio.eq.0) then
        write(*,*)'SS setup: ifpstime',ifpstime 
        write(*,*)'SS setup: ifjfnk  ',ifjfnk
        write(*,*)'SS setup: ifbdff  ',ifbdff
        write(*,*)'SS setup: ifusrodr',ifusrodr
        write(*,*)'SS setup: ifprec  ',ifprec
        write(*,*)'SS setup: npvts   ',(npvts(ic),ic=1,nfldc)
        write(*,*)'SS setup: alpha_nt',alpha_nt
        write(*,*)'SS setup: Jeps_op ',Jeps_op
        write(*,*)'SS setup: maxit_ps',maxit_ps
        write(*,*)'SS setup: maxit_nt',maxit_nt
        write(*,*)'SS setup: maxit_g ',maxit_g
      endif

c      if (ifbdff) then ! ToDo: check if this part has to be done in userdat2
c         param(27) = -1 ! BDF1 only
c         param(93) = 0.0! no projection
c         param(94) = 0.0
c         param(95) = 0.0
c      endif

c     ics
      ! todo
      if (.not.ifbdfF) then
        if (ifflow) then
          call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)
          if (lx1.ne.lxd) call incomprn(vx,vy,vz,pr)
        endif
        if (ifheat) then
          do ifield=2,nfield
            call bcdirsc(t(1,1,1,1,ifield-1))
          enddo
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_set_order
      implicit none
      include 'SIZE' 
      include 'TOTAL' 
      include 'STDY' 

      integer ic
      integer ntmpid,tmp_id(toteq),test_id(toteq*2)

      if (ifusrodr) then ! ToDo: set up a test case using this

        ! frozen ss_ord +  do nothing to cp_ord
        do ic=1,nfld
          ss_ord(ic) = ic
        enddo

        ! check if cp_ord valid
        ! size check
        if (nfldc.gt.toteq*2) then
          if(nio.eq.0) write(*,*)'nfldc > toteq*2',nfldc,toteq
          call exitt
        endif
        ! duplicate check
        do ic=1,toteq*2
          test_id(ic) = 0
        enddo
        do ic=1,nfldc
          if (test_id(cp_ord(ic)).gt.0) then
            if(nio.eq.0) write(*,*)'cp_ord is wrong plz double check'
            call exitt
          endif
          test_id(cp_ord(ic)) = 1
        enddo

      else

        ! Index table
        ! | fld | ss_ind | ifield 
        ! | vx  | 1      | 1
        ! | vy  | 2      | 1
        ! | vz  | ndim   | 1
        ! | pr  | ndim+1 | x
        ! | t(1)| ndim+2 | 2
        ! | t(2)| ndim+3 | 3
        ! Ex: 
        !  2D flow + t: ss_ord = 1 2 4 cp_ord = 1 2 4 3
        !  3D t-only  : ss_ord = 5     cp_ord = 5 1 2 3 4
  
        nfld  = 0
        nfldc = 0
        ntmpid= 0
        if (ifflow) then            ! solve flow, put pr into cp list
          ss_ord(nfld+1) = 1        ! vx
          ss_ord(nfld+2) = 2        ! vy
          ss_ord(nfld+ndim) = ndim  ! vz
  
          nfld = nfld + ndim
          nfldc= nfldc+ ndim + 1
  
          ntmpid = ntmpid + 1
          tmp_id(ntmpid) = ndim+1   ! pr
        else
          tmp_id(ntmpid+1) = 1
          tmp_id(ntmpid+2) = 2
          tmp_id(ntmpid+ndim)   = ndim
          tmp_id(ntmpid+ndim+1) = ndim + 1
  
          nfldc= nfldc+ ndim + 1
          ntmpid = ntmpid + ntmpid+ndim+1
        endif
        if (ifheat) then            ! solve scalar or put into cp list
          nfldc = nfldc + nfield-1
          do ifield = 2,nfield
            if (idpss(ifield-1).eq.0) then      
              ss_ord(nfld+1) = ifield+ndim
              nfld = nfld + 1
            else
              ntmpid = ntmpid + 1
              tmp_id(ntmpid) = ifield+ndim
            endif
          enddo
        endif
        ! build cp list
        do ic=1,nfld              ! cp solved fld  first
          cp_ord(ic) = ss_ord(ic)
        enddo
        do ic=1,ntmpid            ! cp all
          cp_ord(nfld+ic) = tmp_id(ic)
        enddo

      endif

      if (nio.eq.0) then
        write(*,*) 'SSorder, ifusrodr = ',ifusrodr
        write(*,*) 'SSorder, nfld/nfldc/toteq',nfld,nfldc,toteq
        write(*,*) 'SSorder, ss_ord:',(ss_ord(ic),ic=1,nfld)
        write(*,*) 'SSorder, cp_ord:',(cp_ord(ic),ic=1,nfldc)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ps_drive
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON_loc'
      include 'STDY'

      common /jfnk_tmp/ BDFcnt
      integer BDFcnt

      real fnorms(toteq),fnorm
      real gnorms(toteq),gnorm,ginorm,res_nt,ratio
      real res_gmres

      integer ic,i

      real c_in(lvt1,toteq),c_out(lvt1,toteq),timeNT
      integer iter_GMRES,tot_iterg
      real max_GMRESiter
      ! todo
      if (nid.eq.0) write(*,*) 'Hello jfnk.f'
      if (nid.eq.0) write(*,*) 'ps_drive: Pseudo-transient algorithm'

      istepPS = 0
      istepNT = 0
      tot_iterg = 0

      istep = 0
      timeNT = time

      call copy_nek2newton(ck_0) ! from nek to newton

      istep = 1 ! force Nek treating BDF w/o past info.

      call SS_PSchk
      call SS_output(ck_0,timeNT)      ! initial
c      call copy(cn_k,ck_0,lvt1*toteq)  ! including background fields
      call SSop_copy(cn_k,ck_0,nfld,ss_ord,npvts)  ! including background fields

      ! comp f0
      call compute_fk_nt(fn,ck_0,dt,timeNT)

      ! comp f0 g0 norms
      call ss_comp_norm_l2(fnorms,fnorm,fn,'f','nt',istepNT)
      f_now = fnorm

      ifntcf = .false. ! compute f at NT(0)
      ifntcg = .false. ! compute g at NT(0)
c     Newton Solver
      do istepPS=1,nsteps

        ! comp g0
        call compute_gn_nt(gn_n,fn,cn_k,ck_0,dtNT)
        call SSop_chsign(gn_n,nfld,ss_ord,npvts)
        call ss_comp_norm_l2(gnorms,ginorm,gn_n,'g','nt',istepNT)
        if(nid.eq.0) write(6,90)0,0,1.0,ginorm,ginorm,dtNT,dt

          call nt_solver(gn_n,cn_k,ifbdff,fn,ck_0,ifntcf,ifntcg
     $                  ,istepPS,timeNT,dtNT,mults,nfld,ss_ord,npvts
     $                  ,alpha_nt,maxit_nt,atol_nt,rtol_nt
     $                  ,maxit_g,atol_gmres,rtol_gmres
     $                  ,tot_iterg,gnorm,fnorm)

          call ss_comp_norm_l2(fnorms,fnorm,fn,'f','nt',istepNT)
          f_pre = f_now
          f_now = fnorm

          call SS_PSchk
          call SS_output(cn_k,timeNT)
c          call SS_copy(ck_0,cn_k,lvt1*toteq)
          call SSop_copy(ck_0,cn_k,nfld,ss_ord,npvts)

c         call compute_dt(dt)
          call get_dtNT_loc(dtNT,max_dtNT,f_pre,f_now,min_ftol
     $                     ,max_GMRESiter)

      enddo
 90   format('newton iter',2i6,1p5e12.4)

      call SS_output(cn_k,timeNT)

      return
      end
c-----------------------------------------------------------------------
      subroutine nt_drive
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON_loc'
      include 'STDY'

      common /jfnk_tmp/ BDFcnt ! ToDo merge to somewhere
      integer BDFcnt

      real fnorms(toteq),fnorm
      real gnorms(toteq),gnorm,ginorm,res_nt,ratio
      real res_gmres

      real c_in(lvt1,toteq),c_out(lvt1,toteq),timeNT
      integer iter_GMRES,tot_iterg

      if (nid.eq.0) write(*,*) 'Hello jfnk.f'
      if (nid.eq.0) write(*,*) 'nt_drive: Newton-Krylov based methods'
      

      istepPS = 0
      istepNT = 0

      tot_iterg = 0

      istep = 0
      timeNT = time

      call copy_nek2newton(ck_0) ! from nek to newton

      istep = 1 ! force Nek treating BDF w/o past info.

      call SS_NTchk(ck_0,istepNT)
      call SS_output(ck_0,timeNT)      ! initial
      call copy(cn_k,ck_0,lvt1*toteq)  ! including all background fields


      istepPS = 1

c     Newton Solver
      ifntcf = .true. ! compute f at NT(0)
      ifntcg = .true. ! compute g at NT(0)

      call nt_solver(gn_n,cn_k,ifbdff,fn,ck_0,ifntcf,ifntcg
     $              ,istepPS,timeNT,dtNT,mults,nfld,ss_ord,npvts
     $              ,alpha_nt,maxit_nt,atol_nt,rtol_nt
     $              ,maxit_g,atol_gmres,rtol_gmres
     $              ,tot_iterg,gnorm,fnorm)

      ! ToDo: add stat
 90   format('newton iter',2i6,1p5e12.4)

      call SS_output(cn_k,timeNT)

      return
      end
c-----------------------------------------------------------------------
      subroutine nt_solver(gn_n,cn_k,ifbdff,fn,ck_0,ifntcf,ifntcg
     $              ,istepPS,timeNT,dtNT,mults,nfld,ord,npvts
     $              ,alpha_nt,maxit_nt,atolNT,rtolNT
     $              ,maxit_g,atol_gmres,rtol_gmres
     $              ,tot_iterg,gnorm,fnorm)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      logical ifbdff,ifntcf,ifntcg
      real gn_n(lvt1,toteq)  ! up-to-date
     $    ,cn_k(lvt1,toteq)  ! main in/output
     $    ,fn  (lvt1,toteq)  ! up-to-date
     $    ,ck_0(lvt1,toteq)  ! fixed
     $    ,sn_k(lvt1,toteq)  ! work array
     $    ,mults(lvt1,toteq) ! input
      real atol_gmres,rtol_gmres,res_gmres
      real atolNT,rtolNT,alpha_nt
      integer nfld,npvts(toteq),ord(toteq)
     $       ,istepPS,istepNT,maxit_nt,maxit_g
     $       ,iter_nt,tot_iterg,iter_gmres
     $       ,iocomm_nt

      real gnorms(toteq),gnorm,ginorm,ratio
      real fnorms(toteq),fnorm
      real snorms(toteq),snorm
      real timeNT,dtNT

      istepNT = 0
c      iocomm_nt = 5 ! ToDo, merge the iocomm
c      if (ifdbgnt) iocomm_nt = 1
      iocomm_nt = 1

      ! comp f0
      if (ifntcf) call compute_fk_nt(fn,ck_0,dt,timeNT)

      ! comp g0
      if (ifntcg) then
      call compute_gn_nt(gn_n,fn,cn_k,ck_0,dtNT)
      call SSop_chsign(gn_n,nfld,ord,npvts)
      endif

      ! comp f0 g0 norms
      call ss_comp_norm_l2(fnorms,fnorm,fn,'f','nt',istepNT)
      call ss_comp_norm_l2(gnorms,ginorm,gn_n,'g','nt',istepNT)
      if(nid.eq.0) write(6,90)istepPS,0,1.0,ginorm,ginorm,dtNT,dt


      ! Main loop
      do istepNT=1,maxit_nt

        call gmres_newton
     $         (sn_k,gn_n,cn_k,fn,timeNT,mults,nfld,ord,npvts
     $         ,maxit_g,atol_gmres,rtol_gmres,iter_gmres,res_gmres
     $         ,istepPS,istepNT)
        tot_iterg = tot_iterg + iter_gmres

        call SSop_add2s2(cn_k,sn_k,alpha_nt,nfld,ord,npvts)! cn_k = cn_k + alpha*sn_k

        if (.not.ifbdfF) call NT_dir_incomp(cn_k)    ! pr/BC correction

        call SS_NTchk(cn_k,istepNT) ! chk per NT iter 

        ! update f g and norms
        call compute_fk_nt(fn,cn_k,dt,timeNT)
        call compute_gn_nt(gn_n,fn,cn_k,ck_0,dtNT)
        call SSop_chsign(gn_n,nfld,ord,npvts)
        call ss_comp_norm_l2(snorms,snorm,sn_k,'s','nt',istepNT)
        call ss_comp_norm_l2(gnorms,gnorm,gn_n,'g','nt',istepNT)

        ratio=gnorm/ginorm

        if (nio.eq.0.AND.mod(istepNT,iocomm_nt).eq.0) 
     $     write(6,90) istepPS,istepNT,ratio,gnorm,ginorm,dtNT,dt

        ! chk conv.
        if (ratio.lt.rtolNT) goto 900
        if (gnorm.lt.atolNT) goto 900

        timeNT = timeNT + real(dtNT)/maxit_nt
      enddo
 90   format('newton iter',2i6,1p5e12.4)
900   continue

      iter_nt= istepNT

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_comp_norm_l2(norm_i,norm_a,c_in,s1,s2,iter)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'

      real glsc3
      real c_in(lvt1,toteq),norm_i(toteq),norm_a
      integer ic,iter,ii
      character*1 s1
      character*2 s2

      norm_a = 0.0
      do ii=1,nfld
        ic = ss_ord(ii)
        norm_i(ic) = glsc3(c_in(1,ic),c_in(1,ic),mults(1,ic),npvts(ic))
        norm_i(ic) = sqrt(norm_i(ic))
        norm_a = max(norm_a,norm_i(ic))
      enddo

      if(nio.eq.0) write(6,95)s1,s2,iter,
     $                        (ss_ord(ic),norm_i(ss_ord(ic)),ic=1,nfld)
      if(nio.eq.0) write(6,96)s1,s2,iter,norm_a
 95   format('jfnk ',A,'_now ',A,i6,7(i4,1p1e12.4))
 96   format('jfnk ',A,'_now ',A,i6,' all',1p1e12.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_comp_norm_linf(norm_i,norm_a,c_in)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'

      real glamax
      real c_in(lvt1,toteq),norm_i(toteq),norm_a
      integer ic,ii

      norm_a = 0.0
      do ii=1,nfld
        ic = ss_ord(ii)
        norm_i(ic) = glamax(c_in(1,ic),npvts(ic))
        norm_a = max(norm_a,norm_i(ic))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine NT_dir_incomp(c_in)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'
      real c_in(lvt1,toteq)

      if (ifbdff) return ! fixme: outer flags

      call copy_newton2nek(c_in)
      if (ifflow) then
        call bcdirvc(vx,vy,vz,v1mask,v2mask,v3mask)
        if (lx1.eq.lxd) then
          write(*,*)'PNPN is not supported for analytic F'
        else 
          call incomprn(vx,vy,vz,pr)
        endif
      endif
      if (ifheat) then
        do ifield=2,nfield
          call bcdirsc(t(1,1,1,1,ifield-1))
        enddo
      endif

      call copy_nek2newton(c_in)

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_nek2newton(c_out)! ToDo: cp_fld copy
c     user interface
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY' ! nfldc, cp_ord

      integer nnv,nnt,nnpr,ic,ii
      real c_out(lvt1,toteq) !FIXME max lvt1, lvt2, Lan

      nnv = nx1*ny1*nz1*nelv
      nnt = nx1*ny1*nz1*nelt
      nnpr= nx2*ny2*nz2*nelv

      if (ifusrodr) then                ! user specified order

        do ii=1,nfldc
          ic = cp_ord(ii)
          if    (ic.eq.1) then
                                          call copy(c_out(1,ii),vx,nnv)
          elseif(ic.eq.2) then
                                          call copy(c_out(1,ii),vy,nnv)
          elseif(ic.eq.ndim.AND.ndim.eq.3) then
                                          call copy(c_out(1,ii),vz,nnv)
          elseif(ic.eq.ndim+1) then
                                          call copy(c_out(1,ii),pr,nnpr)
          else
                        call copy(c_out(1,ii),t(1,1,1,1,ic-ndim-1),nnt)
          endif
        enddo

      else                              ! deefault order

        call copy(c_out(1,1),vx,nnv)
        call copy(c_out(1,2),vy,nnv)
        if (ndim.eq.3) call copy(c_out(1,ndim),vz,nnv)
  
        call copy(c_out(1,ndim+1),pr,nnpr)
  
        do ic=1,ldimt
        call copy(c_out(1,ndim+1+ic),t(1,1,1,1,ic),nnt)
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_newton2nek(c_in)
c     user interface
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY' ! nfldc, cp_ord

      integer nnv,nnt,nnpr,ic,ii
      real c_in(lvt1,toteq)

      nnv = nx1*ny1*nz1*nelv
      nnt = nx1*ny1*nz1*nelt
      nnpr= nx2*ny2*nz2*nelv

      if (ifusrodr) then                ! user specified order

        do ii=1,nfldc
          ic = cp_ord(ii)
          if    (ic.eq.1) then
                                          call copy(vx,c_in(1,ii),nnv)
          elseif(ic.eq.2) then
                                          call copy(vy,c_in(1,ii),nnv)
          elseif(ic.eq.ndim.AND.ndim.eq.3) then
                                          call copy(vz,c_in(1,ii),nnv)
          elseif(ic.eq.ndim+1) then
                                          call copy(pr,c_in(1,ii),nnpr)
          else
                        call copy(t(1,1,1,1,ic-ndim-1),c_in(1,ii),nnt)
          endif
        enddo

      else                              ! default order

        call copy(vx,c_in(1,1),nnv)
        call copy(vy,c_in(1,2),nnv)
        if (ndim.eq.3) call copy(vz,c_in(1,ndim),nnv)
  
        call copy(pr,c_in(1,ndim+1),nnpr)
        do ic=1,ldimt
        call copy(t(1,1,1,1,ic),c_in(1,ndim+1+ic),nnt)
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ss_set_length !Todo, change name and merged
c     init. local variables
c-----------------------------------------------------------------------
      implicit none 
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'

      integer ic,ii,nnv,nnpr,nnt

      nnv = nx1*ny1*nz1*nelv
      nnt = nx1*ny1*nz1*nelt
      nnpr= nx2*ny2*nz2*nelv


      if (ifusrodr) then    ! usr-provided index

        do ii=1,nfldc
          ic = cp_ord(ii)
          if (ic.le.ndim) then 
            npvts(ii) = nnv
            call copy(mults(1,ii),vmult,nnv)
          elseif (ic.eq.ndim+1) then
            npvts(ii) = nnpr
            call rone(mults(1,ii),nnpr)
          else
            npvts(ii) = nnt
            call copy(mults(1,ii),tmult,nnt)
          endif
        enddo

      else                  ! default index

        npvts(1) = nnv
        npvts(2) = nnv
        npvts(ndim) = nnv
        npvts(ndim+1) = nnpr
        do ic=1,ldimt
        npvts(ndim+1+ic) = nnt
        enddo
  
        call copy(mults(1,1),vmult,nnv)
        call copy(mults(1,2),vmult,nnv)
        call copy(mults(1,ndim),vmult,nnv)
        call rone(mults(1,ndim+1),nnpr)
        do ic=1,ldimt
        call copy(mults(1,ndim+1+ic),tmult,nnt)
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine newton_init
c     init. global variables
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON_loc'
      include 'STDY'

      common /jfnk_tmp/ BDFcnt
      integer BDFcnt

      integer ic

c     Initialize working arrays
      do ic=1,toteq
        call rzero(fn  (1,ic),lvt1)
        call rzero(gn_n(1,ic),lvt1)
        call rzero(ck_0(1,ic),lvt1)
        call rzero(cn_k(1,ic),lvt1)
        call rzero(sn_k(1,ic),lvt1)
      enddo

      BDFcnt = 0 ! number that BDF1 is called

c     TTTTTT
c      maxit_nt  = 200           ! iteration for Newton method
      dtNT      = dtau           ! pseudo-transient time step ! ToDo merge
c      dtau      = dt           ! pseudo-transient time step
      tolNT     = 1e-12        ! tolerance for Newton method
      max_dtNT  = 1E12         ! Max. of dtNT, avoid NaN
      min_ftol  = 1E-12        ! Min. of f norm, avoid NaN

      if (nid.eq.0) then
        write(*,*)'jfnk_param, maxit_nt',maxit_nt
        write(*,*)'jfnk_param, dtNT    ',dtNT
        write(*,*)'jfnk_param, tolNT   ',tolNT
        write(*,*)'jfnk_param, max_dtNT',max_dtNT
        write(*,*)'jfnk_param, min_ftol',min_ftol
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_gn_NT(g_out,f_in,ck,c0,dtNT)
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'STDY'

      integer ic,ii
      real g_out(lx1*ly1*lz1*lelv,toteq)
     $   , f_in (lx1*ly1*lz1*lelv,toteq)
     $   , ck   (lx1*ly1*lz1*lelv,toteq)
     $   , c0   (lx1*ly1*lz1*lelv,toteq)
      real dtNT,beta

      beta=0.0 ! reduced to g = -f
      if (ifpstime) beta=1.0/dtNT

      do ii=1,nfld
        ic = ss_ord(ii)
        call sub3 (g_out(1,ic),ck(1,ic),c0(1,ic),npvts(ic))
c        call cmult(g_out(1,ic),beta,npvts(ic)) 
        call col2c(g_out(1,ic),bm1,beta,npvts(ic))   ! g_out = g_out*bm1*beta
        call sub2 (g_out(1,ic),f_in(1,ic),npvts(ic)) ! g = bm1*(ck-c0)/dtNT - f
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine compute_fk_NT(f_out,c_in,dt_loc,timeNT)
c     couple
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'

      integer ic,ii
      real f_out(lx1*ly1*lz1*lelv,toteq)
     $   , c_in (lx1*ly1*lz1*lelv,toteq)
     $   , c_out(lx1*ly1*lz1*lelv,toteq)
      real dt_loc ! local usage
      real timeNT
      logical ifh_bak,iff_bak

      real rk(lvt1,toteq)

      dt=dt_loc

      time = timeNT
      istep = 1

      if (ifbdff) then
        ! Option 1: bdf F approx [f(u+ep)-f(u)]/eps
        call bdf1_in_newton(c_in,c_out,istepNT,timeNT,nfld,npvts)
  
        do ii=1,nfld
          ic = ss_ord(ii)
          call sub3(f_out(1,ic),c_out(1,ic),c_in(1,ic),npvts(ic))
          call cmult(f_out(1,ic),1.0/dt,npvts(ic))
c          call col2(f_out(1,ic),bm1,npvts(ic))
        enddo

      else
        ! Option 2: analytic f (steady.f)
        call copy_newton2nek(c_in)

        ! ToDo: add interface to support heat and coupled eqns
        call tt_scnv ! conv term

        do ii=1,nfld
          ic = ss_ord(ii)
          if (ic.eq.1) then ! flow
            ifield = 1
            call cfkf(rk(1,ic))
          elseif (ic.gt.ndim+1) then
            ifield = ic-ndim
c            call cfks(rk(1,ic))
            call tt_crks(rk(1,ic))
c            call comp_Fs(rk(1,ic),ifield)
          endif
          call copy (f_out(1,ic),rk(1,ic),npvts(ic))
        enddo
      endif

      return
      end
c----------------------------------------------------------------------- 
      subroutine JacobiMatVec_loc(Jacp,p,uc,f_k0,dt_loc,timeNT)
c     fully coupled
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP' ! ifield
      include 'MASS' ! bm1
      include 'STDY'

      real Jacp (lvt1,toteq)
     $   , p    (lvt1,toteq)
     $   , uc   (lvt1,toteq) ! Jacobian Mat at uc
     $   , f_k0 (lvt1,toteq) 
     $   , uep  (lvt1,toteq) 
     $   , f_out(lvt1,toteq) 
     $   , foi  (lvt1,toteq) ! f_out - f_in
      real dt_loc,timeNT  ! local usage
      real eps,epsinv,beta
      integer ic,ii

      beta = 0.0
      if (ifpstime) beta = 1.0/dtau
      

      if (ifjfnk) then

        call SS_get_Jeps(eps,p,uc)
        epsinv = 1./eps

        call copy(uep,uc,lvt1*toteq) ! copy all to preserve background fields
        do ii=1,nfld
          ic = ss_ord(ii)
          call add3s2(uep(1,ic),uc(1,ic),p(1,ic),1.0,eps,npvts(ic)) ! uep = u + eps*p
        enddo

        call compute_fk_NT(f_out,uep,dt_loc,timeNT)

        do ii=1,nfld
          ic = ss_ord(ii)
          call sub3   (foi(1,ic),f_out(1,ic),f_k0(1,ic),npvts(ic))
        enddo
        do ii=1,nfld
          ic = ss_ord(ii)
          call cmult  (foi(1,ic),epsinv,npvts(ic)) ! foi = (fo-fi)/eps
          call add3s2 (Jacp(1,ic),p(1,ic),foi(1,ic),beta,-1.0,npvts(ic))
                                                   ! Jp = p/dt_newton - (fo-f)/eps
c         call cmult  (Jacp(1,ic),foi(1,ic),-epsinv,npvts(ic)) ! -(fo-fi)/eps
c         call admcol3(Jacp(1,ic),p(1,ic),bm1,beta,npvts(ic))
                                                   ! Jp = - (fo-f)/eps + bm1*p/dt_newton
c         call cmult  (Jacp(1,ic),foi(1,ic),-epsinv,npvts(ic)) ! -(fo-fi)/eps
c         call admcol3(Jacp(1,ic),p(1,ic),tmask,beta,npvts(ic))

c          call col2(Jacp(1,ic),bm1,npvts(ic))
        enddo

      else  ! analytic J

        do ii=1,nfld
          ic = ss_ord(ii)
          if (ic.eq.1) then         ! flow
            ifield = 1
            call cjskf_pr(Jacp(1,ic),p(1,ic)) ! do vx,vy,vz at once
          elseif(ic.gt.ndim+1) then ! heat
            ifield = ic-ndim
            call tt_cjsks(Jacp(1,ic),p(1,ic))
          endif
c          call add2sxy(Jacp(1,ic),-1.0,p,beta,npvts(ic)) ! Jac_g*p = beta*p - Jac_f*P
          call admcol3(Jacp(1,ic),p,bm1,-beta,npvts(ic)) ! Jac_g*p = beta*B*p - Jac_f*P
          call chsign(Jacp(1,ic),npvts(ic))
        enddo

      endif

      return
      end
c----------------------------------------------------------------------- 
      subroutine gmres_newton
     $           (phi_fc,res_fc,uc,f_fc,timeNT,mults,nlfd,ord,npvts
     $           ,maxit_g,atolg,rtolg,iter_out,res_out
     $           ,istepPS,istepNT)
c     fully coupled
c     copy from NekCEM (solving only Ax=b)
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL' !dt is in here, ToDo fix it, Lan
      include 'NMGMRES'

      common /snsflags/ ifpseudo, ifbouss, ifprec ! precond, Kento's code
      logical ifpseudo,ifbouss,ifprec
      logical ifrefc

      integer  npvts(toteq),nlfd,ord(toteq),outer,maxit_g,iter_out
      real     phi_fc(lvt1,toteq),res_fc(lvt1,toteq),mults(lvt1,toteq)
      real     res_out,forcing,tol,alpha,l,temp
      real     eps,uc(lvt1,toteq),f_fc(lvt1,toteq)
      real*8   etime1,dnekclock

      integer  i,j,k,ii,ic,iconv,iter,m,itmax
      real     rinorm,rnorm,tolpss,ratio,rtolg,atolg
      real     glsc3,glamax
      real     timeNT
      integer  istepPS,istepNT

      real     sc(toteq),sr(toteq),tmpx(lvt1,toteq)

      real tmp_r(10),glsc2! Lan-debug

      ifprec = .false.
      if (param(51).gt.0) then
        ifprec=.true. ! follow Kento's setup ! FIXME
        call copy_nek2newton(uc)
c        call tt_set_jac_h ! old
        call tmp_tt_set_jac_h !ToDo: rename later -> ss_set_prec_h ! new
      endif
      ifrefc=.true.
      if (glamax(res_fc,npvts(ord(1))).lt.1E-10) ifrefc = .false.


      iter  = 0
      m     = min(lgmres_NT,maxit_g)

c     FFFF
      itmax = min(lgmres_NT,maxit_g) ! FIXME: restart is broken, Lan
      tolps = atolg
      tolpss= tolps
      iconv = 0

      call rzero(x_NT,lvt1*toteq)
      call rzero(h_NT,m*m)

      call SS_get_sc_norms(uc,ord,sc)
      do ii=1,nlfd
        ic = ord(ii)
        sr(ic) = 1./sc(ic)
      enddo
      if (nid.eq.0) write(*,*)'Lan_>> newtongmres_sc='
     $             ,(sc(ord(ii)),ii=1,nlfd)

      do ii=1,nlfd
        ic = ord(ii)
        call cmult(res_fc(1,ic),sr(ic),npvts(ic))
      enddo

      outer = 0
      do while (iconv.eq.0.and.iter.lt.itmax)
         outer = outer+1
         if(iter.eq.0) then
            do ii=1,nlfd
            ic = ord(ii)
            call copy  (r_NT(1,ic),res_fc(1,ic),npvts(ic))           ! r = res
            call chsign(r_NT(1,ic),npvts(ic))
            enddo
         else
            !update residual
            do ii=1,nlfd
            ic = ord(ii)
            call copy  (r_NT(1,ic),res_fc(1,ic),npvts(ic))          ! r = res
            call chsign(r_NT(1,ic),npvts(ic))
            enddo
            do ii=1,nlfd
              ic = ord(ii)
              call copy (tmpx(1,ic),x_NT(1,ic),npvts(ic))
              call cmult(tmpx(1,ic),sc(ic),npvts(ic))  ! w =     A S_c x
            enddo
            call JacobiMatVec_loc(w_NT,tmpx,uc,f_fc,dt,timeNT)! w =     A x
            do ii=1,nlfd
              ic = ord(ii)
              call cmult(w_NT(1,ic),sr(ic),npvts(ic))  ! w = S_r A x
            enddo
ccc            call incomprn(w_NT,w_NT(1,2),w_NT(1,3),pr)

            do ii=1,nlfd
            ic = ord(ii)
            call add2s2 (r_NT(1,ic),w_NT(1,ic),-1.,npvts(ic)) ! r = r - w
            enddo
         endif

         gamma_NT(1) = 0.0
         do ii=1,nlfd
         ic = ord(ii)
         gamma_NT(1) = gamma_NT(1)
     $               + glsc3(r_NT(1,ic),r_NT(1,ic)
     $                      ,mults(1,ic),npvts(ic))
! gamma  = (r,r)
         enddo
         gamma_NT(1) = sqrt(gamma_NT(1))                 ! gamma  = sqrt{ (r,r) }

         !check for lucky convergence
         rnorm = 0.
         if(gamma_NT(1) .eq. 0.) goto 9000
         temp = 1./gamma_NT(1)
         do ii=1,nlfd
         ic = ord(ii)
         call cmult2(v_NT(1,ic,1),r_NT(1,ic),temp,npvts(ic))       !  v  = r / gamma
         enddo
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            if (j.eq.1.AND.iter.gt.1) then !ToDo: need more test
c               call ss_add2(uc,x_nt,nlfd,npvts) ! tmp_uc = uc + x_nt
               call NT_dir_incomp(uc) ! force correct dir + incomp
c               call copy_newton2nek(tmp_uc)
c               call opadd2(vx,vy,vz,x_nt,x_nt(1,2),x_nt(1,3))
c               call incomprn(vx,vy,vz,pr)
c               call copy_nek2newton(uc)
            endif
            call apply_prec(z_NT(1,1,j),v_NT(1,1,j),uc,ifrefc) ! z = P^{-1} v,  P = Id or P ~ A
            ifrefc=.false.
c            if (.not.ifprec) then
c              do ii=1,nlfd
c              ic = ord(ii)
c              call copy(z_NT(1,ic,j),v_NT(1,ic,j),npvts(ic))
c              enddo
c            endif
cc            if (j.eq.1) then
cc              if (ifprec) then
cc                if (iter.gt.1) then
cc                  ! do nothing
cc                  call copy_newton2nek(uc)
cc                  call opadd2(vx,vy,vz,x_nt,x_nt(1,2),x_nt(1,3))
cc                  call incomprn(vx,vy,vz,pr)
cc                endif
cc                ! copy from fixed to vx
cc                call copy_newton2nek(uc)
ccc                call flprec(z_NT(1,1,j),v_NT(1,1,j),.true.)
cc              endif
cc            else
cc              if (ifprec) then
cc                call copy_newton2nek(uc)
ccc                call flprec(z_NT(1,1,j),v_NT(1,1,j),.false.)
cc              endif
cc            endif


cc         !write(6,*) 'start form m-th krylov subspace'
cc         do j=1,m
cc            iter = iter+1
cc            if (.not.ifprec) then
cc              do ii=1,nlfd
cc              ic = ord(ii)
cc              call copy(z_NT(1,ic,j),v_NT(1,ic,j),npvts(ic))
cc              enddo
cc            endif
cc            if (j.eq.1) then
cc              if (ifprec) then
cc                if (iter.gt.1) then
cc                  ! do nothing
cc                  call copy_newton2nek(uc)
cc                  call opadd2(vx,vy,vz,x_nt,x_nt(1,2),x_nt(1,3))
cc                  call incomprn(vx,vy,vz,pr)
cc                endif
cc                ! copy from fixed to vx
cc                call copy_newton2nek(uc)
ccc                call flprec(z_NT(1,1,j),v_NT(1,1,j),.true.)
cc              endif
cc            else
cc              if (ifprec) then
cc                call copy_newton2nek(uc)
ccc                call flprec(z_NT(1,1,j),v_NT(1,1,j),.false.)
cc              endif
cc            endif

            do ii=1,nlfd
              ic = ord(ii)
              call copy (tmpx(1,ic),z_NT(1,ic,j),npvts(ic))
              call cmult(tmpx(1,ic),sc(ic),npvts(ic))    ! w =     A S_c v
            enddo
            call JacobiMatVec_loc(w_NT,tmpx,uc,f_fc,dt,timeNT)! w = A v
            do ii=1,nlfd
              ic = ord(ii)
              call cmult(w_NT(1,ic),sr(ic),npvts(ic))    ! w = S_r A v
            enddo

            !modified Gram-Schmidt
            do i=1,j
               h_NT(i,j)=0.0
               do ii=1,nlfd
               ic = ord(ii)
               h_NT(i,j)=h_NT(i,j)
     $                  + glsc3(w_NT(1,ic),v_NT(1,ic,i)
     $                         ,mults(1,ic),npvts(ic))        ! h = (w,v )
               enddo                                   ! i,j       i
               do ii=1,nlfd
               ic = ord(ii)
               call add2s2(w_NT(1,ic),v_NT(1,ic,i),-h_NT(i,j)
     $                    ,npvts(ic))! w = w - h    v
               enddo
            enddo                                 !         i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_NT(i,j)
               h_NT(i  ,j)=  c_NT(i)*temp + s_NT(i)*h_NT(i+1,j)
               h_NT(i+1,j)= -s_NT(i)*temp + c_NT(i)*h_NT(i+1,j)
            enddo
            alpha = 0.0
            do ii=1,nlfd                               ! ______
            ic = ord(ii)
            alpha = alpha + glsc3(w_NT(1,ic),w_NT(1,ic)
     $                           ,mults(1,ic),npvts(ic))     ! alpha =  \/ (w,w)
            enddo
            alpha = sqrt(alpha)

            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h_NT(j,j)*h_NT(j,j)+alpha*alpha)
            temp = 1./l
            c_NT(j) = h_NT(j,j) * temp
            s_NT(j) = alpha  * temp
            h_NT(j,j) = l
            gamma_NT(j+1) = -s_NT(j) * gamma_NT(j)
            gamma_NT(j)   =  c_NT(j) * gamma_NT(j)

            rnorm = abs(gamma_NT(j+1))
            if (iter.eq.1) rinorm=rnorm
c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
c   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')
            ratio = rnorm/rinorm
c     OOOOO
            if (rnorm .lt. atolg) goto 900 !converged by absolute tol
            if (ratio .lt. rtolg) goto 900 !converged by relative tol

            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            do ii=1,nlfd
            ic = ord(ii)
            call cmult2(v_NT(1,ic,j+1),w_NT(1,ic),temp,npvts(ic))   ! v    = w / alpha
            enddo                                 !  j+1
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
         res_out = rnorm   ! output |Ax-b|
         if (nid.eq.0) write(6,*) istepPS,istepNT,
     $                 'forcing_term',rnorm/tolpss
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
            do ii=1,nlfd
            ic = ord(ii)
            call add2s2(x_NT(1,ic),z_NT(1,ic,i),c_NT(i),npvts(ic))     ! x = x +
c  z
            enddo
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      do ii=1,nlfd
      ic = ord(ii)
      call copy(phi_fc(1,ic),x_NT(1,ic),npvts(ic))
      enddo
      do ii=1,nlfd
      ic = ord(ii)
      call cmult(phi_fc(1,ic),sc(ic),npvts(ic))
      call chsign(phi_fc(1,ic),npvts(ic))
      enddo

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if (nid.eq.0) then
          write(6,9999) istepPS,istepNT,iter,rtolg,atolg
      endif

 9999 format(' ',' ',i9,i6,'  1gmres_newton_iteration#',i6,1p2e12.4)
      iter_out = iter

      return
      end
c----------------------------------------------------------------------- 
      subroutine apply_prec(z,v,uc,ifrefc)
      ! z = P^{-1} v, where P can be Id or Prec
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'
      real z (lx1*ly1*lz1*lelv,toteq)
     $   , v (lx1*ly1*lz1*lelv,toteq)
     $   , uc(lvt1,toteq)
      logical ifrefc
      integer ic,ii

      if (.not.ifprec) then
        call SSop_copy(z,v,nfld,ss_ord,npvts)
      else !ToDo add couple prec, provide several couple modes
        call copy_newton2nek(uc)
        do ii=1,nfld
        ic = ss_ord(ii)
        if (ic.eq.1) then         ! flow
          ifield = 1
          call flprec(z(1,ic),v(1,ic),ifrefc) ! ToDo, fix ui, update visc myself
        elseif(ic.gt.ndim+1) then ! heat
          ifield = 2
          call scprec(z(1,ic),v(1,ic),ifrefc) ! ToDo, fix ui, update visc myself
        endif
        enddo
      endif

      return
      end
c----------------------------------------------------------------------- 
      subroutine get_dtNT_loc(dtNT,max_dtNT,f_pre,f_now,min_ftol
     $                       ,iterGMRES)
c     SER
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'
      real dtNT,max_dtNT,f_pre,f_now,min_ftol
      real ratio,iterGMRES

c      real tol_ps,SER_min,SER_max,tau_max
c      integer SER_op
      ! ToDo global or all local
      SER_op = 0
      SER_max= 1.2
      SER_min= 0.8
      tau_max= 1.E12


c     GGGG
      ratio= f_pre/ f_now

      ratio = min(ratio,2.0)
      ratio = max(ratio,1.05)

      ratio= 1.0 ! dtNT fixed
      dtNT = dtNT * ratio

      if (dtNT.gt.max_dtNT.or.f_now.lt.min_ftol) dtNT = dtNT

      
c      ifntcf = .false. ! compute f at NT(0)
      ifntcg = .true. ! compute g at NT(0)

      return
      end
c----------------------------------------------------------------------- 
      subroutine get_dt(dt_out)
c     not used for now
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real dt_out

      return
      end
c----------------------------------------------------------------------- 
      subroutine bdf1_in_newton(c_in,c_out,istepNT,timeNT,nfld,npvts)
c     fully couple
c     a clean copy-in-copy-out
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /jfnk_tmp/ BDFcnt
      integer BDFcnt

      common /cgeom/ igeom  ! for breaking nek_advance, Lan
      integer igeom,ntot

      real time_loc
      save time_loc
      data time_loc /0.0/

      integer nfld,ic,npvts(toteq)
      real c_in (lx1*ly1*lz1*lelv,toteq)
     $   , c_out(lx1*ly1*lz1*lelv,toteq)
      real timeNT
      integer istepNT

      istep = 1
      time_loc = time_loc + 1.E-5
      time     = time_loc

      call copy_newton2nek(c_in) ! copy c_in into vx vy vz sro srk

c     Nek5000 BDF1/EXT
      call nek_advance ! main solver

cc     BBBB Backup plan: break BDF1 to make seperation
c       call nekgsync
c       call setup_convect(2) ! Save conv vel
c       if (iftran) call settime
c       if (ifmhd ) call cfl_check
c       call setsolv
c       call comment
c
c       call setprop         ! v,k,w -> uservp -> vdiff,vtran
c      igeom=1
c       call heat(igeom)     ! makeq 
c       call fluid(igeom)    ! makef <- ffx_new
c      igeom=2
cc       call setprop         ! v,k,w -> uservp -> vdiff,vtran
c       call heat(igeom)     ! vdiff,vtrans -> h1 h2
c       call fluid(igeom)    ! solve v

      call userchk_nt
      call copy_nek2newton(c_out) ! copy vx vy vz sro srk into c_out
      time = timeNT
      BDFcnt = BDFcnt + 1

      return
      end
c----------------------------------------------------------------------- 
      subroutine chk_amax_self(s3,a,nfld,ord,nss)
c     check max abs of multiple arrays
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      character*3 s3
      integer lt,nfld,ord(toteq)
      parameter (lt=lx1*ly1*lz1*lelt)
      integer nss(toteq),ic,ii
      real a(lt,toteq),amx,glamax,glmin,glmax

      do ii=1,nfld
        ic = ord(ii)
        amx=glamax(a(1,ic),nss(ic))
        if(nid.eq.0) write(6,*) 'Lan Lcheck amx: ',s3,ic,amx
        amx=glmax(a(1,ic),nss(ic))
        if(nid.eq.0) write(6,*) 'Lan Lcheck max: ',s3,ic,amx
        amx=glmin(a(1,ic),nss(ic))
        if(nid.eq.0) write(6,*) 'Lan Lcheck min: ',s3,ic,amx
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine chk_diff_self(s3,a,b,nfld,ord,nss)
c     check max abs difference of mutiple arrays
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      character*3 s3
      integer lt,nfld,ord(toteq)
      parameter (lt=lx1*ly1*lz1*lelt)
      integer nss(nfld),ic,ii
      real a(lt,nfld),b(lt,nfld),tmp(lt,nfld),amx,glamax,all_mx

      all_mx=0.0
      do ii=1,nfld
        ic = ord(ii)
        call sub3(tmp(1,ic),a(1,ic),b(1,ic),nss(ic))

        amx=glamax(tmp(1,ic),nss(ic))
        if(nid.eq.0) write(6,*) 'Lan Lcheck diff(amx): ',s3,ic,amx
        all_mx = max(amx,all_mx)
      enddo
        if(nid.eq.0) write(6,*) 'Lan Lcheck diff(amx): ',s3,'all',all_mx

      return
      end
c----------------------------------------------------------------------- 
      subroutine chk_l2_self(s3,a,b,nfld,ord,nss,mults)
c     check max abs difference of mutiple arrays
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      character*3 s3
      integer lt,nfld
      parameter (lt=lx1*ly1*lz1*lelt)
      integer nss(nfld),ic,ii,ord(toteq)
      real a(lt,nfld),b(lt,nfld),tmp(lt,nfld),amx,glsc3,all_mx
      real mults(lt,toteq)

      all_mx=0.0
      do ii=1,nfld
        ic = ord(ii)
        call sub3(tmp(1,ic),a(1,ic),b(1,ic),nss(ic))

        amx=glsc3(tmp(1,ic),tmp(1,ic),mults(1,ic),nss(ic))
        if(nid.eq.0) write(6,*) 'Lan Lcheck diff(l2): ',s3,ic,amx
        all_mx = amx+all_mx
      enddo
        if(nid.eq.0) write(6,*) 'Lan Lcheck diff(l2): ',s3,'all',all_mx

      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_NTchk(out_fld,iter)
      ! called by each NT iter
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'NEWTON_loc'
      include 'STDY'
      
      real out_fld(lvt1,toteq)
      integer iter,iocomm_nt

c      iocomm_nt = 5
c      if (ifdbgnt) iocomm_nt = 1

      if (ifdbgnt) then
c       call SS_output(ck_0,istepNT,timeNT,nfld,npvts)
c        call SS_vs_ref_solution(out_fld,iter)
      endif

      istep = iter
      if (mod(istep,20).eq.0)  then
        call copy_newton2nek(out_fld)
        call outpost(vx,vy,vz,pr,t,'ooo')
c        call sns_psi_omega ! kav
      endif
      istep = 1

      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_PSchk
      ! called by each PS iter
      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_output(out_fld,timeNT)
c     write customized output, measure error, etc
c----------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON_loc'
      include 'STDY'
      real timeNT
      real out_fld(lvt1,toteq)

      istep = istepNT ! FIXME
      time = timeNT

      call copy_newton2nek(out_fld)

      if(nid.eq.0)write(*,*)'newton , timeNT',istepNT,timeNT
      if(nid.eq.0)write(*,*)'newton_, timeNT',istepNT,timeNT
      if(nid.eq.0)write(*,*)'f_now, timeNT',istepNT,timeNT
      if(nid.eq.0)write(*,*)'uex, timeNT',istepNT,timeNT


c      if (ifsspso) istep = istepPS
c      if (ifssnto) istep = istepNT


c     Nek5000
      call check_ioinfo
      call set_outfld
      call userchk_nt ! will output initial fields inside
      call sns_psi_omega ! kav
      call prepost (ifoutfld,'his')
      call chk_diff_self('uex',bdf_saved(1,1,1),out_fld
     $                 ,nfld,ss_ord,npvts)

      istep = 1

      return
      end
c----------------------------------------------------------------------- 
      subroutine SSop_copy(c_out,c_in,nfld,ord,npvts)
      implicit none
      include 'SIZE'
      include 'SOLN'
      real c_out(lvt1,toteq),c_in(lvt1,toteq)
      integer nfld,npvts(toteq),ii,ic,ord(toteq)

      do ii=1,nfld
        ic = ord(ii)
        call copy(c_out(1,ic),c_in(1,ic),npvts(ic))
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine SSop_add2(a,b,nfld,ord,npvts)
      implicit none
      include 'SIZE'
      include 'SOLN'
      real a(lvt1,toteq),b(lvt1,toteq)
      integer nfld,npvts(toteq),ii,ic,ord(toteq)

      do ii=1,nfld
        ic = ord(ii)
        call add2(a(1,ic),b(1,ic),npvts(ic))
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine SSop_add3(a,b,c,nfld,ord,npvts)
      implicit none
      include 'SIZE'
      include 'SOLN'
      real a(lvt1,toteq),b(lvt1,toteq),c(lvt1,toteq)
      integer nfld,npvts(toteq),ii,ic,ord(toteq)

      do ii=1,nfld
        ic = ord(ii)
        call add3(a(1,ic),b(1,ic),c(1,ic),npvts(ic))
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine SSop_add2s2(c_out,c_in,alpha,nfld,ord,npvts)
      implicit none
      include 'SIZE'
      include 'SOLN'
      real c_out(lvt1,toteq),c_in(lvt1,toteq),alpha
      integer nfld,npvts(toteq),ii,ic,ord(toteq)

      do ii=1,nfld
        ic = ord(ii)
        call add2s2(c_out(1,ic),c_in(1,ic),alpha,npvts(ic))
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine SSop_chsign(c_in,nfld,ord,npvts)
      implicit none
      include 'SIZE'
      include 'SOLN'
      real c_in(lvt1,toteq)
      integer nfld,npvts(toteq),ii,ic,ord(toteq)

      do ii=1,nfld
        ic = ord(ii)
        call chsign(c_in(1,ic),npvts(ic))
      enddo

      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_get_Jeps(eps,p,uc)
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'STDY'

      real eps,p(lvt1,toteq),uc(lvt1,toteq)
      real eps_max,eps_min
      real pnorm,unorm,glsc3
      integer ii,ic,i,n_tmp

      eps_min = 1.E-6 ! JJJJ 
      eps_max = 1.E-2

      if (Jeps_op.eq.1) then
        pnorm = 0.0
        unorm = 0.0
        do ii=1,nfld
        ic = ss_ord(ii)
          pnorm = pnorm+glsc3(p(1,ic),p(1,ic),mults(1,ic),npvts(ic))
          unorm = unorm+glsc3(uc(1,ic),uc(1,ic),mults(1,ic),npvts(ic))
        enddo
        pnorm = sqrt(pnorm+1.e-16)
        unorm = sqrt(unorm)
        eps = (1+unorm)*1.e-14 
        eps = sqrt(eps)      
        eps = eps/pnorm

      elseif (Jeps_op.eq.2) then
        pnorm = 0.0
        unorm = 0.0
        n_tmp = 0
        do ii=1,nfld
        ic = ss_ord(ii)
        do i=1,npvts(ic)
          unorm = unorm + abs(uc(i,ic))
          pnorm = pnorm + p(i,ic)*p(i,ic)*mults(i,ic)
          n_tmp = n_tmp + npvts(ic) 
        enddo
        enddo
        unorm = unorm / real(n_tmp)
        pnorm = sqrt(pnorm)
        eps = (unorm/pnorm + 1)*1.E-6

      else ! Jeps = 0
        eps = Jeps
        return ! overwrite eps_min/max
      endif

      eps = max(eps,eps_min)
      eps = min(eps,eps_max)

      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_get_sc_norms(c_in,ord,norms)
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'STDY'
      real norms(toteq),c_in(lvt1,toteq),glsc3,tmp,tmpv,tmp1,eps
      integer ic,ii,ord(toteq)
     
      tmp1 = 0.0 
      tmpv = 0.0 
      eps= 1.E-16
      do ii=1,nfld
        ic = ord(ii)
        norms(ic)=sqrt(max(glsc3(c_in,c_in,mults(1,ic),npvts(ic)),eps))
        if (ii.eq.1) tmp1 = tmp1 + norms(ic)
        if (ic.le.ndim) tmpv = tmpv + norms(ic)
      enddo
      tmpv = tmpv / real(ndim)

      tmp = tmp1                      ! normalize by fld 1
      if (ord(1).le.ndim) tmp = tmpv ! fld1 = vel
      tmp = max(tmp,eps)
      tmp = min(tmp,1./eps)
      tmp = max(tmp,1.E-2)
      tmp = min(tmp,1.E2)

      do ii=1,nfld
        ic = ord(ii)
        if (ic.le.ndim) then
          norms(ic) = tmpv/tmp
          norms(ic) = 1.0 !tmpv/tmp
        else
          norms(ic) = norms(ic)/tmp
          norms(ic) = 1.0 !norms(ic)/tmp
        endif
      enddo
      return
      end
c----------------------------------------------------------------------- 
      subroutine SS_vs_ref_solution(out_fld,iter)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'STDY'

      real out_fld(lvt1,toteq)
      real bak_fld(lvt1,toteq)
      integer iter

      character*80 s80
      integer icalld
      data icalld /0/
      save icalld

      real p67
      real ref_fld
      common /NTit_chk/ ref_fld(lvt1,toteq)

      if (icalld.eq.0) then
        s80 = 'ref.fld'

        call copy_nek2newton(bak_fld)

        p67 = param(67)
        param(67) = 6.00
        call chcopy (initc,s80,80)
        call bcast  (initc,80)
        call restart(1)
        param(67) = p67

        call copy_nek2newton(ref_fld)
        call copy_newton2nek(bak_fld)

      icalld = 1
      endif
      
      if(nid.eq.0)write(*,*)'NTe',iter
      call chk_diff_self('NTe',out_fld,ref_fld,nfld,ss_ord,npvts)
      call chk_l2_self('NTe',out_fld,ref_fld,nfld,npvts,ss_ord,mults)

      return
      end
c----------------------------------------------------------------------- 
