c-----------------------------------------------------------------------
c     Inside cem_drift_newton.F we use Jacobi-free Newton Krylov method
c     (JFNK) to solve the steady state of PNP equation
c     By Ping Hsuan 06/27/2017
c     Modified by Yu-Hsiang Lan
c---------------------------------------------------------------------
c     This subroutine is the main core of JFNK method which involves
c     using variables cN and cP in DRIFT
      subroutine cem_drift_newton
c-----------------------------------------------------------------------      
      include 'SIZE'
      include 'TOTAL'
      include 'NGMRES'
      include 'NEWTON'
      include 'CTIMER'
c     include 'RTIMER'

      parameter (lcdim=ldim)

      real alphan,tt
      real tolNT,tolGMRES
      real glmax,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real temp(lx1*ly1*lz1*lelt)

      integer i,j,k,n,ic
      integer maxit
      real glsum
      real f_pre,f_now
      real max_dtNT,min_tol      
      real sk  (lx1*ly1*lz1*lelt,lcdim),  !
     $     fk  (lx1*ly1*lz1*lelt,lcdim),  !
     $     gk  (lx1*ly1*lz1*lelt,lcdim),  !
     $     cn_n(lx1*ly1*lz1*lelt,lcdim),  ! working array in Pseudo time step
     $     cn_k(lx1*ly1*lz1*lelt,lcdim)   ! working array in Newton iteration

      common /ctmp1e/ uxe (lx1,ly1,lz1,lelv), uye (lx1,ly1,lz1,lelv)

C     Parametar settings ToDo: move to user file, Lan

      alphan = 1.            ! Newton relaxation parameter alphan in (0,1]
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxit = 100           ! #iter of Newton method
      jaceps = 1.e-8        ! eps for approx. Jacobian
      epsinv = 1./jaceps
      tolGMRES = 1.e-10     ! tol for GMRES
      tolNT = 1.e-5         ! tol for Newton method

      dtNT = 1.!param(1) ! pseudo-transient time step
                      ! it can grow as istep increase 
      max_dtNT  = 1E8   ! Max. of dtNT, avoid NaN
      min_tol   = 1E-12 ! Min. of tol,  avoid NaN

      dtinv = 1./dt
      dtNTinv = 1./dtNT

      cpu_t = 0.0
      cpu_chk = 0.0

      npts = lx1*ly1*lz1*nelv

C     Initialization
      do ic = 1,lcdim
        call rzero(sk(1,ic),npts)
        call rzero(fk(1,ic),npts)
        call rzero(gk(1,ic),npts)
        call rzero(cn_n(1,ic),npts) ! working array in Pseudo time
        call rzero(cn_k(1,ic),npts) ! working array in NT
      enddo

      call copy(cn_k(1,1),vx,npts)
      call copy(cn_k(1,2),vy,npts)
      if (ldim.eq.3) call copy(cn_k(1,3),vz,npts)

      call copy(cn_n(1,1),vx,npts)
      call copy(cn_n(1,2),vy,npts)
      if (ldim.eq.3) call copy(cn_n(1,3),vz,npts)

      call compute_f_nt(fk,cn_k,lcdim,npts) ! compute f(u_0)
      call compute_gk(gk,fk,cn_k,cn_n,lcdim,npts) ! compute rhs g_k for gmres

      if (nio.eq.0) then
      do i=1,lx1*ly1*lz1*nelv
         write (6,*) 'vx=',fk(i,1)
         write (6,*) 'vy=',fk(i,2)
      enddo
      endif
c     stop

      do ic = 1,lcdim
        call chsign (gk(1,ic),npts)
      enddo
      
      do ic = 1,lcdim
        fnorms(ic) = glsc3(fk(1,ic),fk(1,ic),vmult,npts)
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm = glmax(fnorms,lcdim)
      f_now = fnorm
      call outpost(cn_n(1,1),cn_n(1,2),vz,pr,t,'   ')

C     Start pseudo-time step
      do istepn = 1,nsteps 
         if (nio.eq.0) write (6,*) 'pseudo step',istepn
         if (mod(istepn,iostep).eq.0) then
            istep=istepn
            call outpost(cn_n(1,1),cn_n(1,2),vz,pr,t,'   ')
            istep=1
         endif
         if (nio.eq.0) write (6,*) 'istepn=',istepn

c       cpu_dtime = dclock() !FIXME

C       SER, dt, dtNT varying
        write (6,*) 'f_pre,f_now',f_pre,f_now
        if ((istepn .eq. 1)) then
          dtNT = param(1)
c          dt = dtNT
        else
          if (f_now.lt.min_tol) f_now = min_tol
          dtNT = dtNT * f_pre/f_now
c            dt = dtNT
        endif
        if (dtNT.gt.max_dtNT) dtNT = max_dtNT

        dtinv = 1./dt
        dtNTinv = 1./dtNT
         write (6,*) 'before newton'

C       Start Newton iteration
        do iterNT=1,maxit
            write (6,*) 'iternt=',iternt

          do ic = 1,lcdim ! solving newton sk
            write (6,*) 'ic=',ic

c           call hmh_gmres_newton4(gk(1,ic),vmult,iter)
            
c           call hmh_gmres_newton3(sk(1,ic),gk(1,ic),vmult,iter,
c    $                             cn_k,fk(1,ic),npts,tolgmres,ic,ldim)
c           call hmh_gmres_newton3(gk(1,ic),vmult,iter)
c           call hmh_gmres_newton2(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),npts,
c    $                             tolgmres,ic,vmult,iter)

c           call hmh_gmres_newton(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),vmult,
c    $                            npts,tolgmres,ic)
c           stop

            call drift_hmh_gmres_newton(sk(1,ic),gk(1,ic)
     $           ,cn_k,fk(1,ic),vmult,ldim,npts,tolGMRES,ic)

c           write (6,*) 'gmres iter=',iter
c           stop
          enddo
          write (6,*) 'newton int=',iternt
c           stop

          write (6,*) 'check 1'
          do ic = 1,lcdim ! update Newton iteration
            call add2s2(cn_k(1,ic),sk(1,ic),alphan,npts) ! cp_k=cp_k+alphan*sp_k
          enddo
          write (6,*) 'check 2'

        call compute_f_nt(fk,cn_k,lcdim,npts)
        call compute_gk(gk,fk,cn_k,cn_n,lcdim,npts)  ! checking tol. + next NT iter
        write (6,*) 'check 3'
        do ic = 1,lcdim
           call chsign(gk(1,ic),npts)
        enddo
        write (6,*) 'check 4'

        do ic = 1,lcdim
           rnorms(ic) = glsc3(gk(1,ic),gk(1,ic),vmult,npts)
           rnorms(ic) = sqrt(rnorms(ic))
        enddo

        write (6,*) 'check 5'
        rnorm = glmax(rnorms,lcdim)
        if (iterNT.eq.1) rinorm=rnorm
        ratio=rnorm/rinorm 
        write (6,*) 'check 6'

c       if ((nid.eq.0).and.(mod(istepn,iocomm).eq.0)) 
        write(6,90)
     $     istepn,iterNT,ratio,rnorm,rinorm,dtNT,dt
        write (6,*) 'check 7'

        if (ratio.lt.tolNT) goto 900 ! Newton conv
c       if (rnorm.lt.new_tol) goto 900 ! Newton conv

        enddo
 90     format('newton iter',2i6,1p5e12.4)   
900     continue

        do ic = 1,lcdim
           call copy(cn_n(1,ic),cn_k(1,ic),npts)   ! update cn_n for computing gn
        enddo                                     ! and next pseudo-time step

c       Compute the norm of f
        do ic = 1,lcdim
           fnorms(ic) = glsc3(fk(1,ic),fk(1,ic),vmult,npts)  
           fnorms(ic) = sqrt(fnorms(ic))
        enddo

        fnorm  = glmax(fnorms,lcdim)
        f_pre = f_now  ! set up old fnorm to f_pre
        f_now = fnorm  ! set up new fnorm to f_now

c       Compute the CPU_time
c       cpu_dtime = dclock()-cpu_dtime !FIXME

        cpu_t = cpu_t+cpu_dtime
        cpu_t_step = cpu_t/istep
        cpu_p_t = glsum(cpu_t_step /npts,1)/np

        time = time + dt

c       call errchk(cn_k(1,1),cn_k(1,2))

C       Copy for cem_out and cem_chk
c       do ic = 1,lcdim
c         call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
c       enddo
c       call userchk
c       call cem_out

c       call errorchk(cn_k(1,1),cn_k(1,2))

      enddo
      
c     call cem_end !FIXME

      return
      end
c-----------------------------------------------------------------------
c     This routine computes gk for each Newton iteration
c
c          g^{n}_k = u^{n}_k - \delta t f(u^{n}_k) - u^{n}_0
c
c     Input    ckn
c              c0n denote the initial of each Newton iteration 
c     Input    f
c     Output   gout
      subroutine compute_gk(gout,f,ckn,c0n,nion,n)
c-----------------------------------------------------------------------      
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'
      integer nion,n,ic
      real gout(lx1*ly1*lz1*lelt,nion)
     $   , ckn (lx1*ly1*lz1*lelt,nion)
     $   , c0n (lx1*ly1*lz1*lelt,nion)
     $   , f   (lx1*ly1*lz1*lelt,nion)

      do ic = 1,nion
c        write (6,*) 'ic=',ic
c        call add3s3(gout(1,ic),ckn(1,ic),f(1,ic),c0n(1,ic)
c    $       ,1.0,-1.0*dtNT,-1.0,n)

c        call sub3s2(gout,ckn,f,1.,-dtnt,n)

         call add3s2(gout(1,ic),ckn(1,ic),f(1,ic),1.,-dtnt,n)
         write (6,*) 'after add3s2'

         call sub2(gout(1,ic),c0n(1,ic),n)
         write (6,*) 'after sub2'
         call cmult(gout(1,ic),dtNTinv,n)   ! case with divided by dt_newton
         write (6,*) 'after cmult'
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f_nt(fout,cin,nion,n)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      real fout(lx1*ly1*lz1*lelt,nion) ! f output
      real cin (lx1*ly1*lz1*lelt,nion) ! c input
      real cout(lx1*ly1*lz1*lelt,nion) ! c output (sem_bdf)

      write (6,*) 'check first'

      call cem_drift_sem_bdf1_newton(cin,cout,nion,n,0)

      write (6,*) 'check second'

      do ic = 1,nion
        call sub3(fout(1,ic),cout(1,ic),cin(1,ic),n)
        call cmult(fout(1,ic),dtinv,n)
      enddo

      write (6,*) 'check third'

      return
      end
c-----------------------------------------------------------------------
c     This routine computes Jp where J is the Jacobian matrix and p 
c     is a vector. Note that we do not construct the Jacobian matrix 
c     exactly but use the following formula to compute Jp
c   
c        J_k s_k = s_k - (dt /eps) *  ( f(u_k + eps*s_k) - f(u_k) )
c   
c     where f(u_k) has been store in each Newton iteration
c           f(u_k + eps*s_k) has to be computed in each GMRES iteration
      subroutine jacobimatvec(jacp,p,uc,fi,nion,n,iflag)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      real     jacp(lx1*ly1*lz1*lelt),p(lx1*ly1*lz1*lelt)
      real     uc (lx1*ly1*lz1*lelt,nion)     ! treat other ions as variable
      real     fi (lx1*ly1*lz1*lelt)          ! f_in
      real     fo (lx1*ly1*lz1*lelt,nion)     ! f_out (multi-ions)
      real     foi(lx1*ly1*lz1*lelt)          ! f_out - f_in (restricted on iflag)
      real     uep(lx1*ly1*lz1*lelt,nion)     ! uc + eps p

C     Determine eps with norms of p and u_iflag
      pnorm = glsc3(p,p,vmult,n)
      pnorm = sqrt(pnorm)
      unorm = glsc3(uc(1,iflag),uc(1,iflag),vmult,n)
      unorm = sqrt(unorm)

      eps = (1+unorm)*1e-14   ! formula for varing eps
      eps = sqrt(eps)
      eps = eps/pnorm
      epsinv = 1./eps         

C     Compute Jp

      do ic = 1,nion
        call copy(uep(1,ic),uc(1,ic),n)
      enddo

      call add3s2(uep(1,iflag),uc(1,iflag),p,1.0,eps,n) ! uep = u + eps*p

      call compute_f_nt(fo,uep,nion,n)

      call sub3  (foi,fo(1,iflag),fi,n)
      call cmult (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3  (jacp,p,foi,n)
      call cmult (jacp,dtntinv,n)      ! Jp = p/dt_newton - (fo-f)/eps
                                     ! case with divided by dt_newton
c     call copy(jacp,p,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,uc,f,wt,nion,n,tol,iflag)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
c     include 'FDMH1'
c     include 'GMRES1'
      include 'NGMRES'

      integer  n,nion,outer,iflag
      real     phi(lx1*ly1*lz1*lelt),res(lx1*ly1*lz1*lelt),
     $         wt(lx1*ly1*lz1*lelt)
      real     tol,alphan,l,temp
      real     eps,uc(lx1*ly1*lz1*lelt,nion),f(lx1*ly1*lz1*lelt)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x_gmres,n)
      call rzero(h_gmres,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r_gmres,res,n)                  ! r = res
         else
            !update residual
            call copy   (r_gmres,res,n)                  ! r = res
            call JacobiMatVec(w_gmres,x_gmres,uc,f,nion,n,iflag)! w = A x
            call add2s2 (r_gmres,w_gmres,-1.,n) ! r = r - w
         endif

         gamma_gmres(1) = glsc3(r_gmres,r_gmres,wt,n)                ! gamma_gmres  = (r,r)
         gamma_gmres(1) = sqrt(gamma_gmres(1))                 ! gamma_gmres  = sqrt{ (r,r) }

         tolps = 0.1*gamma_gmres(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable
         tolpsss = tolpss

         write (6,*) 'tolps,tolpss',tolps,tolpss

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,n)             !  v  = r / gamma_gmres
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call JacobiMatVec(w_gmres,v_gmres(1,j),uc,f,nion,n,iflag) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h_gmres(i,j)=glsc3(w_gmres,v_gmres(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $         + s_gmres(i)*h_gmres(i+1,j)
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $         + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                              !            ______
            alphan = sqrt(glsc3(w_gmres,w_gmres,wt,n))     ! alphan =  \/ (w,w)
            if(alphan.eq.0.) goto 900 !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alphan*alphan)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alphan  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            rnorm = abs(gamma_gmres(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c                 write (6,*) rnorm,'newton'
                  write (6,66) iter,tolpsss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolpsss) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alphan
            call cmult2(v_gmres(1,j+1),w_gmres,temp,n)   ! v    = w / alphan
                                             !  j+1
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma_gmres
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,v_gmres(1,i),c_gmres(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x_gmres,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

c     if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,newton_iter,iter,tolpss
c     endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end

c---------------------------------------------------------------------
C     This subroutine compute sem_bdf1, user can provide any(*) 
C     bdf1 time swapper in here
      subroutine cem_drift_sem_bdf1_newton(cinput,coutput,nion,n,iflag)
c---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
c     include 'BCS'

      integer iflag,nion,n,ic
      real cinput(lx1*ly1*lz1*lelt,nion)
      real coutput(lx1*ly1*lz1*lelt,nion)

      real dummy1(lx1,ly1,lz1,lelt)
      real dummy2(lx2,ly2,lz2,lelt)

C     Copy into working array of sem_bdf1

      write (6,*) 'check pre'

      if (ldim.eq.3) call opcopy(vx,vy,vz,
     $                           cinput(1,1),cinput(1,2),cinput(1,3))

      if (ldim.eq.2) call opcopy(vx,vy,vz,
     $                           cinput(1,1),cinput(1,2),dummy1)

      call plan5(2)

      call opcopy(coutput(1,1),coutput(1,2),dummy1,vx,vy,vz)

      write (6,*) 'check post'

C     Copy out working array of sem_bdf1

      return
      end
c-----------------------------------------------------------------------
c     subroutine hmh_gmres_newton(res,h1,h2,wt,iter)
      subroutine hmh_gmres_newton(out,res,uc,f,wt,n,tol,isd)

c           call hmh_gmres_newton(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),vmult,
c    $                            lcdim,npts,tolgmres,ic)

c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.

      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'NGMRES'

      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint

      real             res  (lx1*ly1*lz1*lelv)
c     real             h1   (lx1,ly1,lz1,lelv)
c     real             h2   (lx1,ly1,lz1,lelv)
      real             wt   (lx1,ly1,lz1,lelv)

      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)

      real alpha, l, temp
      integer outer

      logical iflag,if_hyb
      save    iflag,if_hyb
c     data    iflag,if_hyb  /.false. , .true. /
      data    iflag,if_hyb  /.false. , .false. /
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      n = nx1*ny1*nz1*nelv

      etime1  = dnekclock()
      etime_p = 0.
      divex   = 0.
      iter    = 0
      m       = lgmres

      write (6,*) 'gm 1'

      if (.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml_gmres,mu_gmres,bm1,binvm1,
     $                          nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif

      write (6,*) 'gm 2'

c     if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
      do i=1,lx1*ly1*lz1*nelv
         d(i) = 1.
      enddo

      write (6,*) 'gm 3'

c     call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)

c     if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
c    $   tolps = abs(param(21))
c     if (istep.eq.0) tolps = 1.e-4
      tolps=tol
      tolpss = tolps

      iconv = 0
      call rzero(x_gmres,n)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1

         if (iter.eq.0) then                   !      -1
            call col3(r_gmres,ml_gmres,res,n) ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r_gmres,res,n)           ! r = res

            call jacobimatvec(w_gmres,x_gmres,uc,f,ldim,n,isd)
            call copy(w_gmres,x_gmres,n)
c           call ax    (w_gmres,x_gmres,h1,h2,n) ! w = A x
            call add2s2(r_gmres,w_gmres,-1.,n)   ! r = r - w
                                                 !      -1
            call col2(r_gmres,ml_gmres,n)        ! r = L   r
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wt,n)) ! gamma  = \/ (r,r) 
                                                            !      1
         if(iter.eq.0) then
            div0 = gamma_gmres(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if (gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,n) ! v  = r / gamma
                                                  !  1            1

         ifmgrid=.false.

         do j=1,m
            iter = iter+1
                                                       !       -1
            call col3(w_gmres,mu_gmres,v_gmres(1,j),n) ! w  = U   v
                                                       !           j

c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .

            etime2 = dnekclock()

c           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
            if (ifmgrid) then
               call h1mg_solve(z_gmres(1,j),w_gmres,if_hyb) ! z  = M   w
            else                                            !  j
               kfldfdm = ndim+1
               if (param(100).eq.2) then
                   call h1_overlap_2 (z_gmres(1,j),w_gmres,pmask)
               else
                   call fdm_h1
     $               (z_gmres(1,j),w_gmres,d,pmask,vmult,nelv,
     $                ktype(1,1,kfldfdm),wk)
               endif
               call crs_solve_h1 (wk,w_gmres)        ! z  = M   w
               call add2         (z_gmres(1,j),wk,n) !  j        
            endif

            call ortho        (z_gmres(1,j)) ! Orthogonalize wrt null space, if present
            etime_p = etime_p + dnekclock()-etime2
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

            call jacobimatvec(w_gmres,z_gmres(1,j),uc,f,ldim,n,isd)
            call copy(w_gmres,z_gmres(1,j),n)

c           call ax  (w_gmres,z_gmres(1,j),h1,h2,n) ! w = A z
                                                    !        j
     
                                                    !      -1
            call col2(w_gmres,ml_gmres,n)           ! w = L   w

c           !modified Gram-Schmidt

c           do i=1,j
c              h_gmres(i,j)=glsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c                                                            !  i,j       i

c              call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
c           enddo                                                !          i,j  i

c           2-PASS GS, 1st pass:

            do i=1,j
               h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
            enddo                                            !  i,j       i

            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
            enddo                                                !          i,j  i


c           2-PASS GS, 2nd pass:
c
c           do i=1,j
c              wk1(i)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c           enddo                                      !  i,j       i
c                                                      !
c           call gop(wk1,wk2,'+  ',j)                  ! sum over P procs
c
c           do i=1,j
c              call add2s2(w_gmres,v_gmres(1,i),-wk1(i),n)    ! w = w - h    v
c              h_gmres(i,j) = h_gmres(i,j) + wk1(i)           !          i,j  i
c           enddo

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)                   
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                        + s_gmres(i)*h_gmres(i+1,j)  
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                      !            ______
            alpha = sqrt(glsc3(w_gmres,w_gmres,wt,n)) ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence newton')

#ifndef TST_WSCAL
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,n) ! v    = w / alpha
                                                       !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),n) ! x = x + c  z
         enddo                                             !          i  i
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      enddo
 9000 continue

      divex = rnorm
      call copy(out,x_gmres,n)

      call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,iter,divex,div0,tolpss,etime_p,
     &                            etime1,if_hyb
c     call flush_hack
 9999 format(4x,i7,'  newton gmres ',4x,i5,1p5e13.4,1x,l4)

      if (outer.le.2) if_hyb = .false.

      return
      end
c-----------------------------------------------------------------------
c     subroutine hmh_gmres_newton2(out,res,uc,f,n,tol,isd,wt,iter)
      subroutine hmh_gmres_newton2(res,wt,iter)
c     subroutine hmh_gmres_newton(out,res,uc,f,wt,n,tol,isd)

c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.

     
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES1'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
c     real             out  (lx1*ly1*lz1*lelv)
      real             res  (lx1*ly1*lz1*lelv)
      real             wt   (lx1,ly1,lz1,lelv)
c     real             f    (lx1,ly1,lz1,lelv)
c     real             uc   (lx1,ly1,lz1,lelv)

      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      real alpha, l, temp
      integer outer

      logical iflag,if_hyb
      save    iflag,if_hyb
c     data    iflag,if_hyb  /.false. , .true. /
      data    iflag,if_hyb  /.false. , .false. /
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock


      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml_gmres,mu_gmres,bm1,binvm1,
     $                          nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif

c     if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
      call rone(d,lx1*ly1*lz1*nelv)

c     call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)

      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      iconv = 0
      call rzero(x_gmres,n)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1

         if(iter.eq.0) then                   !      -1
            call col3(r_gmres,ml_gmres,res,n) ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r_gmres,res,n)           ! r = res
            call copy  (w_gmres,x_gmres,n)
c           call ax    (w_gmres,x_gmres,h1,h2,n) ! w = A x
c           call jacobimatvec(w_gmres,x_gmres,uc,f,ldim,n,isd)
            call add2s2(r_gmres,w_gmres,-1.,n)   ! r = r - w
                                                 !      -1
            call col2(r_gmres,ml_gmres,n)        ! r = L   r
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wt,n)) ! gamma  = \/ (r,r) 
                                                            !      1
         if(iter.eq.0) then
            div0 = gamma_gmres(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,n) ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                       !       -1
            call col3(w_gmres,mu_gmres,v_gmres(1,j),n) ! w  = U   v
                                                       !           j

c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .

            etime2 = dnekclock()

c           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
            if (ifmgrid) then
               call h1mg_solve(z_gmres(1,j),w_gmres,if_hyb) ! z  = M   w
            else                                            !  j
               kfldfdm = ndim+1
               if (param(100).eq.2) then
                   call h1_overlap_2 (z_gmres(1,j),w_gmres,pmask)
               else
                   call fdm_h1
     $               (z_gmres(1,j),w_gmres,d,pmask,vmult,nelv,
     $                ktype(1,1,kfldfdm),wk)
               endif
               call crs_solve_h1 (wk,w_gmres)        ! z  = M   w
               call add2         (z_gmres(1,j),wk,n) !  j        
            endif


            call ortho        (z_gmres(1,j)) ! Orthogonalize wrt null space, if present
            etime_p = etime_p + dnekclock()-etime2
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

     
            call copy(w_gmres,z_gmres(1,j),n)
c           call ax  (w_gmres,z_gmres(1,j),h1,h2,n) ! w = A z
                                                    !        j
     
                                                    !      -1
            call copy(w_gmres,z_gmres(1,j),n)
c           call jacobimatvec(w_gmres,z_gmres(1,j),uc,f,ldim,n,isd)
            call col2(w_gmres,ml_gmres,n)           ! w = L   w

c           !modified Gram-Schmidt

c           do i=1,j
c              h_gmres(i,j)=glsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c                                                            !  i,j       i

c              call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
c           enddo                                                !          i,j  i

c           2-PASS GS, 1st pass:

            do i=1,j
               h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
            enddo                                            !  i,j       i

            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
            enddo                                                !          i,j  i


c           2-PASS GS, 2nd pass:
c
c           do i=1,j
c              wk1(i)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c           enddo                                      !  i,j       i
c                                                      !
c           call gop(wk1,wk2,'+  ',j)                  ! sum over P procs
c
c           do i=1,j
c              call add2s2(w_gmres,v_gmres(1,i),-wk1(i),n)    ! w = w - h    v
c              h_gmres(i,j) = h_gmres(i,j) + wk1(i)           !          i,j  i
c           enddo

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)                   
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                        + s_gmres(i)*h_gmres(i+1,j)  
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                      !            ______
            alpha = sqrt(glsc3(w_gmres,w_gmres,wt,n)) ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence n2')

#ifndef TST_WSCAL
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,n) ! v    = w / alpha
                                                       !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),n) ! x = x + c  z
         enddo                                             !          i  i
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      enddo
 9000 continue

      divex = rnorm

      call copy(res,x_gmres,n)
      call ortho(res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,iter,divex,div0,tolpss,etime_p,
     &                            etime1,if_hyb
c     call flush_hack
 9999 format(4x,i7,'  n2 gmres ',4x,i5,1p5e13.4,1x,l4)

      if (outer.le.2) if_hyb = .false.

      return
      end
c-----------------------------------------------------------------------
c     subroutine hmh_gmres_newton3(out,res,wt,iter)
      subroutine hmh_gmres_newton3(out,res,wt,iter,
     $                             cn,fk,npts,tol,isd,nion)
c           hmh_gmres_newton3(sk(1,ic),gk(1,ic),vmult,iter,
c                             cn_k,fk(1,ic),npts,tolgmres,ic,ldim)

c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.

     
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES1'

      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint

      real out(lx1*ly1*lz1*lelt)
      real res(lx1*ly1*lz1*lelt)
      real wt (lx1,ly1,lz1,lelt)
      real cn (lx1,ly1,lz1,lelt)
      real fk (lx1,ly1,lz1,lelt)

      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      real alpha, l, temp
      integer outer

      logical iflag,if_hyb
      save    iflag,if_hyb
c     data    iflag,if_hyb  /.false. , .true. /
      data    iflag,if_hyb  /.false. , .false. /
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock


      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml_gmres,mu_gmres,bm1,binvm1,
     $                          nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif

c     if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
      call rone(d,lx1*ly1*lz1*nelv)

c     call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      iconv = 0
      call rzero(x_gmres,n)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1

         if(iter.eq.0) then                   !      -1
            call col3(r_gmres,ml_gmres,res,n) ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r_gmres,res,n)           ! r = res
c           call ax    (w_gmres,x_gmres,h1,h2,n) ! w = A x
            call copy  (w_gmres,x_gmres,n)
c           call jacobimatvec(w_gmres,x_gmres,cn,fk,nion,n,isd)
            call add2s2(r_gmres,w_gmres,-1.,n)   ! r = r - w
                                                 !      -1
            call col2(r_gmres,ml_gmres,n)        ! r = L   r
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wt,n)) ! gamma  = \/ (r,r) 
                                                            !      1
         if(iter.eq.0) then
            div0 = gamma_gmres(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,n) ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                       !       -1
            call col3(w_gmres,mu_gmres,v_gmres(1,j),n) ! w  = U   v
                                                       !           j

c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .

            etime2 = dnekclock()

c           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
            if (ifmgrid) then
               call h1mg_solve(z_gmres(1,j),w_gmres,if_hyb) ! z  = M   w
            else                                            !  j
               kfldfdm = ndim+1
               if (param(100).eq.2) then
                   call h1_overlap_2 (z_gmres(1,j),w_gmres,pmask)
               else
                   call fdm_h1
     $               (z_gmres(1,j),w_gmres,d,pmask,vmult,nelv,
     $                ktype(1,1,kfldfdm),wk)
               endif
               call crs_solve_h1 (wk,w_gmres)        ! z  = M   w
               call add2         (z_gmres(1,j),wk,n) !  j        
            endif


c           call ortho        (z_gmres(1,j)) ! Orthogonalize wrt null space, if present
            etime_p = etime_p + dnekclock()-etime2
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

     
c           call copy(w_gmres,z_gmres(1,j),n)
c           call ax  (w_gmres,z_gmres(1,j),h1,h2,n) ! w = A z
                                                    !        j
     
                                                    !      -1
c           call jacobimatvec(w_gmres,z_gmres,cn,fk,nion,n,isd)
            call col2(w_gmres,ml_gmres,n)           ! w = L   w

c           !modified Gram-Schmidt

c           do i=1,j
c              h_gmres(i,j)=glsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c                                                            !  i,j       i

c              call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
c           enddo                                                !          i,j  i

c           2-PASS GS, 1st pass:

            do i=1,j
               h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
            enddo                                            !  i,j       i

            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
            enddo                                                !          i,j  i


c           2-PASS GS, 2nd pass:
c
c           do i=1,j
c              wk1(i)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c           enddo                                      !  i,j       i
c                                                      !
c           call gop(wk1,wk2,'+  ',j)                  ! sum over P procs
c
c           do i=1,j
c              call add2s2(w_gmres,v_gmres(1,i),-wk1(i),n)    ! w = w - h    v
c              h_gmres(i,j) = h_gmres(i,j) + wk1(i)           !          i,j  i
c           enddo

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)                   
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                        + s_gmres(i)*h_gmres(i+1,j)  
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                      !            ______
            alpha = sqrt(glsc3(w_gmres,w_gmres,wt,n)) ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
c           if (ifprint.and.nio.eq.0) 
               write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence newton')
            write (6,*) 'stop after write statement'
            stop

#ifndef TST_WSCAL
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,n) ! v    = w / alpha
                                                       !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),n) ! x = x + c  z
         enddo                                             !          i  i
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      enddo
 9000 continue

      divex = rnorm
      call copy(out,x_gmres,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,iter,divex,div0,tolpss,etime_p,
     &                            etime1,if_hyb
c     call flush_hack
 9999 format(4x,i7,'  PRES gmres ',4x,i5,1p5e13.4,1x,l4)

      if (outer.le.2) if_hyb = .false.

      return
      end
