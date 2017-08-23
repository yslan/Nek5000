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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'   ! lcdim, cn
      include 'POISSON' ! mult
      include 'NEWTON'
      include 'CTIMER'
      include 'RTIMER'

      real alpha,tt
      real tolNT,tolGMRES
      real glmax,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real temp(lpts1)

      integer i,j,k,n,ic
      integer maxit
      real glsum
      real f_pre,f_now
      real max_dtNT,min_tol      
      real sk  (lpts1,lcdim),  !
     $     fk  (lpts1,lcdim),  !
     $     gk  (lpts1,lcdim),  !
     $     cn_n(lpts1,lcdim),  ! working array in Pseudo time step
     $     cn_k(lpts1,lcdim)   ! working array in Newton iteration

C     Parametar settings ToDo: move to user file, Lan

      alpha = 1            ! Newton relaxation parameter alpha in (0,1]
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxit = 20           ! #iter of Newton method
      jaceps = 1e-5        ! eps for approx. Jacobian
      epsinv = 1./jaceps
      tolGMRES = 1e-10     ! tol for GMRES
      tolNT = 1e-5         ! tol for Newton method

      dtNT = param(1) ! pseudo-transient time step
                      ! it can grow as istep increase 
      max_dtNT  = 1E5   ! Max. of dtNT, avoid NaN
      min_tol   = 1E-12 ! Min. of tol,  avoid NaN

      dtinv = 1./dt
      dtNTinv = 1./dtNT

      cpu_t = 0.0
      cpu_chk = 0.0

C     Initialization
      do ic = 1,lcdim
        call rzero(sk(1,ic),npts)
        call rzero(fk(1,ic),npts)
        call rzero(gk(1,ic),npts)
        call rzero(cn_n(1,ic),npts) ! working array in Pseudo time
        call rzero(cn_k(1,ic),npts) ! working array in NT
      enddo
      
      call compute_f(fk,cn_k,lcdim,npts)
      call compute_gk(gk,fk,cn_k,cn_n,lcdim,npts) ! compute rhs g_k for GMRES
      do ic = 1,lcdim
        call chsign (gk(1,ic),npts)
      enddo

      do ic = 1,lcdim
        fnorms(ic) = glsc3(fk(1,ic),fk(1,ic),mult,npts)
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm = glmax(fnorms,lcdim)
      f_now = fnorm


C     Start pseudo-time step
      do istep = 1,nsteps 

        cpu_dtime = dclock()

C       SER, dt, dtNT varying
        if ((istep .eq. 1)) then
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


C       Start Newton iteration
        do iterNT=1,maxit

          do ic = 1,lcdim ! solving newton sk
            call drift_hmh_gmres_newton(sk(1,ic),gk(1,ic)
     $           ,cn_k,fk(1,ic),mult,lcdim,npts,tolGMRES,ic)
          enddo
          do ic = 1,lcdim ! update Newton iteration
            call add2s2(cn_k(1,ic),sk(1,ic),alpha,npts) ! cp_k=cp_k+alpha*sp_k
          enddo

          call compute_f(fk,cn_k,lcdim,npts)
          call compute_gk(gk,fk,cn_k,cn_n,lcdim,npts)  ! checking tol. + next NT iter
          do ic = 1,lcdim
            call chsign(gk(1,ic),npts)
          enddo

          do ic = 1,lcdim
            rnorms(ic) = glsc3(gk(1,ic),gk(1,ic),mult,npts)
            rnorms(ic) = sqrt(rnorms(ic))
          enddo
          rnorm = glmax(rnorms,lcdim)
          if (iterNT.eq.1) rinorm=rnorm
          ratio=rnorm/rinorm 

          if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) write(6,90)
     $       istep,iterNT,ratio,rnorm,rinorm,dtNT,dt

          if (ratio.lt.tolNT) goto 900 ! Newton conv
c          if (rnorm.lt.new_tol) goto 900 ! Newton conv

        enddo
 90     format('newton iter',2i6,1p5e12.4)   
900     continue

        do ic = 1,lcdim
          call copy(cn_n(1,ic),cn_k(1,ic),npts)   ! update cn_n for computing gn
        enddo                                     ! and next pseudo-time step

c       Compute the norm of f
        do ic = 1,lcdim
          fnorms(ic) = glsc3(fk(1,ic),fk(1,ic),mult,npts)  
          fnorms(ic) = sqrt(fnorms(ic))
        enddo
        fnorm  = glmax(fnorms,lcdim)
        f_pre = f_now  ! set up old fnorm to f_pre
        f_now = fnorm  ! set up new fnorm to f_now

c       Compute the CPU_time
        cpu_dtime = dclock()-cpu_dtime
        cpu_t = cpu_t+cpu_dtime
        cpu_t_step = cpu_t/istep
        cpu_p_t = glsum(cpu_t_step /npts,1)/np
        
        time = time + dt

C       Copy for cem_out and cem_chk
        do ic = 1,lcdim
          call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
        enddo
        call userchk
        call cem_out
      
      enddo
      
      call cem_end

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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'
      integer nion,n,ic
      real gout(lpts1,nion)
     $   , ckn (lpts1,nion)
     $   , c0n (lpts1,nion)
     $   , f   (lpts1,nion)

      do ic = 1,nion
        call add3s3(gout(1,ic),ckn(1,ic),f(1,ic),c0n(1,ic)
     $       ,1.0,-1.0*dtNT,-1.0,n)
        call cmult(gout(1,ic),dtNTinv,n)   ! case with divided by dt_newton
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f(fout,cin,nion,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'
      include 'POISSON' 
      integer nion,n,ic
      real fout(lpts1,nion) ! f output
      real cin (lpts1,nion) ! c input
      real cout(lpts1,nion) ! c output (sem_bdf)

      call cem_drift_sem_bdf1_newton(cin,cout,nion,n,0)

      do ic = 1,nion
        call sub3(fout(1,ic),cout(1,ic),cin(1,ic),n)
        call cmult(fout(1,ic),dtinv,n)
      enddo

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
      subroutine JacobiMatVec(Jp,p,uc,fi,nion,n,iflag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'
      include 'POISSON'

      real     Jp (lpts1),p(lpts1)
      real     uc (lpts1,nion)     ! treat other ions as variable
      real     fi (lpts1)          ! f_in
      real     fo (lpts1,nion)     ! f_out (multi-ions)
      real     foi(lpts1)          ! f_out - f_in (restricted on iflag)
      real     uep(lpts1,nion)     ! uc + eps p
      integer  nion,n,iflag,ic
      real     eps,pnorm,unorm,glsc3 

C     Determine eps with norms of p and u_iflag
      pnorm = glsc3(p,p,mult,n)
      pnorm = sqrt(pnorm)
      unorm = glsc3(uc(1,iflag),uc(1,iflag),mult,n)
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

      call compute_f(fo,uep,nion,n)

      call sub3  (foi,fo(1,iflag),fi,n)
      call cmult (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3  (Jp,p,foi,n)
      call cmult (Jp,dtNTinv,n)      ! Jp = p/dt_newton - (fo-f)/eps
                                     ! case with divided by dt_newton
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
      include 'FDMH1'
      include 'GMRES'

      integer  n,nion,outer,iflag
      real     phi(lpts1),res(lpts1),wt(lpts1)
      real     tol,alpha,l,temp
      real     eps,uc(lpts1,nion),f(lpts1)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            !update residual
            call copy   (r,res,n)                  ! r = res
            call JacobiMatVec(w,x,uc,f,nion,n,iflag)! w = A x
            call add2s2 (r,w,-1.,n) ! r = r - w
         endif

         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }

         tolps = 0.1*gamma(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call JacobiMatVec(w,v(1,j),uc,f,nion,n,iflag) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                              !            ______
            alpha = sqrt(glsc3(w,w,wt,n))     ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,v(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,newton_iter,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end

c---------------------------------------------------------------------
C     This subroutine compute sem_bdf1, user can provide any(*) 
C     bdf1 time swapper in here
      subroutine cem_drift_sem_bdf1_newton(cinput,coutput,nion,n,iflag)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' ! cn
      include 'ZPER'
      include 'BCS'

      integer iflag,nion,n,ic
      real cinput(lpts1,nion)
      real coutput(lpts1,nion)


C     Copy into working array of sem_bdf1
      do ic = 1,lcdim
        call copy(cN(1,ic),cinput(1,ic),n)
      enddo

C     Call sem_bdf1
C     vvvv Put sem_bdf subroutine here
      call cem_drift_op_bdf


C     Copy out working array of sem_bdf1
      do ic = 1,lcdim
        call copy(coutput(1,ic),cN(1,ic),n)
      enddo

      return
      end
c-----------------------------------------------------------------------
