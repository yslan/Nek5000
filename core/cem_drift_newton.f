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
      include 'NEWTON'
      include 'CTIMER' ! nekcem timer, fixme
c     include 'RTIMER'

      parameter (lcdim=ldim) ! number of field, lcdim(vx,vy)=2

      real alphaNT,tt
      real tolNT,tolNTGMRES
      real glmax,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real temp(lx1*ly1*lz1*lelt)

      integer i,j,k,n,ic
      integer maxNTit
      real glsum
      real nfNT_pre,nfNT_now                ! norm of fNT
      real max_dtNT,min_tolNT
      real skNT  (lx1*ly1*lz1*lelt,lcdim),  !
     $     fkNT  (lx1*ly1*lz1*lelt,lcdim),  !
     $     gkNT  (lx1*ly1*lz1*lelt,lcdim),  !
     $     cnNT  (lx1*ly1*lz1*lelt,lcdim),  ! working array in Pseudo time step
     $     ckNT  (lx1*ly1*lz1*lelt,lcdim)   ! working array in Newton iteration
      real dummy1(lx1*ly1*lz1*lelt)        ! dummy field for third dimension

c      common /ctmp1e/ uxe (lx1,ly1,lz1,lelv), uye (lx1,ly1,lz1,lelv)

C     Parametar settings ToDo: move to user file, Lan

      alphaNT = 1.            ! Newton relaxation parameter alphan in (0,1]
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxNTit = 20           ! #iter of Newton method
      jaceps = 1e-8        ! eps for approx. Jacobian
      epsinv = 1./jaceps
      tolNTGMRES = 1e-8     ! tol for GMRES
      tolNT = 1e-5         ! tol for Newton method

      dtNT = 1.!param(1) ! pseudo-transient time step
                      ! it can grow as istep increase 
      max_dtNT  = 1E8   ! Max. of dtNT, avoid NaN
      min_tolNT   = 1E-12 ! Min. of tol,  avoid NaN

      dtinv = 1./dt
      dtNTinv = 1./dtNT

      cpu_t = 0.0
      cpu_chk = 0.0

      npts = lx1*ly1*lz1*nelv

C     Initialization
      do ic = 1,lcdim
        call rzero(skNT(1,ic),npts)
        call rzero(fkNT(1,ic),npts)
        call rzero(gkNT(1,ic),npts)
        call rzero(cnNT(1,ic),npts) ! working array in Pseudo time
        call rzero(cnNT(1,ic),npts) ! working array in NT
      enddo

C     Get initial value
      call copy(ckNT(1,1),vx,npts) ! initial value
      call copy(ckNT(1,2),vy,npts)
      if (ldim.eq.3) call copy(ckNT(1,3),vz,npts)

      call copy(cnNT(1,1),vx,npts)
      call copy(cnNT(1,2),vy,npts)
      if (ldim.eq.3) call copy(cnNT(1,3),vz,npts)

      call compute_fkNT(fkNT,ckNT,lcdim,npts) ! compute f(u_0)
      call compute_gkNT(gkNT,fkNT,ckNT,cnNT,lcdim,npts) ! compute rhs g_k for gmres

c      if (nio.eq.0) then
c      do i=1,lx1*ly1*lz1*nelv
c         write (6,*) 'vx=',fk(i,1)
c         write (6,*) 'vy=',fk(i,2)
c      enddo
c      endif
c      stop

      do ic = 1,lcdim
        call chsign (gkNT(1,ic),npts)
      enddo
      
      do ic = 1,lcdim
        fnorms(ic) = glsc3(fkNT(1,ic),fkNT(1,ic),vmult,npts)
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm = glmax(fnorms,lcdim)
      nfNT_now = fnorm
      call outpost(cnNT(1,1),cnNT(1,2),vz,pr,t,'   ')

C     Start pseudo-time step: update cnNT in each iteration
      do istepNT = 1,nsteps !fixme by lan: istep in Nek5000 is different
         if (nio.eq.0) write (6,*) 'pseudo step',istepNT
         if (mod(istepNT,iostep).eq.0) then
            istep=istepNT
            call outpost(cnNT(1,1),cnNT(1,2),vz,pr,t,'   ')
            istep=1
         endif
         if (nio.eq.0) write (6,*) 'istepnt=',istepNT

c       cpu_dtime = dclock() !FIXME by Lan, NekCEM timer

C       SER, dt, dtNT varying
        if (nio.eq.0) write (6,*) 'nfnt_pre,nfnt_now',nfNT_pre,nfNT_now
        if ((istepNT .eq. 1)) then
          dtNT = 1.0 !param(1) fixme by Lan, initial tau
c          dt = dtNT
        else
          if (nfNT_now.lt.min_tolNT) then
            dtNT = max_dtNT
          else
            dtNT = dtNT * nfNT_pre/nfNT_now
c            dt = dtNT
          endif
        endif
        if (dtNT.gt.max_dtNT) dtNT = max_dtNT
        if (nio.eq.0) write(6,*) istepNT,'dtnt',dtNT
        if (nio.eq.0) write (6,*) 'nfnt_pre,nfnt_now',nfNT_pre,nfNT_now

        dtinv = 1./dt
        dtNTinv = 1./dtNT
        if(nio.eq.0) write (6,*) 'before newton iteration'

C       Start Newton iteration: update c_{k} = c_{k-1} + s_{k}
        do iterNT=1,maxNTit
            if(nio.eq.0) write (6,*) 'before_gmresNT, iternt=',iterNT

          do ic = 1,lcdim ! solving newton sk by GMRES with JacMatVec
            if(nio.eq.0) write (6,*) 'ic=',ic

c           call hmh_gmres_newton4(gk(1,ic),vmult,iter)
            
c           call hmh_gmres_newton3(sk(1,ic),gk(1,ic),vmult,iter,
c    $                             cn_k,fk(1,ic),npts,tolgmres,ic,ldim)
c           call hmh_gmres_newton3(gk(1,ic),vmult,iter)
c           call hmh_gmres_newton2(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),npts,
c    $                             tolgmres,ic,vmult,iter)

c           call hmh_gmres_newton(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),vmult,
c    $                            npts,tolgmres,ic)
c           stop
            if(nio.eq.0) write(6,*) 'tolntgmres',tolNTgmres
            call drift_hmh_gmres_newton(skNT(1,ic),gkNT(1,ic)
     $           ,ckNT,fkNT(1,ic),vmult,ldim,npts,tolNTGMRES,ic)
          enddo
          if(nio.eq.0) write (6,*) 'after__gmresNT=',iterNT

c          write (6,*) 'check 1'
          do ic = 1,lcdim ! update Newton iteration
            call add2s2(ckNT(1,ic),skNT(1,ic),alphaNT,npts) ! cp_k=cp_k+alphan*sp_k
          enddo
c          write (6,*) 'check 2'

c       Update fk, gk, check tol
        call compute_fkNT(fkNT,ckNT,lcdim,npts)
        call compute_gkNT(gkNT,fkNT,ckNT,cnNT,lcdim,npts)  ! checking tol. + next NT iter
c        write (6,*) 'check 3'
        do ic = 1,lcdim
           call chsign(gkNT(1,ic),npts)
        enddo
c        write (6,*) 'check 4'

        do ic = 1,lcdim
           rnorms(ic) = glsc3(gkNT(1,ic),gkNT(1,ic),vmult,npts)
           rnorms(ic) = sqrt(rnorms(ic))
        enddo

c        write (6,*) 'check 5'
        rnorm = glmax(rnorms,lcdim)
        if (nio.eq.0) write(6,*) istepNT,iterNT,'Lan, rnorm',rnorm
        if (iterNT.eq.1) rinorm=rnorm
        ratio=rnorm/rinorm 
c        write (6,*) 'check 6'

c       if ((nid.eq.0).and.(mod(istepn,iocomm).eq.0)) 
        if (nio.eq.0) write(6,90) istepNT,iterNT,
     $      ratio,rnorm,rinorm,dtNT,dt
c        write (6,*) 'check 7'

        if (ratio.lt.tolNT) goto 900 ! Newton conv
c       if (rnorm.lt.new_tol) goto 900 ! Newton conv

        enddo ! END of Newton Loop
 90     format(2i6,'newton iter',1p5e12.4)   
900     continue

        do ic = 1,lcdim
           call copy(cnNT(1,ic),ckNT(1,ic),npts)   ! update cn_n for computing gn
        enddo                                      ! and next pseudo-time step

c       Compute the norm of f
        do ic = 1,lcdim
           fnorms(ic) = glsc3(fkNT(1,ic),fkNT(1,ic),vmult,npts)  
           fnorms(ic) = sqrt(fnorms(ic))
           if (nio.eq.0) write(6,*)'Lan, ic, fnorms',ic,fnorms(ic)
        enddo

        fnorm  = glmax(fnorms,lcdim)
        nfNT_pre = nfNT_now  ! set up old fnorm to f_pre
        nfNT_now = fnorm  ! set up new fnorm to f_now
        if (nio.eq.0) write(6,*)'Lan, fnorm nfnt_pre,nfnt_now',fnorm,
     $   nfNT_pre,nfNT_now
        

c       Compute the CPU_time
c       cpu_dtime = dclock()-cpu_dtime !FIXME by Lan, NekCEM timer
        cpu_t = cpu_t+cpu_dtime
        cpu_t_step = cpu_t/istepnt
        cpu_p_t = glsum(cpu_t_step /npts,1)/np

        time = time + dtNT ! update physical timing

c       call errchk(cn_k(1,1),cn_k(1,2))

C       Copy for cem_out and cem_chk
c       do ic = 1,lcdim
c         call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
c       enddo
        if (ldim.eq.3) call opcopy(vx,vy,vz,
     $                           cnNT(1,1),cnNT(1,2),cnNT(1,3))
        if (ldim.eq.2) call opcopy(vx,vy,vz,
     $                           cnNT(1,1),cnNT(1,2),dummy1)

c       call userchk
       if(nio.eq.0) write(6,*)'before errchk'
       call errchk(vx(1,1,1,1),vy(1,1,1,1)) ! in user file
       if(nio.eq.0) write(6,*)'after  errchk'
c       call cem_out
       if(nio.eq.0) write(6,*)'before outpost' ! fixme, output in Nek5000
c       call outpost(vx,vy,vz,pr,t,'   ')
       if(nio.eq.0) write(6,*)'after  outpost'

c       call errorchk(cn_k(1,1),cn_k(1,2))

      enddo ! end of pseudo-time step Loop

c     call cem_end !FIXME end of everything

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
      subroutine compute_gkNT(gout,f,ckn,c0n,nion,n)
c-----------------------------------------------------------------------      
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'
      integer nion,n,ic
      real gout(lx1*ly1*lz1*lelt,nion)
     $   , ckn (lx1*ly1*lz1*lelt,nion)
     $   , c0n (lx1*ly1*lz1*lelt,nion) ! saved ckn from last pseudo-time step
     $   , f   (lx1*ly1*lz1*lelt,nion) ! precomputed by compute_fkNT

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
      subroutine compute_fkNT(fout,cin,nion,n)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      real fout(lx1*ly1*lz1*lelt,nion) ! f output (of subroutine)
      real cin (lx1*ly1*lz1*lelt,nion) ! c input
      real cout(lx1*ly1*lz1*lelt,nion) ! c output (from sem_bdf)

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
c     is any vector. Instead of constructing the Jacobian matrix, we
c     aproximate Jp by the following formula
c   
c        J_k s_k = s_k - (dt /eps) *  ( f(u_k + eps*s_k) - f(u_k) )
c   
c     where f(u_k) has been store in each Newton iteration
c           f(u_k + eps*s_k) has to be computed in each GMRES iteration
      subroutine JacobiMatVec(jacp,p,uc,fi,nion,n,iflag)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      real     jacp(lx1*ly1*lz1*lelt),p(lx1*ly1*lz1*lelt)
      real     uc (lx1*ly1*lz1*lelt,nion)     ! treat other ck as fix point
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

      call compute_fkNT(fo,uep,nion,n)

      call sub3  (foi,fo(1,iflag),fi,n)
      call cmult (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3  (jacp,p,foi,n)
      call cmult (jacp,dtntinv,n)      ! Jp = p/dt_newton - (fo-f)/eps
                                     ! case with divided by dt_newton
c      call copy(jacp,p,n) ! for testing identity matrix

      return
      end
c-----------------------------------------------------------------------
C     GMRES that solving Newton iteration. Use JacMacVec as Ax
C     This solve linear system for single array, but the other arrays
C     are also invloved as a referenced points
C     from NekCEM (and NekCEM from Nek5000)
      subroutine drift_hmh_gmres_newton
     $           (phi,res,uc,f,wt,nion,n,tol,ic)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
c     include 'FDMH1'
c     include 'GMRES1'
      include 'NTGMRES' ! fixme
      include 'NEWTON'

      integer  n,nion,outer,ic
      real     phi(lx1*ly1*lz1*lelt),res(lx1*ly1*lz1*lelt),
     $         wt(lx1*ly1*lz1*lelt)
      real     tol,alphant,l,temp
      real     eps,uc(lx1*ly1*lz1*lelt,nion),f(lx1*ly1*lz1*lelt)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x_ntgmres,n)
      call rzero(h_ntgmres,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r_ntgmres,res,n)                  ! r = res
         else
            !update residual
            call copy   (r_ntgmres,res,n)                  ! r = res
            call JacobiMatVec(w_ntgmres,x_ntgmres,uc,f,nion,n,ic)! w = A x
            call add2s2 (r_ntgmres,w_ntgmres,-1.,n) ! r = r - w
         endif

         gamma_ntgmres(1) = glsc3(r_ntgmres,r_ntgmres,wt,n)          ! gamma_gmres  = (r,r)
         gamma_ntgmres(1) = sqrt(gamma_ntgmres(1))           ! gamma_gmres  = sqrt{ (r,r) }

         tolps = 0.1*gamma_ntgmres(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable
         tolpsss = tolpss

         write (6,*) 'tolps,tolpss',tolps,tolpss

         !check for lucky convergence
         rnorm = 0.
         if(gamma_ntgmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_ntgmres(1)
         call cmult2(v_ntgmres(1,1),r_ntgmres,temp,n)             !  v  = r / gamma_gmres
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call JacobiMatVec(w_ntgmres,v_ntgmres(1,j),uc,f,nion,n,ic) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h_ntgmres(i,j)=glsc3(w_ntgmres,v_ntgmres(1,i),wt,n)     ! h    = (w,v )
                                                  ! i,j       i
               call add2s2(w_ntgmres,v_ntgmres(1,i),-h_ntgmres(i,j),n) ! w = w - h    v
            enddo                                 !         i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_ntgmres(i,j)
               h_ntgmres(i  ,j)=  c_ntgmres(i)*temp 
     $         + s_ntgmres(i)*h_ntgmres(i+1,j)
               h_ntgmres(i+1,j)= -s_ntgmres(i)*temp 
     $         + c_ntgmres(i)*h_ntgmres(i+1,j)
            enddo
                                              !            ______
            alphant = sqrt(glsc3(w_ntgmres,w_ntgmres,wt,n))     ! alphan =  \/ (w,w)
            if(alphant.eq.0.) goto 900 !converged
            l = sqrt(h_ntgmres(j,j)*h_ntgmres(j,j)+alphant*alphant)
            temp = 1./l
            c_ntgmres(j) = h_ntgmres(j,j) * temp
            s_ntgmres(j) = alphant  * temp
            h_ntgmres(j,j) = l
            gamma_ntgmres(j+1) = -s_ntgmres(j) * gamma_ntgmres(j)
            gamma_ntgmres(j)   =  c_ntgmres(j) * gamma_ntgmres(j)

            rnorm = abs(gamma_ntgmres(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c                 write (6,*) rnorm,'newton'
                  write (6,66) iter,tolpsss,rnorm,istepnt,iternt
   66       format(i5,1p2e12.5,2i8,' gmres_newton rnorm')

            if (rnorm .lt. tolpsss) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alphant
            call cmult2(v_ntgmres(1,j+1),w_ntgmres,temp,n)   ! v    = w / alphan
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
            temp = gamma_ntgmres(k)
            do i=j,k+1,-1
               temp = temp - h_ntgmres(k,i)*c_ntgmres(i)
            enddo
            c_ntgmres(k) = temp/h_ntgmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_ntgmres,v_ntgmres(1,i),c_ntgmres(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x_ntgmres,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

c     if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
      if (nio.eq.0)  write(6,9999) istepNT,iterNT,iter,rnorm,tolpsss
c     endif

 9999 format(' ',' ',i9,i6,'  gmresNT_iteration#',i6,2(1p1e12.4))

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

C     Copy into working arrays of sem_bdf1

      write (6,*) 'check pre'

      if (ldim.eq.3) call opcopy(vx,vy,vz,
     $                           cinput(1,1),cinput(1,2),cinput(1,3))

      if (ldim.eq.2) call opcopy(vx,vy,vz,
     $                           cinput(1,1),cinput(1,2),dummy1)

C     *** main bdf1 solver ***
      call plan5(2)


C     Copy out working arrays from sem_bdf1

      if (ldim.eq.2)
     $      call opcopy(coutput(1,1),coutput(1,2),dummy1,vx,vy,vz)
      if (ldim.eq.3) 
     $      call opcopy(coutput(1,1),coutput(1,2),coutput(1,3),vx,vy,vz)

      write (6,*) 'check post'

      return
      end
c-----------------------------------------------------------------------
cc     subroutine hmh_gmres_newton(res,h1,h2,wt,iter)
c      subroutine hmh_gmres_newton(out,res,uc,f,wt,n,tol,isd)
c
cc           call hmh_gmres_newton(sk(1,ic),gk(1,ic),cn_k,fk(1,ic),vmult,
cc    $                            lcdim,npts,tolgmres,ic)
c
cc     Solve the Helmholtz equation by right-preconditioned 
cc     GMRES iteration.
c
c      include 'SIZE'
c      include 'TOTAL'
c      include 'FDMH1'
cc      include 'NTGMRES'
c      include 'NEWTON' !istepNT
c
c      common  /ctolpr/ divex
c      common  /cprint/ ifprint
c      logical          ifprint
c
c      real             res  (lx1*ly1*lz1*lelv)
cc     real             h1   (lx1,ly1,lz1,lelv)
cc     real             h2   (lx1,ly1,lz1,lelv)
c      real             wt   (lx1,ly1,lz1,lelv)
c
c      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)
c
c      common /cgmres1/ y(lgmres)
c      common /ctmp0/   wk1(lgmres),wk2(lgmres)
c
c      real alpha, l, temp
c      integer outer
c
c      logical iflag,if_hyb
c      save    iflag,if_hyb
cc     data    iflag,if_hyb  /.false. , .true. /
c      data    iflag,if_hyb  /.false. , .false. /
c      real    norm_fac
c      save    norm_fac
c
c      real*8 etime1,dnekclock
c
c      n = nx1*ny1*nz1*nelv
c
c      etime1  = dnekclock()
c      etime_p = 0.
c      divex   = 0.
c      iter    = 0
c      m       = lgmres
c
c      write (6,*) 'gm 1'
c
c      if (.not.iflag) then
c         iflag=.true.
c         call uzawa_gmres_split(ml_gmres,mu_gmres,bm1,binvm1,
c     $                          nx1*ny1*nz1*nelv)
c         norm_fac = 1./sqrt(volvm1)
c      endif
c
c      write (6,*) 'gm 2'
c
cc     if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
c      do i=1,lx1*ly1*lz1*nelv
c         d(i) = 1.
c      enddo
c
c      write (6,*) 'gm 3'
c
cc     call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
c
cc     if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
cc    $   tolps = abs(param(21))
cc     if (istep.eq.0) tolps = 1.e-4
c      tolps=tol
c      tolpss = tolps
c
c      iconv = 0
c      call rzero(x_gmres,n)
c
c      outer = 0
c      do while (iconv.eq.0.and.iter.lt.500)
c         outer = outer+1
c
c         if (iter.eq.0) then                   !      -1
c            call col3(r_gmres,ml_gmres,res,n) ! r = L  res
cc           call copy(r,res,n)
c         else
c            !update residual
c            call copy  (r_gmres,res,n)           ! r = res
c
c            call jacobimatvec(w_gmres,x_gmres,uc,f,ldim,n,isd)
c            call copy(w_gmres,x_gmres,n)
cc           call ax    (w_gmres,x_gmres,h1,h2,n) ! w = A x
c            call add2s2(r_gmres,w_gmres,-1.,n)   ! r = r - w
c                                                 !      -1
c            call col2(r_gmres,ml_gmres,n)        ! r = L   r
c         endif
c                                                            !            ______
c         gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wt,n)) ! gamma  = \/ (r,r) 
c                                                            !      1
c         if(iter.eq.0) then
c            div0 = gamma_gmres(1)*norm_fac
c            if (param(21).lt.0) tolpss=abs(param(21))*div0
c         endif
c
c         !check for lucky convergence
c         rnorm = 0.
c         if (gamma_gmres(1) .eq. 0.) goto 9000
c         temp = 1./gamma_gmres(1)
c         call cmult2(v_gmres(1,1),r_gmres,temp,n) ! v  = r / gamma
c                                                  !  1            1
c
c         ifmgrid=.false.
c
c         do j=1,m
c            iter = iter+1
c                                                       !       -1
c            call col3(w_gmres,mu_gmres,v_gmres(1,j),n) ! w  = U   v
c                                                       !           j
c
cc . . . . . Overlapping Schwarz + coarse-grid . . . . . . .
c
c            etime2 = dnekclock()
c
cc           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergence
c            if (ifmgrid) then
c               call h1mg_solve(z_gmres(1,j),w_gmres,if_hyb) ! z  = M   w
c            else                                            !  j
c               kfldfdm = ndim+1
c               if (param(100).eq.2) then
c                   call h1_overlap_2 (z_gmres(1,j),w_gmres,pmask)
c               else
c                   call fdm_h1
c     $               (z_gmres(1,j),w_gmres,d,pmask,vmult,nelv,
c     $                ktype(1,1,kfldfdm),wk)
c               endif
c               call crs_solve_h1 (wk,w_gmres)        ! z  = M   w
c               call add2         (z_gmres(1,j),wk,n) !  j        
c            endif
c
c            call ortho        (z_gmres(1,j)) ! Orthogonalize wrt null space, if present
c            etime_p = etime_p + dnekclock()-etime2
cc . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
c
c            call jacobimatvec(w_gmres,z_gmres(1,j),uc,f,ldim,n,isd)
c            call copy(w_gmres,z_gmres(1,j),n)
c
cc           call ax  (w_gmres,z_gmres(1,j),h1,h2,n) ! w = A z
c                                                    !        j
c     
c                                                    !      -1
c            call col2(w_gmres,ml_gmres,n)           ! w = L   w
c
cc           !modified Gram-Schmidt
c
cc           do i=1,j
cc              h_gmres(i,j)=glsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
cc                                                            !  i,j       i
c
cc              call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
cc           enddo                                                !          i,j  i
c
cc           2-PASS GS, 1st pass:
c
c            do i=1,j
c               h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
c            enddo                                            !  i,j       i
c
c            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs
c
c            do i=1,j
c               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
c            enddo                                                !          i,j  i
c
c
cc           2-PASS GS, 2nd pass:
cc
cc           do i=1,j
cc              wk1(i)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
cc           enddo                                      !  i,j       i
cc                                                      !
cc           call gop(wk1,wk2,'+  ',j)                  ! sum over P procs
cc
cc           do i=1,j
cc              call add2s2(w_gmres,v_gmres(1,i),-wk1(i),n)    ! w = w - h    v
cc              h_gmres(i,j) = h_gmres(i,j) + wk1(i)           !          i,j  i
cc           enddo
c
c            !apply Givens rotations to new column
c            do i=1,j-1
c               temp = h_gmres(i,j)                   
c               h_gmres(i  ,j)=  c_gmres(i)*temp 
c     $                        + s_gmres(i)*h_gmres(i+1,j)  
c               h_gmres(i+1,j)= -s_gmres(i)*temp 
c     $                        + c_gmres(i)*h_gmres(i+1,j)
c            enddo
c                                                      !            ______
c            alpha = sqrt(glsc3(w_gmres,w_gmres,wt,n)) ! alpha =  \/ (w,w)
c            rnorm = 0.
c            if(alpha.eq.0.) goto 900  !converged
c            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
c            temp = 1./l
c            c_gmres(j) = h_gmres(j,j) * temp
c            s_gmres(j) = alpha  * temp
c            h_gmres(j,j) = l
c            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
c            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)
c
c            rnorm = abs(gamma_gmres(j+1))*norm_fac
c            ratio = rnorm/div0
c            if (ifprint.and.nio.eq.0) 
c     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istepnt
c   66       format(i5,1p4e12.5,i8,' Divergence newton')
c
c#ifndef TST_WSCAL
c            if (rnorm .lt. tolpss) goto 900  !converged
c#else
c            if (iter.gt.param(151)-1) goto 900
c#endif
c            if (j.eq.m) goto 1000 !not converged, restart
c
c            temp = 1./alpha
c            call cmult2(v_gmres(1,j+1),w_gmres,temp,n) ! v    = w / alpha
c                                                       !  j+1            
c         enddo
c  900    iconv = 1
c 1000    continue
c         !back substitution
c         !     -1
c         !c = H   gamma
c         do k=j,1,-1
c            temp = gamma_gmres(k)
c            do i=j,k+1,-1
c               temp = temp - h_gmres(k,i)*c_gmres(i)
c            enddo
c            c_gmres(k) = temp/h_gmres(k,k)
c         enddo
c         !sum up Arnoldi vectors
c         do i=1,j
c            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),n) ! x = x + c  z
c         enddo                                             !          i  i
cc        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
c      enddo
c 9000 continue
c
c      divex = rnorm
c      call copy(out,x_gmres,n)
c
c      call ortho   (res) ! Orthogonalize wrt null space, if present
c
c      etime1 = dnekclock()-etime1
c      if (nio.eq.0) write(6,9999) istepnt,iter,divex,div0,tolpss,
c     $                            etime_p,etime1,if_hyb
cc     call flush_hack
c 9999 format(4x,i7,'  newton gmres ',4x,i5,1p5e13.4,1x,l4)
c
c      if (outer.le.2) if_hyb = .false.
c
c      return
c      end
cc-----------------------------------------------------------------------
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
