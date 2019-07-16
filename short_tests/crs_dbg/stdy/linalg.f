c=======================================================================
c
c     LIBRARY ROUTINES FOR LINEAR ALGEBRA
c
c     June 2018
c
c     For questions, comments or suggestions, please contact:
c
c     Pablo Daniel Brubeck
c     brubeck@protonmail.com
c     
c-----------------------------------------------------------------------
      subroutine display(m,n,A,lda)
      implicit none
      integer m,n,lda
      real A(lda,*)

      integer i,j

      do i=1,m
         write(*,*) ( A(i,j), j=1,n )
      enddo ! i
      return
      end
c-----------------------------------------------------------------------
      subroutine savemat(io,m,n,A,lda)
      implicit none
      integer io
      integer m,n,lda
      real A(lda,*)

      integer i,j

      do i=1,m
         write(io,*) ( A(i,j), j=1,n )
      enddo ! i
      return
      end
c-----------------------------------------------------------------------
      subroutine swapi(a,b)
      integer a,b,temp
      temp = a
      a = b
      b = temp
      return
      end
c-----------------------------------------------------------------------
      subroutine set_ident(eye,ldi,n) ! move from gfdm_op.f, removed in v19-rc1
      real eye(ldi,1)
      do j=1,n
         do i=1,n
            eye(i,j) = 0.
         enddo
         eye(j,j) = 1.
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine xasolv(m,n,X,A)
      implicit none
      integer m,n
      real X(*), A(*)

      real one
      parameter (one=1.0E0)

      call dtrsm('R','L','T','N',m,n,one,A,n,X,m)
      call dtrsm('R','L','N','N',m,n,one,A,n,X,m)
      return
      end
c-----------------------------------------------------------------------
      subroutine normc(m,n,A,nrm)
      implicit none
      integer m,n
      real A(m,*), nrm(*)
      real sca
      integer j

      real zero,one,eps
      parameter (zero=0.0E0, one=1.0E0, eps=1E-14)

      real dnrm2

      do j=1,n
         nrm(j) = dnrm2(m,A(1,j),1)
         if(nrm(j).lt.eps) then
            nrm(j) = one
         else
            sca = one/nrm(j)
            call dscal(m,sca,A(1,j),1)
         endif
      enddo ! j
      return
      end
c-----------------------------------------------------------------------
      subroutine dkron(m1,m2, n1,n2, alpha,A,B, beta,C)
      implicit none
      integer m1,m2, n1,n2
      real alpha, beta
      real A(m1,*), B(m2,*), C(m2,m1, n2,*)
      integer i,j, p,q

      do q=1,n1
      do j=1,n2
      do p=1,m1
      do i=1,m2
         C(i,p,j,q) = beta*C(i,p, j,q) + alpha*A(p,q)*B(i,j)
      enddo ! i
      enddo ! p
      enddo ! j
      enddo ! q
      return
      end
c-----------------------------------------------------------------------
      subroutine dkron3(m1,m2,m3, n1,n2,n3, alpha,A,B,C, beta,D)
      implicit none
      integer m1,m2,m3, n1,n2,n3
      real alpha, beta
      real A(m1,*), B(m2,*), C(m3,*), D(m3,m2,m1, n3,n2,*)
      integer i,j,k, l,m,n

      do n=1,n1
      do m=1,n2
      do l=1,n3
      do k=1,m1
      do j=1,m2
      do i=1,m3
         D(i,j,k,l,m,n) = beta*D(i,j,k,l,m,n)+alpha*A(k,n)*B(j,m)*C(i,l)
      enddo ! i
      enddo ! j
      enddo ! k
      enddo ! l
      enddo ! m
      enddo ! n
      return
      end
c-----------------------------------------------------------------------
      subroutine diag(m,n,A,lda,d)
      implicit none
      integer m,n,lda
      real A(lda,*), d(*)
      integer i

      do i=1,m
         d(i) = A(i,i)
      enddo ! i
      return
      end
c-----------------------------------------------------------------------
      subroutine diagxpy(m,alpha,x,Y,ldy)
      implicit none
      integer m,ldy
      real alpha, x(*), Y(ldy,*)
      integer i

      do i=1,m
         Y(i,i) = Y(i,i) + alpha*x(i)
      enddo ! i
      return 
      end
c-----------------------------------------------------------------------
      subroutine dxm(m,n,d,A, alpha,B)
      implicit none
      integer m,n
      real alpha
      real d(*), A(m,*), B(m,*)
      
      integer i,j
      do j=1,n
      do i=1,m
         B(i,j) = alpha*B(i,j) + d(i)*A(i,j) 
      enddo ! i
      enddo ! j
      return
      end
c-----------------------------------------------------------------------
      subroutine dotdxm(m,n,d1,d2,d3,A,B)
      implicit none
      integer m,n
      real d1(*), d2(*), d3(*), A(m,n,*), B(m,*)

      integer i,j
      do j=1,n
      do i=1,m
         B(i,j) = d1(i)*A(i,j,1)+d2(i)*A(i,j,2)+d3(i)*A(i,j,3)
      enddo ! i
      enddo ! j
      return
      end
c-----------------------------------------------------------------------
      subroutine mxmc(A,m,B,k,C,n)
c     C = A*B
      implicit none
      integer m,k,n
      complex*16 A(m,*),B(k,*),C(m,*)
      complex*16 one,zero
      parameter (one=(1.0E0,0.0E0), zero=(0.0E0,0.0E0))
      
      call zgemm('N','N',m,n,k,one,A,m,B,k,zero,C,m)
      return
      end
c-----------------------------------------------------------------------
      subroutine kron_mxv1(m1, n1,n2,n3, A, X, Y)
c     Y = kron(I,I,A)*X
      implicit none
      integer m1,n1,n2,n3
      real A(m1,*), X(n1,n2,*), Y(m1,n2,*)

      call mxm(A,m1,X,n1,Y,n2*n3)
      return
      end
c-----------------------------------------------------------------------
      subroutine kron_mxv2(m2, n1,n2,n3, At, X, Y)
c     Y = kron(I,A,I)*X
      implicit none
      integer m2,n1,n2,n3
      real At(n2,*), X(n1,n2,*), Y(n1,m2,*)
      integer k

      do k=1,n3
         call mxm(X(1,1,k),n1,At,n2,Y(1,1,k),m2)
      enddo ! k
      return
      end
c-----------------------------------------------------------------------
      subroutine kron_mxv3(m3, n1,n2,n3, At, X, Y)
c     Y = kron(A,I,I)*X
      implicit none
      integer m3,n1,n2,n3
      real At(n3,*), X(n1,n2,*), Y(n1,n2,*)

      call mxm(X,n1*n2,At,n3,Y,m3)
      return
      end
c-----------------------------------------------------------------------
      subroutine symrank1(ndim,n,C1,C2,C3,C4,C5,C6,a,b1,b2,b3)
c     Symmetric rank-1 update to tensor, C = C + a *( b*b' )
      implicit none
      integer ndim, n
      real C1(*), C2(*), C3(*), C4(*), C5(*), C6(*)
      real  a(*), b1(*), b2(*), b3(*)
      call addcol4(C1, a,b1,b1, n)
      call addcol4(C2, a,b2,b2, n)
      call addcol4(C4, a,b1,b2, n)
      if (ndim.gt.2) then
      call addcol4(C3, a,b3,b3, n)
      call addcol4(C5, a,b1,b3, n)
      call addcol4(C6, a,b2,b3, n)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine diag_innerp(m,n,A,lda,B,ldb,C,ldc,x,WORK)
c     x = diag(A'*B*C) = sum(A.*(B*C),1)
c     x : n x 1
c     A : m x n
c     B : m x m
c     C : m x n

      implicit none
      integer m,n,lda,ldb,ldc
      real A(lda,*), B(ldb,*), C(ldc,*), x(*), WORK(m,*)

      integer j
      real one, zero
      parameter (one=1.0E0, zero=0.0E0)
      real ddot

c     WORK = B*C
      call dgemm('N','N',m,n,m,one,B,ldb,C,ldc,zero,WORK,m)

c     x = sum(A.*WORK,1)
      do j=1,n
         x(j) = ddot(m,WORK(1,j),1,A(1,j),1)
      enddo ! j
      return
      end
c-----------------------------------------------------------------------
      subroutine diag_outerp(m,n,A,lda,B,ldb,C,ldc,x,WORK) 
c     x = diag(A*B*C') = sum((A*B).*C,2)
c     x : m x 1
c     A : m x n
c     B : n x n
c     C : m x n

      implicit none
      integer m,n,lda,ldb,ldc
      real A(lda,*), B(ldb,*), C(ldc,*), x(*), WORK(m,*)

      integer i

      real one, zero
      parameter (one=1.0E0, zero=0.0E0)
      real ddot

c     WORK = A*B
      call dgemm('N','N',m,n,n,one,A,lda,B,ldb,zero,WORK,m)

c     v = sum(WORK.*C,2)
      do i=1,m
         x(i) = ddot(n,WORK(i,1),m,C(i,1),ldc)
      enddo ! i
      return
      end
c-----------------------------------------------------------------------
      subroutine innerp_diag(m,n,alpha,A,lda,b,C,ldc,beta,X,ldx,WORK)
c     X = alpha*A'*diag(b)*C + beta*X
c     X : n x n
c     A : m x n
c     b : m x 1
c     C : m x n

      implicit none
      integer m,n,lda,ldc,ldx
      real alpha,beta,A(lda,*),b(*),C(ldc,*),X(ldx,*),WORK(m,*)

      integer i,j
      
c     WORK = diag(b)*A
      do j=1,n
         do i=1,m
            WORK(i,j)=b(i)*C(i,j)
         enddo ! i 
      enddo ! j

c     X = alpha*A'*WORK + beta*X
c     FIXME Transpose flag
      call dgemm('T','N',n,n,m,alpha,A,lda,WORK,m,beta,X,ldx)
      return
      end
c-----------------------------------------------------------------------
      subroutine schwarz1d(n,p,nmat,A,lda,S,lds,bcs)
c     1D assembly for p x p, n x n, p x p diagonal blocks of A
      implicit none
      integer n,p,nmat
      integer lda,lds
      real A(lda,lda,*)
      real S(lds,lds,*)
      integer bcs(*)

      integer i,k, m,pa,ps

      ! TODO move this subroutine to tens.f

c     Initialize to the identity    
      call rzero(S,lds*lds*nmat)
      do k=1,min(nmat,2)
      do i=1,lds
         S(i,i,k) = 1.0E0
      enddo ! i
      enddo ! k

      do k=1,nmat
         m = n
         pa = p+1
         ps = max(p,1)
         if(bcs(2).eq.1) m = m-1
         if(bcs(1).eq.1) then
            pa = pa+1
            ps = ps+1
            m = m-1
         endif
         call dlacpy('A',m,m,A(pa,pa,k),lda,S(ps,ps,k),lds)
         
         if((bcs(1).eq.0).and.(p.gt.0)) then
            if(p.gt.1)call dlacpy('A',p,p,A(1,1,k),lda,S(1,1,k),lds)
            S(ps,ps,k) = A(pa-1,pa-1,k) + A(pa,pa,k)
         endif

         if((bcs(2).eq.0).and.(p.gt.0)) then
            ps = n+p-1
            pa = n+p+1
            if(p.gt.1)call dlacpy('A',p,p,A(pa,pa,k),lda,S(ps,ps,k),lds)
            S(ps,ps,k) = A(pa-1,pa-1,k) + A(pa,pa,k)
         endif
      enddo ! k
      return
      end
c-----------------------------------------------------------------------
      subroutine imass(m,A,ipiv)
      implicit none
      integer m
      real A(m,m,3)
      integer ipiv(*)

      integer info
      external dgetrf,dgetrs

c     M = lu(M)
      call dgetrf(m,m,A(1,1,2),m,ipiv,info)
      if (info.ne.0) then
         write(6,*)'dgetrf broken, info =',info
         call exitt
      endif

c     K = M\K
      call dgetrs('N',m,m,A(1,1,2),m,ipiv,A(1,1,1),m,info)
      if (info.ne.0) then
         write(6,*)'dgetrs broken, info =',info
         call exitt
      endif

c     L = M\L
      call dgetrs('N',m,m,A(1,1,2),m,ipiv,A(1,1,3),m,info)
      if (info.ne.0) then
         write(6,*)'dgetrs broken, info =',info
         call exitt
      endif
      return
      end
c-----------------------------------------------------------------------
      function eigsort(alphar,alphai,beta)
      implicit none
      logical eigsort
      real alphar, alphai, beta
      real zero
      parameter (zero=0.0E0)
      eigsort = (alphai.ne.zero)
      return
      end
c-----------------------------------------------------------------------
      subroutine schbal(TRANS,m,S,Q,Z,alphar,alphai)
c     Balances Schur decomposition
c     Applies similarity transformation to S such that the 2-by-2
c     diagonal blocks have the form [a, -b; b, a]
c     
c     S <- P\S*P
c     Q <- Q*P
c     Z <- Z/P
c     Retains Z'*A*Q = S, Z'*B*Q = I, but Q'*Q = P^2
      implicit none
      character*1 TRANS
      integer m
      real S(m,m), Q(m,m), Z(m,m), alphar(m), alphai(m)

      real a,b
      integer k  

      real zero
      parameter (zero=0.0E0)

      k = m
      do while (k.gt.0)
      if (alphai(k).eq.zero) then
         k = k-1
      else
         a = alphai(k)/S(k,k-1)
         b = S(k,k-1)/alphai(k)
         if (TRANS.eq.'N') then
         call dscal(m,a,Q(1,k-1),1) ! Q =Q*D
         call dscal(m,a,S(1,k-1),1) ! S =S*D
         call dscal(m,b,S(k-1,1),m) ! S =D\S
         call dscal(m,b,Z(k-1,1),m) ! Z'=D\Z'
         elseif (TRANS.eq.'T') then
         call dscal(m,a,Q(k-1,1),m) ! Q'=D*Q'
         call dscal(m,a,S(1,k-1),1) ! S =S*D
         call dscal(m,b,S(k-1,1),m) ! S =D\S
         call dscal(m,b,Z(1,k-1),1) ! Z =Z/D'
         endif
         k = k-2
         endif
      enddo ! k
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact1(TRANS,p,m,S,Q,Z,wr,wi,B,ipiv,WORK,BWORK)
c     Solves Z'*A*Q = S,  Z'*B*Q = I, with Q'*Q = I
c     
c     TRANS : (on entry) 'T' or 'N'
c     p : (on entry) number of blocks of A
c     m : (on entry) size of each block of A
c     S : (on entry) B\A, (on exit) Schur upper quasi-triangular form
c     Q : (on exit) Orthonormal right Schur vectors, transposed if TRANS='T' 
c     Z : (on exit) Left Schur vectors, transposed if TRANS='N'
c     alphar, alphai : (on exit) real and imaginary parts of eigenvalues
c     B : (on entry) LU decomp of B as returned by dgetrf
c     ipiv : (on entry) pivots of B as returned by dgetrf
c     WORK : size max(8,p*m)*(p*m)

      implicit none
      character*1 TRANS
      integer p,m
      real S(*),Q(m,p,m,p),Z(m,p,m,p)
      real wr(*),wi(*)
      real B(*)

      integer ipiv(*)  ! m
      real    WORK(*)  ! max(p*m,8)**2
      logical BWORK(*) ! p*m

      integer n,i,j, sdim,lw,info

      real zero
      parameter (zero=0.0E0)

      external dgetrs, dgees
      logical  eigsort
      external eigsort
      
      n = m*p
      lw = max(n,8)**2

c     S = Q*S*Q'
      call dgees('V','S',eigsort,n,S,n,sdim,wr,wi,Q,n,
     $                   WORK,lw,BWORK,info)
      if (info.ne.0) then
         write(6,*) 'dgees broken, info =',info
         call exitt
      endif


c     Z = B'\Q
      call dcopy(n*n,Q,1,Z,1)
      do i=1,p
      call dgetrs('T',m,n,B,m,ipiv,Z(1,i,1,1),n,info)
      if (info.ne.0) then
         write(6,*)'dgetrs broken, info =',info
         call exitt
      endif
      enddo ! i
      
      if (TRANS.eq.'N') then
c     Q = Q
c     Z = Z'
      call dcopy(n*n,Z,1,WORK,1)
      call transpose(Z,n,WORK,n)
      elseif (TRANS.eq.'T') then
c     Q = Q'
c     Z = Z
      call dcopy(n*n,Q,1,WORK,1)
      call transpose(Q,n,WORK,n)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact22(m1,m2,S1,S2,Q1,Q2,Z1,Z2,A,B,
     $                     WORK,IWORK,BWORK)
      implicit none
      integer m1,m2
      real S1(*),S2(*),Q1(*),Q2(*),Z1(*),Z2(*)
      real A(m1,m1,*),B(m2,m2,*)
      
      real     WORK(*) ! max(8,m1+2,m2+2)**2
      integer IWORK(*) ! max(m1,m2)
      logical BWORK(*) ! max(m1,m2)

      integer p1,p2,pw
      
      p1 = 1
      p2 = p1 + max(m1,m2)
      pw = p2 + max(m1,m2)

      call imass(m1,A,IWORK)
      call dcopy(m1*m1,A(1,1,1),1,S1,1)
      call schfact1('N',1,m1,S1,Q1,Z1,WORK(p1),WORK(p2),
     $     A(1,1,2),IWORK,WORK(pw),BWORK)
      
      call imass(m2,B,IWORK)
      call dcopy(m2*m2,B(1,1,1),1,S2,1)
      call schfact1('T',1,m2,S2,Q2,Z2,WORK(p1),WORK(p2),
     $     B(1,1,2),IWORK,WORK(pw),BWORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine schsolv22(m1,m2,X,S1,S2,Q1,Q2,Z1,Z2,WORK)
      implicit none
      integer m1,m2
      real X(*),S1(*),S2(*),Q1(*),Q2(*),Z1(*),Z2(*),WORK(*)! WORK(m1,m2)

      real scal
      integer info
      integer isgn
      parameter (isgn=1)
      external dtrsyl

c     X(:) = kron(Z2,Z1)'*X(:)
      call kron_mxv1(m1, m1,m2,1, Z1, X, WORK)
      call kron_mxv2(m2, m1,m2,1, Z2, WORK, X)

c     X(:) = (kron(I,S1) + kron(S2,I))\X(:)
      call dtrsyl('N','T',isgn,m1,m2,S1,m1,S2,m2,X,m1,scal,info)
      if (info.ne.0) then
         write(6,*)'dtrsyl broken, info =',info
         call exitt
      endif

c     X(:) = kron(Q2,Q1)*X(:)
      call kron_mxv1(m1, m1,m2,1, Q1, X, WORK)
      call kron_mxv2(m2, m1,m2,1, Q2, WORK, X)

      call dscal(m1*m2,1.0E0/scal,X,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact222(m1,m2,m3,S,Q,Z,ls,A,B,C,WORK,IWORK,BWORK)
      implicit none
      integer m1,m2,m3
      real S(*), Q(*), Z(*) ! S(ls), Q(ls), Z(ls)
      integer ls
      real A(m1,m1,*), B(m2,m2,*), C(m3,m3,*)

      real     WORK(*) ! 4*max(m1+1,m2+1,m3+1,4)**2
      integer IWORK(*) ! m1+m2+m3
      logical BWORK(*) ! max(2*m1,2*m2,m3)
      
c     TODO Proper implementation
      call dcopy(m2*m2,B(1,1,2),1,B(1,1,3),1)
      call schfact232(m1,m2,m3,S,Q,Z,ls,A,B,C,WORK,IWORK,BWORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine schsolv222(m1,m2,m3,X,S,Q,Z,WORK)
c     X (on entry) RHS, (on exit) X
      implicit none
      integer m1,m2,m3
      real X(m1,m2,m3)
      real S(*), Q(*), Z(*), WORK(*)

c     TODO Proper implementation
      call schsolv232(m1,m2,m3,X,S,Q,Z,WORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact232(m1,m2,m3,S,Q,Z,ls,A,B,C,WORK,IWORK,BWORK)
      implicit none
      integer m1,m2,m3
      real S(*), Q(*), Z(*) ! S(ls), Q(ls), Z(ls)
      integer ls
      real A(m1,m1,*), B(m2,m2,*), C(m3,m3,*)

      real     WORK(*) ! 4*max(m1+1,m2+1,m3+1,4)**2
      integer IWORK(*) ! m1+m2+m3
      logical BWORK(*) ! max(m1,2*m2,m3)
      
      real alpha,beta
      real E(2,2)
      integer k
      integer p1,p2,p3
      integer w1,w2,wa,wb,wr
      integer piv1,piv2,piv3

      real one, zero
      parameter ( one=1.0E0, zero=0.0E0 )

      w1 = 1
      w2 = w1 + 2*max(m1,m2)
      wa = w2 + 2*max(m1,m2)
      wb = wa + m3
      wr = wb + m3

      piv1 = 1
      piv2 = piv1 + m1
      piv3 = piv2 + m2

      p3 = 1
      p1 = p3 + m3*m3
      p2 = p1 + m1*m1

      call dcopy(m1*m1,A(1,1,1),1,S(p1),1)
      call dcopy(m1*m1,A(1,1,2),1,Q(p1),1)
      call dcopy(m1*m1,A(1,1,3),1,Z(p1),1)

      call transpose(S(p2),m2,B(1,1,1),m2)
      call transpose(Q(p2),m2,B(1,1,2),m2)
      call transpose(Z(p2),m2,B(1,1,3),m2)

      call imass(m1,A,IWORK(piv1))
      call imass(m2,B,IWORK(piv2))
      call imass(m3,C,IWORK(piv3))

c     S3 balanced
      call dcopy(m3*m3,C(1,1,1),1,S(p3),1)
      call schfact1('T',1,m3,S(p3),Q(p3),Z(p3),WORK(wa),WORK(wb),
     $     C(1,1,2),IWORK(piv3),WORK(wr),BWORK)
      call schbal('T',m3,S(p3),Q(p3),Z(p3),WORK(wa),WORK(wb))

      p1 = p2 + m2*m2
      p2 = p1 + m1*m1

      call dcopy(m1*m1,A(1,1,1),1,S(p1),1)
      call schfact1('N',1,m1,S(p1),Q(p1),Z(p1),WORK(w1),WORK(w2),
     $     A(1,1,2),IWORK(piv1),WORK(wr),BWORK)

      k = m3
      do while(k.gt.0)
      alpha = WORK(wa-1+k)
      beta  = WORK(wb-1+k)
      if(beta.eq.zero) then ! real eigenvalue
         call dcopy(m2*m2,         B(1,1,1),1,S(p2),1)
         call daxpy(m2*m2,alpha,B(1,1,3),1,S(p2),1)
         call schfact1('T',1,m2,S(p2),Q(p2),Z(p2),WORK(w1),WORK(w2),
     $        B(1,1,2),IWORK(piv2),WORK(wr),BWORK)
         p2 = p2 + m2*m2
         k = k-1
      else ! complex eigenvalue pair
         call set_ident(E,2,2)
         call dkron(2,m2,2,m2,one,E,B(1,1,1),zero,S(p2))
         E(1,1)=alpha
         E(2,2)=alpha
         E(1,2)=-beta
         E(2,1)= beta
         call dkron(2,m2,2,m2,one,E,B(1,1,3),one,S(p2))
         call schfact1('T',2,m2,S(p2),Q(p2),Z(p2),WORK(w1),WORK(w2),
     $        B(1,1,2),IWORK(piv2),WORK(wr),BWORK)
         p2 = p2 + 4*m2*m2
         k = k-2
      endif
      enddo ! k
      ls = p2-p3
      return
      end
c-----------------------------------------------------------------------
      subroutine schsolv232(m1,m2,m3,X,S,Q,Z,W)
c     X : (on entry) RHS, (on exit) X
      implicit none
      integer m1,m2,m3
      real X(m1,m2,m3)
      real S(*), Q(*), Z(*)
      real W(m1,m2,*) ! W(m1,m2,m3+2)

      integer i,j,k, p,pk,pw, p1,p2,p3, r1,r2,r3

      real beta
      real zero, one, minus
      parameter ( zero=0.0E0, one=1.0E0, minus=-1.0E0 )
      
      p3 = 1
      r1 = p3 + m3*m3
      r2 = r1 + m1*m1
      p1 = r2 + m2*m2
      p2 = p1 + m1*m1
      pw = m3 + 1

c     W(:) = kron(Z3,I,I)'*X(:)
      call kron_mxv3(m3, m1,m2,m3, Z(p3), X, W)

      k = m3
      do while(k.gt.0)

c     beta = S3(k,k-1)
      if(k.eq.1) then
         beta = zero
      else
         r3 = p3+(k-1)+m3*(k-2)
         beta = S(r3)
      endif

c     p = number of 2d slices in block 
      if(beta.eq.zero) then ! real eigenvalue  
         p = 1
      else ! complex eigenvalue pair
         p = 2
      endif
      pk = k-p+1 ! index of first slice

c     W(:,:,pk)=kron(Q2,Q1)*(kron(I,S1)+kron(S2,I))\kron(Z2,Z1)'*W(:,:,pk)
      call schsolv22(m1,p*m2,W(1,1,pk),
     $     S(p1),S(p2),Q(p1),Q(p2),Z(p1),Z(p2),W(1,1,pw))
      p2 = p2 + (p*m2)**2

c     Column-oriented back substitution Golub-Van Loan p. 109
c     W(:,:,1:pk-1) = W(:,:,1:pk-1) - kron(S3(1:pk-1,pk:k),M12)*W(:,:,pk:k)
      call kron_mxv1(m1, m1,m2,p, Q(r1), W(1,1,pk), W(1,1,pw))
      call kron_mxv2(m2, m1,m2,p, Z(r2), W(1,1,pw), X(1,1,pk))
      do i=1,pk-1
         do j=pk,k
            r3 = p3 + (i-1)+m3*(j-1)
            call daxpy(m1*m2, -S(r3), X(1,1,j),1, W(1,1,i),1)
         enddo ! j
      enddo ! i

      k = k-p
      enddo ! k

c     X(:) = kron(Q3,I,I)*W(:)
      call kron_mxv3(m3, m1,m2,m3, Q(p3), W, X)
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact233(m1,m2,m3,S,Q,Z,ls,A,B,C,WORK,IWORK,BWORK)
      implicit none
      integer m1,m2,m3
      real S(*), Q(*), Z(*) ! S(ls), Q(ls), Z(ls)
      integer ls
      real A(m1,m1,*), B(m2,m2,*), C(m3,m3,*)

      real     WORK(*) ! 4*max(m1+1,m2+1,m3+1,4)**2
      integer IWORK(*) ! m1+m2+m3
      logical BWORK(*) ! max(2*m1,2*m2,m3)
      
      real alpha,beta
      real E(2,2)
      integer k
      integer p1,p2,p3
      integer w1,w2,wa,wb,wr
      integer piv1,piv2,piv3

      real one, zero
      parameter ( one=1.0E0, zero=0.0E0 )

      w1 = 1
      w2 = w1 + 2*max(m1,m2)
      wa = w2 + 2*max(m1,m2)
      wb = wa + m3
      wr = wb + m3

      piv1 = 1
      piv2 = piv1 + m1
      piv3 = piv2 + m2
      
      p3 = 1
      p1 = p3 + m3*m3
      p2 = p1 + m1*m1

      call dcopy(m1*m1,A(1,1,1),1,S(p1),1)
      call dcopy(m1*m1,A(1,1,2),1,Q(p1),1)
      call dcopy(m1*m1,A(1,1,3),1,Z(p1),1)

      call transpose(S(p2),m2,B(1,1,1),m2)
      call transpose(Q(p2),m2,B(1,1,2),m2)
      call transpose(Z(p2),m2,B(1,1,3),m2)

      call imass(m1,A,IWORK(piv1))
      call imass(m2,B,IWORK(piv2))
      call imass(m3,C,IWORK(piv3))

c     S3 balanced
      call dcopy(m3*m3,C(1,1,1),1,S(p3),1)
      call schfact1('T',1,m3,S(p3),Q(p3),Z(p3),WORK(wa),WORK(wb),
     $     C(1,1,2),IWORK(piv3),WORK(wr),BWORK)
      call schbal('T',m3,S(p3),Q(p3),Z(p3),WORK(wa),WORK(wb))

      p1 = p2 + m2*m2
      k = m3
      do while(k.gt.0)
      alpha = WORK(wa-1+k)
      beta  = WORK(wb-1+k)
      if(beta.eq.zero) then ! real eigenvalue
         p2 = p1 + m1*m1
         call dcopy(m1*m1,      A(1,1,1),1,S(p1),1)
         call daxpy(m1*m1,alpha,A(1,1,3),1,S(p1),1)
         call schfact1('N',1,m1,S(p1),Q(p1),Z(p1),WORK(w1),WORK(w2),
     $        A(1,1,2),IWORK(piv1),WORK(wr),BWORK)

         call dcopy(m2*m2,      B(1,1,1),1,S(p2),1)
         call daxpy(m2*m2,alpha,B(1,1,3),1,S(p2),1)
         call schfact1('T',1,m2,S(p2),Q(p2),Z(p2),WORK(w1),WORK(w2),
     $        B(1,1,2),IWORK(piv2),WORK(wr),BWORK)
         p1 = p2 + m2*m2
         k = k-1
      else ! complex eigenvalue pair
         p2 = p1 + 4*m1*m1
         call set_ident(E,2,2)
         call dkron(2,m1,2,m1,one,E,A(1,1,1),zero,S(p1))
         call dkron(2,m2,2,m2,one,E,B(1,1,1),zero,S(p2))
         E(1,1)=alpha
         E(2,2)=alpha
         E(1,2)=-beta
         E(2,1)= beta
         call dkron(2,m1,2,m1,one,E,A(1,1,3),one,S(p1))
         E(1,2)= beta
         E(2,1)=-beta
         call dkron(2,m2,2,m2,one,E,B(1,1,3),one,S(p2))
         call schfact1('N',2,m1,S(p1),Q(p1),Z(p1),WORK(w1),WORK(w2),
     $        A(1,1,2),IWORK(piv1),WORK(wr),BWORK)
         call schfact1('T',2,m2,S(p2),Q(p2),Z(p2),WORK(w1),WORK(w2),
     $        B(1,1,2),IWORK(piv2),WORK(wr),BWORK)   
         p1 = p2 + 4*m2*m2
         k = k-2
      endif
      enddo ! k
      ls = p1-p3
      return
      end
c-----------------------------------------------------------------------
      subroutine schsolv233(m1,m2,m3,X,S,Q,Z,W)
      implicit none
      integer m1,m2,m3
      real X(m1,m2,*) ! X(m1,m2,m3)
      real S(*), Q(*), Z(*)
      real W(m1,m2,*) ! W(m1,m2,m3+8)

      integer i,j,k, j11,j21,j12,j22, p,pk,pw, p1,p2,p3, r1,r2,r3

      real beta
      real zero, one, minus
      parameter ( zero=0.0E0, one=1.0E0, minus=-1.0E0 )
      
      p3 = 1
      r1 = p3 + m3*m3
      r2 = r1 + m1*m1
      p1 = r2 + m2*m2
      pw = m3+1

c     W(:) = kron(Z3,I,I)'*X(:)
      call kron_mxv3(m3, m1,m2,m3, Z(p3), X, W)

      k = m3
      do while(k.gt.0)
c     beta = S3(k,k-1)
      if(k.eq.1) then
         beta = zero
      else
         r3 = p3+(k-1)+m3*(k-2)
         beta = S(r3)
      endif

c     p = number of 2d slices in block
      if(beta.eq.zero) then ! real eigenvalue
c        W(:,:,k)=kron(Q2,Q1)*(kron(I,S1)+kron(S2,I))\kron(Z2,Z1)'*W(:,:,k)
         p = 1
         p2 = p1 + (p*m1)**2
         call schsolv22(p*m1,p*m2,W(1,1,k),
     $        S(p1),S(p2),Q(p1),Q(p2),Z(p1),Z(p2),W(1,1,pw))
         p1 = p2 + (p*m2)**2
      else ! complex eigenvalue pair
         p = 2
c        Y=[W(:,:,k-1), -W(:,:,k); W(:,:,k), W(:,:,k-1)]
         do j=1,m2
            j11 = 2*j-1
            j21 = 2*j
            j12 = 2*(j+m2)-1
            j22 = 2*(j+m2)
            call dcopy(m1, W(1,j,k-1),1, W(1,j11,pw),1)
            call dcopy(m1, W(1,j,k-1),1, W(1,j22,pw),1)
            call dcopy(m1, W(1,j,k  ),1, W(1,j21,pw),1)
            call dcopy(m1, W(1,j,k  ),1, W(1,j12,pw),1)
            call dscal(m1, minus       , W(1,j12,pw),1)
         enddo ! j
c        Y(:)=kron(Q2,Q1)*(kron(I,S1)+kron(S2,I))\kron(Z2,Z1)'*Y(:)
         p2 = p1 + (p*m1)**2
         call schsolv22(p*m1,p*m2,W(1,1,pw),
     $        S(p1),S(p2),Q(p1),Q(p2),Z(p1),Z(p2),W(1,1,pw+4))
         p1 = p2 + (p*m2)**2
c        W(:,:,k-1) = Y(:,1,:,1)
c        W(:,:,k  ) = Y(:,2,:,1)
         do j=1,m2
            j11 = 2*j-1
            j21 = 2*j
            call dcopy(m1, W(1,j11,pw),1, W(1,j,k-1),1)
            call dcopy(m1, W(1,j21,pw),1, W(1,j,k  ),1)
         enddo ! j
      endif
      pk = k-p+1 ! index of first slice

c     Column-oriented back substitution Golub-Van Loan p. 109
c     W(:,:,1:pk-1) = W(:,:,1:pk-1)-kron(S3(1:pk-1,pk:k),M12)*W(:,:,pk:k)
      call kron_mxv1(m1, m1,m2,p, Q(r1), W(1,1,pk),  W(1,1,pw  ))
      call kron_mxv2(m2, m1,m2,p, Z(r2), W(1,1,pw),  X(1,1,pk  ))
      call kron_mxv1(m1, m1,m2,p, Z(r1), W(1,1,pk),  W(1,1,pw  ))
      call kron_mxv2(m2, m1,m2,p, Q(r2), W(1,1,pw),  W(1,1,pw+p))
      call daxpy(m1*m2*p, one, W(1,1,pw+p),1, X(1,1,pk),1)
      do i=1,pk-1
         do j=pk,k
            r3 = p3 + (i-1)+m3*(j-1)
            call daxpy(m1*m2, -S(r3), X(1,1,j),1, W(1,1,i),1)
         enddo ! j
      enddo ! i

      k = k-p
      enddo ! k
c     X(:) = kron(Q3,I,I)*W(:)
      call kron_mxv3(m3, m1,m2,m3, Q(p3), W, X)
      return
      end
c-----------------------------------------------------------------------
      subroutine rank_info(arank,irank,tflag)
      implicit none
      character*3 arank
      integer irank, tflag
c     Total rank
      if (arank.eq.'222') irank=6
      if (arank.eq.'223'.or.arank.eq.'232'.or.arank.eq.'322') irank=7
      if (arank.eq.'233'.or.arank.eq.'323'.or.arank.eq.'332') irank=8
c     Transpose flag
      if (arank.eq.'222'.or.arank.eq.'232'.or.arank.eq.'233') tflag=0
      if (arank.eq.'322'.or.arank.eq.'323') tflag=1
      if (arank.eq.'332') tflag=2
      if (arank.eq.'223') tflag=3
      return
      end
c-----------------------------------------------------------------------
      subroutine schfactc(irank,m1,m2,m3,S,Q,Z,slen,A,B,C,
     $                    WORK,IWORK,BWORK)
c     Factor canonical decompositions only (222, 232, 233)
      implicit none
      integer irank, m1,m2,m3
      real S(*), Q(*), Z(*)
      integer slen
      real A(*), B(*), C(*)

      real     WORK(*) ! 4*max(m1,m2,m3,4)**2
      integer IWORK(*) ! m1+m2+m3
      logical BWORK(*) ! max(m1,m2,m3)

      if    (irank.eq.6) then
      call schfact222(m1,m2,m3,S,Q,Z,slen,A,B,C,WORK,IWORK,BWORK)
      elseif(irank.eq.7) then
      call schfact232(m1,m2,m3,S,Q,Z,slen,A,B,C,WORK,IWORK,BWORK)
      elseif(irank.eq.8) then
      call schfact233(m1,m2,m3,S,Q,Z,slen,A,B,C,WORK,IWORK,BWORK)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine schfact3d(arank,m1,m2,m3,S,Q,Z,slen,A,B,C,
     $                     WORK,IWORK,BWORK)
      implicit none
      character*3 arank
      integer m1,m2,m3
      real S(*), Q(*), Z(*)
      integer slen
      real A(*), B(*), C(*)

      real     WORK(*) ! 4*max(m1,m2,m3,4)**2
      integer IWORK(*) ! m1+m2+m3
      logical BWORK(*) ! max(m1,m2,m3)

      integer irank, tflag
      
      call rank_info(arank,irank,tflag)

      if    (tflag.eq.0) then ! xyz
      call schfactc(irank,m1,m2,m3,S,Q,Z,slen,A,B,C,WORK,IWORK,BWORK)
      elseif(tflag.eq.1) then ! xzy
      call schfactc(irank,m1,m3,m2,S,Q,Z,slen,A,C,B,WORK,IWORK,BWORK)
      elseif(tflag.eq.2) then ! zyx
      call schfactc(irank,m3,m2,m1,S,Q,Z,slen,C,B,A,WORK,IWORK,BWORK)
      elseif(tflag.eq.3) then ! yxz
      call schfactc(irank,m2,m1,m3,S,Q,Z,slen,B,A,C,WORK,IWORK,BWORK)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine schsolv3d(arank,m1,m2,m3,X,S,Q,Z,F,WORK)
      implicit none
      character*3 arank
      integer m1,m2,m3
      real X(m1,m2,m3), F(m1,m2,m3)
      real S(*), Q(*), Z(*), WORK(*) ! WORK(m1,m2,2*m3+8)
      
      integer irank, tflag
      integer i,j,k, idx,pw
      integer s1,s2,s3, n1,n2,n3
      
      call rank_info(arank,irank,tflag)

      pw = m1*m2*m3+1
      s1 = 1
      s2 = m1
      s3 = m1*m2
      n1 = m1
      n2 = m2
      n3 = m3

c     Stride
      if    (tflag.eq.1) then ! xzy
      s3 = m1
      s2 = m1*m3
      n2 = m3
      n3 = m2
      elseif(tflag.eq.2) then ! zyx
      s3 = 1
      s1 = m3*m2
      n1 = m3
      n3 = m1
      elseif(tflag.eq.3) then ! yxz
      s2 = 1
      s1 = m2
      n1 = m2
      n2 = m1
      endif

      do k=1,m3
      do j=1,m2
      do i=1,m1
         idx = 1 + s1*(i-1) + s2*(j-1) + s3*(k-1)
         WORK(idx) = F(i,j,k)
      enddo ! i
      enddo ! j
      enddo ! k

      if    (irank.eq.6) then
      call schsolv222(n1,n2,n3,WORK,S,Q,Z,WORK(pw))
      elseif(irank.eq.7) then
      call schsolv232(n1,n2,n3,WORK,S,Q,Z,WORK(pw))
      elseif(irank.eq.8) then
      call schsolv233(n1,n2,n3,WORK,S,Q,Z,WORK(pw))
      endif

      do k=1,m3
      do j=1,m2
      do i=1,m1
         idx = 1 + s1*(i-1) + s2*(j-1) + s3*(k-1)
         X(i,j,k) = WORK(idx)
      enddo ! i
      enddo ! j
      enddo ! k

      return
      end 
c-----------------------------------------------------------------------
      subroutine fdmfact1(TRANS,m,lambda,WR,WL,A,B)
c     Setup for 1D generalized eignevalue problem
c
c     Returns:
c
c     Eigenvalues lamdba
c     Right eigenvectors A* WR = B *WR*diag(lambda)
c     Left eigenvectors  A'*WL = B'*WL*diag(lambda)
c     

      implicit none

      character*1 TRANS
      integer m
      complex*16 lambda(m),WL(m,m),WR(m,m)
      real A(m,m),B(m,m)

      complex*16 nrm
c     FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
      real ACP(m,m),BCP(m,m),ZR(m),ZI(m),ZD(m),VL(m,m),VR(m,m),
     $     WORK(8*m*m)
      real atemp,btemp
      integer i,j,k,lwork,info
      complex*16 one,zero
      parameter (one=(1.0E0,0.0E0), zero=(0.0E0,0.0E0))
      external dggev
      
      call dcopy(m*m,A,1,ACP,1)
      call dcopy(m*m,B,1,BCP,1)

      lwork=8*m*m
      call dggev('V', 'V', m, ACP, m, BCP, m, 
     $           ZR, ZI, ZD, VL, m, VR, m, WORK, lwork, info)

      do j=1,m
         lambda(j)=one
         do i=1,m
            WL(i,j)=zero
            WR(i,j)=zero
         enddo
      enddo

      if (TRANS.eq.'N') then
c     Right eigenvectors are stored normally
c     Left eigenvectors are stored conjugated transposed
      do j=1,m
         if (ZI(j).eq.0.0E0) then
            do i=1,m
               WL(j,i)=cmplx(VL(i,j),0.0E0)
               WR(i,j)=cmplx(VR(i,j),0.0E0)
            enddo
         elseif (ZI(j).gt.0.0E0) then
            do i=1,m
               WL(j,i)=cmplx(VL(i,j),-VL(i,j+1))
               WR(i,j)=cmplx(VR(i,j), VR(i,j+1))
            enddo
         else
            do i=1,m
               WL(j,i)=cmplx(VL(i,j-1), VL(i,j))
               WR(i,j)=cmplx(VR(i,j-1),-VR(i,j))
            enddo
         endif
      enddo

c     FIXME
      do k=1,m
         nrm=zero
         do j=1,m
         do i=1,m
            nrm=nrm+WL(k,i)*B(i,j)*WR(j,k)
         enddo
         enddo
         nrm = one/sqrt(nrm)
         call zscal(m,nrm,WR(1,k),1)
         call zscal(m,nrm,WL(k,1),m)
         lambda(k)=cmplx(ZR(k)/ZD(k),ZI(k)/ZD(k))
      enddo
c
      elseif (TRANS.eq.'T') then
c     Right eigenvectors are stored transposed
c     Left eigenvectors are stored complex conjugated
      do j=1,m
         if (ZI(j).eq.0.0E0) then
            do i=1,m
               WL(i,j)=cmplx(VL(i,j),0.0E0)
               WR(j,i)=cmplx(VR(i,j),0.0E0)
            enddo
         elseif (ZI(j).gt.0.0E0) then
            do i=1,m
               WL(i,j)=cmplx(VL(i,j),-VL(i,j+1))
               WR(j,i)=cmplx(VR(i,j), VR(i,j+1))
            enddo
         else
            do i=1,m
               WL(i,j)=cmplx(VL(i,j-1), VL(i,j))
               WR(j,i)=cmplx(VR(i,j-1),-VR(i,j))
            enddo
         endif
      enddo

      do k=1,m
         nrm=zero
         do j=1,m
         do i=1,m
            nrm=nrm+WL(i,k)*B(i,j)*WR(k,j)
         enddo
         enddo
         nrm = one/sqrt(nrm)
         call zscal(m,nrm,WR(k,1),m)
         call zscal(m,nrm,WL(1,k),1)
         lambda(k)=cmplx(ZR(k)/ZD(k),ZI(k)/ZD(k))
      enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfact22(m1,m2,L,V1,V2,W1,W2,A,B)
      implicit none 
      integer m1,m2
      
      complex*16 L(m1,m2),
     $        V1(m1,m1),W1(m1,m1),
     $        V2(m2,m2),W2(m2,m2)
      real A(m1,m1,2),B(m2,m2,2)

c     FIXME FIXME
      complex*16 L1(m1),L2(m2)
      integer i,j

      call fdmfact1('N',m1,L1,V1,W1,A(1,1,1),A(1,1,2))
      call fdmfact1('T',m2,L2,V2,W2,B(1,1,1),B(1,1,2))

      do j=1,m2
      do i=1,m1
         L(i,j)=L1(i)+L2(j)
      enddo
      enddo
       
      return
      end
c-----------------------------------------------------------------------
      subroutine zfdmfact1(TRANS,m,lambda,WR,WL,z,A)
c     
c     Solves generalized complex-shifted eigenproblem
c
c     (A(:,:,1)+z*A(:,:,3)) *WR = A(:,:,2) *WR*diag(lambda)
c     (A(:,:,1)+z*A(:,:,3))'*WL = A(:,:,2)'*WL*diag(lambda)
c     
c     With normalization
c     
c     WL'*A(:,:,2)*WR = I

      implicit none
      character*1 TRANS
      integer m
      complex*16 lambda(m),WR(m,m),WL(m,m),z
      real A(m,m,3)

c     FIXME FIXME FIXME FIXME FIXME
      real RWORK(8*m)
      complex*16 nrm,alpha(m),beta(m),AZ(m,m),BZ(m,m),WORK(m,m*8)
      complex*16 atemp
      integer i,j,k, lwork,info

      complex*16 one,zero
      parameter ( one=(1.0E0,0.0E0), zero=(0.0E0,0.0E0) )

      external zggev

      do j=1,m
      do i=1,m
         AZ(i,j)=A(i,j,1)+z*A(i,j,3)
         BZ(i,j)=A(i,j,2)
         WR(i,j)=zero
         WL(i,j)=zero
      enddo ! i
      enddo ! j


      lwork=8*m*m
      call zggev('V','V',m,AZ,m,BZ,m,alpha,beta,
     $           WL,m,WR,m,WORK,lwork,RWORK,info)

      lambda(1)=one
      lambda(m)=one
      do i=1,m
         lambda(i)=alpha(i)/beta(i)
      enddo

c     FIXME
c     Left eigenvectors always complex conjugated
      do j=1,m
      do i=1,m
         WL(i,j)=conjg(WL(i,j))
      enddo ! i
      enddo ! j

c     FIXME
c     Normalization
      do k=1,m
         nrm=zero
         do j=1,m
         do i=1,m
         nrm=nrm+WL(i,k)*A(i,j,2)*WR(j,k)
         enddo ! i
         enddo ! j
         nrm=one/sqrt(nrm)
         call zscal(m,nrm,WL(1,k),1)
         call zscal(m,nrm,WR(1,k),1)
      enddo ! k

      if    (TRANS.eq.'N') then
c     Left eigenvectors transposed
         call zcopy(m*m,WL,1,WORK,1)
         do i=1,m
            call zcopy(m,WORK(i,1),m,WL(1,i),1)
         enddo ! i

      elseif(TRANS.eq.'T') then
c     Right eigenvectors transposed
         call zcopy(m*m,WR,1,WORK,1)
         do i=1,m
            call zcopy(m,WORK(i,1),m,WR(1,i),1)
         enddo ! i

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmsolv22(m1,m2,X,L,V1,V2,W1,W2,F)
      implicit none

      integer m1,m2
      real X(m1,m2), F(m1,m2)
      complex*16 L(m1,m2),
     $        V1(m1,m1),W1(m1,m1),
     $        V2(m2,m2),W2(m2,m2) 

      complex*16 R(m1,m2),Z(m1,m2)
      integer i,j
      
      do j=1,m2
      do i=1,m1
         R(i,j)=cmplx(F(i,j),0.0E0)
      enddo ! i
      enddo ! j

c     R(:)=kron(W2,W1)'*R(:)      
c     R=conj(W1).'*R*conj(W2)      
      call mxmc(W1,m1,R,m1,Z,m2) 
      call mxmc(Z,m1,W2,m2,R,m2) 

c     R=R./L
      do j=1,m2
      do i=1,m1
         R(i,j)=R(i,j)/L(i,j)
      enddo ! i
      enddo ! j

c     R(:)=kron(V2,V1)*R(:)      
c     R=V1*R*V2.'
      call mxmc(V1,m1,R,m1,Z,m2) 
      call mxmc(Z,m1,V2,m2,R,m2) 

      do j=1,m2
      do i=1,m1
         X(i,j)=real(R(i,j))
      enddo ! i
      enddo ! j

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfact222(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      implicit none 
      integer m1,m2,m3
      
      complex*16 L(m1,m2,m3),
     $        V1(m1,m1),W1(m1,m1),
     $        V2(m2,m2),W2(m2,m2),
     $        V3(m3,m3),W3(m3,m3)
      real A(m1,m1,3),B(m2,m2,3),C(m3,m3,3)

      complex*16 L1(m1),L2(m2),L3(m3)
      integer i,j,k

      call fdmfact1('N',m1,L1,V1,W1,A(1,1,1),A(1,1,2))
      call fdmfact1('T',m2,L2,V2,W2,B(1,1,1),B(1,1,2))
      call fdmfact1('T',m3,L3,V3,W3,C(1,1,1),C(1,1,2))

      do k=1,m3
      do j=1,m2
      do i=1,m1
         L(i,j,k)=L1(i)+L2(j)+L3(k)
      enddo ! i
      enddo ! j
      enddo ! k

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfact232(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      implicit none
      integer m1,m2,m3
      
      complex*16 L(m1,m2,m3),
     $        L1(m1),L2(m2),L3(m3),
     $        V1(m1,m1),   W1(m1,m1),
     $        V2(m2,m2,m3),W2(m2,m2,m3),
     $        V3(m3,m3),   W3(m3,m3)
      real A(m1,m1,3),B(m2,m2,3),C(m3,m3,3)

      integer i,j,k
      complex*16 one
      parameter ( one=(1.0E0,0.0E0) )


      do j=1,m2
      do i=1,m1
         L(i,j,1 ) = one
         L(i,j,m3) = one
      enddo ! i
      enddo ! j

      call fdmfact1('N',m1,L1,V1,W1,A(1,1,1),A(1,1,2))
      call fdmfact1('T',m3,L3,V3,W3,C(1,1,1),C(1,1,2))
      do k=1,m3
      call zfdmfact1('T',m2,L2,V2(1,1,k),W2(1,1,k),L3(k),B)
      do j=1,m2
      do i=1,m1
         L(i,j,k) = L1(i)+L2(j)
      enddo ! i
      enddo ! j
      enddo ! k

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfact233(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      implicit none
      integer m1,m2,m3
      
      complex*16 L(m1,m2,m3),
     $        L1(m1),L2(m2),L3(m3),
     $        V1(m1,m1,m3),W1(m1,m1,m3),
     $        V2(m2,m2,m3),W2(m2,m2,m3),
     $        V3(m3,m3),   W3(m3,m3)
      real A(m1,m1,3),B(m2,m2,2),C(m3,m3,2)

      integer i,j,k
      complex*16 one
      parameter ( one=(1.0E0,0.0E0) )


      do j=1,m2
      do i=1,m1
         L(i,j,1 ) = one
         L(i,j,m3) = one
      enddo ! i
      enddo ! j

      call fdmfact1('T',m3,L3,V3,W3,C(1,1,1),C(1,1,2))
      do k=1,m3
      call zfdmfact1('N',m1,L1,V1(1,1,k),W1(1,1,k),L3(k),A)
      call zfdmfact1('T',m2,L2,V2(1,1,k),W2(1,1,k),L3(k),B)
      do j=1,m2
      do i=1,m1
         L(i,j,k) = L1(i)+L2(j)
      enddo ! i
      enddo ! j
      enddo ! k

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmsolv222(m1,m2,m3,Z,L,V1,V2,V3,W1,W2,W3)
      implicit none
      integer m1,m2,m3

      complex*16 Z(m1,m2,m3),L(m1,m2,m3),
     $        V1(m1,m1),   W1(m1,m1),
     $        V2(m2,m2),   W2(m2,m2),
     $        V3(m3,m3),   W3(m3,m3)

c     FIXME
      complex*16 R(m1,m2,m3) 
      integer i,j,k

      call mxmc(W1,m1,Z,m1,R,m2*m3) 
      call mxmc(R,m1*m2,W3,m3,Z,m3)
      do k=1,m3
      call mxmc(Z(1,1,k),m1,W2,m2,R(1,1,k),m2) 
      do j=1,m2
      do i=1,m1
         R(i,j,k)=R(i,j,k)/L(i,j,k)
      enddo ! i
      enddo ! j
      call mxmc(R(1,1,k),m1,V2,m2,Z(1,1,k),m2)
      enddo ! k
      call mxmc(Z,m1*m2,V3,m3,R,m3)
      call mxmc(V1,m1,R,m1,Z,m2*m3) 
      return
      end
c-----------------------------------------------------------------------
      subroutine fdmsolv232(m1,m2,m3,Z,L,V1,V2,V3,W1,W2,W3)
      implicit none
      integer m1,m2,m3

      complex*16 Z(m1,m2,m3),L(m1,m2,m3),
     $        V1(m1,m1),   W1(m1,m1),
     $        V2(m2,m2,m3),W2(m2,m2,m3),
     $        V3(m3,m3),   W3(m3,m3)

c     FIXME
      complex*16 R(m1,m2,m3)
      integer i,j,k

      call mxmc(W1,m1,Z,m1,R,m2*m3) 
      call mxmc(R,m1*m2,W3,m3,Z,m3)
      do k=1,m3
      call mxmc(Z(1,1,k),m1,W2(1,1,k),m2,R(1,1,k),m2) 
      do j=1,m2
      do i=1,m1
         R(i,j,k)=R(i,j,k)/L(i,j,k)
      enddo ! i
      enddo ! j
      call mxmc(R(1,1,k),m1,V2(1,1,k),m2,Z(1,1,k),m2)
      enddo ! k
      call mxmc(Z,m1*m2,V3,m3,R,m3)
      call mxmc(V1,m1,R,m1,Z,m2*m3) 
      return
      end
c-----------------------------------------------------------------------
      subroutine fdmsolv233(m1,m2,m3,Z,L,V1,V2,V3,W1,W2,W3)
      implicit none
      integer m1,m2,m3

      complex*16 Z(m1,m2,m3),L(m1,m2,m3),
     $        V1(m1,m1,m3),W1(m1,m1,m3),
     $        V2(m2,m2,m3),W2(m2,m2,m3),
     $        V3(m3,m3),   W3(m3,m3)

c     FIXME
      complex*16 R(m1,m2,m3)
      integer i,j,k

      call mxmc(Z,m1*m2,W3,m3,R,m3)
      do k=1,m3
      call mxmc(R(1,1,k),m1,W2(1,1,k),m2,Z(1,1,k),m2) 
      call mxmc(W1(1,1,k),m1,Z(1,1,k),m1,R(1,1,k),m2) 
      do j=1,m2
      do i=1,m1
         R(i,j,k)=R(i,j,k)/L(i,j,k)
      enddo ! i
      enddo ! j
      call mxmc(V1(1,1,k),m1,R(1,1,k),m1,Z(1,1,k),m2) 
      call mxmc(Z(1,1,k),m1,V2(1,1,k),m2,R(1,1,k),m2)
      enddo ! k
      call mxmc(R,m1*m2,V3,m3,Z,m3)
      return
      end
c-----------------------------------------------------------------------
      subroutine rank_len(arank,m1,m2,m3,vlen,p1,p2,p3,tflag,irank)
      implicit none
      character*3 arank
      integer m1,m2,m3

      integer vlen, tflag, p1,p2,p3, irank
      integer r1,r2,r3

      read(arank,'(I1,I1,I1)')r3,r2,r1
c     Get total rank
      irank = r1 + r2 + r3
c     Size of V
      vlen = m1**r1 + m2**r2 + m3**r3
c     Transpose flag
      if (arank.eq.'222'.or.arank.eq.'232'.or.arank.eq.'233') tflag=0
      if (arank.eq.'322'.or.arank.eq.'323') tflag=1
      if (arank.eq.'332') tflag=2
      if (arank.eq.'223') tflag=3
      if    (tflag.eq.1) then ! xzy
      call swapi(m2,m3)
      call swapi(r2,r3)
      elseif(tflag.eq.2) then ! zyx
      call swapi(m1,m3)
      call swapi(r1,r3)
      elseif(tflag.eq.3) then ! yxz
      call swapi(m1,m2)
      call swapi(r1,r2)
      endif
c     Start of V1, V2, V3
      p1 = 1
      p2 = 1 + m1**r1
      p3 = 1 + m1**r1 + m2**r2
      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfactc(irnk,m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
c     Canonical factorization
      implicit none
      integer irnk,m1,m2,m3
      real A(m1,m1,3),B(m2,m2,3),C(m3,m3,3)

      complex*16 L(m1*m2*m3),
     $ V1(m1,m1,*),V2(m2,m2,*),V3(m3,m3),
     $ W1(m1,m1,*),W2(m2,m2,*),W3(m3,m3)
      if    (irnk.eq.6) then
      call fdmfact222(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      elseif(irnk.eq.7) then
      call fdmfact232(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      elseif(irnk.eq.8) then
      call fdmfact233(m1,m2,m3,L,V1,V2,V3,W1,W2,W3,A,B,C)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine fdmfact3d(arank,m1,m2,m3,L,V,W,vlen,A,B,C)
c     Setup driver for Fast Diagonalization Method
      implicit none
      character*3 arank
      integer m1,m2,m3,vlen
      real A(m1,m1,3),B(m2,m2,3),C(m3,m3,3)

      integer p1,p2,p3, tflag, irank
      complex*16 L(m1*m2*m3),V(*),W(*)

      call rank_len(arank,m1,m2,m3,vlen,p1,p2,p3,tflag,irank)

      if    (tflag.eq.1) then
      call fdmfactc(irank,m1,m2,m3,L,
     $     V(p1),V(p2),V(p3),W(p1),W(p2),W(p3),A,C,B)
      elseif(tflag.eq.2) then
      call fdmfactc(irank,m1,m2,m3,L,
     $     V(p1),V(p2),V(p3),W(p1),W(p2),W(p3),C,B,A)
      elseif(tflag.eq.3) then
      call fdmfactc(irank,m1,m2,m3,L,
     $     V(p1),V(p2),V(p3),W(p1),W(p2),W(p3),B,A,C)
      else
      call fdmfactc(irank,m1,m2,m3,L,
     $     V(p1),V(p2),V(p3),W(p1),W(p2),W(p3),A,B,C)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fdmsolv3d(arank,m1,m2,m3,X,L,V,W,F)
c     Solver driver for Fast Diagonalization Method
      implicit none
      character*3 arank
      integer m1,m2,m3
      complex*16 L(m1*m2*m3),V(*),W(*)
      real X(m1,m2,m3), F(m1,m2,m3)

      integer vlen,tflag,irank
      integer n1,n2,n3  ! sizes
      integer s1,s2,s3  ! strides
      integer p1,p2,p3  ! pointers

      complex*16 Z(m1*m2*m3)
      integer i,j,k, idx

      n1 = m1
      n2 = m2
      n3 = m3
      call rank_len(arank,n1,n2,n3,vlen,p1,p2,p3,tflag,irank)
      
      s1 = m3*m2
      s2 = m1*m3
      s3 = m1*m2
      if    (tflag.eq.1) then
         s1 = 1
         s3 = m1
      elseif(tflag.eq.2) then
         s3 = 1
         s2 = m3
      elseif(tflag.eq.3) then
         s2 = 1
         s1 = m2
      else 
         s1 = 1
         s2 = m1
      endif

      do k=1,m3
      do j=1,m2
      do i=1,m1
         idx = 1 + (i-1)*s1 + (j-1)*s2 + (k-1)*s3
         Z(idx) = cmplx(F(i,j,k),0.0E0)
      enddo ! i
      enddo ! j
      enddo ! k

      if    (irank.eq.6) then
      call fdmsolv222(n1,n2,n3,Z,L,V(p1),V(p2),V(p3),W(p1),W(p2),W(p3))
      elseif(irank.eq.7) then
      call fdmsolv232(n1,n2,n3,Z,L,V(p1),V(p2),V(p3),W(p1),W(p2),W(p3))
      elseif(irank.eq.8) then
      call fdmsolv233(n1,n2,n3,Z,L,V(p1),V(p2),V(p3),W(p1),W(p2),W(p3))
      endif

      do k=1,m3
      do j=1,m2
      do i=1,m1
         idx = 1 + (i-1)*s1 + (j-1)*s2 + (k-1)*s3
         X(i,j,k) = real(Z(idx))
      enddo ! i
      enddo ! j
      enddo ! k

      return
      end
c-----------------------------------------------------------------------
