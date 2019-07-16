c-----------------------------------------------------------------------
      subroutine axhelm_new(au,u,h1,h2,h3)
c     Computes anisotropic Helmholtz operator on mesh 1
c     au = au + D^T ( h1 G + h3 B k k^T ) D u + h2 B u 
c     Assumes all elements are deformed

      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ANISO'

      common /helm_work/ w1,w2
      real w1(lx1,ly1,lz1,lelt)
      real w2(lx1,ly1,lz1,lelt)

      common /helm_grad/ ur,us,ut
      real ur(lx1,ly1,lz1,lelt)
      real us(lx1,ly1,lz1,lelt)
      real ut(lx1,ly1,lz1,lelt)
      
      real au(lx1,ly1,lz1,*), u(lx1,ly1,lz1,*), h1(*), h2(*), h3(*)

      logical ifh2,ifh3
      integer i, nx,ny,nz,ne,nt

      nx = nx1
      ny = ny1
      nz = nz1
      ne = nelt
      if(imesh.eq.1) ne=nelv
      nt = nx*ny*nz*ne

      ifh2 = .true.
      ifh3 = .true.
     
      if (nz.eq.1) then ! 2D

c     ur = Dr u, us = Ds u
      call mxm(DXM1, nx, u, nx, ur, ny*ne)
      do i=1,ne
         call mxm(u(1,1,1,i), nx, DYTM1, ny, us(1,1,1,i), ny)
      enddo ! i

      if(ifh3) then
c        w1 = h3( kr ur + ks us )
         call vdot2(w1, kxd,kyd, ur,us, nt)
         call col2(w1,h3,nt)
      endif


c     au = au + Dr'*(h1(Grr ur + Grs us ) + B kr w1)
      call vdot2(w2, g1m1,g4m1, ur,us, nt)
      call col2(w2,h1,nt)
      if(ifh3) call addcol4(w2, bm1,kxd,w1, nt)
      call mxma(DXTM1, nx, w2, nx, au, ny*ne)


c     au = au + Ds'*(h1(Gsr ur + Gss us ) + B ks w1)
      call vdot2(w2, g4m1,g2m1, ur,us, nt)
      call col2(w2,h1,nt)
      if(ifh3) call addcol4(w2, bm1,kyd,w1, nt)
      do i=1,ne
         call mxma(w2(1,1,1,i), nx, DYM1, ny, au(1,1,1,i), ny)
      enddo ! i


      else ! 3D

c     ur = Dr u, us = Ds u, ut = Dt u 
      call mxm(DXM1, nx, u, nx, ur, ny*nz*ne)
      do i=1,nz*ne
         call mxm(u(1,1,i,1), nx, DYTM1, ny, us(1,1,i,1), ny)
      enddo ! i
      do i=1,ne
         call mxm(u(1,1,1,i), nx*ny, DZTM1, nz, ut(1,1,1,i), nz)
      enddo ! i

      if(ifh3) then
c        w1 = h3( kr ur + ks us + kt ut )
         call vdot3(w1, kxd,kyd,kzd, ur,us,ut, nt)
         call col2(w1,h3,nt)
      endif

c     au = au + Dr'*(h1(Grr ur + Grs us + Grt ut) + B kr w1)
      call vdot3(w2, g1m1,g4m1,g5m1, ur,us,ut, nt)
      call col2(w2,h1,nt)
      if(ifh3) call addcol4(w2, bm1,kxd,w1, nt)
      call mxma(DXTM1, nx, w2, nx, au, ny*nz*ne)

c     au = au + Ds'*(h1(Gsr ur + Gss us + Gst ut) + B ks w1)
      call vdot3(w2, g4m1,g2m1,g6m1, ur,us,ut, nt)
      call col2(w2,h1,nt)
      if(ifh3) call addcol4(w2, bm1,kyd,w1, nt)
      do i=1,nz*ne
         call mxma(w2(1,1,i,1), nx, DYM1, ny, au(1,1,i,1), ny)
      enddo ! i

c     au = au + Dt'*(h1(Gtr ur + Gts us + Gtt ut) + B kt w1)
      call vdot3(w2, g5m1,g6m1,g3m1, ur,us,ut, nt)
      call col2(w2,h1,nt)
      if(ifh3) call addcol4(w2, bm1,kzd,w1, nt)
      do i=1,ne
         call mxma(w2(1,1,1,i), nx*ny, DZM1, nz, au(1,1,1,i), nz)
      enddo ! i

      endif

      if(ifh2) call addcol4(au, bm1,h2,u, nt)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_anisotr_new(axd,ayd,azd,ax,ay,az)
      implicit none
      include 'SIZE'
      include 'GEOM'

      real axd(*),ayd(*),azd(*)
      real  ax(*), ay(*), az(*)

      real w, wx,wy,wz, ux,uy,uz
      integer i, nt

      real eps
      parameter (eps=1.0E-14)

      nt = lx1*ly1*lz1*nelt


      if (lz1.eq.1) then ! 2D

      do i=1,nt
         wx = ax(i)
         wy = ay(i)
         w  = jacm1(i,1,1,1)*sqrt( wx*wx + wy*wy )
         if( w.gt.eps ) w = 1.0E0/w
         wx = wx*w
         wy = wy*w
         axd(i) = rxm1(i,1,1,1)*wx + rym1(i,1,1,1)*wy
         ayd(i) = sxm1(i,1,1,1)*wx + sym1(i,1,1,1)*wy
      enddo ! i

      else ! 3D

      do i=1,nt
         wx = ax(i)
         wy = ay(i)
         wz = az(i)
         w  = jacm1(i,1,1,1)*sqrt( wx*wx + wy*wy + wz*wz )
         if( w.gt.eps ) w = 1.0E0/w
         wx = wx*w
         wy = wy*w
         wz = wz*w
         axd(i) = rxm1(i,1,1,1)*wx + rym1(i,1,1,1)*wy + rzm1(i,1,1,1)*wz
         ayd(i) = sxm1(i,1,1,1)*wx + sym1(i,1,1,1)*wy + szm1(i,1,1,1)*wz
         azd(i) = txm1(i,1,1,1)*wx + tym1(i,1,1,1)*wy + tzm1(i,1,1,1)*wz
      enddo ! i

      endif
      return
      end
c-----------------------------------------------------------------------
