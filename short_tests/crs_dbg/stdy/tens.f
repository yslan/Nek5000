c=======================================================================
c
c     LIBRARY ROUTINES FOR FAST TENSOR PRODUCTS
c
c     June 2018
c
c     For questions, comments or suggestions, please contact:
c
c     Pablo Daniel Brubeck
c     brubeck@protonmail.com
c     
c-----------------------------------------------------------------------
      subroutine tens_init(h1,h2,bcs)
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'TOTAL'

      real h1(*), h2(*)
      integer bcs(*)

      call tens_gen_crs(lxc)
      call tens_gen_dealias(.true.)
      call tens_gen_ratio(h1,bcs)
      return
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_coef(e, h1,h2,h3, ifconv,ifanis)
c     Initialize advection-diffusion coefficients on lx1 x ly1 x lz1
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'ANISO'

      integer e
      real h1(lx1,ly1,lz1,*), h2(lx1,ly1,lz1,*), h3(lx1,ly1,lz1,*)
      logical ifconv, ifanis

      common /tens_coef/B00,C00,C11,C22,C33,C12,C13,C23,C01,C02,C03
      real B00(lx1,ly1,lz1), C00(lx1,ly1,lz1), 
     $     C11(lx1,ly1,lz1), C22(lx1,ly1,lz1), C33(lx1,ly1,lz1),
     $     C12(lx1,ly1,lz1), C13(lx1,ly1,lz1), C23(lx1,ly1,lz1),
     $     C01(lxd,lyd,lzd), C02(lxd,lyd,lzd), C03(lxd,lyd,lzd)

      common /tens_flag/ifdealias
      logical ifdealias

      integer n
 
      n = lx1*ly1*lz1
      call col3(C00, BM1(1,1,1,e),h2(1,1,1,e),n)
      call col3(C11,G1M1(1,1,1,e),h1(1,1,1,e),n)
      call col3(C22,G2M1(1,1,1,e),h1(1,1,1,e),n)
      call col3(C33,G3M1(1,1,1,e),h1(1,1,1,e),n)
      call col3(C12,G4M1(1,1,1,e),h1(1,1,1,e),n)
      call col3(C13,G5M1(1,1,1,e),h1(1,1,1,e),n)
      call col3(C23,G6M1(1,1,1,e),h1(1,1,1,e),n)
      
      if (ifanis) then
         call col3(B00, h3(1,1,1,e), BM1(1,1,1,e), n)
         call symrank1(ldim,n,C11,C22,C33,C12,C13,C23, 
     $     B00,kxd(1,1,1,e),kyd(1,1,1,e),kzd(1,1,1,e))
      endif
      call copy(B00, BM1(1,1,1,e),n)

      n = lxd*lyd*lzd
      if (ifconv) then
         if(ifdealias) then
            call copy(C01,vxd(1,1,1,e),n)
            call copy(C02,vyd(1,1,1,e),n)
            call copy(C03,vzd(1,1,1,e),n)
         else
            call set_convect_old(e,C01,C02,C03,vx,vy,vz)
         endif
      else
         call rzero(C01,n)
         call rzero(C02,n)
         call rzero(C03,n)
      endif 

      return   
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_dealias(ifd)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      logical ifd

      common /tens_flag/ifdealias
      logical ifdealias

      common /tens_dealias/
     $ DXMD,DYMD,DZMD,DXMDT,DYMDT,DZMDT,wxmd,wymd,wzmd   
      real
     $ DXMD (lxd,lxd), DYMD (lyd,lyd), DZMD (lzd,lzd),
     $ DXMDT(lxd,lxd), DYMDT(lyd,lyd), DZMDT(lzd,lzd),
     $ wxmd (lxd),     wymd (lyd),     wzmd (lzd)

      real WK(lxd*lxd)

      ifdealias = ifd

      if(ifdealias) then
         call zwgl(WK,wxmd,lxd)
         call gen_dgl(DXMD,DXMDT,lxd,lxd,WK)

         call zwgl(WK,wymd,lyd)
         call gen_dgl(DYMD,DYMDT,lyd,lyd,WK)

         if(nz1.gt.1) then
         call zwgl(WK,wzmd,lzd)
         call gen_dgl(DZMD,DZMDT,lzd,lzd,WK)
         endif
      else
         call dcopy(nx1,wxm1     ,1,wxmd ,1)
         call dcopy(nx1*nx1,DXM1 ,1,DXMD ,1)
         call dcopy(nx1*nx1,DXTM1,1,DXMDT,1)

         call dcopy(ny1,wym1     ,1,wymd ,1)
         call dcopy(ny1*ny1,DYM1 ,1,DYMD ,1)
         call dcopy(ny1*ny1,DYTM1,1,DYMDT,1)

         if(nz1.gt.1) then
         call dcopy(nz1,wzm1     ,1,wzmd ,1)
         call dcopy(nz1*nz1,DZM1 ,1,DZMD ,1)
         call dcopy(nz1*nz1,DZTM1,1,DZMDT,1)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_basis(mx1,my1,mz1)
c     Initialize GLL/GL interpolation and derivative operators
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'HSMGL'
      
      integer mx1,my1,mz1

      common /tens_dealias/
     $ DXMD,DYMD,DZMD,DXMDT,DYMDT,DZMDT,wxmd,wymd,wzmd   
      real
     $ DXMD (lxd,lxd), DYMD (lyd,lyd), DZMD (lzd,lzd),
     $ DXMDT(lxd,lxd), DYMDT(lyd,lyd), DZMDT(lzd,lzd),
     $ wxmd (lxd),     wymd (lyd),     wzmd (lzd)

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      real zgx(lx1,2), zgy(ly1,2), zgz(lz1,2), WORK(lxd,lxd)

      integer px,py,pz
      integer i,j

      px = 1+(lxb-lx1)/2
      call zwgll(zgx(1,1),WORK,mx1)
      call zwgll(zgx(1,2),WORK,lx1)
      call igllm(JX1(1,px),WORK,zgx(1,1),zgx(1,2),mx1,lx1,lx1,lx1)
      if (lx1.ne.lxd) then
      call gen_int(JXD(1,px),WORK(1,1),lxd,mx1,WORK(1,2))
      else
      call dcopy(lx1*lx1,JX1(1,px),1,JXD(1,px),1)
      endif
      call mxm(DXM1,lx1,JX1(1,px),lx1,DX1(1,px),mx1)
      call mxm(DXMD,lxd,JXD(1,px),lxd,DXD(1,px),mx1)

      py = 1+(lyb-ly1)/2
      call zwgll(zgy(1,1),WORK,my1)
      call zwgll(zgy(1,2),WORK,ly1)
      call igllm(JY1(1,py),WORK,zgy(1,1),zgy(1,2),my1,ly1,ly1,ly1)
      if (ly1.ne.lyd) then
      call gen_int(JYD(1,py),WORK(1,1),lyd,my1,WORK(1,2))
      else
      call dcopy(ly1*ly1,JY1(1,py),1,JYD(1,py),1)
      endif
      call mxm(DYM1,ly1,JY1(1,py),ly1,DY1(1,py),my1)
      call mxm(DYMD,lyd,JYD(1,py),lyd,DYD(1,py),my1)

      if (mz1.gt.1) then
      pz = 1+(lzb-lz1)/2
      call zwgll(zgz(1,1),WORK,mz1)
      call zwgll(zgz(1,2),WORK,lz1)
      call igllm(JZ1(1,pz),WORK,zgz(1,1),zgz(1,2),mz1,lz1,lz1,lz1)
      if (lz1.ne.lzd) then
      call gen_int(JZD(1,pz),WORK(1,1),lzd,mz1,WORK(1,2))
      else
      call dcopy(lz1*lz1,JZ1(1,pz),1,JZD(1,pz),1)
      endif
      call mxm(DZM1,lz1,JZ1(1,pz),lz1,DZ1(1,pz),mz1)
      call mxm(DZMD,lzd,JZD(1,pz),lzd,DZD(1,pz),mz1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_ratio(h1,bcs)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /tens_ratio/ hnorm,vol,volvisc
      real hnorm(6,lelt), vol(6,lelt), volvisc(6,lelt)

      common /adfswapl/ l
      real l(lx1,ly1,lz1,lelt)

      common /adflf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt


      real h1(lx1,ly1,lz1,*)
      integer bcs(6,*)
      integer n,e,f

      real vp(lelt)
      real temp(lelt)
      real vlsum, vlsc2

      n = lx1*ly1*lz1

c     Ratio of length along face normal
      call adf_swap_lengths
      call adf_swap_lengths_fix
      do e=1,nelt
         hnorm(1,e) = llr(e)/lmr(e)
         hnorm(2,e) = lrr(e)/lmr(e)
         hnorm(3,e) = lls(e)/lms(e)
         hnorm(4,e) = lrs(e)/lms(e)
         hnorm(5,e) = llt(e)/lmt(e)
         hnorm(6,e) = lrt(e)/lmt(e)
      enddo ! e

c     Ratio of (int dV)
      do e=1,nelt
         vp(e) = vlsum(BM1(1,1,1,e),n)
         call cfill(l(1,1,1,e),vp(e),n)
      enddo ! e
      call dssum(l,lx1,ly1,lz1)
      do e=1,nelt
         vol(1,e) = max(l(  1,2,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         vol(2,e) = max(l(lx1,2,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         vol(3,e) = max(l(2,  1,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         vol(4,e) = max(l(2,ly1,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         vol(5,e) = max(l(2,2,  1,e)/vp(e)-1.0E0, 0.0E0)
         vol(6,e) = max(l(2,2,lz1,e)/vp(e)-1.0E0, 0.0E0)
         do f=1,2*ldim
            if(bcs(f,e).ne.0) vol(f,e) = 0.0E0 ! don't assemble
         enddo ! f
      enddo ! e

c     Ratio of (int visc dV)
      do e=1,nelt
         vp(e) = vlsc2(h1(1,1,1,e),BM1(1,1,1,e),n)
         call cfill(l(1,1,1,e),vp(e),n)
      enddo ! e
      call dssum(l,lx1,ly1,lz1)
      do e=1,nelt
         volvisc(1,e) = max(l(  1,2,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         volvisc(2,e) = max(l(lx1,2,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         volvisc(3,e) = max(l(2,  1,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         volvisc(4,e) = max(l(2,ly1,ldim-1,e)/vp(e)-1.0E0, 0.0E0)
         volvisc(5,e) = max(l(2,2,  1,e)/vp(e)-1.0E0, 0.0E0)
         volvisc(6,e) = max(l(2,2,lz1,e)/vp(e)-1.0E0, 0.0E0)
         do f=1,2*ldim
            if(bcs(f,e).ne.0) volvisc(f,e) = 0.0E0 ! don't assemble
         enddo ! f
      enddo ! e

c     Debug communication
      if (.false.) then
      ifxyo = .true.
      do e=1,nelt
         call cfill(l(1,1,1,e),vp(e),n)
      enddo ! e
      call outpost(l,l,l,pr,l,'vis')
      do f=1,2*ldim
         do e=1,nelt
            call cfill(l(1,1,1,e),volvisc(f,e)*vp(e),n)
         enddo ! e
         call outpost(l,l,l,pr,l,'vis')
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_omegabar(e,mx1,my1,mz1)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      integer e,mx1,my1,mz1

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      common /tens_ratio/ hnorm,vol,volvisc
      real hnorm(6,lelt), vol(6,lelt), volvisc(6,lelt)

      integer ns
      real scl, scr

      ns = (lxb-lx1)/2
      scl = sqrt(volvisc(1,e))
      scr = sqrt(volvisc(2,e))
      call omegabar_basis(mx1,ns,JX1,lx1,scl,scr)
      scl = scl/hnorm(1,e)
      scr = scr/hnorm(2,e)
      call omegabar_basis(mx1,ns,DX1,lx1,scl,scr)
      scl = sqrt(vol(1,e))
      scr = sqrt(vol(2,e))
      call omegabar_basis(mx1,ns,JXD,lxd,scl,scr)
      scl = scl/hnorm(1,e)
      scr = scr/hnorm(2,e)
      call omegabar_basis(mx1,ns,DXD,lxd,scl,scr)

      ns = (lyb-ly1)/2
      scl = sqrt(volvisc(3,e))
      scr = sqrt(volvisc(4,e))
      call omegabar_basis(my1,ns,JY1,ly1,scl,scr)
      scl = scl/hnorm(3,e)
      scr = scr/hnorm(4,e)
      call omegabar_basis(my1,ns,DY1,ly1,scl,scr)
      scl = sqrt(vol(3,e))
      scr = sqrt(vol(4,e))
      call omegabar_basis(my1,ns,JYD,lyd,scl,scr)
      scl = scl/hnorm(3,e)
      scr = scr/hnorm(4,e)
      call omegabar_basis(my1,ns,DYD,lyd,scl,scr)

      if (lz1.gt.1) then
      ns = (lzb-lz1)/2
      scl = sqrt(volvisc(5,e))
      scr = sqrt(volvisc(6,e))
      call omegabar_basis(mz1,ns,JZ1,lz1,scl,scr)
      scl = scl/hnorm(5,e)
      scr = scr/hnorm(6,e)
      call omegabar_basis(mz1,ns,DZ1,lz1,scl,scr)
      scl = sqrt(vol(5,e))
      scr = sqrt(vol(6,e))
      call omegabar_basis(mz1,ns,JZD,lzd,scl,scr)
      scl = scl/hnorm(5,e)
      scr = scr/hnorm(6,e)
      call omegabar_basis(mz1,ns,DZD,lzd,scl,scr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine omegabar_basis(n,ns,A,lda,scl,scr)
      implicit none
      integer n, ns, lda
      real A(lda,*)
      real scl, scr

      if(ns.gt.0) then
         call cmult2(A(1,1     ), A(1,n+1 ), scl, ns*lda)
         call cmult2(A(1,n+ns+1), A(1,ns+1), scr, ns*lda)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine tens_gen_crs(nxc)
c     Generate stiffness and mass on coarse grid
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'TOTAL'

      common /tens_crs/ acrs,bcrs 
      real acrs(lcr,lcr,lelt),bcrs(lcr,lcr,lelt)

      real k0(lx1,ly1,lz1,lelt), k1(lx1,ly1,lz1,lelt)
      real w1(lx1,ly1,lz1,lelt), w2(lx1,ly1,lz1,lelt)

      integer nxc, ncr, n
      logical ifconv

      ncr = nxc**ldim
      n = lx1*ly1*lz1*nelt

      ifconv = .false.
      call rzero(k0,n)
      call rone (k1,n)
      call get_local_crs_galerkin_t(acrs,ncr,nxc,k1,k0,k0,ifconv,w1,w2)
      call get_local_crs_galerkin_t(bcrs,ncr,nxc,k0,k1,k0,ifconv,w1,w2)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_old(e,cr,cs,ct,vx,vy,vz)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'WZ'
      
      integer e
      real cr(lx1,ly1,lz1), 
     $     cs(lx1,ly1,lz1), 
     $     ct(lx1,ly1,lz1),
     $     vx(lx1,ly1,lz1,*),
     $     vy(lx1,ly1,lz1,*),
     $     vz(lx1,ly1,lz1,*)

      integer i,j,k
      real w,wx,wy,wz,wyz

      if (nz1.gt.1) then

      do k=1,nz1
      do j=1,ny1
      wyz = wym1(i)*wzm1(k)
      do i=1,nx1
        w = wxm1(i)*wyz
        wx = vx(i,j,k,e)
        wy = vy(i,j,k,e)
        wz = vz(i,j,k,e)
        cr(i,j,k)=w*(rxm1(i,j,k,e)*wx+rym1(i,j,k,e)*wy+rzm1(i,j,k,e)*wz)
        cs(i,j,k)=w*(sxm1(i,j,k,e)*wx+sym1(i,j,k,e)*wy+szm1(i,j,k,e)*wz)
        ct(i,j,k)=w*(txm1(i,j,k,e)*wx+tym1(i,j,k,e)*wy+tzm1(i,j,k,e)*wz)
      enddo ! i
      enddo ! j
      enddo ! k

      else ! 2D

      do j=1,ny1
      wyz = wym1(i)
      do i=1,nx1
        w = wxm1(i)*wyz
        wx = vx(i,j,1,e)
        wy = vy(i,j,1,e)
        cr(i,j,1)=w*(rxm1(i,j,1,e)*wx+rym1(i,j,1,e)*wy)
        cs(i,j,1)=w*(sxm1(i,j,1,e)*wx+sym1(i,j,1,e)*wy)
      enddo ! i
      enddo ! j

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_shuff2d(d, my,mx, Y,X,S, WORK)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      common /tens_coef/B00,C00,C11,C22,C33,C12,C13,C23,C01,C02,C03
      real B00(lx1,ly1,lz1), C00(lx1,ly1,lz1), 
     $     C11(lx1,ly1,lz1), C22(lx1,ly1,lz1), C33(lx1,ly1,lz1),
     $     C12(lx1,ly1,lz1), C13(lx1,ly1,lz1), C23(lx1,ly1,lz1),
     $     C01(lxd,lyd,lzd), C02(lxd,lyd,lzd), C03(lxd,lyd,lzd)

      common /tens_flag/ifdealias
      logical ifdealias

      integer d, my,mx
      real Y(*), X(*), S, WORK(*)
      integer ny,nx

      nx = lx1
      ny = ly1
c     Symmetric terms (lx1 x ly1)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, JY1,JX1,C00,JY1,JX1, WORK)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, JY1,DX1,C11,JY1,DX1, WORK)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, DY1,JX1,C22,DY1,JX1, WORK)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, DY1,JX1,C12,JY1,DX1, WORK)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, JY1,DX1,C12,DY1,JX1, WORK)

      if(ifdealias) then
         nx = lxd
         ny = lyd
      endif
c     Skew-symmetric terms (lxd x lyd)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, JYD,JXD,C01,JYD,DXD, WORK)
      call shuff2d(d,my,mx,ny,nx, Y,X,S, JYD,JXD,C02,DYD,JXD, WORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_shuff(d, mz,my,mx, Z,Y,X,S, WORK)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      common /tens_coef/B00,C00,C11,C22,C33,C12,C13,C23,C01,C02,C03
      real B00(lx1,ly1,lz1), C00(lx1,ly1,lz1), 
     $     C11(lx1,ly1,lz1), C22(lx1,ly1,lz1), C33(lx1,ly1,lz1),
     $     C12(lx1,ly1,lz1), C13(lx1,ly1,lz1), C23(lx1,ly1,lz1),
     $     C01(lxd,lyd,lzd), C02(lxd,lyd,lzd), C03(lxd,lyd,lzd)

      common /tens_flag/ifdealias
      logical ifdealias

      integer d, mz,my,mx
      real Z(*), Y(*), X(*), S, WORK(*)
      integer nz,ny,nx

      nx = lx1
      ny = ly1
      nz = lz1
c     Symmetric terms (lx1 x ly1 x lz1)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S, 
     $     JZ1,JY1,JX1,C00,JZ1,JY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S, 
     $     JZ1,JY1,DX1,C11,JZ1,JY1,DX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,JY1,DX1,C12,JZ1,DY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,JY1,DX1,C13,DZ1,JY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,DY1,JX1,C12,JZ1,JY1,DX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,DY1,JX1,C22,JZ1,DY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,DY1,JX1,C23,DZ1,JY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     DZ1,JY1,JX1,C13,JZ1,JY1,DX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     DZ1,JY1,JX1,C23,JZ1,DY1,JX1, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     DZ1,JY1,JX1,C33,DZ1,JY1,JX1, WORK)

      if(ifdealias) then
         nx = lxd
         ny = lyd
         nz = lzd
      endif
c     Skew-symmetric terms (lxd x lyd x lzd)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZD,JYD,JXD,C01,JZD,JYD,DXD, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S, 
     $     JZD,JYD,JXD,C02,JZD,DYD,JXD, WORK)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S, 
     $     JZD,JYD,JXD,C03,DZD,JYD,JXD, WORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine mass_shuff(d, mz,my,mx, Z,Y,X,S, WORK)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      common /tens_coef/B00,C00,C11,C22,C33,C12,C13,C23,C01,C02,C03
      real B00(lx1,ly1,lz1), C00(lx1,ly1,lz1), 
     $     C11(lx1,ly1,lz1), C22(lx1,ly1,lz1), C33(lx1,ly1,lz1),
     $     C12(lx1,ly1,lz1), C13(lx1,ly1,lz1), C23(lx1,ly1,lz1),
     $     C01(lxd,lyd,lzd), C02(lxd,lyd,lzd), C03(lxd,lyd,lzd)

      integer d, mz,my,mx
      real Z(*), Y(*), X(*), S, WORK(*)
      integer nz,ny,nx

      nx = lx1
      ny = ly1
      nz = lz1
c     Mass term (lx1 x ly1 x lz1)
      call shuff(d, mz,my,mx, nz,ny,nx, Z,Y,X,S,
     $     JZ1,JY1,JX1,B00,JZ1,JY1,JX1, WORK)
      return
      end
c-----------------------------------------------------------------------
      subroutine shuff2d(d,my,mx,ny,nx, Y,X,S, B2,B1,C,A2,A1, WORK)
      implicit none
      integer d, my,mx, ny,nx 

      real Y(*), X(*), S ! X(m,m)
      real B2(*), B1(*)  ! B(n,m)
      real C(*)          ! C(n,n)
      real A2(*), A1(*)  ! A(n,m)
      real WORK(*) ! WORK(nx+1,ny+1)

      integer p1,p2,pw

      real one
      parameter (one=1.0E0)

      p1 = 1
      p2 = p1 + max(mx,nx)
      pw = p2 + max(my,ny)

      if (d.eq.1) then
c     X = X + S*B1'*diag(C*diag(B2*Y*A2'))*A1
c     w2 = diag(B2*Y*A2')
      call diag_outerp(ny,my,B2,ny,Y,my,A2,ny,WORK(p2),WORK(pw))
c     w1 = C*w2
      call mxm(C, nx, WORK(p2), ny, WORK(p1), 1)
c     X = X + S*B1'*diag(w1)*A1
      call innerp_diag(nx,mx,S,B1,nx,WORK(p1),A1,nx,one,X,mx,WORK(pw))
      elseif (d.eq.2) then
c     Y = Y + S*B2'*diag(C'*diag(B1*X*A1'))*A2
c     w1 = diag(B1*X*A1')
      call diag_outerp(nx,mx,B1,nx,X,mx,A1,nx,WORK(p1),WORK(pw))
c     w2' = w1'*C
      call mxm(WORK(p1), 1, C, nx, WORK(p2), ny)
c     Y = Y + S*B2'*diag(w2)*A2
      call innerp_diag(ny,my,S,B2,ny,WORK(p2),A2,ny,one,Y,my,WORK(pw))
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine shuff(d, mz,my,mx,    nz,ny,nx, Z,Y,X,S,
     $                    B3,B2,B1, C, A3,A2,A1, WORK)
      implicit none
      integer d, mz,my,mx, nz,ny,nx

      real X(*), Y(*), Z(*), S ! X(m,m)
      real B3(*), B2(*), B1(*) ! B(n,m)
      real C(*)                ! C(n,n)
      real A3(*), A2(*), A1(*) ! A(n,m)
      real WORK(*) ! WORK(nx+1,ny+1,nz+1)
      
      integer p1,p2,p3,pw,pu
      
      real one
      parameter (one=1.0E0)

      p1 = 1
      p2 = p1 + max(mx,nx)
      p3 = p2 + max(my,ny)
      pw = p3 + max(mz,nz)
      pu = pw + max(my,ny)*max(mz,nz)

      if (d.eq.0) then
c     w1 = diag(B1*X*A1')
c     w2 = diag(B2*Y*A2')
c     w3 = diag(B3*Z*A3')
      call diag_outerp(nx,mx,B1,nx,X,mx,A1,nx,WORK(p1),WORK(pw))
      call diag_outerp(ny,my,B2,ny,Y,my,A2,ny,WORK(p2),WORK(pw))
      call diag_outerp(nz,mz,B3,nz,Z,mz,A3,nz,WORK(p3),WORK(pw))
c     S = S + kron(w3,w2,w1)'*C(:)
c     S = S + w1'*( kron(I,w2,I)'*C(:) )*w3
      call kron_mxv1(1, nx,ny,nz, WORK(p1), C,        WORK(pw))
      call kron_mxv2(1, 1 ,ny,nz, WORK(p2), WORK(pw), WORK(pu))
      call kron_mxv3(1, 1 ,1 ,nz, WORK(p3), WORK(pu), WORK(pw))
      S = S + WORK(pw)
      elseif (d.eq.1) then
c     w2 = diag(B2*Y*A2')
c     w3 = diag(B3*Z*A3')
      call diag_outerp(ny,my,B2,ny,Y,my,A2,ny,WORK(p2),WORK(pw))
      call diag_outerp(nz,mz,B3,nz,Z,mz,A3,nz,WORK(p3),WORK(pw))
c     w1 = kron(w3,w2,I)'*C(:)
      call kron_mxv2(1, nx,ny,nz, WORK(p2), C,        WORK(pw))
      call kron_mxv3(1, nx,1 ,nz, WORK(p3), WORK(pw), WORK(p1))
c     X = X + S*B1'*diag(w1)*A1
      call innerp_diag(nx,mx,S,B1,nx,WORK(p1),A1,nx,one,X,mx,WORK(pw))
      elseif (d.eq.2) then
c     w1 = diag(B1*X*A1')
c     w3 = diag(B3*Z*A3')
      call diag_outerp(nx,mx,B1,nx,X,mx,A1,nx,WORK(p1),WORK(pw))
      call diag_outerp(nz,mz,B3,nz,Z,mz,A3,nz,WORK(p3),WORK(pw))
c     w2 = kron(w3,I,w1)'*C(:)
      call kron_mxv1(1, nx,ny,nz, WORK(p1), C,        WORK(pw))
      call kron_mxv3(1, 1, ny,nz, WORK(p3), WORK(pw), WORK(p2))
c     Y = Y + S*B2'*diag(w2)*A2
      call innerp_diag(ny,my,S,B2,ny,WORK(p2),A2,ny,one,Y,my,WORK(pw))
      elseif (d.eq.3) then
c     w1 = diag(B1*X*A1')
c     w2 = diag(B2*Y*A2')
      call diag_outerp(nx,mx,B1,nx,X,mx,A1,nx,WORK(p1),WORK(pw))
      call diag_outerp(ny,my,B2,ny,Y,my,A2,ny,WORK(p2),WORK(pw))
c     v1 = kron(w3,I,w1)'*C(:)
      call kron_mxv1(1, nx,ny,nz, WORK(p1), C,        WORK(pw))
      call kron_mxv2(1, 1 ,ny,nz, WORK(p2), WORK(pw), WORK(p3))
c     Z = Z + S*B3'*diag(w3)*A3
      call innerp_diag(nz,mz,S,B3,nz,WORK(p3),A3,nz,one,Z,mz,WORK(pw))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tens_set_ptr()
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)
      integer qtr(3,4,11)
      data qtr/ 1,1,1, 0,0,0, 0,0,0, 0,0,0 ! 111 mass  
     $        , 2,1,1, 1,2,1, 0,0,0, 0,0,0 ! 122 mass  
     $        , 2,1,1, 1,1,2, 0,0,0, 0,0,0 ! 212 mass  
     $        , 1,2,1, 1,1,2, 0,0,0, 0,0,0 ! 221 mass  
     $        , 2,2,1, 2,1,2, 1,2,2, 0,0,0 ! 222 stiff
     $        , 3,2,1, 2,1,2, 1,2,2, 0,0,0 ! 223 stiff
     $        , 2,3,1, 2,1,2, 1,2,2, 0,0,0 ! 232 stiff
     $        , 2,2,1, 1,2,2, 2,1,3, 0,0,0 ! 322 stiff
     $        , 3,2,1, 2,3,1, 2,1,2, 1,2,2 ! 233 stiff
     $        , 2,2,1, 3,1,2, 1,2,2, 2,1,3 ! 323 stiff
     $        , 2,2,1, 2,1,2, 1,3,2, 1,2,3 ! 332 stiff
     $        /
      
      character*3 zyx
      call tens_get_rank(zyx)
      if ( (zyx.eq.'111').or.(zyx.eq.'112').or.
     $     (zyx.eq.'121').or.(zyx.eq.'211').or. 
     $     (zyx.eq.'113').or.(zyx.eq.'131').or. 
     $     (zyx.eq.'311') ) then ! Degenerate cases
         xyzrank(1) = 1
         xyzrank(2) = 1
         xyzrank(3) = 1
         terms = 1
         call icopy(ptr,qtr(1,1,1),12)
      elseif (zyx.eq.'122') then
         terms = 2
         call icopy(ptr,qtr(1,1,2),12)
      elseif (zyx.eq.'212') then
         terms = 2
         call icopy(ptr,qtr(1,1,3),12)
      elseif (zyx.eq.'221') then
         terms = 2
         call icopy(ptr,qtr(1,1,4),12)
      elseif (zyx.eq.'222') then
         terms = 3
         call icopy(ptr,qtr(1,1,5),12)
      elseif (zyx.eq.'223') then
         terms = 3
         call icopy(ptr,qtr(1,1,6),12)
      elseif (zyx.eq.'232') then
         terms = 3
         call icopy(ptr,qtr(1,1,7),12)
      elseif (zyx.eq.'322') then
         terms = 3
         call icopy(ptr,qtr(1,1,8),12)
      elseif (zyx.eq.'233') then
         terms = 4
         call icopy(ptr,qtr(1,1,9),12)
      elseif (zyx.eq.'323') then
         terms = 4
         call icopy(ptr,qtr(1,1,10),12)
      elseif (zyx.eq.'332') then
         terms = 4
         call icopy(ptr,qtr(1,1,11),12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tens_set_rank(z,y,x)
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)
      integer x,y,z
      xyzrank(1) = x
      xyzrank(2) = y
      xyzrank(3) = z
      call tens_set_ptr()
      return
      end
c-----------------------------------------------------------------------
      subroutine tens_get_rank(arank)
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)
            
      character*3 arank
      integer x,y,z
      x = xyzrank(1)
      y = xyzrank(2)
      z = xyzrank(3)
      write(arank,'(I1,I1,I1)')z,y,x
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_als(e,arank,mx,my,mz,A,lda,B,ldb,C,ldc,
     $                           h1,h2,h3,bcs,W)
c     Driver routine for tensor decomposition
c     mx,my,mz = grid sizes
c     lda,ldb,ldc = basis functions on omegabar

c     Initial guess recovers exactly separable cases in 1 ALS iteration.
c     
c     From geometric arguments, 2nd and 3rd ranks are present in mass.
c     The shuffled opertor acting on the tensor complement of the 1st
c     rank results in a linear combination of all the ranks.

      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'HSMGL'
      include 'TOTAL'

      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      common /tens_basis/JX1,DX1,JXD,DXD,JY1,DY1,JYD,DYD,JZ1,DZ1,JZD,DZD
      real
     $ JX1(lx1,lxb),DX1(lx1,lxb),JXD(lxd,lxb),DXD(lxd,lxb),
     $ JY1(ly1,lyb),DY1(ly1,lyb),JYD(lyd,lyb),DYD(lyd,lyb),
     $ JZ1(lz1,lzb),DZ1(lz1,lzb),JZD(lzd,lzb),DZD(lzd,lzb)

      common /tens_coef/B00,C00,C11,C22,C33,C12,C13,C23,C01,C02,C03
      real B00(lx1,ly1,lz1), C00(lx1,ly1,lz1), 
     $     C11(lx1,ly1,lz1), C22(lx1,ly1,lz1), C33(lx1,ly1,lz1),
     $     C12(lx1,ly1,lz1), C13(lx1,ly1,lz1), C23(lx1,ly1,lz1),
     $     C01(lxd,lyd,lzd), C02(lxd,lyd,lzd), C03(lxd,lyd,lzd)
      
      common /tens_crs/ acrs,bcrs
      real acrs(lcr,lcr,lelt), bcrs(lcr,lcr,lelt)

      integer e, mx,my,mz, lda,ldb,ldc
      character*3 arank,brank
      real A(lda*lda,*), B(ldb*ldb,*), C(ldc*ldc,*)
      real h1(*), h2(*), h3(*), W(*)
      integer bcs(*)
      real S(4)

      real tol
      integer maxit, d,i,j,k, nc,nx,ny,nz
      integer mxs,mys,mzs

      real G1(lx1),G2(ly1),G3(lz1)

      logical ifconv, ifanis
      integer n1

      real zero,one
      parameter (zero=0.0E0, one=1.0E0)

c     ALS parameters
      maxit = 20
      tol = 1.0E-16

      nc = lxc
      nx = lda*lda
      ny = ldb*ldb
      nz = ldc*ldc

      call rone(S,4)

c     Element setup
      ifconv = .true.
      ifanis = .true.
      n1 = lx1*ly1*lz1
      call tens_gen_coef(e, h1,h2,h3, ifconv,ifanis)
      call tens_gen_omegabar(e,mx,my,mz)

c     Compute mass rank
      call crs_rank(2,nc,bcrs(1,1,e))
      call tens_get_rank(brank)

c     Naive initial guess for 2nd rank
      call rzero(G1,lx1)
      call rzero(G2,ly1)
      call rzero(G3,lz1)
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         G1(i) = G1(i) + 0.25*G1M1(i,j,k,e)
         G2(j) = G2(j) + 0.25*G2M1(i,j,k,e)
         G3(k) = G3(k) + 0.25*G3M1(i,j,k,e)
      enddo ! i
      enddo ! j
      enddo ! k
      call innerp_diag(lx1,lda,one,JX1,lx1,G1,JX1,lx1,zero,A(1,2),lda,W)
      call innerp_diag(ly1,ldb,one,JY1,ly1,G2,JY1,ly1,zero,B(1,2),ldb,W)
      call innerp_diag(lz1,ldc,one,JZ1,lz1,G3,JZ1,lz1,zero,C(1,2),ldc,W)

c     Uncoupled mass ALS guess for 3rd rank
      call rzero(A(1,3),nx)
      call rzero(B(1,3),ny)
      call rzero(C(1,3),nz)
      do d=1,3
      do k=1,terms
      if ((1+ptr(d,k)).eq.3) then
         call mass_shuff(d,mz,my,mx,
     $   C(1,1+ptr(3,k)),B(1,1+ptr(2,k)),A(1,1+ptr(1,k)),S(k),W)
      endif
      enddo ! k
      enddo ! d      

c     Perform ALS on mass
      if(xyzrank(1).eq.1) call copy(A(1,1),A(1,2),nx)
      if(xyzrank(2).eq.1) call copy(B(1,1),B(1,2),ny)
      if(xyzrank(3).eq.1) call copy(C(1,1),C(1,2),nz)
      call mass_als(lda,ldb,ldc,A(1,2),B(1,2),C(1,2),maxit,tol,W)
      if(xyzrank(1).eq.1) call copy(A(1,3),A(1,1),nx)
      if(xyzrank(2).eq.1) call copy(B(1,3),B(1,1),ny)
      if(xyzrank(3).eq.1) call copy(C(1,3),C(1,1),nz)

c     Compute stiff rank
      call crs_rank(3,nc,acrs(1,1,e))
      call tens_get_rank(arank)

c     Uncoupled stiff ALS guess for 1st rank
      call rzero(A(1,1),nx)
      call rzero(B(1,1),ny)
      call rzero(C(1,1),nz)
      do d=1,3
      do k=1,terms
      if (ptr(d,k).eq.1) then
         call adf_shuff(d,lda,ldb,ldc,
     $   C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),W)
      endif
      enddo ! k
      enddo ! d

c     Perform ALS on stiff
      call stiff_als(lda,ldb,ldc,A,B,C,maxit,tol,W)

c     Assemble 1D tensors
      mxs = (lda-mx)/2
      mys = (ldb-my)/2
      mzs = (ldc-mz)/2
      call copy(A(1,4), A(1,1), 3*nx)
      call copy(B(1,4), B(1,1), 3*ny)
      call copy(C(1,4), C(1,1), 3*nz)
      call schwarz1d(mx,mxs,3, A(1,4),lda, A(1,1),max(mx,lda-2),bcs(1))
      call schwarz1d(my,mys,3, B(1,4),ldb, B(1,1),max(my,ldb-2),bcs(3))
      call schwarz1d(mz,mzs,3, C(1,4),ldc, C(1,1),max(mz,ldc-2),bcs(5))
      return
      end
c-----------------------------------------------------------------------
      subroutine mass_als(mx,my,mz,A,B,C,maxit,tol,WORK)
c     Alternating Least Squares Procedure for separable
c     Mass matrix (maximum rank is 2)
      implicit none
      include 'SIZE'
      include 'TOTAL'
   
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer m
      parameter (m=2)

      integer mx,my,mz, maxit
      real tol
      real A(mx*mx,*), B(my*my,*), C(mz*mz,*)
      real WORK(*)

      character*3 arank
      integer nx,ny,nz, i,j,k, it, info
      real fn1,fn2, rel
      real S(4), ST(4), nrm(4)
      real F1(m,m,3), F2(m,m,3), IFRO(m,m)

      real one, zero
      parameter (one=1.0E0, zero=0.0E0)

      real normfro

      nx = mx*mx
      ny = my*my
      nz = mz*mz

c     Normalize tensors A,B,C
      call rone(S,terms)
      call normc(nx,m,A,nrm)
      call normc(ny,m,B,nrm)
      call normc(nz,m,C,nrm)

c     Transpose is not redudant, it also stores previous iteration
      call transpose(A(1,m+1),m,A,nx)
      call transpose(B(1,m+1),m,B,ny)
      call transpose(C(1,m+1),m,C,nz)
c     Compute Frobenius inner products
      call mxm(A(1,m+1),m,A,nx,F1(1,1,1),m)
      call mxm(B(1,m+1),m,B,ny,F1(1,1,2),m)
      call mxm(C(1,m+1),m,C,nz,F1(1,1,3),m)
      fn1 = normfro(m,S,S,F1)

      it = 0
      rel = one
      do while ( (rel.gt.tol).and.(it.lt.maxit) )
c     ALS iteration
      call copy(ST,S,terms)

c     Minimize on C
      call rzero(C,nz*m)
      call normc(terms,1,S,nrm)
      call compfro(3,m,IFRO,F1,S)
      do k=1,terms
      call mass_shuff(3,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(nz,m,C,IFRO)
      call normc(nz,m,C,nrm)
      call vxptr(3,S,nrm)
      call mxm(C(1,m+1),m,C,nz,F2(1,1,3),m)
      call transpose(C(1,m+1),m,C,nz)
      call mxm(C(1,m+1),m,C,nz,F1(1,1,3),m)

c     Minimize on A
      call rzero(A,nx*m)
      call normc(terms,1,S,nrm)
      call compfro(1,m,IFRO,F1,S)
      do k=1,terms
      call mass_shuff(1,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(nx,m,A,IFRO)
      call normc(nx,m,A,nrm)
      call vxptr(1,S,nrm)
      call mxm(A(1,m+1),m,A,nx,F2(1,1,1),m)
      call transpose(A(1,m+1),m,A,nx)
      call mxm(A(1,m+1),m,A,nx,F1(1,1,1),m)
      
c     Minimize on B
      call rzero(B,ny*m)
      call normc(terms,1,S,nrm)
      call compfro(2,m,IFRO,F1,S)
      do k=1,terms
      call mass_shuff(2,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(ny,m,B,IFRO)
      call normc(ny,m,B,nrm)
      call vxptr(2,S,nrm)
      call mxm(B(1,m+1),m,B,ny,F2(1,1,2),m)
      call transpose(B(1,m+1),m,B,nz)
      call mxm(B(1,m+1),m,B,ny,F1(1,1,2),m)

c     Compute relative difference
      fn2 = fn1
      fn1 = normfro(m,S,S,F1)
      rel = (fn1-2*normfro(m,ST,S,F2)+fn2)/fn1
      it = it+1
      enddo ! it

c     if(nio.eq.0) then
c        call tens_get_rank(arank)
c        write(6,*) it,rel,'   ',arank,'   alsb'
c     endif
      call renorm(nx,ny,nz, S, A,B,C)
      return
      end
c-----------------------------------------------------------------------
      subroutine stiff_als(mx,my,mz,A,B,C,maxit,tol,WORK)
c     Alternating Least Squares Procedure for separable
c     Advection-difussion stiffness matrix (maximum rank is 3)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer m
      parameter (m=3)

      integer mx,my,mz, maxit
      real tol
      real A(mx*mx,*), B(my*my,*), C(mz*mz,*)
      real WORK(*)

      character*3 arank
      integer nx,ny,nz, i,j,k, it, info
      real fn1,fn2, rel
      real S(4), ST(4), nrm(4)
      real F1(m,m,3), F2(m,m,3), IFRO(m,m)

      real one, zero
      parameter (one=1.0E0, zero=0.0E0)

      real normfro

      nx = mx*mx
      ny = my*my
      nz = mz*mz

c     Normalize tensors A,B,C
      call rone(S,terms)
      call normc(nx,m,A,nrm)
      call normc(ny,m,B,nrm)
      call normc(nz,m,C,nrm)

c     Transpose is not redudant, it also stores previous iteration
      call transpose(A(1,m+1),m,A,nx)
      call transpose(B(1,m+1),m,B,ny)
      call transpose(C(1,m+1),m,C,nz)
c     Compute Frobenius inner products
      call mxm(A(1,m+1),m,A,nx,F1(1,1,1),m)
      call mxm(B(1,m+1),m,B,ny,F1(1,1,2),m)
      call mxm(C(1,m+1),m,C,nz,F1(1,1,3),m)
      fn1 = normfro(m,S,S,F1)

      it = 0
      rel = one
      do while ( (rel.gt.tol).and.(it.lt.maxit) )
c     ALS iteration
      call copy(ST,S,terms)

c     Minimize on C
      call rzero(C,nz*m)
      call normc(terms,1,S,nrm)
      call compfro(3,m,IFRO,F1,S)
      do k=1,terms
      call adf_shuff(3,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(nz,m,C,IFRO)
      call normc(nz,m,C,nrm)
      call vxptr(3,S,nrm)
      call mxm(C(1,m+1),m,C,nz,F2(1,1,3),m)
      call transpose(C(1,m+1),m,C,nz)
      call mxm(C(1,m+1),m,C,nz,F1(1,1,3),m)

c     Minimize on A
      call rzero(A,nx*m)
      call normc(terms,1,S,nrm)
      call compfro(1,m,IFRO,F1,S)
      do k=1,terms
      call adf_shuff(1,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(nx,m,A,IFRO)
      call normc(nx,m,A,nrm)
      call vxptr(1,S,nrm)
      call mxm(A(1,m+1),m,A,nx,F2(1,1,1),m)
      call transpose(A(1,m+1),m,A,nx)
      call mxm(A(1,m+1),m,A,nx,F1(1,1,1),m)
      
c     Minimize on B
      call rzero(B,ny*m)
      call normc(terms,1,S,nrm)
      call compfro(2,m,IFRO,F1,S)
      do k=1,terms
      call adf_shuff(2,mz,my,mx,
     $     C(1,ptr(3,k)),B(1,ptr(2,k)),A(1,ptr(1,k)),S(k),WORK)
      enddo ! k
      call xasolv(ny,m,B,IFRO)
      call normc(ny,m,B,nrm)
      call vxptr(2,S,nrm)
      call mxm(B(1,m+1),m,B,ny,F2(1,1,2),m)
      call transpose(B(1,m+1),m,B,nz)
      call mxm(B(1,m+1),m,B,ny,F1(1,1,2),m)

c     Compute relative difference
      fn2 = fn1
      fn1 = normfro(m,S,S,F1)
      rel = (fn1-2*normfro(m,ST,S,F2)+fn2)/fn1
      it = it+1
      enddo ! it

c     if(nio.eq.0) then
c        call tens_get_rank(arank)
c        write(6,*) it,rel,'   ',arank,'   alsa'
c     endif
      call renorm(nx,ny,nz, S, A,B,C)
      return
      end
c-----------------------------------------------------------------------
      subroutine vxptr(d,a,b)
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer d
      real a(*), b(*)
      integer k

      do k=1,terms
         a(k) = a(k)*b(ptr(d,k))
      enddo ! k
      return
      end
c-----------------------------------------------------------------------
      subroutine compfro(d,m,IFRO,FRO,S)
c     Cholesky factorization of the normal equation Gramian
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer d, m
      real IFRO(m,m), FRO(m,m,3), S(*)
      integer i,j, d1,d2, r, info
      
      r = 1
      d1 = mod(d,3)+1
      d2 = mod(d+1,3)+1
      call rzero(IFRO,m*m)
      do j=1,terms
      do i=1,terms
         IFRO(ptr(d,i),ptr(d,j)) = IFRO(ptr(d,i),ptr(d,j)) +
     $ S(i)*S(j)*FRO(ptr(d1,i),ptr(d1,j),d1)*FRO(ptr(d2,i),ptr(d2,j),d2)
      enddo ! i
      r = max(r,ptr(d,j))
      enddo ! j
      do j=r+1,m
         IFRO(j,j) = 1.0E0
      enddo ! j
      call dpotrf('L',m,IFRO,m,info)
      return
      end
c-----------------------------------------------------------------------
      function normfro(m,S1,S2,F)
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      real normfro
      integer m
      real S1(*),S2(*),F(m,m,3)
      integer i,j

      normfro = 0.0E0
      do j=1,terms
      do i=1,terms
          normfro = normfro +  S1(i)*S2(j)*F(ptr(1,i),ptr(1,j),1)*
     $              F(ptr(2,i),ptr(2,j),2)*F(ptr(3,i),ptr(3,j),3)
      enddo ! i
      enddo ! j
      return
      end
c-----------------------------------------------------------------------
      subroutine renorm(nx,ny,nz, S, A,B,C)
      implicit none
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer nx,ny,nz
      real S(*), A(nx,*), B(ny,*), C(nz,*)

      real sca
      logical btr(3,4)
      integer i,j,k,n

c     Check if ptr(k,i) is unique inside a row
c     i.e. find which A(:,i), B(:,j), C(:,k) are unique
      do k=1,3
      do i=1,terms
         btr(k,i) = .true.
         do j=1,terms
            if((ptr(k,i).eq.ptr(k,j)).and.(i.ne.j)) btr(k,i)=.false.
         enddo ! j
      enddo ! i
      enddo ! k

c     Let unique factors absorb normalization equitatively
      do i=1,terms
         n = 0
         do k=1,3
            if(btr(k,i)) n = n+1
         enddo ! k

         sca = 1.0E0
         if (n.eq.1) sca = S(i)
         if (n.eq.2) sca = sqrt(S(i))
         if (n.gt.2) sca = S(i)**(1.0/n)

         if(btr(1,i)) call dscal(nx,sca,A(1,ptr(1,i)),1)
         if(btr(2,i)) call dscal(ny,sca,B(1,ptr(2,i)),1)
         if(btr(3,i)) call dscal(nz,sca,C(1,ptr(3,i)),1)
      enddo ! i
      return
      end
c-----------------------------------------------------------------------
      subroutine crs_rank(ranmax,nxc,acrs)
      implicit none
      include 'SIZE'
      include 'DOMAIN'
      include 'PARALLEL'
      
      integer m,n,lw
      parameter (m=lxc**2, n=lxc**4)
      parameter (lw=8*m*n)

      integer ranmax, nxc
      real acrs(nxc**ldim,*)

      real ashuff(lcr**2)
      real sig(m),U(1),V(1),W(lw)
      real kap(ldim), kapmax,kapmin, tol
      integer d,ran(3),info

      tol = 1.0E-2
      kapmax = 0.0E0
      kapmin = 2.0E0

      do d=1,ldim
         call full_shuff(d,nxc,nxc,nxc,ashuff,acrs)
         call dgesvd('N','N',m,n,ashuff,m,sig,U,1,V,1,W,lw,info)
         kap(d) = sig(ranmax)
         if(kap(d).gt.kapmax) kapmax = kap(d)
         if(kap(d).le.kapmin) kapmin = kap(d)
      enddo ! d

      if(kapmax.lt.1.E-10) kapmax=1.0E0

      do d=1,ldim
         kap(d)=(kap(d)-kapmin)/(kapmax-kapmin)
         ran(d) = ranmax-1
         if(kap(d).ge.tol) ran(d) = ranmax
      enddo ! d

      call tens_set_rank(ran(3),ran(2),ran(1))
c     if(nio.le.0) write(6,*) info,'irank',kap(1),kap(2),kap(3)
      return
      end
c-----------------------------------------------------------------------
      subroutine full_shuff(d, nx,ny,nz, ashuff, afull)
      implicit none
      integer d,nx,ny,nz
      real ashuff((nx*ny*nz)**2), afull(nx,ny,nz, nx,ny,nz)

      integer i,j,k, l,m,n, p

      p = 1
      if     (d.eq.1) then
         do n=1,nz
         do m=1,ny
         do k=1,nz
         do j=1,ny
         do l=1,nx
         do i=1,nx
            ashuff(p) = afull(i,j,k,l,m,n)
            p = p+1
         enddo ! i
         enddo ! l
         enddo ! j
         enddo ! k
         enddo ! m
         enddo ! n
      elseif (d.eq.2) then
         do n=1,nz
         do l=1,nx
         do k=1,nz
         do i=1,nx
         do m=1,ny
         do j=1,ny
            ashuff(p) = afull(i,j,k,l,m,n)
            p = p+1
         enddo ! j
         enddo ! m
         enddo ! i
         enddo ! k
         enddo ! l
         enddo ! n
      elseif (d.eq.3) then
         do m=1,ny
         do l=1,nx
         do j=1,ny
         do i=1,nx
         do n=1,nz
         do k=1,nz
            ashuff(p) = afull(i,j,k,l,m,n)
            p = p+1
         enddo ! k
         enddo ! n
         enddo ! i
         enddo ! j
         enddo ! l
         enddo ! m
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_svd(e,mx,my,A,lda,B,ldb,h1,h2,h3,bcs,WORK)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      integer e, mx,my, lda,ldb
      real A(lda*lda,*), B(ldb*ldb,*)
      real h1(*), h2(*), h3(*), WORK(*)
      integer bcs(*)

      logical ifconv
      logical ifanis

      integer mxs, mys

      ifconv = .true.
      ifanis = .true.

      mxs = (lxb-lx1)/2
      mys = (lyb-ly1)/2

c     Element setup
      call tens_gen_coef(e, h1,h2,h3, ifconv,ifanis)
      call tens_gen_omegabar(e,mx,my,1)

      call adf_svd_arpack(lda,ldb,A(1,4),B(1,4),WORK)

      call schwarz1d(mx,mxs,2, A(1,4),lda, A(1,1),max(mx,lda-2),bcs(1))
      call schwarz1d(my,mys,2, B(1,4),ldb, B(1,1),max(my,ldb-2),bcs(3))
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_svd_arpack(mx,my,A,B,works)
      implicit none
      include 'SIZE'
      include 'HSMGL'

      integer mx,my
      real A(mx*mx,*), B(my*my,*)
      real works(*)
      real scal

c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | It is assumed that A is M by N with M .ge. N.        |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | M .le. MAXM  and  N .le. MAXN                        |
c     |                                                      |
c     | The NEV right singular vectors will be computed in   |
c     | the N by NCV array V.                                |
c     |                                                      |
c     | The NEV left singular vectors will be computed in    |
c     | the M by NEV array U.                                |
c     |                                                      |
c     | NEV is the number of singular values requested.      |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXM:   Maximum number of rows of the A allowed.     |
c     | MAXN:   Maximum number of columns of the A allowed.  |
c     | MAXNEV: Maximum NEV allowed                          |
c     | MAXNCV: Maximum NCV allowed                          |
c     %------------------------------------------------------%
c
      integer maxm,maxn, maxnev,maxncv, ldu,ldv
      parameter ( maxm=lyb**2, maxn=lxb**2,
     $            maxnev=2, maxncv=25, ldu=maxm, ldv=maxn)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      real
     & v(ldv, maxncv), u(ldu, maxnev), 
     & workl(maxncv*(maxncv+8)), workd(3*maxn), 
     & s(maxncv,2), resid(maxn), ax(maxm)
      logical select(maxncv)
      integer iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character bmat*1, which*2
      integer ido, m, n, nev, ncv, lworkl, info, ierr,
     &        j,k, ishfts, maxitr, mode1, nconv
      logical rvec
      real tol, sigma, temp
c
c     %------------%
c     | Parameters |
c     %------------%

      real one, zero
      parameter (one = 1.0E0, zero = 0.0E0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      real dnrm2
      external dnrm2, daxpy, dcopy, dscal
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
 
c     include 'debug.h'
c     ndigit = -3
c     logfil = 6
c     msgets = 0
c     msaitr = 0 
c     msapps = 0
c     msaupd = 1
c     msaup2 = 0
c     mseigt = 0
c     mseupd = 0
 
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      m = my*my
      n = mx*mx
      scal = one
c
c     %------------------------------------------------%
c     | Specifications for ARPACK usage are set        | 
c     | below:                                         |
c     |                                                |
c     |    1) NEV = 4 asks for 4 singular values to be |  
c     |       computed.                                | 
c     |                                                |
c     |    2) NCV = 20 sets the length of the Arnoldi  |
c     |       factorization                            |
c     |                                                |
c     |    3) This is a standard problem               |
c     |         (indicated by bmat  = 'I')             |
c     |                                                |
c     |    4) Ask for the NEV singular values of       |
c     |       largest magnitude                        |
c     |         (indicated by which = 'LM')            |
c     |       See documentation in DSAUPD for the      |
c     |       other options SM, BE.                    | 
c     |                                                |
c     | Note: NEV and NCV must satisfy the following   |
c     |       conditions:                              |
c     |                 NEV <= MAXNEV,                 |
c     |             NEV + 1 <= NCV <= MAXNCV           |
c     %------------------------------------------------%
c
      nev   = 2
      ncv   = 10 
      bmat  = 'I'
      which = 'LM'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SVD: N is greater than MAXN '
         go to 9000
      else if ( m .gt. maxm ) then
         print *, ' ERROR with _SVD: M is greater than MAXM '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SVD: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SVD: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |              (machine precision) is used.           |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1.)             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = n
      mode1 = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
            call rzero(ax, m)
            call rzero(workd(ipntr(2)), n)
            call adf_shuff2d(2,my,mx,ax,workd(ipntr(1)),scal,works)
            call adf_shuff2d(1,my,mx,ax,workd(ipntr(2)),scal,works)
            go to 10
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %--------------------------------------------%
c        | No fatal errors occurred.                  |
c        | Post-Process using DSEUPD.                 |
c        |                                            |
c        | Computed singular values may be extracted. |  
c        |                                            |
c        | Singular vectors may also be computed now  |
c        | if desired.  (indicated by rvec = .true.)  | 
c        |                                            |
c        | The routine DSEUPD now called to do this   |
c        | post processing                            | 
c        %--------------------------------------------%
c           
         rvec = .true.
c
         call dseupd ( rvec, 'All', select, s, v, ldv, sigma, 
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
c
c        %-----------------------------------------------%
c        | Singular values are returned in the first     |
c        | column of the two dimensional array S         |
c        | and the corresponding right singular vectors  | 
c        | are returned in the first NEV columns of the  |
c        | two dimensional array V as requested here.    |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of DSEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' '
c
         else
c
            nconv =  iparam(5)
            do 20 j=1, nconv
c
               s(j,1) = sqrt(s(j,1))
c
c              %-----------------------------%
c              | Compute the left singular   |
c              | vectors from the formula    |
c              |                             |
c              |     u = Av/sigma            |
c              |                             |
c              | u should have norm 1 so     |
c              | divide by norm(Av) instead. |
c              %-----------------------------%
c
               call rzero(ax,m)
               call adf_shuff2d(2,my,mx,ax,v(1,j),scal,works) 
               call dcopy(m, ax, 1, u(1,j), 1)
               temp = one/dnrm2(m, u(1,j), 1)
               call dscal(m, temp, u(1,j), 1)
c
c              %---------------------------%
c              |                           |
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*v - sigma*u ||    |
c              |                           |
c              | for the NCONV accurately  |
c              | computed singular values  |
c              | and vectors.  (iparam(5)  |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance).               |
c              | Store the result in 2nd   |
c              | column of array S.        |
c              %---------------------------%
c
               call daxpy(m, -s(j,1), u(1,j), 1, ax, 1)
               s(j,2) = dnrm2(m, ax, 1)
c
 20         continue
c
c           %-------------------------------%
c           | Display computed residuals    |
c           %-------------------------------%
c
c           call dmout(6, nconv, 2, s, maxncv, -6,
c    &                'Singular values and direct residuals')
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit',
     &               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      
 
c        print *, ' '
c        print *, ' _SVD '
c        print *, ' ==== '
c        print *, ' '
c        print *, ' Size of the matrix is ', n
c        print *, ' The number of Ritz values requested is ', nev
c        print *, ' The number of Arnoldi vectors generated',
c    &            ' (NCV) is ', ncv
c        print *, ' What portion of the spectrum: ', which
c        print *, ' The number of converged Ritz values is ', 
c    &              nconv 
c        print *, ' The number of Implicit Arnoldi update',
c    &            ' iterations taken is ', iparam(3)
c        print *, ' The number of OP*x is ', iparam(9)
c        print *, ' The convergence criterion is ', tol
c        print *, ' '
 
      end if

c     P = kron(B2,A1) + kron(B1,A2)
      do k=1,2
      do j=1,m
         B(j,3-k)=sqrt(s(k,1))*u(j,k)
      enddo ! j
      do j=1,n
         A(j,k)=sqrt(s(k,1))*v(j,k)
      enddo ! j
      enddo ! k
 9000 continue
      end
c-----------------------------------------------------------------------
c     gen_fast_spacing, from fast3d.f, for pressure grid
c-----------------------------------------------------------------------
      subroutine adf_swap_lengths

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'
      common /adfswapl/ l
      real l(lx1,ly1,lz1,lelt) ! Pablo lelv -> lelt

      common /adflf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt

      integer e

      n2 = lx1-1
      nz0 = 1
      nzn = 1
      nx  = lx1-2
      if (if3d) then
         nz0 = 0
         nzn = n2
      endif
      call adf_plane_space(lmr,lms,lmt,0,n2,wxm1,xm1,ym1,zm1,
     $                     nx,n2,nz0,nzn)

      n=n2+1
      if (if3d) then
         do e=1,nelt ! really should be nelt here !!
         do j=2,n2
         do k=2,n2
            l(1,k,j,e) = lmr(e)
            l(n,k,j,e) = lmr(e)
            l(k,1,j,e) = lms(e)
            l(k,n,j,e) = lms(e)
            l(k,j,1,e) = lmt(e)
            l(k,j,n,e) = lmt(e)
         enddo
         enddo
         enddo
         call dssum(l,lx1,ly1,lz1)
         do e=1,nelt
            llr(e) = l(1,2,2,e)-lmr(e)
            lrr(e) = l(n,2,2,e)-lmr(e)
            lls(e) = l(2,1,2,e)-lms(e)
            lrs(e) = l(2,n,2,e)-lms(e)
            llt(e) = l(2,2,1,e)-lmt(e)
            lrt(e) = l(2,2,n,e)-lmt(e)
         enddo
      else
         do e=1,nelt
         do j=2,n2
            l(1,j,1,e) = lmr(e)
            l(n,j,1,e) = lmr(e)
            l(j,1,1,e) = lms(e)
            l(j,n,1,e) = lms(e)
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
         enddo
         enddo
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         call dssum(l,lx1,ly1,lz1)
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         do e=1,nelt
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
            llr(e) = l(1,2,1,e)-lmr(e)
            lrr(e) = l(n,2,1,e)-lmr(e)
            lls(e) = l(2,1,1,e)-lms(e)
            lrs(e) = l(2,n,1,e)-lms(e)
         enddo
      endif

c     do ie=1,nelt
c     write(6,*) llr(ie),lmr(ie),lrr(ie),lls(ie),lms(ie),lrs(ie),ie
c    $  ,nelv,nelt, 'inside'
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_swap_lengths_fix

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'
      common /adfswapl/ l(lx1,ly1,lz1,lelt) ! Pablo lelv -> lelv
      real l
      common /adflf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt

      integer e

      !  hack so that none of the lengths is 0
      do ie=1,nelt
        if(llr(ie).le.1.e-7) then
          llr(ie) = lmr(ie)
        endif
        if(lrr(ie).le.1.e-7) then
          lrr(ie) = lmr(ie)
        endif
        if(lls(ie).le.1.e-7) then
          lls(ie) = lms(ie)
        endif
        if(lrs(ie).le.1.e-7) then
          lrs(ie) = lms(ie)
        endif
        if(lmr(ie).le.1.e-7) then
          write(6,*) 'wrong lmr value',ie,lmr(ie)
          call exitt
        endif
        if(lms(ie).le.1.e-7) then
          write(6,*) 'wrong lms value',ie,lms(ie)
          call exitt
        endif
        if(ldim.eq.3) then
          if(llt(ie).le.1.e-7) then
            llt(ie) = lmt(ie)
          endif
          if(lrt(ie).le.1.e-7) then
            lrt(ie) = lmt(ie)
          endif
          if(lmt(ie).le.1.e-7) then
            write(6,*) 'wrong lmt value',ie,lmt(ie)
            call exitt
          endif
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_gen_fast_spacing(x,y,z)
c
c     Generate fast diagonalization matrices for each element
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'WZ'
c
      parameter(lxx=lx1*lx1)
c
      common /adflf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
      real x(lx1,ly1,lz1,nelt)
      real y(lx1,ly1,lz1,nelt)
      real z(lx1,ly1,lz1,nelt)
      real axwt(lx2)

      ierr = 0

c     if (param(44).eq.1) then ! for FEM?
c                                       __ __ __
c        Now, for each element, compute lr,ls,lt between specified planes
c
         n1 = lx2
         n2 = lx2+1

         nz0 = 1
         nzn = 1
         if (if3d) then
            nz0= 0
            nzn=n2
         endif
         eps = 1.e-7
         if (wdsize.eq.8)  eps = 1.e-14

         ! why passing in xm1, ym1, zm1, but use wxm2?
c
c        Find mean spacing between "left-most" planes
         call adf_plane_space2(llr,lls,llt, 0,wxm2,x,y,z,n1,n2,nz0,nzn)
c
c        Find mean spacing between "middle" planes
         call adf_plane_space (lmr,lms,lmt, 1,n1,wxm2,x,y,z
     $                        ,n1,n2,nz0,nzn)
c
c        Find mean spacing between "right-most" planes
         call adf_plane_space2(lrr,lrs,lrt,n2,wxm2,x,y,z,n1,n2,nz0,nzn)
c
c     else
c        call load_semhat_weighted    !   Fills the SEMHAT arrays
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_plane_space(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
c
c     Here, spacing is based on harmonic mean.  pff 2/10/07
c
c
      include 'SIZE'
      include 'INPUT'
c
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
      j2 = i2
      k2 = i2
c
      do ie=1,nelt
c
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
c              lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
c    $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
c    $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
c    $                     *   weight
               lr2  = lr2  +   weight /
     $                       ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
     $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
     $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
c              ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
c    $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
c    $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
c    $                     *   weight
               ls2  = ls2  +   weight /
     $                       ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
     $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
     $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
c              lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
c    $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
c    $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
c    $                     *   weight
               lt2  = lt2  +   weight /
     $                       ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
     $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
     $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = 1./sqrt(lt2)
c
         else              ! 2D
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
c              lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
c    $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
c    $                     *   weight
               lr2  = lr2  + weight /
     $                       ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
     $                       + (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
c              ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
c    $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
c    $                     *   weight
               ls2  = ls2  + weight /
     $                       ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
     $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine adf_plane_space2(lr,ls,lt,i1,w,x,y,z,nx,nxn,nz0,nzn)
c
c     Here, the local spacing is already given in the surface term.
c     This addition made to simplify the periodic bdry treatment.
c
c
      include 'SIZE'
      include 'INPUT'
c
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
c
      do ie=1,nelt
c
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
               lr2  = lr2  + ( (x(i1,j,k,ie))**2
     $                     +   (y(i1,j,k,ie))**2
     $                     +   (z(i1,j,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
               ls2  = ls2  + ( (x(i,j1,k,ie))**2
     $                     +   (y(i,j1,k,ie))**2
     $                     +   (z(i,j1,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
               lt2  = lt2  + ( (x(i,j,k1,ie))**2
     $                     +   (y(i,j,k1,ie))**2
     $                     +   (z(i,j,k1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = sqrt(lt2)
c           write(6,1) 'lrlslt',ie,lr(ie),ls(ie),lt(ie)
    1       format(a6,i5,1p3e12.4)
c
         else
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
               lr2  = lr2  + ( (x(i1,j,1,ie))**2
     $                     +   (y(i1,j,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
               ls2  = ls2  + ( (x(i,j1,1,ie))**2
     $                     +   (y(i,j1,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
            write(6,*) 'lrls',ie,lr(ie),ls(ie),lt(ie)
         endif
      enddo
      return
      end
c----------------------------------------------------------------------
c     2d tests      
c----------------------------------------------------------------------
      function test_schur(nx,ny,S1,S2,V1,V2,W1,W2,A,B,WORK)
      implicit none
      real test_schur
      integer nx,ny
      real S1(*),S2(*),V1(*),V2(*),W1(*),W2(*)
      real A(nx,nx,*),B(ny,ny,*)
      real WORK(*)
      
      integer i,j,m,n
      real X(nx,ny), R(nx,ny)
      real dnrm2

      call rzero(R,nx*ny)

      do j=2,ny-1
      do i=2,nx-1

      do n=1,ny
      do m=1,nx
      X(m,n) = B(n,j,2)*A(m,i,1) + B(n,j,1)*A(m,i,2)
      enddo ! m
      enddo ! n

      call schsolv22(nx,ny,X,S1,S2,V1,V2,W1,W2,WORK)

      X(i,j) = X(i,j)-1.0E0
      R(i,j) = dnrm2(nx*ny,X,1)

      enddo ! i
      enddo ! j

      test_schur = dnrm2(nx*ny,R,1)
      return
      end
c-----------------------------------------------------------------------
      function test_fdm(nx,ny,LL,V1,V2,W1,W2,A,B)
      implicit none
      real test_fdm
      integer nx,ny
      complex*16 LL(nx,ny), V1(nx,nx),V2(ny,ny),W1(nx,nx),W2(ny,ny)
      real A(nx,nx,2),B(ny,ny,2)
      
      integer i,j,m,n
      real X(nx,ny), Y(nx,ny), R(nx,ny)
      real dnrm2

      call rzero(R,nx*ny)

      do j=2,ny-1
      do i=2,nx-1

      do n=1,ny
      do m=1,nx
      X(m,n) = B(n,j,1)*A(m,i,1) + B(n,j,2)*A(m,i,2)
      enddo ! m
      enddo ! n

      call fdmsolv22(nx,ny,Y,LL,V1,V2,W1,W2,X)

      Y(i,j) = Y(i,j)-1.0E0
      R(i,j) = dnrm2(nx*ny,Y,1)

      enddo ! i
      enddo ! j

      test_fdm = dnrm2(nx*ny,R,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine fullprec(nx,ny,nz,A,B,C)
      implicit none
      include 'SIZE'
      common /tens_rank/ xyzrank,terms,ptr
      integer xyzrank(3),terms,ptr(3,4)

      integer nx,ny,nz
      real A(nx*nx,*), B(ny*ny,*), C(nz*nz,*)
      real P(1)!((lx1*ly1*lz1)**2)
      integer n, i,j

      real one
      parameter (one=1.0E0)
      
      n = nx*ny*nz

      call rzero(P,n*n)
      do j=1,terms
         call dkron3(nz,ny,nx, nz,ny,nx, 
     $        one,C(1,ptr(3,j)),B(1,ptr(2,j)),A(1,ptr(1,j)), one,P)
      enddo ! j
      call savemat(nx,n,n,P,n)
      return
      end
c-----------------------------------------------------------------------
