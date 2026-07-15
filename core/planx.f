      SUBROUTINE PLAN3 (IGEOM)
      use scrns_mod
      use scrvh_mod
C-----------------------------------------------------------------------
C
C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
C
      real, pointer :: RESV1(:,:,:,:), RESV2(:,:,:,:), RESV3(:,:,:,:)
     $               , DV1(:,:,:,:), DV2(:,:,:,:), DV3(:,:,:,:)
      real, pointer :: H1(:,:,:,:), H2(:,:,:,:)
C
      H1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrvh(0*lx1*ly1*lz1*lelv+1 : 1*lx1*ly1*lz1*lelv)
      H2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrvh(1*lx1*ly1*lz1*lelv+1 : 2*lx1*ly1*lz1*lelv)
      RESV1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(0*lx1*ly1*lz1*lelv+1 : 1*lx1*ly1*lz1*lelv)
      RESV2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(1*lx1*ly1*lz1*lelv+1 : 2*lx1*ly1*lz1*lelv)
      RESV3(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(2*lx1*ly1*lz1*lelv+1 : 3*lx1*ly1*lz1*lelv)
      DV1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(3*lx1*ly1*lz1*lelv+1 : 4*lx1*ly1*lz1*lelv)
      DV2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(4*lx1*ly1*lz1*lelv+1 : 5*lx1*ly1*lz1*lelv)
      DV3(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scrns(5*lx1*ly1*lz1*lelv+1 : 6*lx1*ly1*lz1*lelv)
C
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL MAKEF
C
      ELSE
C
C        New geometry, new b.c.
C
         intype = -1
         call sethlm  (h1,h2,intype)
         call cresvif (resv1,resv2,resv3,h1,h2)

         call ophinv  (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2  (vx,vy,vz,dv1,dv2,dv3)
c
         call incomprn(vx,vy,vz,pr)
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE LAGPRES 
C--------------------------------------------------------------------
C
C     Keep old pressure values
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      common /cgeom/ igeom

      IF (NBDINP.EQ.3.and.igeom.le.2) THEN
         NTOT2 = lx2*ly2*lz2*NELV
         CALL COPY (PRLAG,PR,NTOT2)
      ENDIF
      RETURN
      END
C
      subroutine cresvif (resv1,resv2,resv3,h1,h2)
      use scruz_mod
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the velocity solver
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      REAL           RESV1 (LX1,LY1,LZ1,1)
      REAL           RESV2 (LX1,LY1,LZ1,1)
      REAL           RESV3 (LX1,LY1,LZ1,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      real, pointer :: W1(:,:,:,:), W2(:,:,:,:), W3(:,:,:,:)

      common /cgeom/ igeom

      W1(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scruz(0*lx1*ly1*lz1*lelv+1 : 1*lx1*ly1*lz1*lelv)
      W2(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scruz(1*lx1*ly1*lz1*lelv+1 : 2*lx1*ly1*lz1*lelv)
      W3(1:lx1,1:ly1,1:lz1,1:lelv) =>
     $   cb_scruz(2*lx1*ly1*lz1*lelv+1 : 3*lx1*ly1*lz1*lelv)

      NTOT1 = lx1*ly1*lz1*NELV
      NTOT2 = lx2*ly2*lz2*NELV
      if (igeom.eq.2) CALL LAGVEL 
      CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      CALL BCNEUTR
C
      call extrapp (pr,prlag)
      call opgradt (resv1,resv2,resv3,pr)
      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)
      CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)
C
      RETURN
      END
c-----------------------------------------------------------------------
