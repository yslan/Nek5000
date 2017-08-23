c-----------------------------------------------------------------------
      subroutine newton

      include 'SIZE'
      include 'TOTAL'

      real cn(lx1,ly1,lz1,lelt),
     $     sk(lx1,ly1,lz1,lelt),
     $     gk(lx1,ly1,lz1,lelt)

      ! set parameters

c     npseudo=...
c     nnewton=...
      nel = lx1*ly1*lz1*nelv

      ! zeroing

      ! pseudo-time-step iteration
      do i=1,npseudo
         ! newton iteration
c        dtnt = 10.
         do j=1,nnewton
c           call hmh_gmres(...)
c           call add2s2(cn,sk,alpha,nel)
         enddo
c        call compute_f(...)
c        call compute_gk(...)
c        call chsign(gk,npts)

c        rnorm=sqrt(glsc3(gk,gk,mult,nel))
         if (i.eq.1) rinorm=rnorm
         ratio=rnrom/rinorm
      enddo

      return
      end
c-----------------------------------------------------------------------
