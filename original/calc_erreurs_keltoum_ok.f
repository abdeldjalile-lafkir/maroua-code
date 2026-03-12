c      Etat de modification : fini 
      subroutine prodsc_L2(gamma1, gamma2, u, v, prod,
     & 	                  prodi, prodf)
!       this subroutine computes the scalar product in L2 
!       of u and v 
!       (u, v)_L2 = int(u*v)dx     
	  implicit none
      include 'fondements_keltoum.h'
	  include  'constants_keltoum.h'

          

 	  double precision u(*), v(*), prod, gamma1, gamma2, ress,
     &                   prodi, prodf, res1
	  integer i, j, nb, s, k, i1, j1, l,  nbd

	  double precision res, som, ww0(10, nbtetra)

       
c  	   write(*, *) '  calcul des produits d''ordre 0...'
	   call ordre0(gamma1, gamma2,  N, Ng, nbtetra,
     &	         tetra, xin, ww0)


      res = 0.D0
	  prodi = 0.D0
	  prodf = 0.D0
	  do nb = 1, nbtetra
	    ress = 0.D0
        do i = 1, 4
	     s = tetra(i, nb)
	     res1 = 0.D0
	     if (s.ne.1) then
	       do j = 1, 4
	        k = tetra(j, nb)
	        if (k.ne.1) then
	         j1 = max(i,j)
	         i1 = min(i,j)
	         l  = (j1-1)*j1/2 + i1
             res1 = res1 + ww0(l,nb)*v(k-1)
	        endif
	       enddo
	       res = res + res1*u(s-1)
	       ress =  ress + res1*u(s-1)
          endif
	    enddo
        nbd = domtet(nb)
        if (nbd.lt.5) then
         prodi = prodi +  ress 
	    else
         prodf = prodf + ress
        endif
       enddo
	   prod  = res
	  return
	  end
c-----------------------------------------------------------
c      Etat de modification : fini
      subroutine prodsc_L2_grad(u, v, prod, prodi, prodf)
	  implicit none
      include 'fondements_keltoum.h'
	  include 'constants_keltoum.h'
!       this subroutine computes the scalar product in L2 
!       of grad_u and grad_v

 	  double precision u(*), v(*), prod, res1,
     &                    prodi, prodf, ress
	  integer i, j, nb, s, k, i1, j1, l, nbd
	  double precision res

      res = 0.D0
	  prodi = 0.D0
	  prodf = 0.D0
	  do  nb = 1, nbtetra
	    ress = 0.D0
        do i = 1, 4
         s = tetra(i, nb)
	     if (s.ne.1) then
	      res1 = 0.D0
	      do j = 1, 4
	       k = tetra(j, nb)
	       if (k.ne.1) then
	        j1 = max(i, j)
	        i1 = min(i, j)
	        l  = (j1 - 1)*j1/2 + i1
            res1  = res1 + grad_grad(l,nb)*v(k-1) 
	       endif
	      enddo	    
	      res = res +  res1*u(s-1) 
	      ress = ress + res1*u(s-1) 
	     endif
	   enddo	
       nbd = domtet(nb)
       if (nbd.lt.5) then
         prodi = prodi +  ress 
	   else
         prodf = prodf + ress
	   endif
      enddo
	  prod = res
      return
	  end       
c-------------------------------------------------------------
c       Etat d'avancement: fini (aucun changement)
        subroutine norm_L2(gamma, u, nrm2, nrm2i, nrm2f)
!       this subroutine computes the result
!        norm_L2( u )= (int(u*u)dx)^(1/2)	    
	    implicit none
	    double precision u(*), nrm2, res, gamma, nrm2i, nrm2f               
               
        call prodsc_L2(gamma, gamma,  u, u, res, nrm2i, nrm2f)
	    nrm2 = dsqrt(res)
  	    nrm2i = dsqrt(nrm2i)
        nrm2f = dsqrt(nrm2f)		   
        return
	    end
c-------------------------------------------------------------
c       Etat d'avancement: fini (aucun changement)
        subroutine norm_grad_L2(u, nrm2, nrm2i, nrm2f)
!       this subroutine computes the result
!        norm_L2(grad_u)= (int(grad_u*grad_u))^(1/2)	    
	    implicit none 
	    double precision u(*), nrm2, res, nrm2i, nrm2f
               
        call prodsc_L2_grad(u, u, res, nrm2i, nrm2f)
	    nrm2 = dsqrt(res)
	    nrm2i = dsqrt(nrm2i)
        nrm2f = dsqrt(nrm2f)
        return
	    end
c-------------------------------------------------------------
c       Etat d'avancement: fini (aucun changement)
        subroutine erreur_L2( gamma, u, v, erreur, 
     &                       erreuri, erreurf)
!       computing  norm_L2(u - v)
	    implicit none
        include 'fondements_keltoum.h'
	   
	    double precision u(*), v(*), erreur, gamma,
     &                     erreuri, erreurf
		integer i
        double precision w(nbnodes-1)

        do i = 1, nbnodes - 1
          w(i) = u(i) - v(i)
	    enddo

        call norm_L2(gamma, w, erreur, erreuri, erreurf)
	    return
        end
c-------------------------------------------------------------
c       Etat d'avancement: fini (aucun changement)

        subroutine erreur_grad_L2(gamma, u, v, erreur, 
     &   erreuri, erreurf)
     
!       comuting norm_L2(grad_u - grad_v)     
	    implicit none
        include 'fondements_keltoum.h'
        double precision u(nbnodes-1), v(nbnodes-1),
     &           erreur, gamma, erreuri, erreurf
          
		
		integer i
        double precision w(nbnodes-1)

        do i = 1, nbnodes - 1
         w(i) = u(i) - v(i)
	    enddo

        call norm_grad_L2(w, erreur, erreuri, erreurf)
	    return
        end 
c------------------------------------------------------------
c      Etat d'avancement: fini (aucun changement)
       subroutine norms_erreurs_w10(gamma, u, v, nrm0, nrm1, 
     &            nrm0i, nrm0f, nrm1i, nrm1f,
     &            erreurL2, errGrd, erreurL2i, errGrdi,
     &            erreurL2f, errGrdf)
     
          
	   implicit none
  	   include 'fondements_keltoum.h'
       double precision u(*), v(*),
     &         nrm0i, nrm0f, nrm1i, nrm1f,
     &         erreurL2, errGrd, erreurL2i, errGrdi,
     &         erreurL2f, errGrdf,
     &         gamma, nrm0, nrm1	
       double precision u1(nbnodes-1), v1(nbnodes-1)   
   
       double precision gammau1, gammav1, ss
       
        ss = -1.0

        call weighting(ss, gamma, u, gammau1, u1)
        call weighting(ss, gamma, v, gammav1, v1)
        call erreur_grad_L2(gamma, u, v, errGrd,
     &    errGrdi, errGrdf)
        call norm_grad_L2(u,  nrm1, nrm1i, nrm1f)
              

               
       call norm_L2(gammau1, u1,  nrm0, nrm0i, nrm0f)
       call erreur_L2(gammau1, u1, v1, erreurL2, 
     &                   erreurL2i, erreurL2f)
	   return 
	   end
c
c----------------------------------------------------------
c      Etat de modification : fini

       subroutine weighting(alpha, gammau, u, gammav, v)
	   implicit none
       include 'fondements_keltoum.h'
	   double precision gammau, gammav, alpha, u(*), v(*)

	   integer i, dom, k, m
       double precision x(3), vx(3), rho, rr, res

	   external rr

	   do i = 2, nbnodes
	    dom = domnod(i)
	    do k = 1, 3
         x(k) = xin(k, i)
	    enddo

        if (dom.lt.5.and.dom.gt.0) then
	      call fphi(dom, x, vx)
	      res =  rr(dom,  vx)
	      rho = 1.D0
	      do m = 1, 3
           rho  = rho + vx(m)*vx(m)
	      enddo
          rho = dsqrt(rho)
          v(i-1) = u(i-1)*(rho**alpha)/(res**alpha)
	    else if (dom.eq.5) then
	     rho = 1.D0
	     do m = 1, 3
          rho  = rho + x(m)*x(m)
	     enddo
	     rho = dsqrt(rho)   
	     v(i-1) = u(i-1)*(rho**alpha)    
	    else 
	     write(*, *) 'probleme'
		 read(*, *)            
	    endif
	   enddo
       gammav = gammau - alpha
       return
       end







