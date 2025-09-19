c------------------------------------------------------------------
      subroutine grad_grad_wei(exmp, iso, gamma, N, Ng, 
     &                 nbtetra, tetra, xin, grad_grad)
	  implicit none	         
	  include 'simplexes_keltoum.h'
	  include 'constants_keltoum.h'


	  integer  N, nbtetra,  tetra(4, *), Ng, exmp, iso 
      double precision xin(3,*), grad_grad(10,nbtetra),
     &                   gamma

      integer l, i, j, s,  nm1, nm2, nb, pow_sum, nprisme, nbtet
	  external pow_sum 
	  double precision BT(3,3),IBT(3,3),IBT_T(3,3), IBKT(3,3),
     & 	               BK(3,3), IBK(3, 3), avrvalue_a
	  double precision  IBK_T(3,3), MTK(3,3), IBKT_T(3,3),
     &                  MAT(3, 3), JACO(3, 3), VBK(3), 
     &                  vertices(4, 3), Avrgmat(3,3)

 	  double precision gradref(3,4), delt, mu(3,4), AA(3,3)
	  double precision  res, res_g, res40, res41, res42, res31, 
     &         ress(4), res1, res32, res21, res22, xh(3),
     &         delt0, delt1, delt2, delt3, idelt, verct(3)	  
      double precision xx(3,4), vc0(3), vec1(3), vec2(3), xvec(3), 
     &   veca(3), detBT, ncc, nmu(3, nq), mu_n(3,3,2), dimin, dimax	   
	  double precision  va(3, 4), val0(4), Aveca(3), Avec2(3),
     &	               det, rval0(4), absdete, prd, prd1
	  double precision dgau(Ng), rdgau(Ng), rwgau(Ng), wgau(Ng)
	  integer  nbprisme, k, nm, m, ss, nnb, kk, nbr, ni, ii,
     &          i1, k1, j1       
	  external det
 
 

	  call grad_ref(gradref)                
	  do i = 2, 4
	    rval0(i) = 0.D0
	  enddo
	 
	  rval0(1) = 1.D0

	  nbprisme = pow_sum(2,N-1) + pow_sum(1,N-1)
	  nbtet    = 3*nbprisme + N
     
	  call gauss_lobatto(rdgau, rwgau, Ng) 
	  do s = 1, 4   ! index of the big subdomain
	   call bt_big_simp(s, BT)
       call inversion(BT, IBT)
	   detBT = det(BT) 
	   call transpos(3, 3, IBT, IBT_T)
       call matv3(IBT_T, cc, veca)   

       do nnb = 1, nbtet 
          
          nb = nnb + (s-1)*nbtet

	      do l = 1, 10
             grad_grad(l, nb) = 0.D0
	      enddo
	      kk = 1
          if (nnb.eq.1) kk = 2 
          nm1 = tetra(1, nb)
	      do k = 2, 4
             nm =  tetra(k, nb)
              do i = 1, 3 
               bk(i,k-1) = xin(i,nm) - xin(i, nm1)   
             enddo
          enddo
          
		  do i = 1, 3
           vbk(i) = xin(i, nm1)   ! a vertex of the tetrahedra
	      enddo
	      
          call inversion(BK, IBK)
!         call matmat(3, 3, 3, IBT, BK, MTK)
          call matmat(3, 3, 3, IBK, BT, IBKT)
	      call transpos(3, 3, IBKT, IBKT_T) 
   
          call matv3(IBK, VBK, vc0)   

    
         do k = 1, 4  ! vertices

  	       call nmat_v3(k, IBKT_T, gradref, va)  ! OK 
           res = 0.D0
	       do m = 1, 3
	        res = res + gradref(m,k)*vc0(m)
	       enddo
           val0(k) = rval0(k) - res
!          here we compute the coordinates of 
!          the antecedents of the vertices of K (the nb-tetrahedra)
!          by inverse of F_T. The result is saved in xx
 		   nm =  tetra(k, nb)
           do i = 1, 3
            xh(i) = xin(i, nm)  
           enddo
   	       call matv3_bis(k, IBT, xh, xx)
         enddo



	    call delta_extremas(xx,  dimin, dimax)
	    call dec_gauss_lobatto(dimin, dimax, Ng, rwgau, 
     &                        rdgau,  wgau, dgau)
!       here we compute the integrals by a quadraure formula
!     
	   do ss = 1, Ng
         delt = dgau(ss)
         delt0 = 0.D0
	     delt1 = 0.D0
         delt2 = 0.D0
	     if (dabs(delt).gt.zero) then
           delt0 = delt**(2.D0*gamma - 2.D0)
           delt1 = delt**(2.D0*gamma - 1.D0)
           delt2 = delt**(2.D0*gamma)
           idelt = 1./delt
         else 
           idelt = 0.D0
         endif
	     if ((dabs(gamma-1.D0).lt.zero)) delt0 = 1.D0
	     if ((dabs(gamma-0.5D0).lt.zero)) delt1 = 1.D0
	     if ((dabs(gamma).lt.zero)) delt2 = 1.D0

   
 
 
!        Computing the intersection of plane x+y+z=delt
!        with the tetrahedra whose vertices are given by xx
        
         call intersection(xx, delt, mu, nbr)
 
         ni = 0 
          if (nbr.eq.3) then 
!          The intersection is a triangle
	       ni = 1
           do k1 = 1, 3 
	        do i1 = 1, 3
              mu_n(i1, k1, 1) = mu(i1, k1)
		    enddo 
	       enddo
		  else if (nbr.eq.4) then
!          The intersection is a quadrilateral		   
		     ni =  2
	         call  div_quadrilatere(mu, mu_n)
 	      else if (nbr.lt.3.and.nbr.ge.0) then
	        res_g = 0.D0
            ni    = 0
          else
          	write(*, *) 'A bad intersection between a plane 
     &	       and a tetrahedra.'
	        write(*,*)  '(What hapens here is necessarily a mistake)' 
	        write(*, *) 'The program is stopped. Please correct.'
            stop
	     endif	
	     
	      prd1 = wgau(ss)*dabs(detBT)
	      do i1 = 1, ni
	       call triangles(i1, mu_n, absdete, nmu)

!          here we compute the integrals in the 
!          triangles of intersection
!          we used a quadrature formula at seven points	
       
           do i = 1, nq
	        call jacob_phi(i, nmu, jaco)
            call mat_v3(i, bt, nmu, xvec)  ! a verifier
            do ii = 1, 3
               xvec(ii) = xvec(ii)*idelt
            enddo
            call A_MAT(exmp, iso, 0,  xvec, AA)
            
            call matv3(AA, veca, Aveca)
            call prodsc(veca, Aveca, ncc)       
	        ncc = gamma*gamma*ncc
            
		    prd = weiq(i)*absdete*prd1
		    
		    do k = 1, 4
	         res = 0.D0
		     do m = 1, 3
	         res = res + nmu(m,i)*va(m, k)
	         enddo
	         ress(k) = res
		    enddo
		    
		    call matmat(3, 3, 3, ibt_t, jaco, mat)	 
            do k = kk, 4
             do  j = k, 4
!             (alpha(k).xs)*(alpha(j).xs)*delt2             
	          res42 = ress(k)*ress(j)*delt2  
	          
!             (alpha(k)*a(j)+alpha(j)*a(k))*delta1      
              res41 = ress(k)*val0(j) + ress(j)*val0(k)  
	          res41 = res41*delt1
	          
!             a(j)*a(k)*delt0	          
              res40 = val0(j)*val0(k)*delt0
              
!             IBT_T*(I-2xs*ct)*alpha(k)              
 	          call mat_v3(k, mat, va, vec1) 
 	          
!             IBT*(I-2xs*ct)*alpha(j) 	          
		      call mat_v3(j, mat, va, vec2)
		      
!             Avec2 = AA()*IBT*(I-2xs*ct)*alpha(j)		      
		      call matv3(AA, vec2, Avec2) 
		      
!             alpha(k)^T *(I-2 xs c^t) IBT AA() IBT*(I-2 c xs^t)*alpha(j)		      
	          call prodsc(vec1, Avec2, res1) 
              res1  = res1*delt2
              
              
!            Produit des termes d'ordre 1
	          call prodsc(vec1, Aveca, res)
	          res21 = res*val0(j)*delt1
	          res22 = res*ress(j)*delt2

	          call prodsc(vec2, Aveca, res)
	          
		      res31 = res*val0(k)*delt1
		      res32 = res*ress(k)*delt2

	          res_g = ncc*(res40 + res41 + res42)
     &           + res1 - gamma*(res21+res22+res31+res32)

  	          res_g =  res_g*prd
                l = ((j-1)*j)/2 + k
	          grad_grad(l, nb) = grad_grad(l, nb) + res_g
	         enddo
	        enddo
	       enddo
	      enddo 	     
        enddo

	   enddo 
       enddo

!     
!         les elements finis
!
           s = 5
	       do nb = (s-1)*nbtet + 1, s*nbtet

	        nm1 = tetra(1, nb)
	        
 	        do j = 2, 4
	         nm2 = tetra(j, nb)
		     do i = 1, 3
	         bk(i, j-1) = xin(i,nm2) - xin(i, nm1)
		     enddo  
	        enddo 
	        absdete = det(BK)
	        absdete = dabs(absdete)

	        do i = 1, 3
	         verct(i) = xin(i, nm1)
	        enddo   
   	        call inversion(BK, IBK)
	        absdete = det(BK)
	        absdete = dabs(absdete)
	        call transpos(3, 3, IBK, IBK_T)
            do k = 1, 4
  	         call nmat_v3(k, IBK_T, gradref, va)
            enddo

            if (iso.eq.1) then
             call average_Asca_elem(exmp,  BK, verct,
     &                              avrvalue_a)
             do k = 1, 4
	          do j = k, 4
	           res = 0.D0
               do m = 1, 3               
                res =  res + va(m,k)*va(m,j)               
               enddo
	           l = ((j-1)*j)/2 + k
	           grad_grad(l, nb) = res*absdete*avrvalue_a 
	          enddo
             enddo
            else

             call average_Amat_elem(exmp, iso,  BK, verct,
     &                              AvrgMat)  
          
             do k = 1, 4
	         do j = k, 4
	          res = 0.D0
              do m = 1, 3
               do n = 1, 3
                res =  res + va(n,k)*AvrgMat(n,m)*va(m,j)
               enddo
              enddo
	          l = ((j-1)*j)/2 + k
	          grad_grad(l, nb) = res*absdete
	         enddo
            enddo
            endif    
       
           enddo
	       return
          end
            

!--------------------------------------------------------------------------
            subroutine average_Amat_elem(exmp, iso,  BK, verct,
     &                              AvrgMat)
         
!            This subroutine computes the integral of a matrix
!            on a tetrahedra 
!            The result is saved in AvrgMat
            implicit none
            include "constants_keltoum.h"
            
            double precision BK(3, 3), Avrgmat(3,3), AA(3, 3)
            double precision verct(3), xvec(3), svec(3)
            integer exmp, i, j, k, iso
            


            do i = 1, 3
             do j = 1, 3
              Avrgmat(i,j) = 0.D0
             enddo 
            enddo 

            do k = 1, NQT
             call mat_v3(k, BK, xyzqt, xvec)
             call som_vec(3, 1.D0, xvec, 1.D0, verct, svec)
             call A_MAT(exmp, iso, 0, svec, AA)
             do i = 1, 3
              do j = 1, 3
                Avrgmat(i,j) = Avrgmat(i,j)+AA(i, j)*weiqt(k)
              enddo 
             enddo 
            enddo
           end
!--------------------------------------------------------------------------
            subroutine average_Asca_elem(exmp, BK, verct,
     &                              avrvalue_a)
         
!            This subroutine computes the integral of the coefficient
!            "a(.)" on a tetrahedra with the coordinates of
!            vertices given (vertices).
!            The result is saved in avrvalue_a
            implicit none
            include "constants_keltoum.h"
            
            double precision avrvalue_a, a_sca, val
            double precision verct(3), xvec(3), svec(3), BK(3, 3)
            integer exmp, k
            external a_sca

            avrvalue_a = 0.D0
            
            do k = 1, NQT
             call mat_v3(k, BK, xyzqt, xvec)
             call som_vec(3, 1.D0, xvec, 1.D0, verct, svec)
             val = a_sca(exmp, 0, xvec)
             avrvalue_a = avrvalue_a + val*weiqt(k)
            enddo
           end
