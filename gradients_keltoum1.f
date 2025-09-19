c------------------------------------------------------------------
      subroutine grad_grad_wei(exmp, iso, gamma, N, Ng, 
     &                 nbtetra, tetra, xin, grad_grad)
	  implicit none	         
	  include 'simplexes_keltoum.h'
	  include 'constants_keltoum.h'


	  integer  N, nbtetra,  tetra(4, *), Ng, exmp, iso 
      double precision xin(3,*), grad_grad(10,nbtetra),
     &                   gamma

      integer l, i, j, s,  nm1, nm2, nb, pow_sum, nbtet
	  external pow_sum 
	  double precision BT(3,3),IBT(3,3),IBT_T(3,3), IBKT(3,3),
     & 	               BK(3,3), IBK(3, 3), avrvalue_a
	  double precision  IBK_T(3,3),  IBKT_T(3,3), Avrgmat(3,3),
     &  MAT1(3, 3), MAT2(3, 3),MAT(3, 3), JACO(3, 3), VTX0(3) 


 	  double precision gradref(3,4), delt, mu(3,4), AA(3,3)
	  double precision  res, res_g,  xh(3),
     &         delt0, delt1, idelt, verct(3)	  
      double precision xx(3,4), vc0(3),  xvec(3), 
     &   detBT, nmu(3, nq), mu_n(3,3,2), dimin, dimax	   
	  double precision va(3, 4), val0(4), Agwbase(3, 4), 
     &  det, rval0(4), absdete, prd, prd1,
     &   wbase(4), gwbase(3, 4), vec(3),
     &    dgau(Ng), rdgau(Ng), rwgau(Ng), wgau(Ng)
	  integer  nbprisme, k, nm, m, ss, nnb, kk, nbr, ni, ii,
     &          i1, k1 
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

       do nnb = 1, nbtet 
          
          nb = nnb + (s-1)*nbtet

	      do l = 1, 10
             grad_grad(l, nb) = 0.D0
	      enddo
          nm1 = tetra(1, nb)	      
	      do i = 1, 3
		      VTX0(i) = xin(i, nm1)  
	      enddo
	      
	      kk = 1
          if (nnb.eq.1) kk = 2 

!         a vertex of the tetrahedra          

	      do k = 2, 4
             nm =  tetra(k, nb)
              do i = 1, 3 
               bk(i,k-1) = xin(i,nm) - VTX0(i) 
             enddo
          enddo


          call inversion(BK, IBK)
          call matmat(3, 3, 3, IBK, BT, IBKT) 
	      call transpos(3, 3, IBKT, IBKT_T) 
          call matv3(IBK, VTX0, vc0)   
    
         do k = 1, 4  ! vertices

  	       call nmat_v3(k, IBKT_T, gradref, va)  
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

	     if ((dabs(gamma).lt.zero)) delt0 = 1.D0
	     if ((dabs(gamma-1.D0).lt.zero)) delt1 = 1.D0
	     
	     if (dabs(delt).gt.zero) then
            if ((dabs(gamma).gt.zero)) delt0 = delt**gamma
            if ((dabs(gamma-1.D0).gt.zero)) delt1 = delt**(gamma-1.D0)
            idelt = 1./delt
         else 
            idelt = 0.D0
         endif
         
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
            call mat_v3(i, bt, nmu, xvec) 
            do ii = 1, 3
               xvec(ii) = xvec(ii)*idelt
            enddo

!           the matrix A(x)            
            call A_MAT(exmp, iso, 0,  xvec, mat1)  
            
!           the matrix B_T^{-T}(I-2xs*ct)             
		    call matmat(3, 3, 3, ibt_t, jaco, mat2)  
		    
!           the matrix A B_T^{-T}  (I-2xs*ct) 
		    call matmat(3, 3, 3, mat1, mat2, mat)
		    
!           the matrix (I-2 ct*xs^t) B^{-T} A B_T^{-T}(I-2xs*ct)         	 
	        call transpos(3, 3, mat2, mat1)
		    call matmat(3, 3, 3, mat1, mat, AA)     	        		    
            
		    prd = weiq(i)*absdete*prd1 
		    
		    do k = 1, 4
	         call prodcolumns(i, k,  nmu, va, res)
 	         wbase(k) = (val0(k) + res*delt)*delt0
 	         do m = 1, 3
 	          vec(m) = delt0*va(m,k)
     &       + gamma*delt1*(val0(k)+res*delt)
 	          gwbase(m, k) = vec(m)
 	         enddo
	         call matv3_bis(k, AA, vec, Agwbase) 	     	           
		    enddo
		    
	 
            do k = kk, 4
             do  j = k, 4
              call prodcolumns(k, j, Agwbase, gwbase, res_g)
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
   	        call inversion(BK, IBK)
	        call transpos(3, 3, IBK, IBK_T)
	        
	        do i = 1, 3
	          VTX0(i) = xin(i, nm1)
	        enddo   

            do k = 1, 4
  	         call nmat_v3(k, IBK_T, gradref, va)
            enddo

            if (iso.eq.1) then
             call average_Asca_elem(exmp,  BK,  VTX0,
     &                              avrvalue_a) 
             do k = 1, 4
	          do j = k, 4
               call prodcolumns(k, j, va, va, res)
	           l = ((j-1)*j)/2 + k
	           grad_grad(l, nb) = res*absdete*avrvalue_a/6.D0 
	          enddo
             enddo
            else
             call average_Amat_elem(exmp, BK, verct,
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
	          grad_grad(l, nb) = res*absdete/6.D0
	         enddo
            enddo
            endif    
       
           enddo

	       return
          end
            

!--------------------------------------------------------------------------
            subroutine average_Amat_elem(exmp, BK, verct,
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
