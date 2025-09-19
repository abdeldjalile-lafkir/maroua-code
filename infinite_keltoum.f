! dossier:  /Users/TZB/Desktop/TOTAL_COMPTE_OCT2013/TAHAR/MES_RECHERCHES/CODES/code3D_kaliche
! commande : gfortran -Wunused  

        program tester
	    implicit none
        include 'fondements_keltoum.h'
	    include 'constants_keltoum.h'
 	    include 'simplexes_keltoum.h'



         double precision gamma1, gamma2, som, erreurL2, errGrd,
     &   erreurL2i, errGrdi, erreurL2f, errGrdf,
     &   nrm0i, nrm0f, nrm1i, nrm1f, mu, tol, alphm, alphm1

        double precision u(nbnodes-1),v_rhs(nbnodes-1),rhs(nbnodes-1),
     &           v(nbnodes-1), svol(5)


        integer  pow_sum, i
        double precision  integ_rho, finteg_rho
	    double precision pform, form_rate, hh, reso1, reso2, reso, 
     &                   reso4, reso3
	    character*72, fichierres
	    external pow_sum, integ_rho, finteg_rho
 
        integer info, iter, exmp, iso, k
	    double precision resid, nrm0,nrm1,nrm0_exa, nrm1_exa, tab(12)

	    external psolve,  matvec
        

        write(*, *) '--> Initilization of constants...'
        call initialise_constants()

	 
        fichierres = 'res_form_fin_0.4_N60.dat'
	    open (unit = 24, file = fichierres, status = 'unknown')



        tab(1) = -0.9D0
        tab(2) = -0.5D0
        tab(3) = 0.0D0
        tab(4) = 0.5D0
        tab(5) = 0.7D0
        tab(6) = 1.D0
        tab(7) = 1.5D0
        tab(8) = 2.D0
        tab(9) = 2.5D0
        tab(10) = 3.D0
        tab(11) = 3.5D0
        tab(12) = 4.0D0

        reso1 = integ_rho(4)
    	reso2 = integ_rho(6)
	    reso3 = integ_rho(8)
        reso4 = integ_rho(10)
	    reso = integ_rho(12)
        nrm0_exa = 2.D0*pi*(reso1-3.D0*reso2+3.D0*reso3-reso4)/5.D0
	    nrm1_exa = 8.D0*pi*(reso1 -2.D0*reso2 + reso3)/3.D0 
     &             - 16.D0*nrm0_exa
     &             + 32.D0*pi*(reso1 -4.D0*reso2 + 6.D0*reso3 
     &             - 4.D0*reso4 + reso)/5.D0 
        nrm0_exa = dsqrt(nrm0_exa)
        nrm1_exa = dsqrt(nrm1_exa)
        gamma1 = 1.D0
	    gamma2 = 3.6D0 ! must be greater than 3.D0
	    exmp = 1
	    iso = 1	    
	    pform  = 2.D0
        mu = 0.5D0	   
        write(*, *) 'Ng                                  : ', Ng
        write(*, *) 'N                                   : ', N
        write(*, *) 'Number of tetraedra                 : ', nbtetra
        write(*, *) 'Number of vertices                  : ', nbnodes
!        write(*, *) 'Norme W0_{-1}  de la solution (exa): ', nrm0_exa	   
!        write(*, *) 'Norme Grad de la solution     (exa): ', nrm1_exa       
 	    write(*, *) 'Gamma                               : ', gamma1
  	    write(*, *) 'R                                   : ', pform
        write(*, *) 'mu                                  : ', mu


        write(*, *) '--> Meshing...'
        
        call mesh(N, pform, mu, xin, tetra, domnod, domtet)
        call calc_pas(nbtetra, tetra, xin, domtet, hmax,
     &                        form_rate,  5)


        hh = 0.D0
	    do i = 1, 3
	     hh = max(hh, hmax(i))
	    enddo 
	    write(*, *) ' The size h (in the unbounded domain): ', hh 
	    write(*, *) ' The size hs (in the bounded domain):', hmax(5) 

        write(*, *) '... done ...'
        write(*, *) '--> Validity of the mesh: checking volumes'   
        tol = 1.D-4             
        call test_volumes(tol, info)
        if (info.eq.1) then
           write(*, *) '... successful ...  '      
        else
          write(*, *) ' Warning (a mesh problem): '
          write(*, *) '  the sum of volumes of elements is different'
          write(*, *) '  from the total volume of the big subdomains !'
        endif 

        write(*, *) '--> Computing the exact solution'  
        write(*, *) '    and the right-hand-side (by interpolation)'
        

	    call example_and_rhs(exmp,iso, gamma1, gamma2,u,rhs)   
!         call calc_rhs(gamma1, gamma2, rhs, v_rhs)
!         write(*, *) '... done ...'
         	    
!	    som =  gamma1  + gamma2	    	    
!        juste pour tester
!	    call example_and_rhs(exmp,iso,gamma2,gamma2,v,rhs)9
!        rhs = v	  
!--------------------------------------------


        write(*, *) '--> Computing the exact solution'  
        write(*, *) '    and the right-hand-side (directly)'
         call  direct_rhs(exmp, iso, gamma1, N, Ng, nbtetra,
     &                nbnodes, tetra, xin, v_rhs)

         write(*, *) '... done ...'

!        juste pour tester       
!          reso = 0.D0
!           do i=1, nbnodes-1
!            reso = reso + v_rhs(i)*u(i)
!          enddo
!          write(*, *) reso, 4*pi*(integ_rho(6)-integ_rho(8))
   
        write(*, *) '--> Computing some useful quantities...'
        
        call grad_grad_wei(exmp, iso, gamma1, N, Ng, 
     &                 nbtetra, tetra, xin, grad_grad)
        write(*, *) '... done...'

     	do i = 1, nbnodes - 1
          v(i) = 0.D0
	    enddo
	    
        resid = 1.D-6
	    iter = 30000

	    write(*, *) '--> Solving the linear system...'
	    
!        juste pour tester	   
!	    call matvec(1.D0, u, 0.D0, v_rhs)
!	    reso=0.D0
!         do i = 1, nbnodes-1
!          reso = reso+u(i)*v_rhs(i)
!        enddo
!        alphm = 4.D0
!        alphm1 = (alphm*alphm-2.*alphm)*2/3.D0
!        
!        reso1 = (2+alphm1)*integ_rho(6) -
!     &    (2+2*alphm1+alphm*alphm*2/3.)*integ_rho(8)
!     &   + (alphm1 + alphm*alphm*4/3.)*integ_rho(10)
!     &  - alphm*alphm*(2/3.D0)*integ_rho(12)
!        reso1 = reso1*2*pi
!        write(*, *) reso, reso1
!        write(*, *) reso, 64.D0*PI*(integ_rho(8)-
!      &     2.*integ_rho(10) + integ_rho(12))       
!        write(*, *) reso, 4*PI*(integ_rho(6)-integ_rho(8))
!
!        stop
!!----------------------------------------        
        call CG(nbnodes-1, v_rhs, v,  iter, resid,
     &                 matvec, psolve, info) 
        call interpol(exmp, iso, gamma1, v)       
        write(*, *)  '... done ...'
        write(24, *) 'Fin de l''inversion: Info =', info   
        write(24, *) 'Nombre d''iterations      =', iter 
        write(*, *)  '--> Computing errors...'   
       
	    reso = 0.D0
	    reso1 = 0.D0
	    do i = 1, nbnodes - 1
          reso = reso + u(i)*u(i)
	      reso1 = reso1 + (u(i)-v(i))*(u(i)-v(i))
c 	      write(*,*) u(i), v(i)
	    enddo 
	  
        call norms_erreurs_w10(gamma1, u, v,  
     &  nrm0, nrm1, nrm0i, nrm0f, nrm1i, nrm1f,
     &  erreurL2, errGrd, erreurL2i, errGrdi,
     &                  erreurL2f, errGrdf)

	    
 	     write(*, *)  'L''erreur relative euclidienne        :', 
     &                 dsqrt(reso1/reso)     
	      write(*, *) 'L''erreur ponderee L2 relative        :', 
     &                erreurL2/nrm0_exa
	      write(*, *) 'L''erreur ponderee H1 relative        :',
     &                errGrd/nrm1_exa
 	      write(*, *) 'L''erreur ponderee L2 relative (IFEM) :', 
     &                erreurL2i/nrm0i
 	      write(*, *) 'L''erreur ponderee H1 relative (IFEM) :',
     &                errGrdi/nrm1i
 	      write(*, *) 'L''erreur ponderee L2 relative (MEF)  :', 
     &                erreurL2f/nrm0f
 	      write(*, *) 'L''erreur ponderee H1 relative (MEF)  :',
     &                errGrdf/nrm1
	    close(24)
       end

c--------------------------------------------------------------------
c          Etat d'avancement: fini
         subroutine initialise_constants()
	     implicit none	    
	     include 'constants_keltoum.h'


         integer i, j
 
	     do i = 1, 3
            cc(i) = 1.D0
	     enddo



	     do i = 1, nq
           xq(i) = 0.D0
	       yq(i) = 0.D0
	     enddo

	     xq(2) = 0.5D0
         xq(3) = 1.D0
	     xq(4) = 1.D0/3.D0
	     
	     yq(4) = 1.D0/3.D0
	     yq(5) = 0.5D0

	     xq(6) = 0.5D0
	     yq(6) = 0.5D0
 	     
	     yq(7) = 1.D0


	     weiq(1) = 1.D0/40.D0
         weiq(2) = 1.D0/15.D0
	     weiq(3) = 1.D0/40.D0
	     weiq(4) = 9.D0/40.D0
	     weiq(5) = 1.D0/15.D0
	     weiq(6) = 1.D0/15.D0
	     weiq(7) = 1.D0/40.D0 



         do j= 1, 10
          do i = 1, 3
           xyzqt(i, j)= 0.D0
          enddo 
         enddo
          
         xyzqt(1, 1)= 0.5684305841968444D0
         xyzqt(2, 1)= 0.1438564719343852D0
         xyzqt(3, 1)= 0.1438564719343852D0
          
         xyzqt(1, 2)= 0.1438564719343852D0
         xyzqt(2, 2)= 0.1438564719343852D0
         xyzqt(3, 2)= 0.1438564719343852D0
          
         xyzqt(1, 3)= 0.1438564719343852D0
         xyzqt(2, 3)= 0.1438564719343852D0
         xyzqt(3, 3)= 0.5684305841968444D0
          
         xyzqt(1, 4)= 0.1438564719343852D0
         xyzqt(2, 4)= 0.5684305841968444D0
         xyzqt(3, 4)= 0.1438564719343852D0
          
         xyzqt(2, 5)= 0.5D0
         xyzqt(3, 5)= 0.5D0
          
         xyzqt(1, 6)= 0.5D0 
         xyzqt(3, 6)= 0.5D0
          
         xyzqt(1, 7)= 0.5D0
         xyzqt(2, 7)= 0.5D0
          
         xyzqt(1, 8)= 0.5D0
          
         xyzqt(2, 9)= 0.5D0
          
         xyzqt(3, 10)= 0.5D0
             
         do j=1, 4
           weiqt(j)= 0.2177650698804054D0
         enddo
         do j=5, 10
           weiqt(j)= 0.0214899534130631D0
         enddo
        return
	   end
c--------------------------------------------------------------------
         subroutine matvec(alpha, u, beta, v)
             
c        This subroutine computes the product matrix*vector
c        Given alpha, beta, two reals, and two vectors u,  v
c        it makes v = alpha A u  + beta v
            
         implicit none               
         include 'fondements_keltoum.h'
	         

	     double precision alpha, beta, u(*), v(*)
	     integer i, j, s, k, nb,  k1, j1, l
		 double precision res


         do i = 1, nbnodes - 1
            v(i) = beta*v(i)
   	     enddo
	    open (unit = 25, file = 'essai.dat', status = 'unknown')

         do nb = 1, nbtetra
          do j = 1, 4
           i = tetra(j, nb)
	       if (i.ne.1) then
	        res = 0.D0
            do k = 1, 4
	         s = tetra(k, nb)
             if (s.ne.1) then
             j1 = max(k, j)
	         k1 = min(k, j)
	         l = (j1*(j1-1))/2 + k1
		     res = res + grad_grad(l,nb)*u(s-1)
		     endif
            enddo
             v(i-1) = v(i-1) + alpha*res

	       endif
	      enddo
	     enddo
	     close(25)
	    return
        end

        subroutine psolve(v,  b)
c       Preconditionning subroutine. Here the matrix
c       of preconditionning is the identity

        implicit none
	    include 'fondements_keltoum.h'
	    double precision  b(*), v(*)
	    integer i
        do i = 1, nbnodes-1
          v(i) = b(i)
	    enddo
	    end
c----------------------------------------------------------------------
        subroutine psolve1(v,  b)
        implicit none
        include 'fondements_keltoum.h'

	    double precision  b(*), v(*)


	    integer i, j, s, k, nb,  k1, j1, l
	    double precision reso(nbnodes), reso1(nbnodes), res


        do i = 1, nbnodes-1
         reso(i) = 0.D0
	     reso1(i) = 0.D0
	    enddo

        do nb = 1, nbtetra
         do j = 1, 4
          i = tetra(j, nb)
          if (i.ne.1) then
	      res = 0.D0
          do k = 1, 4
	      s = tetra(k, nb)
           if (s.ne.1) then
            j1 = max(k, j)
	        k1 = min(k, j)
	        l = (j1*(j1-1))/2 + k1
			res = res + grad_grad(l, nb)
           endif
          enddo
           reso(i-1) = reso(i-1) + res
	      endif
	      enddo
	      enddo

          do i = 1, nbnodes-1
          v(i) = b(i) / reso(i)
	      reso(i) = 0.D0
	      enddo
         return
  
          do nb = 1, nbtetra
           do j = 1, 4
           i = tetra(j, nb)
	        if (i.ne.1) then
	         res = 0.D0
            do k = 1, 4
	          s = tetra(k, nb)
                  if (s.ne.1.and.s.ne.i) then
                    j1 = max(k, j)
	               k1 = min(k, j)
	               l = (j1*(j1-1))/2 + k1
				  res = res + grad_grad(l,nb)*v(s-1)
                  endif
                 enddo
	             reso(i-1) = reso(i-1) + res
                  l = (j*(j-1))/2 + j
	            reso1(i-1) = reso1(i-1) + grad_grad(l, nb)
	         endif
	     enddo
	 enddo
        do i = 1, nbnodes-1
          v(i) = (b(i) - reso(i)) / reso1(i)
	  enddo

	return

	end


c--------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)
           subroutine  calc_rhs(gam1, gam2, rhs, v_rhs)
c          Compute the right hand side of the linear system
c
c          som : it is the sum of gamma1 (corresponding to test functions)
c               and gamma2 (corresponding to f)
c 
	       implicit none
           include 'fondements_keltoum.h'
	         
          
	       double precision rhs(*), v_rhs(*), gam1, gam2

	       integer i, j, s, k, nb,  k1, j1, l
	       double precision res
           double precision ww0(10, nbtetra)

    
           call ordre0(gam1, gam2,  N, Ng,  nbtetra,
     &                   tetra, xin, ww0)


           do i = 1, nbnodes - 1
              v_rhs(i) = 0.D0
	       enddo

           do nb = 1, nbtetra
            do j = 1, 4
             i = tetra(j, nb)
	         if (i.ne.1) then
	         res = 0.D0
              do k = 1, 4
	           s = tetra(k, nb)
               if (s.ne.1) then
               j1 = max(k, j)
	           k1 = min(k, j)
	           l = (j1*(j1-1))/2 + k1
	           res = res + ww0(l, nb)*rhs(s-1)  
               endif
              enddo
              v_rhs(i-1) = v_rhs(i-1) + res
	         endif
	        enddo
	       enddo
	      return
         end


c--------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)


          subroutine  grad_tetra(nbt, tetra, grad_tet, xin)
        

c        This subroutine computes the L2 product of the gradients
c        of the barycentric coordinates into a give tetrahedra
c        nbt     : an integer. The reference of the tetrahedra
c        tetra   : an array containing the references of 
c                  vertices of the all the tetrahedra in the mesh
c        xin     : an array containing the coordinates of the vertices
c        grad_tet: an array having 10 entries and containing the results 
c                  int_K (grad lamda_i).(grad lambda_j) dx = grad_tet(s)
c                  with s = ((j-1)*j)/2 + i
c                  For this calculus, the gradient are given in terms 
c                  of gradients of the barycentric coordinates in the reference
c                   tetrahedra

          implicit none

	   double precision grad_tet(10),grad_r(3, 4),xin(3,*)
	   integer i, j, s, kk, tetra(4, *), nb1, nbt, nb

	   double precision bk(3,3), ck(3,3), bbk(3,3), 
     &   vec(3), detK, res, det
	   external det

       call grad_ref(grad_r)

	   nb1 = tetra(1, nbt)
	   do j = 1, 3
	     nb = tetra(j+1, nbt)
	    do i = 1, 3
	     bk(i, j) = xin(i,nb) - xin(i, nb1)
	    enddo  
	   enddo    

        detK = det(bk)
	   call inversion(bk, bbk)
	   call transpos(3, 3, bbk, ck)
       call produitAAT(ck, bbk)	
			 	    
  	   do j = 1, 4
	     call mat_v3(j, bbk, grad_r, vec)
	     do i = 1, j 
	        res = 0.D0
	        do kk = 1, 3
               res = res + grad_r(kk, i)*vec(kk)
	        enddo
	          s = ((j-1)*j)/2 + i
              grad_tet(s) = res*detK/6.D0         
	     enddo
       enddo
	   
	   return
        end
c--------------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)


           subroutine grad_ref(grad)

c             This subroutine computes the gradients
c             of the barycentric coordinates
c             lambda1,...,lambda4 (affine functions)
c             in the reference tetrahedra
c             


		   implicit none 
 	       double precision grad(3, 4)
		                
		     grad(1, 1) = - 1.D0
             grad(2, 1) = - 1.D0
             grad(3, 1) = - 1.D0 
          
	         grad(1, 2) = 1.D0
             grad(2, 2) = 0.D0
             grad(3, 2) = 0.D0

	         grad(1, 3) = 0.D0
             grad(2, 3) = 1.D0
             grad(3, 3) = 0.D0

	         grad(1, 4) = 0.D0
             grad(2, 4) = 0.D0
             grad(3, 4) = 1.D0
             return
	      end       
  


c     ------------------------------------------------------------
c       Etat d'avancement: fini
         subroutine jacob_phi(k, nmu, jaco)
 !         computes the matrix I - 2*c*nmu^t
 !         at point nmu(.,k)
	     implicit none
	     double precision nmu(3, *), jaco(3, 3)
		 integer i, j, k


		 do i = 1, 3
	      do j= 1, 3
	       jaco(i, j) =  -2.D0*nmu(j, k)
          enddo
           jaco(i, i) = 1.D0 - 2.D0*nmu(i, k)
		 enddo   
         return
         end
c     ------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)


         subroutine bt_big_simp(dom, bt)

c          This subroutine computes the
c          matrix BT of the a domain "dom"
c          considered as a tetrahedra

	     implicit none 
         include 'simplexes_keltoum.h'
	     
         integer          dom, j, i
	     double precision bt(3,3)

   	     do j = 1, 3
	      do i = 1, 3 
  	       BT(i,j) = vertex(i,j+1,dom)-vertex(i,1, dom) 
  	      enddo
	     enddo
         return
	     end


c----------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)
         subroutine triangles(ni, mu, absdete, nmu)
	      
         implicit none
	     include 'constants_keltoum.h'

         double precision  mu(3,3,2), absdete, nmu(3,*)            
         integer i, ni
         double precision alpha(2), beta(2), res
       
          alpha(1) = mu(1, 2, ni) - mu(1, 1, ni)
		  alpha(2) = mu(2, 2, ni) - mu(2, 1, ni)
		  beta(1) =  mu(1, 3, ni) - mu(1, 1, ni)
		  beta(2) =  mu(2, 3, ni) - mu(2, 1, ni)
           
 	      res = alpha(1)*beta(2) - alpha(2)*beta(1)
	      absdete = dabs(res)            
          do i =  1, nq  
		   nmu(1, i) = alpha(1)*xq(i)+beta(1)*yq(i)+mu(1,1, ni) 
		   nmu(2, i) = alpha(2)*xq(i)+beta(2)*yq(i)+mu(2,1, ni)
		   nmu(3, i) = 1.D0 - (nmu(1,i) + nmu(2, i))		                       
          enddo
	      return
	 end
c-------------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)


        subroutine delta_extremas(x,  dimin, dimax)
        
c        This subroutine computes the minimal and the maximal values
c        of the function delta(x) = x_1 + x_2 + x_3 on the tetrahedra
c        x (x denotes the matrix of coordinates of the vertices) 



         implicit none 
	     double precision dimin, dimax, x(3, *)

  	     integer i, j
	     double precision delt
	   

	     dimax = 0.D0
	     dimin = 100.D0
	     do j = 1, 4
	      delt = 0.D0
	      do i = 1, 3
           delt  =  delt + x(i, j) 
          enddo
          dimin = min(dimin, delt)
          dimax = max(dimax, delt)
	     enddo 
         return
       end
c-------------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)

          subroutine intersection(xx, delt, mu, nbr)
c           This subroutines computes the intersection
c           of the plane x+y+z=delt with the tetrahedra
c           xx (xx is an array containing the coordinates
c           of the vertices). 
c           INPUT:
C              xx   : a matrix containing the coordinates of the vertices
c                     of the tetrahedra
c              delt : a real number 
c            OUTPUT       
c              nbr  : the number of the vertices of the intersection
c              (0: void 1: a point, 2: a segment  3: a triangle  4: a quadrilateral) 
c               mu   : a matrix containing the values of x/delta at the vertices of 
c                      the intersection    
          implicit none
          include 'constants_keltoum.h'

     	  integer nbr
          double precision delt, mu(3, *), xx(3,*)
          
	      integer i, m, j, it
	      double precision x1(3), x2(3), delt1, delt2, desc, tau
          

   	      if (dabs(delt).lt.zero) then
           nbr = 0
	      return
	      endif
	      it = 0
	      do i = 1, 4
	       delt1 = 0.D0
	       do m = 1, 3
            x1(m) = xx(m, i) 
 	        delt1  = delt1 + x1(m)
	       enddo
           if (dabs(delt1-delt).lt.zero) then
             it = it + 1
             do m = 1, 3
               mu(m, it) =  x1(m)/delt
	         enddo
	       endif
	      enddo
	      do i = 1, 3
	       delt1 = 0.D0
	       do m = 1, 3
            x1(m) = xx(m, i) 
 	        delt1  = delt1 + x1(m)
	       enddo
	      do j =  i+1, 4
	        delt2 = 0.D0
		    do m = 1, 3
             x2(m) = xx(m, j)
		     delt2 = delt2 + x2(m) 
	        enddo
            desc = (delt - delt1)*(delt-delt2) 
	        if (desc.le.(-zero)) then
	         tau = (delt - delt2)/(delt1 -delt2)
	         it = it + 1
             do m = 1, 3
               mu(m, it) = (tau*x1(m) + (1.D0-tau)*x2(m))/delt
	         enddo
	        endif
	     enddo
	   enddo
       nbr = it	
	   return
 	  end
c-------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)

          subroutine div_quadrilatere(mu, mu_n)
c         This subroutine divides a quadrilateral
c         into two triangles. 
c         INPUT :
c         mu   :  a matrix containing the values of
c               x/(x_1+x_2+x_3) at the four vertices     
c         mu_n : an array containing the values of
c               x/(x_1+x_2+x_3) at the three vertices
c               of each triangle
      
	      implicit none
	      double precision  mu(3, *), mu_n(3, 3, 2)
	      integer i, j, k, l, rep

	      i = 1
          j = 2
	      k = 3
	      l = 4
	      call div_quadri(i, j, k, l, rep, mu, mu_n)
          if (rep.eq.-1) call div_quadri(l, i, j, k, rep, mu, mu_n)
          if (rep.eq.-1) call div_quadri(k, l, i, j, rep, mu, mu_n)
	      if (rep.eq.-1) call div_quadri(j, k, l, i, rep, mu, mu_n)
	      if (rep.eq.-1) then
	      write(*, *) 'Impossible subdvision of a quadrilateral!'
	     stop
	    endif
	   return
	  end
c-------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)

         subroutine div_quadri(i, j, k, l, rep, mu, mu_n)  
c        This subroutine subdivide a quadrilateral not necessarily convex
c        into two triangles if the point i is not inside the triangle (j, k, l)
c        INPUT
c           i, j, k, l: four integers. {i, j, k, l} = {1, 2, 3, 4}
c           mu   : a matrix containing the values of x/delta at the vertices of 
c                      the quadrilateral
c         OUTPUT
c           rep  : an integer. rep=0 if the subdivision is done 
c                              rep = -1 if the vertex i is inside the triangle
c                              (j, k, l)
 
c           mu_n  : a matrix containing the values of x/delta at the vertices of 
c                      the quadrilateral
         implicit none
         include 'constants_keltoum.h'

         integer i, j, k, l, rep, s
  	     double precision mu(3, 4), mu_n(3, 3, 2)
        

         double precision b(2, 2), y(2), detb, ib(2, 2), 
     &                   res1, res2, res3, prod
 
 
         rep = - 1
	     do s = 1, 2
          b(s, 1) = mu(s, k) - mu(s, j)
		  b(s, 2) = mu(s, l) - mu(s, j) 
	      y(s)    = mu(s, i) - mu(s, j) 
	     enddo
         
	     detb = b(1,1)*b(2,2) - b(1,2)*b(2, 1)

         if (dabs(detb).lt.zero) return

	     ib(1, 1) = b(2,2)/detb
	     ib(2, 2) = b(1,1)/detb
	     ib(1, 2) = - b(1,2)/detb
	     ib(2, 1) = - b(2,1)/detb 

         res1 = ib(1,1)*y(1) + ib(1,2)*y(2) 
	     res2 = ib(2,1)*y(1) + ib(2,2)*y(2)
         res3 = 1.D0 - (res1 + res2)

         prod = res1*res2*res3

         if  (prod.lt.(-zero)) then
	     if (res1.lt.0D0) then
	      rep = 0
	      do s = 1, 3
             mu_n(s, 1, 1) = mu(s, i)
             mu_n(s, 2, 1) = mu(s, j)
             mu_n(s, 3, 1) = mu(s, k)

             mu_n(s, 1, 2) = mu(s, i)
             mu_n(s, 2, 2) = mu(s, l)
             mu_n(s, 3, 2) = mu(s, k)
	       enddo			   
		 else if (res2.lt.0D0) then
	      rep = 0
	      do s = 1, 3
              mu_n(s, 1, 1) = mu(s, i)
              mu_n(s, 2, 1) = mu(s, k)
              mu_n(s, 3, 1) = mu(s, l)

              mu_n(s, 1, 2) = mu(s, i)
              mu_n(s, 2, 2) = mu(s, j)
              mu_n(s, 3, 2) = mu(s, l)
	      enddo				 
		 else if (res3.lt.0D0) then
	       rep = 0
		   do s = 1, 3
	        mu_n(s, 1, 1) = mu(s, i)
            mu_n(s, 2, 1) = mu(s, k)
            mu_n(s, 3, 1) = mu(s, j)
            mu_n(s, 1, 2) = mu(s, i)
            mu_n(s, 2, 2) = mu(s, l)
            mu_n(s, 3, 2) = mu(s, j)
	       enddo
          endif
         endif
	   return	   
	  end   
c
c-------------------------------------------------------------------
c          Etat d'avancement: fini (aucun changement)

       subroutine div_quadrilateral_cvx(mu, mu_n)
c       This subroutine divides a convex quadrilateral
c       into two triangles. 
c       INPUT :
c       mu   :  a matrix containing the values of
c              x/(x_1+x_2+x_3) at the four vertices     
c       mu_n : an array containing the values of
c              x/(x_1+x_2+x_3) at the three vertices
c              of each triangle
      
	    implicit none
	    double precision  mu(3, *), mu_n(3, 3, 2)
	    integer s

	      
	     do s = 1, 3
             mu_n(s, 1, 1) = mu(s, 1)
             mu_n(s, 2, 1) = mu(s, 2)
             mu_n(s, 3, 1) = mu(s, 3)

             mu_n(s, 1, 2) = mu(s, 1)
             mu_n(s, 2, 2) = mu(s, 4)
             mu_n(s, 3, 2) = mu(s, 3)
	      enddo	
	      return
	     end
    
c----------------------------------------------------------------------
            subroutine example_and_rhs(exmp,iso, gamma1,gamma2,u,rhs)
            
            implicit none
            include 'simplexes_keltoum.h'
	        include 'fondements_keltoum.h'

            integer exmp, iso
		    double precision gamma1, gamma2,  u(*), rhs(*)
              
            integer i,dom, k
            double precision  rr, x(3), vx(3), res, val_u, val_rhs
	        external rr

 
		    do i = 2, nbnodes
			 dom = domnod(i)
	         do k = 1, 3
                x(k) = xin(k, i)
	         enddo
	         
	         if (dom.lt.5.and.dom.gt.0) then
	           call fphi(dom, x, vx)
	           res =  rr(dom,  vx)
               call ex_function(exmp, iso, vx, val_u, val_rhs)
               u(i-1) = val_u*(res**gamma1) 
		       rhs(i-1) = val_rhs*(res**gamma2) 
	         else if (dom.eq.5) then
              call ex_function(exmp, iso, x, val_u, val_rhs)
              u(i-1) = val_u
              rhs(i-1) = val_rhs                       
	         endif
	        enddo
	        return
	       end 	    
	    	    
c---------------------------------------------------------------------
c      Etat d'avancement: fini (aucun changement)		     
         double precision function rr(dom, x)
		   implicit none	        
	       include 'simplexes_keltoum.h'


		   integer dom, i
		   double precision x(3), res
		   
		   res = 0.D0
	       do  i = 1, 3
		    res = res + hei(i,dom)*x(i) 
		   enddo
		   rr = res/(nrm_h(dom)*nrm_h(dom))
		   return
         end	
c-------------------------------------------------------------------
c       Etat d'avancement: fini (aucun changement)
        subroutine fphi(dom,  x, vx)

        implicit none
        include 'simplexes_keltoum.h'

	    integer dom, i
        double precision vx(3), x(3), rr, res
        external rr
          
        res = rr(dom, x)

    	 do i = 1,  3
	      vx(i) = x(i)/(res*res)
         enddo
        return 
	    end	


c-------------------------------------------------------------------
          double precision function integ_rho(m)
          implicit none
c          cette subroutine calcule les integrales de la forme
c          I_m = int_0^+infini rho^{-m} avec rho = sqrt(1 + r*r)   
c          On a : m*I_{m+2} = (m-1)*I_m
           integer          m, k
	       double precision res,  pi, tab(m)
	       parameter (pi = 3.14159265358D0)
		  		   
		  		   
		   tab(2) = PI/2.D0
		   tab(3) = 1.D0
		   
		   do k = 4, m
		     tab(k) = (k-3)*tab(k-2)/(k-2)
		   enddo		   
		   
		   integ_rho = tab(m)
		   return 
		   select case(m)
		   case(2)
	        res = PI/2.D0
		   case(3)
		    res = 1.D0
           case(4)
              res = PI/4.D0
		   case(5)
		    res =  2.D0/3.D0
		   case(6) 
		    res = 3.*PI/16.D0
		   case(7)
		    res = 8.D0/15.D0
		   case(8)
		    res = 5*PI/32.D0
		   case(9)
		    res = 16.D0/35.D0
	       case(10)
              res = 35*PI/(32*8.D0)
           case(11)
              res = 8*16.D0/(9*35.D0)
	       case(12)
		    res = 9*35.D0*PI/(320.D0*8)		  
		   end select
	       integ_rho = res
 		  return
	     end 

		   
		    

