c     Etat de modification : fini
       subroutine ordre0(gam1, gam2, N, Ng,  nbtetra,
     &                   tetra, xin, ww0)
	   implicit none	        
	   include 'simplexes_keltoum.h'
	   include 'constants_keltoum.h'

       integer  N, nbtetra,  tetra(4, *), Ng   
       double precision xin(3,*), ww0(10, *), gam1, gam2, agam
       integer l, i, j, s,  nm1, nm2, nb, pow_sum,  nbtet
	   external pow_sum
	   double precision BT(3,3),IBT(3,3),IBT_T(3,3), IBKT(3,3),
     &	   IBKT_T(3,3), VTX0(3), BK(3,3), IBK(3, 3)


 	   double precision gradref(3, 4), delt, mu(3, 4)
	   double precision  res, res_g, res40, res41, res42,
     &   ress(4), xh(3), delt0, delt1, delt2
       double precision xx(3,4),vc0(3),  wbase(3), 
     &  detBT, nmu(3,nq), mu_n(3,3, 2),dimin,dimax	  
	   double precision  va(3, 4), val0(4), idelt, 
     &	               det, rval0(4), absdete, prd, prd1
	   double precision rwgau(Ng),  rdgau(Ng), wgau(Ng), dgau(Ng),
     &                  w0ref(4, 4)
	   integer nbprisme, k, nm, m, ss, nnb, kk, nbr, ni, i1, k1
	   external det
              
	   do i = 2, 4
	    rval0(i) = 0.D0
	   enddo
	
	   rval0(1) = 1.D0
	   
	   nbprisme = pow_sum(2,N-1) + pow_sum(1,N-1)
	   nbtet  = 3*nbprisme + N
	   
       agam = 0.5D0*(gam1+gam2)
       
       call grad_ref(gradref)
	   call ordre0_ref(w0ref)
	   
	   call gauss_lobatto(rdgau, rwgau, Ng)
	   
	   do s = 1, 4
	    call bt_big_simp(s, BT)
        call inversion(BT, IBT)
	    detBT = det(BT)
	    call transpos(3, 3, IBT, IBT_T)
	    
        do nnb = 1, nbtet
         nb = nnb + (s-1)*nbtet
         do l = 1, 10
          ww0(l, nb) = 0.D0
	     enddo
!        One vertex of the tetrhedra

	     nm1 = tetra(1, nb)
	     do i = 1, 3
           VTX0(i) = xin(i, nm1)
	     enddo

!        The matrix of the tetrahedra   
	     do k = 2, 4
           nm =  tetra(k, nb)
           do i = 1, 3
            bk(i, k-1) = xin(i,nm) - VTX0(i)
           enddo
         enddo

   
         call inversion(BK, IBK)       
         call matmat(3, 3, 3, IBK, BT, IBKT)
	     call transpos(3, 3, IBKT, IBKT_T)    
	     call matv3(IBK, VTX0, vc0)


	     kk = 1
         if (nnb.eq.1) kk = 2		     
         do k = 1, 4
  	       call nmat_v3(k, IBKT_T, gradref, va) 
           res = rval0(k)
	       do m = 1, 3
	        res = res - gradref(m, k)*vc0(m)
	       enddo
           val0(k) =   res
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
	     call dec_gauss_lobatto(dimin, dimax, Ng, rwgau, rdgau,
     &                         wgau, dgau)     
	     do ss = 1, Ng
          delt = dgau(ss)
          delt0 = 0.D0
	      if ((dabs(agam-2.D0).lt.zero)) delt0 = 1.D0
	     
	      if (dabs(delt).gt.zero) then
            if ((dabs(agam-2.D0).gt.zero)) delt0 = delt**(agam-2.D0)
            idelt = 1./delt
          else 
            delt0 = 0.D0
            idelt = 0.D0
          endif

c        Computing the intersection of plane x+y+z=delt
c        with the tetrahedra whose vertices are given by xx

          call intersection(xx, delt, mu, nbr)
	      ni = 0
          if (nbr.eq.3) then
c         The intersection is a triangle           
	      ni = 1
          do k1 = 1, 3
	       do i1 = 1, 3
             mu_n(i1, k1, 1) = mu(i1, k1)
		   enddo
	      enddo
		  else if (nbr.eq.4) then
c         The intersection is a quadrilateral			 
		  ni =  2
	      call  div_quadrilatere(mu, mu_n)
 	      else if (nbr.lt.3.and.nbr.ge.0) then
	      res_g = 0.D0
          ni = 0
	      else   
	        write(*, *) ' A bad intersection between a
     &	               plane and a tetrahedra.'
	        write(*,*)  ' (What hapens here is 
     &                necessarily a mistake)' 
	        write(*, *) ' The program is stopped 
     &                 Please correct.'
            stop
	      endif	
	      
	      prd1 = wgau(ss)*dabs(detBT)	
          do i1 = 1, ni
	       call triangles(i1, mu_n, absdete, nmu)
           do i = 1, nq
		    prd = weiq(i)*absdete*prd1
		   
		    do k = 1, 4
	         call prodcolumns(i, k, nmu, va, res)
 	         wbase(k) = (val0(k) + res*delt)*delt0     	          
		    enddo
		    
            do k = kk, 4
            do  j = k, 4
	         res_g = prd*wbase(k)*wbase(j)
             l = ((j-1)*j)/2 + k
	         ww0(l,nb) = ww0(l,nb) + res_g
	        enddo
	       enddo
	      enddo
         enddo

      enddo
	  enddo
      enddo
c    
c      les elements finis
c
      s = 5
	  do nb = (s-1)*nbtet + 1, s*nbtet
	   nm1 = tetra(1, nb)
 	    do j = 2, 4
	     nm2 = tetra(j, nb)
		 do i = 1, 3
	      bk(i, j-1) = xin(i,nm2) - xin(i, nm1)
		 enddo 
	    enddo  
	    absdete = det(Bk)
	    absdete  = dabs(absdete)
        do k = 1, 4
	     do j = k, 4
	      l = ((j-1)*j)/2 + k
	      ww0(l, nb) = w0ref(k,j)*absdete
	     enddo
        enddo
        enddo
	   return
      end
c---------------------------------------------------------------

            subroutine ordre0_ref(w0ref)
            implicit none
            double precision w0ref(4, 4)
            integer i, j

	        do i = 1, 4
             w0ref(i, i) = 1.D0/60.D0
	        enddo

            do i = 1, 3
              do j = i+1, 4
	           w0ref(i, j) = 1.D0/120.D0
	           w0ref(j, i) = 1.D0/120.D0
              enddo
	         enddo
            return
	       end
	       
	       
	       
	       	    
c---------------------------------------------------------------	    
       subroutine direct_rhs(exmp, iso, gam1, N, Ng, nbtetra,
     &                nbnodes, tetra, xin, v_rhs)
	   implicit none	        
	   include 'simplexes_keltoum.h'
	   include 'constants_keltoum.h'

       integer  N, nbtetra, tetra(4,*), Ng, exmp, iso, nbnodes
       double precision xin(3,*), v_rhs(*), gam1, gam2, agam
       integer l, i, j, s,  nm1, nm2, nb, pow_sum,  nbtet
	   external pow_sum
	   double precision BT(3,3),IBT(3,3),IBT_T(3,3), IBKT(3,3),
     &	 IBKT_T(3,3), VTX0(3), BK(3,3), IBK(3, 3), IBK_T(3, 3),
     &   xp(3), vxp(3), xxp(3)


 	   double precision gradref(3, 4), delt, mu(3, 4)
	   double precision  res0,  res, res_g, res40, res41, res42,
     &   ress(4), xh(3), delt0, delt1, delt2, absdet_bis
       double precision xx(3,4), vc0(3),  
     &  detBT, nmu(3,nq), mu_n(3,3, 2),dimin,dimax	  
	   double precision  va(3, 4), val0(4), idelt, 
     &	               det, rval0(4), absdete, prd, prd1
	   double precision rwgau(Ng),  rdgau(Ng), wgau(Ng), dgau(Ng),
     &                  w0ref(4, 4), valu, valrhs, rr
	   integer nbprisme, k, nm, m, ss, nnb, kk, nbr,ni,i1, k1
	   external det, rr

	   rval0(1) = 1.D0              
	   do i = 2, 4
	    rval0(i) = 0.D0
	   enddo
	   
	   do i = 1, nbnodes - 1
         v_rhs(i) = 0.D0
	   enddo

	   
	   nbprisme = pow_sum(2,N-1) + pow_sum(1,N-1)
	   nbtet  = 3*nbprisme + N

       
       call grad_ref(gradref)
	   call ordre0_ref(w0ref)
	   
	   call gauss_lobatto(rdgau, rwgau, Ng)

	   do s = 1, 4
	    call bt_big_simp(s, BT)
        call inversion(BT, IBT)
	    detBT = det(BT)
	    call transpos(3, 3, IBT, IBT_T)
	    
        do nnb = 1, nbtet
         nb = nnb + (s-1)*nbtet

!        One vertex of the tetrhedra

	     nm1 = tetra(1, nb)
	     do i = 1, 3
           VTX0(i) = xin(i, nm1)
	     enddo

!        The matrix of the tetrahedra   
	     do k = 2, 4
           nm =  tetra(k, nb)
           do i = 1, 3
            bk(i, k-1) = xin(i,nm) - VTX0(i)
           enddo
         enddo

   
         call inversion(BK, IBK)       
         call matmat(3, 3, 3, IBK, BT, IBKT)
	     call transpos(3, 3, IBKT, IBKT_T)    
	     call matv3(IBK, VTX0, vc0)




	     kk = 1
         if (nnb.eq.1) kk = 2		     
         do k = 1, 4
  	       call nmat_v3(k, IBKT_T, gradref, va) 
           res = rval0(k)
	       do m = 1, 3
	        res = res - gradref(m, k)*vc0(m)
	       enddo
           val0(k) =   res
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
	     call dec_gauss_lobatto(dimin, dimax, Ng, rwgau, rdgau,
     &                         wgau, dgau)     
	     do ss = 1, Ng
          delt = dgau(ss)
          delt0 = 0.D0
	      if ((dabs(gam1-4.D0).lt.zero)) delt0 = 1.D0
	     
	      if (dabs(delt).gt.zero) then
          if ((dabs(gam1-4.D0).gt.zero)) delt0 = delt**(gam1-4.D0)
            idelt = 1./delt
          else 
            delt0 = 0.D0
            idelt = 0.D0
          endif

c        Computing the intersection of plane x+y+z=delt
c        with the tetrahedra whose vertices are given by xx

          call intersection(xx, delt, mu, nbr)
	      ni = 0
          if (nbr.eq.3) then
c         The intersection is a triangle           
	      ni = 1
          do k1 = 1, 3
	       do i1 = 1, 3
             mu_n(i1, k1, 1) = mu(i1, k1)
		   enddo
	      enddo
		  else if (nbr.eq.4) then
c         The intersection is a quadrilateral			 
		  ni =  2
	      call  div_quadrilatere(mu, mu_n)
 	      else if (nbr.lt.3.and.nbr.ge.0) then
	      res_g = 0.D0
          ni = 0
	      else   
	        write(*, *) ' A bad intersection between a
     &	               plane and a tetrahedra.'
	        write(*,*)  ' (What hapens here is 
     &                necessarily a mistake)' 
	        write(*, *) ' The program is stopped 
     &                 Please correct.'
            stop
	      endif	
	      
	      prd1 = wgau(ss)*dabs(detBT)	
          do i1 = 1, ni
	       call triangles(i1, mu_n, absdete, nmu)
           do i = 1, nq
		    prd = weiq(i)*absdete*prd1
		    do m = 1, 3
              xp(m) = delt*nmu(m, i) 
            enddo 		
            call matv3(BT, xp, xxp)   
 		    call fphi(s, xxp, vxp)
	        res =  rr(s,  vxp)
		    call ex_function(exmp, iso, vxp, valu, valrhs)  
		   
            do k = kk, 4
             call prodcolumns(i, k, nmu, va, res)
             res0 = (val0(k) + res*delt)*delt0
             m = tetra(k, nb)
             v_rhs(m-1) =(v_rhs(m-1) + prd*res0*valrhs)
	       enddo
	      enddo
         enddo

      enddo
	  enddo
      enddo
c    
c      les elements finis
c
          
       s = 5
       
       do k1 = 1, 3
	       do i1 = 1, 3
             mu_n(i1, k1, 1) = 0.D0
		   enddo
	   enddo
	   mu_n(1, 1, 1) = 1.D0    	
       mu_n(2, 2, 1) = 1.D0
       mu_n(3, 3, 1) = 1.D0

	   call triangles(1, mu_n, absdet_bis, nmu)
 
       call dec_gauss_lobatto(0.D0, 1.D0, Ng, rwgau, rdgau,
     &                         wgau, dgau)  
       
	  do nb = (s-1)*nbtet + 1, s*nbtet
	   nm1 = tetra(1, nb)
	   	     
	     do i = 1, 3
           VTX0(i) = xin(i, nm1)
	     enddo
	   
 	    do j = 2, 4
	     nm2 = tetra(j, nb)
		 do i = 1, 3
	      bk(i, j-1) = xin(i,nm2) - VTX0(i) 
		 enddo 
	    enddo  
	    absdete = det(Bk)
	    absdete  = dabs(absdete)
   
        call inversion(BK, IBK)       
	    call transpos(3, 3, IBK, IBK_T) 
	       	    
	    call matv3(IBK, VTX0, vc0)

	     do k = 1, 4
 	       call nmat_v3(k, IBK_T, gradref, va) 
           res = rval0(k)
	       do m = 1, 3
	        res = res - gradref(m, k)*vc0(m)
	       enddo
           val0(k) = res
	     enddo
	    
	     do ss = 1, Ng
          delt = dgau(ss)
	      prd1 = wgau(ss)*absdete*delt*delt	
          do i = 1, nq
		    prd = weiq(i)*absdet_bis*prd1
		    do m = 1, 3
              xp(m) = delt*nmu(m, i) 
            enddo 		
            
            call matv3(BK, xp, vxp)   
            vxp = vxp + VTX0
		    call ex_function(exmp, iso, vxp, valu, valrhs)  
		   
            do k = 1, 4
             res = rval0(k)
	         do m = 1, 3
	          res = res + gradref(m, k)*xp(m)
	         enddo    
!             call prodcolumns(i, k, nmu, va, res)
!             res0 = val0(k) + res*delt
             m = tetra(k, nb)
             v_rhs(m-1) = v_rhs(m-1) + prd*res*valrhs
	        enddo
	      enddo
         enddo
	    
	  enddo 

	   return
      end

	    
	    
	    

	       
	       
	       
	       
	       
