        subroutine interpol(exmp, iso, gamma, v)
!       this subroutine stocks the value of exact 
!       and approach solution in filles in order
!       to visualize them          
        implicit none 
        include 'fondements_keltoum.h'
        include 'simplexes_keltoum.h'  
        include 'constants_keltoum.h'
           
        double precision  vec(3), verx(3)
        double precision valu, valrhs, res1, ress, val_v
        double precision v(*), loc_interp_ifem, pas
        double precision res2, res, gamma, loc_interp_fem
        integer nbr, dom, k, i, j, nm
        integer l,exmp, iso, d, nn
        external rr
        character*72, fichier2, fichier3
        character*72, fichierx2, fichierx3

        d = 20
        nn = 3
        pas = 0.2D0
        
        call zero_vect(3, vec)
!         vec(nn) = -2.5D0
!       file contient values of the exact 
!       solution (second composente) 
!       at evry value of vec
        fichierx2 = 'exact2.dat'
        open (unit = 53, file = fichierx2, status = 'unknown')
         do k = 1, d
         vec(nn) = vec(nn) + pas
         call ex_function(exmp, iso, vec, valu, valrhs )
         write(53,*) vec(nn), valu
         enddo
        close(53)
!       file contient values of the approach  
!       solution (second composente) 
!       at evry value of vec        
        fichier2 = 'interpol2.dat'
        open (unit = 43, file = fichier2, status = 'unknown')
        call zero_vect(3, vec)
!        vec(nn) = -2.5D0
        do k = 1, d
         vec(nn) = vec(nn) + pas
         call findx(vec, dom, nbr)
         if (dom.lt.5) then
          val_v = loc_interp_ifem(gamma, dom, nbr, v, vec)
          write(43,*) vec(nn), val_v, nbr, dom         
          
         else if(dom.eq.5) then         
          val_v = loc_interp_fem(nbr, v, vec)
          write(43,*) vec(nn), val_v, nbr, dom
         endif         
         enddo          
         close(43)
         
        return
        end
     
               
c------------------------------------------------ 
        subroutine findx(x, dom, nbr)
!       This subroutine finds the the number 
!       of domain and tetrahedra in which
!       x belongs         
        include 'fondements_keltoum.h'
        include 'simplexes_keltoum.h'        
        double precision x(3)
        integer nbr, s, nbt, dom, i 
               
        call find_big_domain(vertex, x, s)
        call find_num_tetr(N,xin,x,s,tetra,nbr)
        dom = s

        return
        end
c------------------------------------------------------------
         subroutine find_num_tetr(N,xin,x,num,tetra,nbt)
!        this subroutine finds the number of tetrahedra
!        in which x belongs
         include 'constants_keltoum.h' 
         double precision x(3), ibk(3,3)
         double precision vec(3), vec1(3), bk(3,3),rr
         double precision res, vx(3),xin(3,*)
         integer num, i,j,k, nb, nbt, nbtet,pow_sum,nm1
         integer  nm, tetra(4,*), nbnodes, N
         external pow_sum, rr
        
        
	     nbprisme = pow_sum(2,N-1) + pow_sum(1,N-1)
	     nbtet    = 3*nbprisme + N        
!        unbounded domain
         if (num.lt.5) then
	      call fphi(num, x, vx)         
          do nb = (num-1)*nbtet + 1, num*nbtet
           nm1 = tetra(1, nb)  
!	       do i = 1, 3
!		     VTX0(i) = xin(i, nm1) 
!	       enddo  
	             
	       do j = 2, 4
            nm =  tetra(j, nb)
            do i = 1, 3 
             bk(i,j-1) = xin(i,nm) - xin(i, nm1) 
            enddo
           enddo
           call inversion(bk, ibk)
           do i = 1, 3
            vec(i) = vx(i) - xin(i, nm1)
           enddo
           call matv3(ibk, vec, vec1)
           res= 0.D0
           do i = 1, 3
            res = res +vec1(i)
           enddo           
           
           if ((vec1(1).ge.-zero).and.(vec1(2).ge.-zero).and.
     &       (vec1(3).ge.-zero).and.(res.le.1.D0+zero) ) then
             nbt = nb                       
            return
           endif                  
          enddo
!        bounded domain          
         else if(num.eq.5)then
          do nb = (num-1)*nbtet + 1, num*nbtet 
           nm1 = tetra(1, nb)	      
!	       do i = 1, 3
!		    VTX0(i) = xin(i, nm1)  
!	       enddo        
	       do j = 2, 4
            nm =  tetra(j, nb)
            do i = 1, 3 
             bk(i,j-1) = xin(i,nm) - xin(i, nm1)           
            enddo
           enddo
           call inversion(bk, ibk)
           do i = 1, 3
            vec(i) = x(i) - xin(i, nm1)
           enddo
           call matv3(ibk, vec, vec1)
           res= 0.D0
           do i = 1, 3
            res = res +vec1(i)
           enddo
           if ((vec1(1).ge.-zero).and.(vec1(2).ge.-zero).and.
     &       (vec1(3).ge.-zero).and.(res.le.1.D0+zero) ) then
            nbt = nb
            return
           endif                  
          enddo
          endif
          return
          end     
c------------------------------------------------------------------
         subroutine find_big_domain(vertex, x, nm)
!        this subroutine finds the domain nm 
!        (nm =1, 2, 3, 4, 5) in which
!        the point x belongs to

         double precision vertex(3,4,5), x(3), IBT(3,3)
         double precision vec(3), vec1(3), BT(3,3),detBT
         double precision det, res, vx(3)
         integer i, j, k, nm
         external det

         
         do k= 1, 4
          call fphi(k, x, vx)
          call bt_big_simp(k, BT)
          call inversion(BT, IBT)

          
          do i = 1, 3
           vec(i) = vx(i) - vertex(i,1,k)
          enddo
    
          call matv3(IBT, vec, vec1)

              
          res= 0.D0
          do i = 1, 3
           res = res +vec1(i)
          enddo
         
          if ((vec1(1).ge.0.D0).and.(vec1(2).ge.0.D0).and.
     &       (vec1(3).ge.0.D0).and.(res.le.1.D0) ) then
           nm = k
           return
          endif
         enddo
         k = 5
          call bt_big_simp(k, BT)
          call inversion(BT, IBT)
          
          do i = 1, 3
           vec(i) = x(i) - vertex(i,1,k)
          enddo
    
          call matv3(IBT, vec, vec1)
              
          res= 0.D0
          do i = 1, 3
           res = res +vec1(i)
          enddo
         
          if ((vec1(1).ge.0.D0).and.(vec1(2).ge.0.D0).and.
     &       (vec1(3).ge.0.D0).and.(res.le.1.D0) ) then
           nm = k
           return
          endif

         return
        end
        
        
c----------------------------------------------------------
       function loc_interp_fem(nbr, v, x)
!       This function computes the local interpolation
!       in bounded domain (fem) at x 
        implicit none
        
        include 'fondements_keltoum.h'
        include 'simplexes_keltoum.h'          
        include 'constants_keltoum.h'
        
        double precision x(3), v(*), loc_interp_fem
        integer nbr, nm1, nm, j, i
        
        double precision bary_lambda(4), vec1(3), BK(3,3),
     &                   IBK(3, 3), value, sum, vec(3)
                
        

           nm1 = tetra(1, nbr)	      
	       do j = 2, 4
            nm =  tetra(j, nbr)
            do i = 1, 3 
             bk(i,j-1) = xin(i,nm) -  xin(i, nm1)  
            enddo
           enddo
           
           call inversion(bk, ibk)
                  
           do i = 1, 3
            vec(i) = x(i) - xin(i, nm1)  
           enddo
           call matv3(ibk, vec, vec1)
           
           sum = 0.D0
           
           do i = 2, 4
            bary_lambda(i) = vec1(i-1)
            sum = sum + vec1(i-1)
c            write(*,*) vec1(i-1)
           enddo
           
           bary_lambda(1) = 1.D0 - sum            
c           Checking and confirming that the point x
c           is inside the element number nbr
           
           if ((vec1(1).lt.(-zero)).or.(vec1(2).lt.(-zero)).or.
     &       (vec1(3).lt.(-zero)).or.(sum.gt.1.D0+zero) ) then
              
             write(*,*) 'Error in fem local interpolation'
             stop
           endif
          
          value = 0.D0 

          do i = 1, 4
            nm = tetra(i, nbr)
            value = value + bary_lambda(i)*v(nm)
            
          enddo
          loc_interp_fem = value
          return           
         end 
       
c----------------------------------------------------------
        function loc_interp_ifem(gamma, dom, nbr, v, x)
!       This function computes the local interpolation
!       in unbounded domain (ifem) at x
        implicit none
        
        include 'fondements_keltoum.h'
        include 'simplexes_keltoum.h'          
        include 'constants_keltoum.h'
        
        double precision x(3), v(*), loc_interp_ifem
        integer nbr, dom, i, j, nm1, nm
        
        double precision bary_lambda(4), vec1(3), BK(3,3),
     &                   IBK(3, 3), value, res, gamma,
     &                   vec(3), vx(3), sum, rr
        external rr        
        

           nm1 = tetra(1, nbr)	      
	       do j = 2, 4
            nm =  tetra(j, nbr)
            do i = 1, 3 
             bk(i,j-1) = xin(i,nm) -  xin(i, nm1)  
            enddo
           enddo
           
           call inversion(bk, ibk)
	       call fphi(dom, x, vx)                 
           do i = 1, 3
            vec(i) = vx(i) - xin(i, nm1)  
           enddo
           call matv3(ibk, vec, vec1)
           
           sum = 0.D0
          
           do i = 2, 4
            bary_lambda(i) = vec1(i-1)
            sum = sum + vec1(i-1)
           enddo
           bary_lambda(1) = 1.D0 - sum            
c           Checking and confirming that the point x
c           is inside the element number nbr
             
           if ((vec1(1).lt.(-zero)).or.(vec1(2).lt.(-zero)).or.
     &       (vec1(3).lt.(-zero)).or.(sum.gt.1.D0+zero) ) then
             
             write(*,*) 'Error in ifem local interpolation'
             stop
           endif
          value = 0.D0 
          do i = 1, 4
            nm = tetra(i, nbr)
            
            value = value + bary_lambda(i)*v(nm)

          enddo

          res =  rr(dom,  vx)
          value = value*(res**gamma)

          loc_interp_ifem = value
          return           
       end           
        