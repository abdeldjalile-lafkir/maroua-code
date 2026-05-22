                           
!             Etat d'avancement : fini              
              subroutine A_MAT(exmp, iso, drv, xvec, res)
!             This subroutine computes the entries of matrix 
!             A = (a_{i,j}) in the equation (coefficients of the second
!             order terms)
!             exmp is the example index
!             iso = 1 if A = a(x) I for some scalar function "a" 
!                   0 else. 
!             res(i,j) = a_{i, j} for all i, j if drv = 0
!                      = pt_i a_{i,j}          if drv = 1
!             The matrix A
!             iso = 1 if the matrix is of the form a(.) I
!             iso = 0 in general case

             
              implicit none
              double precision xvec(*), sres
              double precision res(3, 3), a_sca
              integer i, iso, exmp, drv
              
              
              call zero_mat(3, 3, res) 
              if (iso.eq.1) then
               if (drv.eq.0) then 
                 sres = a_sca(exmp, 0, xvec)
                 do i = 1, 3
                  res(i, i) = sres
                 enddo                
               else if (drv.eq.1) then
                do i = 1, 3
                 sres = a_sca(exmp, i, xvec)
                 res(i, i) = sres
                enddo
               endif
              else if (iso.eq.0) then 
               select case(exmp)
                case (1)
!                a modifier ici: mettre un vrai exemple
                  res(1, 1) = 1.D0
                case (2)
!                a modifier ici: mettre un vrai exemple                 
                  res(1, 1) = 1.D0
                 case (3)
!                a modifier ici: mettre un vrai exemple                 
                 res(1, 1) = 0.D0
               end select
              endif
             return
            end
!----------------------------------------------------------------------
           function a_sca(exmp, i, xvec)
!          This function contains to coefficient "a" 
!          of the second order term 
!          when equation writes -div(a grad u)+...
           
           implicit none
           double precision  res, a_sca, xvec(*)
           double precision AA, BB, OMG, rr
           integer exmp, i
            
            res = 1.D0 
             
            select case (exmp)
              case (1)
                 if (i.eq.0) then 
                  res = 1.D0
                 else 
                  res = 0.D0
                 endif
              case (2)
                 AA = 1.D0
                 BB = 0.2D0
                 OMG = 20
                 rr = xvec(1)*xvec(1) + xvec(2)*xvec(2)
     &              + xvec(3)*xvec(3)
                 if (i.eq.0) then 
                  res = AA + BB*sin(OMG*rr)
                 else if ((i.ge.1).and.(i.le.3)) then 
                  res = 2*BB*xvec(i)*cos(OMG*rr)
                 endif
              case (3)
!                a modifier ici: mettre un vrai exemple
                 if (i.eq.0) then 
                  res = 1.D0 - xvec(1)*xvec(1) -xvec(2)*xvec(2)
     &              - xvec(3)*xvec(3)
                 else if ((i.ge.1).and.(i.le.3)) then 
                  res = -2.D0*xvec(i)
                 endif
              case (4)
!                a modifier ici: mettre un vrai exemple
                res = 1.D0 
            end select
            
            a_sca = res
            return
          end
!--------------------------------------------------------------
          subroutine b_vect(exmp, xvec, vect)
!          this subroutine gives value of the vector
!          b
          implicit none
          double precision xvec(*), vect(3)

          integer exmp, m
          
                   
            select case(exmp)
             case (1)
!                a modifier ici: mettre un vrai exemple
             do m = 1, 3 
              vect(m) = 0.D0
             enddo         
             case (2)
!                a modifier ici: mettre un vrai exemple                 
             do m = 1, 3 
              vect(m) = 0.D0
             enddo 
             case (3)
!                a modifier ici: mettre un vrai exemple                 
             do m = 1, 3 
              vect(m) = 0.D0
             enddo 
             end select           
            return
            end
c-----------------------------------------------------------------------------------
           function value_c(exmp, xvec)
!          this function gives the value of c(x) 

           implicit none
           double precision value_c, xvec(*)
           double precision res
           integer exmp 
           
           
             res = 0.D0    
           select case(exmp) 
            case(1)                   
!            a modifier ici: mettre un vrai exemple
             res = 0.D0
            case(2)
!            a modifier ici: mettre un vrai exemple
             res = 0.D0
            case(3)
!            a modifier ici: mettre un vrai exemple
             res = 0.D0
           end select
           value_c = res  
           return                      
           end        
c----------------------------------------------------------------------------------
          subroutine ex_function(exmp, iso, x, val_u, val_rhs)
          implicit none 
!         this subroutine computes the value of the rhs  
!         and of the exact solution         
          double precision x(3), val_rhs
          double precision deriv_fun, val_u, MAA(3, 3), vect(3)       
          double precision grad_u(3), grad_AA(3, 3), 
     &        res0, res11, res12, som12,  res2, res3, res4, value_c     
                          
          integer exmp, iso, j, i, i1, j1, k1
          
           call zero_mat(3, 3, grad_AA)     
           SELECT CASE(exmp)
            CASE(1, 2)
   	         val_u = deriv_fun(exmp, 0, 0, 0, x)	         
	         do i = 1, 3 
	          i1 = (i-2)*(i-3)/2
	          j1 = - (i-1)*(i-3)
	          k1 = (i-1)*(i-2)/2 
	           
	          grad_u(i) = deriv_fun(exmp, i1, j1, k1,  x)
	         enddo	       
	         call A_MAT(exmp, iso, 0, x, MAA) 
	         call A_MAT(exmp, iso, 1, x, grad_AA)
	         call b_vect(exmp, x, vect)	
	        	        
	         res11 = 0.D0
	         res12 = 0.D0
	         res2  = 0.D0
	         do j = 1, 3
	          som12 = 0.D0
	          do i = 1, 3
	           i1 = ((i-2)*(i-3) + (j-2)*(j-3))/2
	           j1 = - (i-1)*(i-3) - (j-1)*(j-3)
	           k1 = ((i-1)*(i-2) + (j-1)*(j-2))/2  
	           res0 = deriv_fun(exmp, i1, j1, k1, x)
	           res11 = res11 + MAA(i,j)*res0
	           som12 = som12 + grad_AA(i, j)
	          enddo 
	          
	          res12 = res12 + som12*grad_u(j)   
	          res2 = res2 + vect(j)*grad_u(j)   
	         enddo
           	 res3 = value_c(exmp, x)*val_u
	         res4 = -res11 - res12 + res2 + res3           
          end select
          
            val_rhs = res4
          return
         end
c------------------------------------------------------------------------
         function deriv_fun(exmp, i, j, k, x)

!        this function computes the partial devivatives
!        of function u with respect to the multi-index
!        (i, j, k). 1 <= i <= 3, 1 <=j <=3, 1 <=k <=3

         implicit none 
         double precision deriv_fun, alpha, rho
         double precision x(3), res
         integer drv, exmp, i, j, k
         
         
         alpha = 4.D0
         drv = i + j + k

         res = 0.
         select case (exmp)    
         case (1, 2)
          rho = 1.D0 + x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
	      rho = dsqrt(rho)
          select case (drv)
          case(0)        
            res = x(3)*x(3)/(rho**alpha)    
c             res = 1.D0/(rho**alpha)            
          case(1)     
           if (i.eq.1) then
            res = - alpha*x(3)*x(3)*x(1)/rho**(alpha+2.D0)
           else if (j.eq.1) then 
            res = - alpha*x(3)*x(3)*x(2)/rho**(alpha+2.D0)
           else
            res =  - alpha*x(3)**3.D0/rho**(alpha+2.D0)
     &         + 2.D0*x(3)/rho**alpha  
           endif  
          case(2)
           if (i.eq.2) then          
            res = alpha*(alpha+2.D0)*x(1)*x(1)*x(3)*x(3)/rho**
     &          (alpha+4.D0)- alpha*x(3)*x(3)/rho**(alpha+2.D0)
           else if (j.eq.2) then          
            res = alpha*(alpha+2.D0)*x(2)*x(2)*x(3)*x(3)/rho**
     &         (alpha+4.D0)- alpha*x(3)*x(3)/rho**(alpha+2.D0)        
           else if (k.eq.2) then        
            res = alpha*(alpha+2.D0)*x(3)**4.D0/rho**(alpha+4.D0)
     &         - 5.D0*alpha*x(3)*x(3)/rho**(alpha+2.D0)
     &         + 2.D0/rho**alpha    
           else if ((i.eq.1).and.(j.eq.1)) then  
            res = alpha*(alpha+2.D0)*x(1)*x(2)*x(3)*x(3)/rho**
     &         (alpha+4.D0)  
           else if ((i.eq.1).and.(k.eq.1)) then          
            res =alpha*(alpha+2.D0)*x(1)*x(3)**3/rho**
     &         (alpha+4.D0)- alpha*2.D0*x(1)*x(3)/rho**(alpha+2.D0)   
           else if ((j.eq.1).and.(k.eq.1)) then  
            res = alpha*(alpha+2.D0)*x(2)*x(3)**3/rho**
     &         (alpha+4.D0)- alpha*2.D0*x(2)*x(3)/rho**(alpha+2.D0)                              
           endif   
          end select
!         case (1, 2)
!          rho = 1.D0 + x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
!	      rho = dsqrt(rho)
!          select case (drv)
!          case(0)        
!            res = 0.D0    
!c             res = 1.D0/(rho**alpha)            
!          case(1)     
!           if (i.eq.1) then
!            res = 0.D0
!           else if (j.eq.1) then 
!            res = 0.D0
!           else
!            res = 0.D0  
!           endif  
!          case(2)
!           if (i.eq.2) then          
!            res =0.D0
!           else if (j.eq.2) then          
!            res = 0.D0        
!           else if (k.eq.2) then        
!            res =0.D0   
!           else if ((i.eq.1).and.(j.eq.1)) then  
!            res =0.D0  
!           else if ((i.eq.1).and.(k.eq.1)) then          
!            res =0.D0  
!           else if ((j.eq.1).and.(k.eq.1)) then  
!            res = 0.D0                            
!           endif   
!          end select
!
        end select
           deriv_fun = res
         return
        end                   