c----------------------------------------------------------
c  Etat de modification: fini (aucun changement)

         function  det(a)
c        computes the determinant of 3*3 matrix

	     implicit none
	     double precision a(3,3), det, res
	  
	     res =   (a(2,2)*a(3,3) - a(2,3)*a(3,2))*a(1,1)
     &         -   (a(1,2)*a(3,3) - a(3,2)*a(1,3))*a(2,1) 
     &         +   (a(1,2)*a(2,3) - a(2,2)*a(1,3))*a(3,1) 
         det = res	  

		 return 
	    end
c----------------------------------------------------------
c  Etat de modification: fini (aucun changement)

          subroutine prodsc(vec1, vec2, resu)
c         computes the inner product of 
c         two vector of size 3

 
		   implicit none
	       double precision vec1(*), vec2(*), resu, res
	       integer i

            res = 0.D0
	       do i = 1, 3
               res = res + vec1(i)*vec2(i) 
	       enddo
             resu = res
	       return
		 end   
c----------------------------------------------------------
c  Etat de modification: fini (aucun changement)

          subroutine prodcolumns(k, j,  mat1, mat2, resu)
c         computes the inner product of 
c         the kth and jth colomns  of matrices mat1, mat2 res

 
		   implicit none
	       double precision mat1(3, *), mat2(3, *), resu, res
	       integer i, k, j

           res = 0.D0
	       do i = 1, 3
               res = res + mat1(i,k)*mat2(i, j) 
	       enddo
            resu = res

	       return
		 end   		   
c -----------------------------------------------------
c     Etat de modification: fini (aucun changement)

           subroutine copy_mat(n, p, mat1, mat2)

c            This subroutine copy a rectangular matrix mat1 into 
c            asecond matrix of the same size mat2

		   implicit none
	       integer n, p
	       double precision mat1(n, p), mat2(n, p)

           integer i, j
	       do  i = 1, n
	        do j = 1, p
              mat2(i, j) = mat1(i, j)
            enddo
           enddo 
	       return
          end
c -----------------------------------------------------
c     Etat de modification: fini (aucun changement)

           subroutine copy_vect(n, vect1, vect2)

c            This subroutine copy a vector vect1
c            into a vector vect2

		   implicit none
	       integer n
	       double precision vect1(*), vect2(*)

           integer i
           
	       do  i = 1, n
              vect2(i) = vect1(i)
           enddo 
	       return
          end
c-------------------------------------------------------------------
c        Etat d'avancement: fini
         function finteg_rho(m, x)
	     implicit none
c        cette subroutine calcule les integrales de la forme
c        int_0^x rho^alpha avec rho = sqrt(1 + r*r)     
         integer          m
	     double precision res, finteg_rho, pi, x, rho
	     parameter (pi = 3.14159265358D0)
		  		   
	     rho = 1.D0 + x*x 
       	 rho = dsqrt(rho)
		 select case(m)
		 case(2)
	        res = atan(x) 
         case(4)
              res = 0.5*x*(rho**(-2.D0)) + 0.5*atan(x) 
		 case(6) 
		    res = 0.25*x*(rho**(-4.D0)) + 0.75D0*
     &	         (0.5*x*(rho**(-2.D0)) + 0.5*atan(x))	 
		 case(8)
            res = x*(rho**(-6.D0))/6.D0 + 5.*
     &       (0.25*x*(rho**(-4.D0)) + 0.75D0*
     &	   (0.5*x*(rho**(-2.D0)) + 0.5*atan(x))	)/6.D0  
		 end select
	     finteg_rho = res
		 return
	    end 
		   
c----------------------------------------------------------------	    
c      Etat d'avancement : fini
       function	test(i, j, x, y)
	   double precision i, j, x, y, test
		 
		 test = x**i*y**j
	  return
	  end	   
		
		   
c----------------------------------------------------------------
c  Etat de modification: fini (aucun changement)

          subroutine zero_mat(n, p, mat)

c         This subroutine makes a rectangular matrix 
c         equal zero 

		  implicit none
	      integer n, p
	      double precision mat(n, p)

          integer i, j
	      do  j = 1, p
	       do i = 1, n
             mat(i, j) = 0.D0
           enddo
          enddo 
	      return
          end
c----------------------------------------------------------------
c  Etat de modification: fini (aucun changement)

          subroutine zero_vect(n, vec)

c         This subroutine makes a vector 
c         equal zero 

		  implicit none
	      integer n
	      double precision vec(n)

          integer i
	      
	       do i = 1, n
             vec(i) = 0.D0
           enddo
          
	      return
          end




c------------------------------------------------------------------------
c  Etat de modification: fini (aucun changement)

           subroutine  produitAAT(ck, bbk)	
	   
c           Computes the matrix product C*C^T
c           where C is a 3*3 matrix
c           The result is put in BBK

    	   implicit none
	       double precision ck(3, 3), bbk(3,3), res
	       integer i, j, k

             do i = 1, 3
	          do j = 1, 3
                res = 0.D0
                do k = 1, 3
                 res = res + ck(i, k)*ck(j, k)
                enddo
	           bbk(i,j) = res
	          enddo
		     enddo
 	       return
           end


c----------------------------------------------------------------------
c  Etat de modification: fini (aucun changement)

 	         subroutine matv3(ak, vec1, vec2)  

c            This subroutine computes the product 
c            matrix vector AK*VEC1, where AK is a 3*3 matrix
c            while VEC1 is a column vector. The result is VEC2

		    implicit none
              
	        double precision ak(3, 3), vec1(3), vec2(3), res
	        integer i, j
             do i = 1, 3
	          res = 0.D0
               do j = 1, 3
		        res = res + ak(i, j)*vec1(j)
			   enddo              
	          vec2(i) = res
             enddo
		   end 
c     ------------------------------------------------------
c  Etat de modification: fini (aucun changement)

 	         subroutine mat_v3(k, ak, mat, vec2) 

c            This subroutine computes the product 
c            matrix vector AK*CL, where AK is a 3*3 matrix
c            while CL is the kth column vector of a 3*3 matrix mat. 
c            The result is VEC2

 
 	         implicit none             
	         double precision ak(3, 3), mat(3, *), vec2(3), res
	         integer i, j, k
		   
             do i = 1, 3
	          res = 0.D0
               do j = 1, 3
			    res =  res + ak(i, j)*mat(j, k)
	           enddo              
	          vec2(i) = res
             enddo
	         return
		     end 
c     -------------------------------------------------------
c       Etat de modification: fini (aucun changement)

 	         subroutine nmat_v3(k, ak, mat1, mat2)  
c            This subroutine computes the product 
c            matrix vector AK*CL, where AK is a 3*3 matrix
c            while CL is the kth column vector of a 3*3 matrix mat1. 
c            The result is put in the kth  column vector  of mat2


 	         implicit none             
	         double precision ak(3, 3), mat1(3, *), mat2(3, *), res
		   
	         integer i, j, k
		   
             do i = 1, 3
	          res = 0.D0
               do j = 1, 3
			    res =  res + ak(i, j)*mat1(j, k)
	           enddo              
	          mat2(i, k) = res
             enddo
	         return
		     end 
c------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)

            subroutine matv3_bis(k, ak, vec, mat)

c            This subroutine computes the product 
c            matrix vector AK*VEC, where AK is a 3*3 matrix
c            while VEC is a column vector 
c            The result is put in the kth  column vector  of MAT

             implicit none
             integer k
	         double precision ak(3, 3), vec(3), mat(3, *)

	         integer i, j
             double precision res

	          do i = 1, 3
               res = 0.D0
	           do j = 1, 3
                 res = res + ak(i, j)*vec(j)
               enddo 
               mat(i, k) = res
	          enddo
	         return
             end    
c--------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)


             subroutine  matmat(n, p, m, a, b, c)

c            Computes the product c = a*b of
c            two rectangular matrices a and b

 
	         implicit none
	         integer n, p, m
	         double precision a(n, p), b(p, m), c(n, m)

	         integer  i, j, k
	         double precision s
         
	         do i  = 1, n
              do j = 1, m
               s = 0.D0
	            do k = 1, p
                s  =  s +  a(i, k)*b(k, j)
	            enddo	   
	           c(i, j) = s
	          enddo
             enddo

	         end
c--------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)

             subroutine  som_mat(n, p, alpha, a, beta, b, c)

c            Computes the linear combination 
c            C = alpha*A + beta*B
c            of two matrices A and B

	         implicit none
	         integer n, p
	         double precision a(n, *), b(n, *), c(n, *), alpha, beta 
  
	         integer  i, j

         
	         do i  = 1, n
              do j = 1, p
	           c(i, j) =  alpha*a(i, j) + beta*b(i, j)
	          enddo
             enddo
             return
	         end
c--------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)

             subroutine transpos(n, p, A, AT)

c            computes the transpose of a rectangular matrix A

	         implicit none
	         integer n, p
	         double precision a(n, p), at(p, n)

             integer i, j

	         do i = 1, n
	          do j = 1, p
		       at(i, j) = a(j, i)
		      enddo
	         enddo	  

	        return
            end


c----------------------------------------------------------------------
c           Etat de modification: fini (aucun changement)

            subroutine som_vec(n, alpha, u, beta, v, w)

c           Computes the linear combination 
c           w = alpha*u + beta*v

	        implicit none
	        integer n
	        double precision alpha, beta, u(*), v(*), w(n)         
	        integer i

	        do i = 1, n
	         w(i) = alpha*u(i) + beta*v(i)
            enddo
            return
	        end
c-------------------------------------------------------------------------
c           Etat de modification: fini (aucun changement)

            subroutine somvec(n, alpha, u, beta, v, w)

c           Computes the linear combination 
c           w = alpha*u + beta*v

	        implicit none
	        integer n
	        double precision alpha, beta, u(*), v(*), w(n, 1)         
	        integer i

	        do i = 1, n
	         w(i, 1) = alpha*u(i) + beta*v(i)
            enddo
            return
	        end
c--------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)

             subroutine inversion(a, ia)
c            computes the inverse of a 3*3 matrix a

	         implicit none
             double precision a(3, 3), ia(3, 3)
             double precision deta, det
	         external det
	         integer i, j
	     
              deta = det(a)



             ia(1, 1) = a(2, 2)*a(3, 3) - a(3, 2)*a(2, 3)
	         ia(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
		     ia(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
			
		     ia(1, 2) = a(3, 2)*a(1, 3) - a(3, 3)*a(1, 2)
		     ia(2, 1) = a(3, 1)*a(2, 3) - a(2, 1)*a(3, 3)

	         ia(1, 3) = a(1, 2)*a(2, 3) - a(2, 2)*a(1, 3)
	         ia(3, 1) = a(2, 1)*a(3, 2) - a(3, 1)*a(2, 2)
			
		     ia(3, 2) = a(3, 1)*a(1, 2) - a(1, 1)*a(3, 2)
		     ia(2, 3) = a(2, 1)*a(1, 3) - a(1, 1)*a(2, 3)
			
			
	         do j = 1, 3
   	          do i = 1, 3
	           ia(i, j) = ia(i, j)/deta
		      enddo
	         enddo 		 
             return
	         end	 
c-------------------------------------------------------------------------
c       Etat de modification: fini (aucun changement)

             subroutine inversion_det(a, ia, deta)

c            Computes the inverse and the determinant of a matrix A

	          implicit none
              double precision a(3, 3), ia(3, 3)
              double precision deta, det
	          external det
	          integer i, j
	     
               deta = det(a)



              ia(1, 1) = a(2, 2)*a(3, 3) - a(3, 2)*a(2, 3)
	          ia(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
		      ia(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
			
		      ia(1, 2) = a(3, 2)*a(1, 3) - a(3, 3)*a(1, 2)
		      ia(2, 1) = a(3, 1)*a(2, 3) - a(2, 1)*a(3, 3)

	          ia(1, 3) = a(1, 2)*a(2, 3) - a(2, 2)*a(1, 3)
	          ia(3, 1) = a(2, 1)*a(3, 2) - a(3, 1)*a(2, 2)
			
		      ia(3, 2) = a(3, 1)*a(1, 2) - a(1, 1)*a(3, 2)
		      ia(2, 3) = a(2, 1)*a(1, 3) - a(1, 1)*a(2, 3)
			
			
	          do j = 1, 3
		       do i = 1, 3
	            ia(i, j) = ia(i, j)/deta
		       enddo
	          enddo 		 
              return
	          end	 
c----------------------------------------------------------------
 	         subroutine prodsc_v3(k, mat, vec, resu) 

c            This subroutine computes the inner product of 
c            the kth column vector of the matrix mat and 
c            the vector vec 

 	         implicit none             
	         double precision  mat(3, *), vec(3), res, resu
	         integer  j, k
		   

	          res = 0.D0
               do j = 1, 3
			    res =  res + mat(j, k)*vec(j)
	           enddo 
	          resu = res              
	         return
		     end 
c--------------------------------------------------------------
 	         subroutine prod_cl(k, j, mat1, mat2, resu) 

c            This subroutine computes the inner product of 
c            the kth column vector of the matrix mat1 and 
c            the jth column vector of the matrix mate 

 	         implicit none             
	         double precision  mat1(3, *), mat2(3, *), res, resu
	         integer  j, k, i
		   

	          res = 0.D0
               do i = 1, 3
			    res =  res + mat1(i, k)*mat2(i, j)
	           enddo 
	          resu = res              
	         return
		     end 

