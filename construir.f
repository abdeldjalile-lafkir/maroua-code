       subroutine construir(gm,nb, xi,tetra)  
       
!      This subroutine gives the nodes in the fem reion 
!      input parameters
!      gm : the number of nodes
!      output parameters:
!      xi : an array contine the coordinates of nodes       
       
       implicit none
       integer i,j, k, gm, nb
       double precision  M1(nb), M2(nb), M3(nb), M4(nb), M5(nb)
       double precision N1(gm), N2(gm), N3(gm), N4(gm)
       double precision tab(nb, 5), tab2(gm,4), xi(3,gm)
       double precision tetra(4,nb)
!      character (72):: fichiers
       
       
     
!       fichier = 'teatredre.dat'
	   open (unit = 24, file = 'tetreades.dat', status = 'old')
	   do j=1,nb	   
	    read(24,*) M1(j), M2(j),M3(j),M4(j)
!	    write(*,*) M1(j), M2(j),M3(j),M4(j),M5(j)
	   enddo	   
	   close(24)
	   open (unit = 20, file = 'nodess.dat', status = 'old')
	   do j=1,gm	   
	    read(20,*) N1(j), N2(j),N3(j)
	   enddo	   
	   close(20)	   
	   do i =1,nb
	    do j =1,5
	     tab(i,j) = 0.D0	     
	    enddo
	   enddo
	   do i =1,gm
	    do j =1,4
	     tab2(i,j) = 0.D0	     
	    enddo
	   enddo	   
	   
	   do i= 1,nb
	     tab(i,1) = M1(i)
	     tab(i,2) = M2(i)
	     tab(i,3) = M3(i)
	     tab(i,4) = M4(i)
c	     tab(i,5) = M5(i)
       enddo
       do i= 1,gm
	     tab2(i,1) = N1(i)
	     tab2(i,2) = N2(i)
	     tab2(i,3) = N3(i)
c	     tab2(i,4) = N4(i)
       enddo
       do i= 1,gm
        do k = 1,3
           xi(k,i)= tab2(i, k)
         enddo
       enddo
       do i = 1, nb
        do k = 1, 4
          tetra(k,i) = tab(i,k)
        enddo
       enddo
!       write(*,*) tetra(1,5)
       end