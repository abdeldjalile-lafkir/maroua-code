       subroutine construir(gm, nt, xi, tetra)
      implicit none

C**** This subroutine gives the nodes in the fem region
C**** input parameters
C**** gm : the number of nodes
C**** output parameters:
C**** xi : an array contine the coordinates of nodes

      integer i, j, k, gm, nt, iostat
      double precision M1(nt), M2(nt), M3(nt), M4(nt)
      double precision N1(gm), N2(gm), N3(gm)
      double precision tab(nt,4), tab2(gm,3)
      double precision xi(3,gm), tetra(4,nt)

C**** Read tetrahedron data from file
      open(unit=24, file='tetra.dat', status='old',
     &     iostat=iostat)
      do j = 1, nt
         read(24, *, iostat=iostat) M1(j), M2(j), M3(j), M4(j)
         if (iostat .ne. 0) exit
      enddo
      close(24)

C**** Read coordinates data from file
      open(unit=20, file='coordinates.dat', status='old',
     &     iostat=iostat)
      do j = 1, gm
         read(20, *, iostat=iostat) N1(j), N2(j), N3(j)
         if (iostat .ne. 0) exit
      enddo
      close(20)

C**** Initialize arrays to zero
      do i = 1, nt
         do j = 1, 4
            tab(i, j) = 0.0d0
         enddo
      enddo

      do i = 1, gm
         do j = 1, 3
            tab2(i, j) = 0.0d0
         enddo
      enddo

C**** Fill arrays with values
      do i = 1, nt
         tab(i, 1) = M1(i)
         tab(i, 2) = M2(i)
         tab(i, 3) = M3(i)
         tab(i, 4) = M4(i)
      enddo

      do i = 1, gm
         tab2(i, 1) = N1(i)
         tab2(i, 2) = N2(i)
         tab2(i, 3) = N3(i)
      enddo

C**** Extract xi
      do i = 1, gm
         do k = 1, 3
            xi(k, i) = tab2(i, k)
  !          write(*,*) 'xi(', k, ',', i, ') = ', xi(k, i)
         enddo
      enddo

C**** Extract tetra
      do i = 1, nt
         do k = 1, 4
            tetra(k, i) = tab(i, k)
    !        write(*,*) 'tetra(', k, ',', i, ') = ', tetra(k, i)
         enddo
      enddo

C**** Write the results to files
      write(*,*) 'The nodes in the fem region have been constructed'
      open(unit=30, file='xi.dat', status='replace',
     &     iostat=iostat)
      do i = 1, gm
         write(30, *, iostat=iostat) xi(1, i), xi(2, i), xi(3, i)
         if (iostat .ne. 0) exit
      enddo
      close(30)
      open(unit=40, file='fortran-tetra.dat', status='replace',
     &     iostat=iostat)
      do i = 1, nt
         write(40, *, iostat=iostat) (tetra(k, i), k=1,4)
         if (iostat .ne. 0) exit
      enddo
      close(40)


      end
