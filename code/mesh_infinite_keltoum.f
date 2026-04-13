!-------------------------------------------------------------------------
! Etat de modification: fini

      subroutine mesh(N, pform, mu, xin1, tetra1, domnod, domtet, gm,
     &                nt,nbtetra,nbtetra1,nbnodes,nbnodes1)


!       This  subroutine generates the mesh of the whole space after
!       sharing it into the union of one big tetrahedra K (the FEM region) and
!       4 infinite  simplices T1, T2, T3 and T4.
!       The finite element region K  is obtained from a fixed reference
!       tetrahedra RK using a homothety centered at the origin and
!       of ratio "pform". Thus, the parameter "pform" characterizes the size
!       of the finite element region. It characterizes also the distance of the
!       frontier between the FEM region and the IFEM region to the origin.
!
!       INPUT PARAMETERS
!        N           : This is an integer. This parameter N and the parameter
!                      mu suffice to generate the "mesh" in both the FEM region and  the virtual
!                      domain (obtained after inverting the infinite simplices)
!                      The greater is N the greater is the number of
!                      tetrahedra
!                      More precisely, N denotes the number of transverse levels i
!                      in the  mesh of one big tetrahedra
!        mu         : This is  the graduation parameter which characterizes the
!                       graduation of the mesh near to origin.
!                       mu must be between 0 and 1.
!                       It could be better to choose  mu greater than 0.5
!                       when mu = 1. the mesh is the virtual region is not graded
!        pform       :  is a real characterizing the  size
!                      of the finite element region. It characterizes also the distance of the
!                      bounbdary between the FEM regions and the IFEM region to the origin.
!                      The user can chose pform = 1 or 2.5 for example.
!      OUTPUT PARAMETERS
!        xin(3, *)   : an array containing the coordinates of the nodes
!        tetra(4, *) : an array containing the vertices numbers for all the tetrahedra
!                      tetra(i, k) returns the global vertex number of the i-th vertex
!                      in the k-th tetrahedra
!        domnod      : an array. for each vertex i, domnod(i) denotes the  number of the domain
!                      to which the vertex belongs. Domnod(i) is an integer
!                      between 1 and 5 (five domains)
!        domtet      : an array. for each tetrahera i, domtet(i) denotes the number of the domain
!                      in which the tetrahera belongs. Domtet(i) is an integer
!                      between 1 and 5 (five domains)
!       THE SAME IN FRENCH
!       Cette subroutine va generer le maillage du demi-espace apres
!       sa decomposition en quatre simplexes infinis et un simplexe fini.
!       Parametres d'entree
!        N           : Un entier qui indique le  nombre des niveaux
!                      du maillage dans chaque simplexe obtenu par inversion.
!        mu          : un parametre reel pour graduer le maillage pr s de l'origine.
!        pform       : un parametre de forme du simplexe element fini
!                      on peut la valeur prendre 1.
!       Parametres de sortie
!        tetra(4, *) : un tableau  contenant les numeros des  4 sommets de chacun
!                      des  tetraedres.
!        domnod      : un tableau qui indique pour chaque noeud le numero
!                      d'un des sous-domaines (grand tetra) auquels il appartient.
!        domtet      : un tableau qui indique pour chaque tetraedre le   numero du sous-domaine (grand tetra) auquel il appartient.
!-----
       implicit none
       integer N,nbnodes1,nbnodes,nbtetra,nbtetra1
	   double precision mu, pform

      integer nrep, i,j,m,k, nbnodesTet, dom, gm, nt,kk,nbprisme
      double precision xin(3,nbnodes), tab(4,nt), xx(3,gm)
	  double precision xin1(3,nbnodes1)
      integer tetra(4, nbtetra), domnod(nbnodes1),  domtet(nbtetra1),
     &    tetra1(4, nbtetra1)

        double precision  xi(3, (5*(N+1)*(N+2)*(N+3))/6)
        integer catg((5*(N+1)*(N+2)*(N+3))/6),
     &   num_catg(10*(N+1)*(N+2)),
     &   rnum((5*(N+1)*(N+2)*(N+3))/6)


	     nbnodesTet = ((N+1)*(N+2)*(N+3))/6
	     nbprisme = ((N-1)*N)/2 + (N*(N-1)*(2*N-1))/6
c	   gm = 284
c	   nt = 996
          
         call first_enumeration_whsp(N, pform, mu, xi)
         call categories_nodes_whsp(N, catg, num_catg, nrep)
         call new_number_whsp(N,catg,nrep, num_catg, rnum, xi)
		
	     do dom = 1, 4
          do k = 1, nbnodesTet
            i = (dom-1)*nbnodesTet + k
 	        j = rnum(i)
            do m = 1, 3
              xin1(m, j) = xi(m, i)
	        enddo
	        domnod(j) = dom
c			write(*,*)  domnod(j)
          enddo
	     enddo
		  
	     call construir(gm,nt,xx,tab)
	     do k = 1, gm
	        i = nbnodes + k
            do m = 1, 3
              xin1(m, i) = xx(m, k)
	        enddo
	        domnod(i) = 5
          enddo
c	        write(*,*) 'k=', domnod(4*nbnodesTet+1)

         call make_tetra_whsp(N, tetra, rnum,  domtet)
         do dom = 1, 4
          do k = 1, 3*nbprisme + N
           i = (dom-1)*(3*nbprisme + N) + k
! 	        j = rnum(j)
           do m = 1, 4! 3
             tetra1(m, i) = tetra(m, i)
           enddo
           domtet(i) = dom
          enddo
         enddo
         do i = 1, nt
            j = 4*(3*nbprisme + N) + i
            do k = 1, 4
              tetra1(k,j) = nbnodes + int(tab(k, i))
            enddo
            domtet(j) = 5
         enddo
c		 do kk= 4*nbnodesTet+1, 4*nbnodesTet+gm
c			 write(*,*) domnod(kk)
c			 if(domnod(kk).eq.5)then
c				 write(*,*)'hello'
c				 write(*,*) xin(1, kk),xin(2, kk),xin(3, kk)
c				 write(*,*) domnod(kk)
c			 endif
c		enddo
		   write(*,*) 'kk=', domnod(nbnodes+1)

   	   return
       end
!-------------------------------------------------------------------------
c Etat de modification: fini
c
        subroutine first_enumeration_whsp(N, pform, mu, xi)
c       This subroutine gives a first list of the vertices
c       of the mesh.
c       Common vertices can appear several times
c          N (input parameter): the number of levels in subdividing
c                              the big domains. The final number of tetrahedra
c                              and nodes depend only on N
c         pform : the size of the fem region (a big tetrahedra)
c         mu : the graduation parameter of the mesh
c         xi : coordinates of the vertices

         implicit none
         integer N
	     double precision mu

         integer deca, nbnodesTet, nd
	     double precision xi(3,*), pform

	     nbnodesTet =  (N+1)*(N+2)*(N+3)/6
         call big_simp_vertices_whsp(pform)
         do nd = 1, 4
  	      deca = (nd-1)*nbnodesTet
  	      call mesh_simp(deca, N, mu, nd, xi)
  	     enddo
         deca = 4*nbnodesTet
	     call mesh_simp(deca, N, 1.D0, 5, xi)

         return
        end
c-------------------------------------------------------------------------
c  Etat de modification : fini

         subroutine categories_nodes_whsp(N, catg, num_catg,
     &                                     nrep)

c         This subroutines identify multiple nodes
c         which are common to more than one domain
c         These nodes have multiple local numbers


c         INPUT  :
c              N       :  an integer, the number of levels in the mesh
c
c         OUTPUT
c              catg    : an array. Its dimension is equal to the total number
c                        of nodes (with multiple numbering for some of them).
c                        If a node i is multiple, then catg(i) = 1, else catg(i) = 0.
c              nrep    : the number of multiple nodes
c              num_catg: an array containing the numbers of multiple nodes
c                      (in a global numbering in which some nodes could
c                          have several numbers).
         implicit none
	     integer N, nbnodesTet, nbnodes,  nb1nodes, nbmult
	     integer catg(*), num_catg(*), nrep, num
	     external num

	     integer i, j, k, l, it, s, deca, nb


	     nbnodesTet = (N+1)*(N+2)*(N+3)/6
	     nbmult     = 5*(N*N + N + 1)
	     nb1nodes   = 5*nbnodesTet
 	     nbnodes    = nb1nodes - nbmult

	     do i = 1, nb1nodes
           catg(i) = 0
         enddo

         it = 0
         do s  = 1, 4
	     deca = (s-1)*nbnodesTet
	     nb = deca + 1
	     catg(nb) =  1
	     it = it + 1
	     num_catg(it) = nb
	     do i = 1, N-1
           k = 0
	       it = it + 1
	       nb = deca + num(i, k, 1)
		   catg(nb) =  1
           num_catg(it) = nb
           do k = 1, i-1
  	        do j = 1, 2
	         l = k*(j-1) + 1
	         it = it + 1
		     nb = deca + num(i, k, l)
		     catg(nb) =  1
             num_catg(it) = nb
            enddo
           enddo
   	       k  = i
  	       do l = 1, k+1
	        it = it + 1
	        nb = deca + num(i, k, l)
	        catg(nb) =   1
            num_catg(it) = nb
           enddo
	      enddo

	     i =  N
          do k = 0, i
  	       do l = 1, k+1
	        it = it + 1
		    nb = deca + num(i, k, l)
	 	    catg(nb) =  1
            num_catg(it) = nb
           enddo
	      enddo
	     enddo

        s  = 5
	    deca = (s-1)*nbnodesTet
	    nb = deca + 1
	    catg(nb) =   1
	    it = it + 1
	    num_catg(it) = nb
	    do i = 1, N-1
          k = 0
	      it = it + 1
	      nb = deca + num(i, k, 1)
	      catg(nb) =  1
          num_catg(it) = nb

          do k = 1, i-1
  	       do j = 1, 2
	         l = k*(j-1) + 1
	         it = it + 1
		     nb = deca + num(i, k, l)
  	         catg(nb) =  1
             num_catg(it) = nb
           enddo
          enddo

   	      k = i
          do l = 1, k+1
	       it = it + 1
	       nb = deca + num(i, k, l)
     	   catg(nb) =  1
           num_catg(it) = nb
          enddo
	     enddo
          i =  N
          do k = 0, i
  	       do l = 1, k+1
	        it = it + 1
		    nb = deca + num(i, k, l)
	 	    catg(nb) =  1
            num_catg(it) = nb
           enddo
	      enddo


          nrep = it

	     return
	   end
c--------------------------------------------------------------------
c      Etat d'avancement: fini

       subroutine new_number_whsp(N, catg, nrep,num_catg,rnum, xi)

c        Cette subroutine part de l'ancienne numerotation ou certains
c        noeuds sont multinumerotes pour construire un nouveau tableau
c        de numeros ou chaque noeud a un et un seul numero.
c        INPUT  :
c              nb1nodes:  the total number of nodes (with multiple numbering)
c              catg    : an array. Its dimension is equal to the total number
c                        of nodes (with multiple numbering for some of them).
c                        If a node i is multiple, then catg(i) = 1, else catg(i) = 0.
c              nrep    : the number of multiple nodes
c              num_catg: an array containing the numbers of multiple nodes
c                      (in a global numbering in which some nodes could have several numbers).
c              xi      : coordinates of nodes
c        OUTPUT
c              rnum    : an array of dimension "nb1nodes". It associates the
c                        old numbering of nodes (with multiplicity of nodes) to the new
c                        one (with multiplicity of nodes).
c              catg    : this array is modified here.
c                        If a node i is dropped (since multiple), then catg(i) = 2, else catg(i) = 0.
c
        implicit none
	    integer nb1nodes, catg(*), num_catg(*), rnum(*), N
		double precision dist, zero, xi(3, *)
	    parameter (zero = 1.D-9)
	    integer nb1, nb2, nrep,  it, i, j, nbn


        nb1nodes = 5*(N+1)*(N+2)*(N+3)/6

c       For the moment, the old and the new numerotations are the same
	     do i = 1, nb1nodes
	      rnum(i) = i
	     enddo
	     do it = 1, nrep
	      nb1 = num_catg(it)
		  if(catg(nb1).eq.1) then
  	       do j =it+1, nrep
	        nb2 = num_catg(j)
            dist = (xi(1,nb1) - xi(1,nb2))*(xi(1,nb1) - xi(1,nb2))
     &           + (xi(2,nb1) - xi(2,nb2))*(xi(2,nb1) - xi(2,nb2))
     &           + (xi(3,nb1) - xi(3,nb2))*(xi(3,nb1) - xi(3,nb2))
	        dist = dsqrt(dist)
	        if ((dist.lt.zero).and.catg(nb2).eq.1) then
             catg(nb2) = 2
             rnum(nb2) = nb1
  	        endif
	       enddo
            catg(nb1) = 0
          endif
		 enddo


         nbn = 0
         do i = 1, nb1nodes
		  if (catg(i).eq.0) then
	       nbn  = nbn + 1
           rnum(i) = nbn
		  endif
	     enddo

         do i = 1, nb1nodes
		  if (catg(i).eq.2) then
	       nbn = rnum(i)
           rnum(i) = rnum(nbn)
		  endif
	     enddo

	    end
c-------------------------------------------------------------------
c      Etat de modification :fini (aucun changement)

        subroutine mesh_simp(deca, N,  mu, kk, xi)

c          This subroutines meshes a big simplex whose
c          vertices are given
c          Input parameters
c             N       : an integer, the numbers of levels in the mesh
c             deca    : The starting shift in  the enumeration of the nodes
c             mu      : the mesh parameter
c             kk      : the index of the big vertex.
c          Output parameters
c             xi     : the coordinates of the nodes of the mesh.
c             The final  number of nodes is (N+1)*(N+2)*(N+3)/6



	     implicit none
         include 'simplexes_keltoum.h'
	     integer N, deca,  nbnodesTet, kk, j, m, nb, i
	     double precision mu, bt(3, 3), res

         double precision xi(3,*), xi_ref(3, (N+1)*(N+2)*(N+3)/6)

		 do  j = 1, 3
		  do m = 1, 3
            bt(m,j)=vertex(m,j+1,kk)-vertex(m, 1,kk)
		  enddo
		 enddo

	     nbnodesTet = (N+1)*(N+2)*(N+3)/6

         call mesh_ref_simp(N,   mu,  xi_ref)

	     do nb =  1, nbnodesTet
	      do i = 1, 3
   		   res = 0.D0
	       do j = 1, 3
			 res = res + bt(i,j)*xi_ref(j, nb)
		   enddo
		   xi(i, nb + deca) = res + vertex(i,1, kk)
		  enddo
	     enddo


	    return
	   end
c-------------------------------------------------------------------
c      Etat de modification :fini (aucun changement)

        subroutine mesh_ref_simp(N,    mu,  xi_ref)
c
c          This subroutines gives a graded meshes of  a reference simplex whose
c          vertices are (0, 0, 0), (1,0, 0), (0, 1, 0), (0, 0, 1)
c          Input parameters
c             N          : an integer, the number of levels in the mesh
c             nbnodesTet : the number of nodes, namely
c                          nbnodesTet = (N+1)*(N+2)*(N+3)/6
c             mu         : the mesh parameter

c          Output parameters
c             nxi_ref    : the coordinates of the nodes of the mesh.


	     implicit none
	     integer N
	     double precision mu, xi_ref(3, *)


c          Local variables
         integer i, j, k, l, m,  nb, num
         double precision xsom(3,3), xlim(3,2), di(N), htot, cf
         external radiale, num



		 htot = 1.D0/dsqrt(1.D0*3)

      	 call radiale (N, mu, htot, di)

		 do m = 1, 3
		   xi_ref(m, 1) = 0.D0
		 enddo

         do j = 1, 3
		  do m = 1, 3
            xsom(m, j) = 0.D0
	      enddo
	     enddo

	    do i = 1, N
          do j = 1, 3
		   xsom(j, j) = di(i)/htot
	      enddo

          k = 0
          nb = num(i, k, 1)
	      do m = 1, 3
           xi_ref(m, nb) = xsom(m, 1)
          enddo

          do k = 1, i
	       cf =  1.D0 - k*1.D0/i
	       do m = 1, 3
            xlim(m,1) = cf*xsom(m,1) + (1-cf)*xsom(m,2)
            xlim(m,2) = cf*xsom(m,1) + (1-cf)*xsom(m,3)
	       enddo

  	       do l = 1, k+1
	        cf =  1.D0 - 1.D0*(l - 1.D0)/k
		    nb = num(i, k, l)
	     	 do m = 1, 3
  		      xi_ref(m,nb) = cf*xlim(m,1) + (1-cf)*xlim(m,2)
	         enddo
	       enddo
		  enddo
	     enddo
	     return
	    end
c-------------------------------------------------------------------
c      Etat de modification :fini (aucun changement)

        function num(i, k, l)
c        This function returns  the local number of
c        the node (i, j, k) in the graded  mesh of the reference
c        element.
c        T. Boulmezaoud June 2003
        implicit none
        integer i, k, l, num, pow_sum
	    external pow_sum

	    num = i + (3*pow_sum(1,i-1)+pow_sum(2,i-1))/2+pow_sum(1,k)+l
	   return
       end
c--------------------------------------------------------------------
c    Etat de modification: fini
c
          subroutine make_tetra_whsp(N, tetra, rnum, domtet)
c         This subroutine constructs all the tetrahedra of the mesh
c          N (input parameter): the number of levels in subdividing
c                              the big domains. The final number of tetrahedra
c                              and nodes depend only on N
c          tetra   : the 2D array containing the definition of tetrahedra
c                    tetra(k, i) denotes the number of the k-vertice (k<=4)
c                     of the i-th tetrahedra.
c          rnum    : (input parameter) this array makes the correspondance
c                    between the old numbers of nodes (domain by domain)
c                    in which some nodes have multiple numbers
c                    and the new one (in which each node has one and only one number)
c          domtet : this table contains the number of the big domain to which
c                    each tetrahedra belong  (output)
c          prisme  : input array) containing the numbers of vertices of the prism



           implicit none
	       integer N, tetra(4, *), rnum(*), domtet(*)

           integer deca,  nbnodesTet,  nbnodes
		   integer prisme(6), num
           integer i, j, k, s, l, nb, nbtet
           external num

	       nbnodesTet = (N+1)*(N+2)*(N+3)/6
	       nbnodes =  5*nbnodesTet - 5*(N*N+N+1)


	       nbtet = 0
		   do s = 1, 4
		    deca = (s-1)*nbnodesTet
		    nbtet = nbtet + 1
            do j = 1, 4
	         nb = deca + j
	         l = rnum(nb)
             tetra(j, nbtet) = l
	        enddo
 	        domtet(nbtet) = s
	        do i = 1, N-1
		     do k = 0, i-1
	          do l = 1, k
                  prisme(1) = deca + num(i, k, l)
	              prisme(2) = deca + num(i, k+1, l)
                  prisme(3) = deca + num(i, k+1, l+1)
	              prisme(4) = deca + num(i+1, k, l)
                  prisme(5) = deca + num(i+1, k+1, l)
	              prisme(6) = deca + num(i+1, k+1, l+1)
                  call divide_prisme(s, nbtet, tetra,
     &	    	                   rnum, domtet,  prisme)

                  prisme(3) = deca + num(i, k, l)
	              prisme(2) = deca + num(i, k, l+1)
                  prisme(1) = deca + num(i, k+1, l+1)
	              prisme(6) = deca + num(i+1, k, l)
                  prisme(5) = deca + num(i+1, k, l+1)
	              prisme(4) = deca + num(i+1, k+1, l+1)
                  call divide_prisme(s, nbtet, tetra,
     &	  	                       rnum, domtet, prisme)
                enddo


                l = k + 1
                prisme(1) = deca + num(i, k, l)
	            prisme(2) = deca + num(i, k+1, l)
                prisme(3) = deca + num(i, k+1, l+1)
	            prisme(4) = deca + num(i+1, k, l)
                prisme(5) = deca + num(i+1, k+1, l)
	            prisme(6) = deca + num(i+1, k+1, l+1)
	            call divide_prisme(s, nbtet, tetra,
     &	                         rnum, domtet, prisme)
              enddo
              k = i
	        do l = 1, k
                prisme(1) = deca + num(i+1, k+1, l)
	            prisme(2) = deca + num(i+1, k, l)
                prisme(3) = deca + num(i, k, l)
	            prisme(4) = deca + num(i+1, k+1, l+1)
                prisme(5) = deca + num(i+1, k, l+1)
	            prisme(6) = deca + num(i, k, l+1)
	          call divide_prisme(s, nbtet,tetra,
     &	                         rnum, domtet, prisme)
              enddo

              l = k+1
              nbtet = nbtet + 1
              domtet(nbtet) = s
	          nb = deca + num(i, k, l)
		      j = rnum(nb)
              tetra(1, nbtet) = j

	          nb = deca + num(i+1, k, l)
			  j = rnum(nb)
              tetra(2, nbtet) = j

	          nb = deca + num(i+1, k+1, l)
		      j = rnum(nb)
              tetra(3, nbtet) = j

	          nb = deca + num(i+1, k+1, l+1)
		      j = rnum(nb)
              tetra(4, nbtet) = j
 	       enddo
  		  enddo

          return
         end
c---------------------------------------------------------------------
c  Etat de modification: fini
c
           subroutine divide_prisme(s, nbtet,tetra,
     &	  			             rnum, domtet, prisme)

c         This subroutine makes some tetrahedra from
c          a given prism.
c          Parameters
c          s       : the number of the big domain to which the prism belongs
c                     (this is an input parameter)
c          nbtet   : the number of each tetrahedra
c                     At input  : the number  of the last constructed tetrahedra
c                     At output : the number  of the last constructed tetrahedra
c                                 constructed here
c          tetra   : the 2D array containing the definition of tetrahedra
c                   tetra(k, i) denotes the number of the k-vertice (k<=4)
c                     of the i-th tetrahedra.
c          rnum    : (input parameter) this array makes the correspondance
c                    between the old numbers of nodes (domain by domain)
c                    in which some nodes have multiple numbers
c                    and the new one (in which each node has one and only one number)
c          domtet  : this table contains the number of the big domain to which
c                    each tetrahedra belong  (output)
c          prisme  : input array) containing the numbers of vertices of the prism
            implicit none
            integer s, nbtet,tetra(4,*),
     &              rnum(*), prisme(6), domtet(*), j


c          -------------
	       nbtet = nbtet + 1
           domtet(nbtet) = s

	       j = rnum(prisme(1))
           tetra(1, nbtet) = j

	       j = rnum(prisme(2))
           tetra(2, nbtet) = j

	       j  = rnum(prisme(3))
           tetra(3, nbtet)  = j

	       j = rnum(prisme(6))
           tetra(4, nbtet) = j

c          -------------
	       nbtet = nbtet + 1
           domtet(nbtet) = s

	       j = rnum(prisme(1))
           tetra(1, nbtet) = j

	       j = rnum(prisme(2))
           tetra(2, nbtet) = j

	       j = rnum(prisme(5))
           tetra(3, nbtet) = j

	       j = rnum(prisme(6))
           tetra(4, nbtet) = j
c          -------------
	       nbtet = nbtet + 1
           domtet(nbtet) = s

	       j = rnum(prisme(1))
           tetra(1, nbtet) = j

	       j = rnum(prisme(4))
           tetra(2, nbtet) = j

	       j = rnum(prisme(5))
           tetra(3, nbtet) = j

	       j = rnum( prisme(6) )
           tetra(4, nbtet) = j
           return
	    end
c---------------------------------------------------------------------
c      Etat de modification :fini (aucun changement)

       subroutine radiale (N, mu, Rtot, di)
       implicit none
c       This subroutine computes the altitudes
c       of the  plans of the meshes

	   integer     N, i
	   double precision di(N), alpha(N), Rtot, mu, h, alp

        alpha(1) = 1
	   do i = 1, N-1
	    alp = alpha(i) + alpha(i)**(1.0 - mu)
        alpha(i+1) = min((i+1)**(1./mu), alp)
	   enddo
	   h = (Rtot/alpha(N))**mu
       do i = 1, N
        di(i) = alpha(i)*(h**(1.D0/mu))
 	   enddo
	  return
      end

c-------------------------------------------------------------------
c    Etat de modification: fini (aucun changement)
c
	     function  pow_sum(k, n)
	     implicit none
c          This function computes the sums
c          1^k + 2^k +....+n^k
c          for k = 1, 2 only!
c          T. Boulmezaoud June 2003

	     integer pow_sum, k, n

		 pow_sum = (n*(n+1)*(2*(k-1)*n + 1))/(4*k - 2)
		 return
	    end
c--------------------------------------------------------------------
c    Etat de modification: fini  (aucun changement)
c   modification fait
         subroutine calc_pas(nbtetra1, tetra1, xin1, domtet, hmax,
     &                        form_rate,  ndom)

c        This subroutine computes the maximal size of tetrahedra in
c        each big subdomain

         implicit none
         integer nbtetra1, tetra1(4, *), domtet(*)
         double precision xin1(3, *), hmax(*), form_rate

         integer  dom, nb, ndom, k, n1, j, n2
	     double precision res1, res2, res3, res4



	     do dom = 1, ndom
            hmax(dom) = 0.D0
	     enddo
         form_rate = 1.D6

	     do nb = 1, nbtetra1
	     res1 = 0.D0
	     res4  = 1.D6
	     dom = domtet(nb)
	      do k = 1, 3
	      n1 = tetra1(k, nb)
           do j = k+1, 4
            n2 = tetra1(j, nb)
            res3 =  (xin1(1, n1)-xin1(1, n2))*(xin1(1, n1)-xin1(1, n2))
     &            + (xin1(2, n1)-xin1(2, n2))*(xin1(2, n1)-xin1(2, n2))
     &            + (xin1(3, n1)-xin1(3, n2))*(xin1(3, n1)-xin1(3, n2))
            res3 = dsqrt(res3)
	        res1 = max(res1, res3)
	        res4 = min(res4, res3)
	       enddo
	      enddo
	      res4 = res4/res1
		  res2 = hmax(dom)
	     form_rate = min(form_rate, res4)
		 hmax(dom) = max(res2, res1)
	    enddo

	   end
c---------------------------------------------------------------------------
c       Etat de modification: fini
        subroutine big_simp_vertices_whsp(pform)

c       This suboutines gives a subdivision of the domain
c       into a finite number of infinite simplices
c       and in some cases (the whole space, the half-space,...etc)
c       a FEM region in the form of a big tetrahedra...
c       pform    : the size of the FEM-region
c       nbsimp   : the number of the subdomains (infinite simplices + bounded domain
c                  which is considered as a big tetrahedra.
c       vertex   : the vertices of the big tetrehedra
c       hei      : the height vectors of the infinite simplices
c       nrm_h    : the norm of the height vector
c       hmax     : the maximal size of elements
c                  of the mesh (or the graded mesh) in each domain

	    include 'simplexes_keltoum.h'

        double precision bt(3, 3), ibt(3, 3), ibt_t(3, 3), cc(3)
	    double precision hc(3), res, pform, res1, det
	    integer j, i, k



        do k = 1, 5
	     do j = 1, 4
	      do i = 1, 3
            vertex(i, j, k) = 0.D0
	      enddo
	     enddo
	    enddo
        vertex(1, 1 ,5) = 2*sqrt(2.D0)/3
        vertex(2, 1 ,5) = 0.0D0
        vertex(3, 1 ,5) = -1.D0/3.0

        vertex(1, 2 ,5) = - sqrt(2.D0)/3
        vertex(2, 2 ,5) = sqrt(2.D0/3.D0)
        vertex(3, 2 ,5) = -1.D0/3.0

        vertex(1, 3 ,5) = - sqrt(2.D0)/3
        vertex(2, 3 ,5) = - sqrt(2.D0/3.D0)
        vertex(3, 3 ,5) = -1.D0/3.0

        vertex(1, 4 ,5) = 0.0D0
        vertex(2, 4 ,5) = 0.0D0
        vertex(3, 4 ,5) = 1.D0
        do j = 1, 4
	      do i = 1, 3
            vertex(i,j,5) = vertex(i,j,5)*pform
	      enddo
	    enddo

c       Infinite simplices
c       The first vertex is the origin

        vertex(:, 2, 1) = vertex(:, 1, 5)
        vertex(:, 3, 1) = vertex(:, 2, 5)
        vertex(:, 4, 1) = vertex(:, 4, 5)

        vertex(:, 2, 2) = vertex(:, 2, 5)
        vertex(:, 3, 2) = vertex(:, 3, 5)
        vertex(:, 4, 2) = vertex(:, 4, 5)


        vertex(:, 2, 3) = vertex(:, 3, 5)
        vertex(:, 3, 3) = vertex(:, 1, 5)
        vertex(:, 4, 3) = vertex(:, 4, 5)


        vertex(:, 2, 4) = vertex(:, 1, 5)
        vertex(:, 3, 4) = vertex(:, 2, 5)
        vertex(:, 4, 4) = vertex(:, 3, 5)


         do i = 1, 3
	      cc(i) = 1.D0
	     enddo


        do k = 1, 5
         call bt_big_simp(k, bt)
         volume(k) = dabs(det(bt))
         call inversion(bt, ibt)
	     call transpos(3, 3, ibt, ibt_t)
  	     call  matv3(ibt_t, cc, hc)
	     call prodsc(hc, hc, res)
	     res1 = 0.D0
	     do i = 1, 3
            hei(i, k) = hc(i)/res
            res1 = res1 + hei(i, k)*hei(i, k)
	     enddo
	     nrm_h(k) = dsqrt(res1)
	    enddo

		return
      end

c------------------------------------------------------------------
      subroutine test_volumes(tol, info)
c      modification fait
	  implicit none
	  include 'fondements_keltoum.h'
	  include 'simplexes_keltoum.h'
	  include 'constants_keltoum.h'

	  integer  info
      double precision svol(5), tol

	  double precision  BK(3,3), det, res

	  integer  nnb, nbtet, nb, nm1, nm, k, i, s
	  external det


	  nbtet    = 3*nbprisme + N

	  do s = 1, 4   ! index of the big subdomain
       svol(s) = 0.0D0
       do nnb = 1, nbtet
          nb = nnb + (s-1)*nbtet
          nm1 = tetra1(1, nb)
!         the matrix of the small  tetrahedra
	      do k = 2, 4
             nm =  tetra1(k, nb)
             do i = 1, 3
               bk(i,k-1) = xin1(i,nm) - xin1(i, nm1)
             enddo
          enddo
          res = det(bk)
          svol(s) = svol(s)  + dabs(res)
       enddo
      enddo

c	  do s = 1, 5   ! index of the big subdomain
	  s = 5
       svol(s) = 0.0D0
       do nnb = 1, nt
          nb = nnb + (s-1)*nbtet
          nm1 = tetra1(1, nb)
!         the matrix of the small  tetrahedra
	      do k = 2, 4
             nm =  tetra1(k, nb)
             do i = 1, 3
               bk(i,k-1) = xin1(i,nm) - xin1(i, nm1)
             enddo
          enddo
          res = det(bk)
          svol(s) = svol(s)  + dabs(res)
       enddo
c      enddo

      info = 1
      res = 0.D0
	  do s = 1, 5   ! index of the big subdomain
       if (dabs((svol(s) - volume(s))/volume(s)).gt.tol) info = 0
      enddo
	  return
      end

