c       Etat de modification : fini
        Integer  N, nbnodes, nbtetra, Ng, nbprisme
        parameter (N = 20, Ng = 10,
     &     nbprisme = ((N-1)*N)/2 + (N*(N-1)*(2*N-1))/6,
     &     nbtetra  = 5*(3*nbprisme + N),
     &     nbnodes =  5*((N+1)*(N+2)*(N+3)/6-(N*N+N+1)))

        integer  tetra(4,nbtetra), domnod(nbnodes), domtet(nbtetra)
        common/datamailla/tetra, domnod, domtet       

        double precision  grad_grad(10, nbtetra), xin(3,nbnodes)
	    common/basicdata/grad_grad,  xin


  
