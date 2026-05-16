!       Etat de modification : fini
        Integer  N, nbnodes, nbtetra, Ng, nbprisme,gm, nt,
     &  nbnodes1, nbtetra1
        parameter (N = 20, Ng = 10, gm= 92, nt= 310,
     &     nbprisme = ((N-1)*N)/2 + (N*(N-1)*(2*N-1))/6,
     &     nbtetra  = 4*(3*nbprisme + N),
     &     nbtetra1  = 4*(3*nbprisme + N)+ nt,
     &     nbnodes = 5*(((N+1)*(N+2)*(N+3))/6-(N*N+N+1)),
     &     nbnodes1 = 5*(((N+1)*(N+2)*(N+3))/6-(N*N+N+1))+gm)

        integer  tetra(4,nbtetra), domnod(nbnodes1),
     &   domtet(nbtetra1),tetra1(4,nbtetra1)

        common/datamailla/tetra, domnod, domtet, tetra1

        double precision  grad_grad(10, nbtetra1),
     &    xin(3,nbnodes), xin1(3,nbnodes1)
	    common/basicdata/grad_grad, xin,  xin1
