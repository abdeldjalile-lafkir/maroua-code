c         This file contains a declaration of the	
c         following common variables

c         nbsimp : the number of infinite simplices (a parameter)
c         vertex : an  array containing the coordinates
c                  of the vertices of the infinite simplices
c                  vertex(i, j, k) denotes the i-coordinates of
c                                 the j-th vertice of the k-th simplex
c         hei    : contains the altitude vector of the infinite simplices    



         integer nbsimp
         parameter (nbsimp = 5)
	     double precision vertex(3, 4, nbsimp), hei(3, nbsimp)
	     double precision nrm_h(nbsimp), hmax(nbsimp), volume(nbsimp)
	     common/datasimp/vertex, hei, nrm_h, hmax, volume

