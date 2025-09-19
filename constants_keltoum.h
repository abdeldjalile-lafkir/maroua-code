c      This file contains the declaration of some 
c      constants         

c     NQ           : denotes the number of gauss-lobatto points 
c                    used in quadrature formula in a 2D triangle 
c                    (mainly involved in inverted elements)
c                    this parameter can be changed by the user
c     xq, yq, weiq : the coordinates of the quadrature points in 
c                     the reference triangle 
c     weig         : an array containing the corresponding weights
c
c     zero        : in some  "if" tests, a real x consider equal to 0
c                    if abs(x) < zero. The parameter "zero" can be changed
c                     by the user, but it must remain small 
c                     (smaller than 1.E-8 for example)        
c     cc           : the constant vector (1, 1, 1)
c
c             
c     NQTET        : denotes the number of points of a quadrature formula
c                    in a tetrahedra. Here NQTET = 10 (ten points)
c     XYZQT        : a matrix containing the coordinates of the quadrature  
c                    points in the reference tetrahedra
c     WEIQT        : the corresponding weights               
       
       integer nq
       parameter (nq=7) 
       double precision zero, pi, xq(nq), yq(nq), weiq(nq),cc(3)
       parameter (zero = 1.D-10, PI = 3.14159265358D0)
       common/quadraTrg/xq, yq, weiq, cc


       integer nqt
       parameter (nqt=10) 
       double precision xyzqt(3, nqt), weiqt(nqt)   
       common/quadraTetra/xyzqt, weiqt
