 
       subroutine dec_gauss_lobatto(a, b, Ng, rwgau, rdgau, 
     &                         wgau, dgau)
 

	   integer Ng
	   double precision a, b, rdgau(*), rwgau(*), wgau(*), dgau(*) 
	   

	   integer i
         

	   do i = 1, Ng
	      dgau(i) =   (b - a)*(rdgau(i)+1.D0)*0.5D0 + a
	      wgau(i) =    rwgau(i)*(b-a)*0.5D0
	   enddo
	   return
	  end

c---------------------------------------------------------------------

      SUBROUTINE gauss_lobatto (Z,W,NP)                                         
C------------------------------------------------------------
C                                                                       
C     Generate NP Gauss-Lobatto Legendre points (Z) and weights (W)     
C     associated with Jacobi polynomial P(N)(alpha=0,beta=0).           
C     The polynomial degree N=NP-1.                                     
C     Z and W are in single precision, but all the arithmetic           
C     operations are done in double precision.                          
C                                                                       
C--------------------------------------------------------------------
      integer NP   
      DOUBLE PRECISION Z(NP),W(NP),ALPHA,BETA
      ALPHA = 0.D0
      BETA  = 0.D0
      CALL ZWGLJD(Z, W, NP, ALPHA, BETA)

      RETURN
      END

C     
C     
      SUBROUTINE ZWGLJD (Z,W,NP,ALPHA,BETA)
C--------------------------------------------------------------------
C     
C     Generate NP GAUSS LOBATTO JACOBI points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
C     The polynomial degree N=NP-1.
C     Double precision version.
C     
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Z(NP),W(NP),ALPHA,BETA
      INTEGER N,NM1
C     
      N     = NP-1
      NM1   = N-1
      ONE   = 1.D0
      TWO   = 2.D0
C     
      IF (NP.LE.1) THEN
         WRITE (6,*) 'Minimum number of Gauss-Lobatto points is 2'
         STOP
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'Alpha and Beta must be greater than -1'
         STOP
      ENDIF
C     
      IF (NM1.GT.0) THEN
         ALPG  = ALPHA+ONE
         BETG  = BETA+ONE
         CALL ZWGJD (Z(2),W(2),NM1,ALPG,BETG)
      ENDIF
      Z(1)  = -ONE
      Z(NP) =  ONE
      DO 100 I=2,NP-1
         W(I) = W(I)/(ONE-Z(I)**2)
 100  CONTINUE
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(1))
      W(1)  = ENDW1 (N,ALPHA,BETA)/(TWO*PD)
      CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,N,ALPHA,BETA,Z(NP))
      W(NP) = ENDW2 (N,ALPHA,BETA)/(TWO*PD)
C     
      RETURN
      END
C     
      DOUBLE PRECISION FUNCTION ENDW1 (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ALPHA,BETA
      ZERO  = 0.D0
      ONE   = 1.D0
      TWO   = 2.D0
      THREE = 3.D0
      FOUR  = 4.D0
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW1 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW1 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+ONE)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (-TWO*(BETA+TWO)*FINT1 + (APB+FOUR)*FINT2)
     $     * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW1 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = DBLE(I-1)
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   = -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW1  = F3
      RETURN
      END
C     
      DOUBLE PRECISION FUNCTION ENDW2 (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ALPHA,BETA
      ZERO  = 0.D0
      ONE   = 1.D0
      TWO   = 2.D0
      THREE = 3.D0
      FOUR  = 4.D0
      APB   = ALPHA+BETA
      IF (N.EQ.0) THEN
         ENDW2 = ZERO
         RETURN
      ENDIF
      F1   = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      F1   = F1*(APB+TWO)*TWO**(APB+TWO)/TWO
      IF (N.EQ.1) THEN
         ENDW2 = F1
         RETURN
      ENDIF
      FINT1 = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+TWO)/GAMMAF(APB+THREE)
      FINT1 = FINT1*TWO**(APB+TWO)
      FINT2 = GAMMAF(ALPHA+TWO)*GAMMAF(BETA+TWO)/GAMMAF(APB+FOUR)
      FINT2 = FINT2*TWO**(APB+THREE)
      F2    = (TWO*(ALPHA+TWO)*FINT1 - (APB+FOUR)*FINT2)
     $     * (APB+THREE)/FOUR
      IF (N.EQ.2) THEN
         ENDW2 = F2
         RETURN
      ENDIF
      DO 100 I=3,N
         DI   = DBLE(I-1)
         ABN  = ALPHA+BETA+DI
         ABNN = ABN+DI
         A1   =  -(TWO*(DI+ALPHA)*(DI+BETA))/(ABN*ABNN*(ABNN+ONE))
         A2   =  (TWO*(ALPHA-BETA))/(ABNN*(ABNN+TWO))
         A3   =  (TWO*(ABN+ONE))/((ABNN+TWO)*(ABNN+ONE))
         F3   =  -(A2*F2+A1*F1)/A3
         F1   = F2
         F2   = F3
 100  CONTINUE
      ENDW2  = F3
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION GAMMAF (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X
      ZERO = 0.0D0
      HALF = 0.5D0
      ONE  = 1.0D0
      TWO  = 2.0D0
      FOUR = 4.0D0
      PI   = FOUR*DATAN(ONE)
      GAMMAF = ONE
      IF (X.EQ.-HALF) GAMMAF = -TWO*DSQRT(PI)
      IF (X.EQ. HALF) GAMMAF =  DSQRT(PI)
      IF (X.EQ. ONE ) GAMMAF =  ONE
      IF (X.EQ. TWO ) GAMMAF =  ONE
      IF (X.EQ. 1.5D0) GAMMAF =  DSQRT(PI)/2.D0
      IF (X.EQ. 2.5D0) GAMMAF =  1.5D0*DSQRT(PI)/2.D0
      IF (X.EQ. 3.5D0) GAMMAF =  2.5D0*1.5D0*DSQRT(PI)/2.D0
      IF (X.EQ. 3.D0 ) GAMMAF =  2.D0
      IF (X.EQ. 4.D0 ) GAMMAF = 6.D0
      IF (X.EQ. 5.D0 ) GAMMAF = 24.D0
      IF (X.EQ. 6.D0 ) GAMMAF = 120.D0
      RETURN
      END
C
      SUBROUTINE JACOBF (POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,
     $                   N,ALP,BET,X)
C--------------------------------------------------------------------
C
C     Computes the Jacobi polynomial (POLY) and its derivative (PDER)
C     of degree N at X.
C
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      APB  = ALP+BET
      POLY = 1.D0
      PDER = 0.D0
      IF (N .EQ. 0) RETURN
      POLYL = POLY
      PDERL = PDER
      POLY  = (ALP-BET+(APB+2.D0)*X)/2.D0
      PDER  = (APB+2.D0)/2.D0
      IF (N .EQ. 1) RETURN
      DO 20 K=2,N
         DK = DFLOAT(K)
         A1 = 2.D0*DK*(DK+APB)*(2.D0*DK+APB-2.D0)
         A2 = (2.D0*DK+APB-1.D0)*(ALP**2-BET**2)
         B3 = (2.D0*DK+APB-2.D0)
         A3 = B3*(B3+1.D0)*(B3+2.D0)
         A4 = 2.D0*(DK+ALP-1.D0)*(DK+BET-1.D0)*(2.D0*DK+APB)
         POLYN  = ((A2+A3*X)*POLY-A4*POLYL)/A1
         PDERN  = ((A2+A3*X)*PDER-A4*PDERL+A3*POLY)/A1
         PSAVE  = POLYL
         PDSAVE = PDERL
         POLYL  = POLY
         POLY   = POLYN
         PDERL  = PDER
         PDER   = PDERN
 20   CONTINUE
      POLYM1 = POLYL
      PDERM1 = PDERL
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      RETURN
      END
C
      SUBROUTINE DGLL (D,DT,Z,NZ)
C-----------------------------------------------------------------
C
C     Compute the derivative matrix D and its transpose DT
C     associated with the Nth order Lagrangian interpolants
C     through the NZ Gauss-Lobatto Legendre points Z.
C     Note: D and DT are square matrices.
C
C-----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Z(NZ),D(NZ,NZ),DT(NZ,NZ)
      N  = NZ-1
      IF (NZ .GT. 65) THEN
         WRITE (6,*) 'Subroutine DGLL'
         WRITE (6,*) 'Maximum polynomial degree = 64'
         WRITE (6,*) 'Polynomial degree         = ',N
      ENDIF
      IF (NZ .EQ. 1) THEN
         D(1,1) = 0.D0
         RETURN
      ENDIF
      FN = DBLE (N)
      D0 = FN*(FN+1.D0)/4.D0
      DO 200 I=1,NZ
      DO 200 J=1,NZ
         D(I,J) = 0.D0
         IF  (I.NE.J) D(I,J) = PNLEG(Z(I),N)/
     $                        (PNLEG(Z(J),N)*(Z(I)-Z(J)))
         IF ((I.EQ.J).AND.(I.EQ.1))  D(I,J) = -D0
         IF ((I.EQ.J).AND.(I.EQ.NZ)) D(I,J) =  D0
         DT(J,I) = D(I,J)
 200  CONTINUE
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PNLEG (Z,N)
C---------------------------------------------------------------------
C
C     Compute the value of the Nth order Legendre polynomial at Z.
C     (Simpler than JACOBF)
C     Based on the recursion formula for the Legendre polynomials.
C
C---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      P1   = 1.D0
      P2   = Z
      P3   = P2
      DO 10 K = 1, N-1
         FK  = DBLE (K)
         P3  = ((2.D0*FK+1.D0)*Z*P2 - FK*P1)/(FK+1.D0)
         P1  = P2
         P2  = P3
 10   CONTINUE
      PNLEG = P3
      RETURN
      END
      SUBROUTINE ZWGJD (Z,W,NP,ALPHA,BETA)
C--------------------------------------------------------------------
C
C     Generate NP GAUSS JACOBI points (Z) and weights (W)
C     associated with Jacobi polynomial P(N)(alpha>-1,beta>-1).
C     The polynomial degree N=NP-1.
C     Double precision version.
C
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Z(NP),W(NP),ALPHA,BETA
C
      N     = NP-1
      DN    = DBLE(N)
      ONE   = 1.D0
      TWO   = 2.D0
      APB   = ALPHA+BETA
C
      IF (NP.LE.0) THEN
         WRITE (6,*) 'Minimum number of Gauss points is 1'
         STOP
      ENDIF
      IF ((ALPHA.LE.-ONE).OR.(BETA.LE.-ONE)) THEN
         WRITE (6,*) 'Alpha and Beta must be greater than -1'
         STOP
      ENDIF
C
      IF (NP.EQ.1) THEN
         Z(1) = (BETA-ALPHA)/(APB+TWO)
         W(1) = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)/GAMMAF(APB+TWO)
     $          * TWO**(APB+ONE)
         RETURN
      ENDIF
C
      CALL JACG (Z,NP,ALPHA,BETA)
C
      NP1   = N+1
      NP2   = N+2
      DNP1  = DBLE(NP1)
      DNP2  = DBLE(NP2)
      FAC1  = DNP1+ALPHA+BETA+ONE
      FAC2  = FAC1+DNP1
      FAC3  = FAC2+ONE
      FNORM = PNORMJ(NP1,ALPHA,BETA)
      RCOEF = (FNORM*FAC2*FAC3)/(TWO*FAC1*DNP2)
      DO 100 I=1,NP
         CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP2,ALPHA,BETA,Z(I))
         W(I) = -RCOEF/(P*PDM1)
 100  CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION PNORMJ (N,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ALPHA,BETA
      ONE   = 1.D0
      TWO   = 2.D0
      DN    = DBLE(N)
      CONST = ALPHA+BETA+ONE
      IF (N.LE.1) THEN
         PROD   = GAMMAF(DN+ALPHA)*GAMMAF(DN+BETA)
         PROD   = PROD/(GAMMAF(DN)*GAMMAF(DN+ALPHA+BETA))
         PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
         RETURN
      ENDIF
      PROD  = GAMMAF(ALPHA+ONE)*GAMMAF(BETA+ONE)
      PROD  = PROD/(TWO*(ONE+CONST)*GAMMAF(CONST+ONE))
      PROD  = PROD*(ONE+ALPHA)*(TWO+ALPHA)
      PROD  = PROD*(ONE+BETA)*(TWO+BETA)
      DO 100 I=3,N
         DINDX = DBLE(I)
         FRAC  = (DINDX+ALPHA)*(DINDX+BETA)/(DINDX*(DINDX+ALPHA+BETA))
         PROD  = PROD*FRAC
 100  CONTINUE
      PNORMJ = PROD * TWO**CONST/(TWO*DN+CONST)
      RETURN
      END
C
      SUBROUTINE JACG (XJAC,NP,ALPHA,BETA)
C--------------------------------------------------------------------
C
C     Compute NP Gauss points XJAC, which are the zeros of the
C     Jacobi polynomial J(NP) with parameters ALPHA and BETA.
C     ALPHA and BETA determines the specific type of Gauss points.
C     Examples:
C     ALPHA = BETA =  0.0  ->  Legendre points
C     ALPHA = BETA = -0.5  ->  Chebyshev points
C
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XJAC(NP)
      DATA KSTOP/10/
      DATA EPS/1.0D-12/
      N   = NP-1
      DTH = 4.D0*DATAN(1.D0)/(2.D0*DFLOAT(N)+2.D0)
      DO 40 J=1,NP
         IF (J.EQ.1) THEN
            X = DCOS((2.D0*(DFLOAT(J)-1.D0)+1.D0)*DTH)
         ELSE
            X1 = DCOS((2.D0*(DFLOAT(J)-1.D0)+1.D0)*DTH)
            X2 = XLAST
            X  = (X1+X2)/2.D0
         ENDIF
         DO 30 K=1,KSTOP
            CALL JACOBF (P,PD,PM1,PDM1,PM2,PDM2,NP,ALPHA,BETA,X)
            RECSUM = 0.D0
            JM = J-1
            DO 29 I=1,JM
               RECSUM = RECSUM+1.D0/(X-XJAC(NP-I+1))
 29         CONTINUE
            DELX = -P/(PD-RECSUM*P)
            X    = X+DELX
            IF (ABS(DELX) .LT. EPS) GOTO 31
 30      CONTINUE
 31      CONTINUE
         XJAC(NP-J+1) = X
         XLAST        = X
 40   CONTINUE
      DO 200 I=1,NP
         XMIN = 2.D0
         DO 100 J=I,NP
            IF (XJAC(J).LT.XMIN) THEN
               XMIN = XJAC(J)
               JMIN = J
            ENDIF
 100     CONTINUE
         IF (JMIN.NE.I) THEN
            SWAP = XJAC(I)
            XJAC(I) = XJAC(JMIN)
            XJAC(JMIN) = SWAP
         ENDIF
 200  CONTINUE
      RETURN
      END	   
