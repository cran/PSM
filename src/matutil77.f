!-----------------------------------------------------------------------!
!                                                                       !
!     V A R I O U S    M A T H E M A T I C A L    U T I L I T I E S     !
!                                                                        !
!                       FORTRAN 77 PROCEDURES                           !
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
      SUBROUTINE DGPADM( IDEG,M,T,H,LDH,WSP,LWSP,IPIV,IEXPH,NS,IFLAG )

      IMPLICIT NONE
      INTEGER IDEG, M, LDH, LWSP, IEXPH, NS, IFLAG, IPIV(M)
      DOUBLE PRECISION T, H(LDH,M), WSP(LWSP)

!-----PURPOSE-----------------------------------------------------------!
!
!     COMPUTES EXP(T*H), THE MATRIX EXPONENTIAL OF A GENERAL MATRIX IN
!     FULL, USING THE IRREDUCIBLE RATIONAL PADE APPROXIMATION TO THE 
!     EXPONENTIAL FUNCTION EXP(X) = R(X) = (+/-)( I + 2*(Q(X)/P(X)) ),
!     COMBINED WITH SCALING-AND-SQUARING.
!
!-----ARGUMENTS---------------------------------------------------------!
!
!     IDEG      : (INPUT) THE DEGRE OF THE DIAGONAL PADE TO BE USED.
!                 A VALUE OF 6 IS GENERALLY SATISFACTORY.
!
!     M         : (INPUT) ORDER OF H.
!
!     H(LDH,M)  : (INPUT) ARGUMENT MATRIX.
!
!     T         : (INPUT) TIME-SCALE (CAN BE < 0).
!                  
!     WSP(LWSP) : (WORKSPACE/OUTPUT) LWSP .GE. 4*M*M+IDEG+1.
!
!     IPIV(M)   : (WORKSPACE)
!
!     IEXPH     : (OUTPUT) NUMBER SUCH THAT WSP(IEXPH) POINTS TO EXP(TH)
!                 I.E., EXP(TH) IS LOCATED AT WSP(IEXPH ... IEXPH+M*M-1)
!                                                                       
!                 NOTE: IF THE ROUTINE WAS CALLED WITH WSP(IPTR), 
!                       THEN EXP(TH) WILL START AT WSP(IPTR+IEXPH-1).
!
!     NS        : (OUTPUT) NUMBER OF SCALING-SQUARING USED.
!
!     IFLAG     : (OUTPUT) EXIT FLAG.
!                      0 - NO PROBLEM
!                     <0 - PROBLEM
!
!-----------------------------------------------------------------------!
!     ROGER B. SIDJE (RBS@MATHS.UQ.EDU.AU)
!     EXPOKIT: SOFTWARE PACKAGE FOR COMPUTING MATRIX EXPONENTIALS.
!     ACM - TRANSACTIONS ON MATHEMATICAL SOFTWARE, 24(1):130-156, 1998
!-----------------------------------------------------------------------!
!     MODIFIED TO RETURN A FLAG INSTEAD OF TERMINATING, WHEN TRYING
!     TO COMPUTE THE EXPONENTIAL OF A MATRIX WITH TOO LARGE ELEMENTS.
!
!     1) NIELS RODE KRISTENSEN, TECHNICAL UNIVERSITY OF DENMARK, 2000
!     2) ANDREAS S. CHRISTENSEN,TECHNICAL UNIVERSITY OF DENMARK, 2006 
!-----------------------------------------------------------------------!
!
      INTEGER MM,I,J,K,IH2,IP,IQ,IUSED,IFREE,IODD,ICOEF,IPUT,IGET
      DOUBLE PRECISION HNORM,SCALE,SCALE2,CP,CQ

      INTRINSIC INT,ABS,DBLE,LOG,MAX

!---  CHECK RESTRICTIONS ON INPUT PARAMETERS ...
      MM = M*M
      IFLAG = 0
      IF ( LDH.LT.M ) IFLAG = -1
      IF ( LWSP.LT.4*MM+IDEG+1 ) IFLAG = -2
      IF ( IFLAG.NE.0 ) RETURN
!
!---  INITIALISE POINTERS ...
!
      ICOEF = 1
      IH2 = ICOEF + (IDEG+1)
      IP  = IH2 + MM
      IQ  = IP + MM
      IFREE = IQ + MM
!
!---  SCALING: SEEK NS SUCH THAT ||T*H/2^NS|| < 1/2; 
!     AND SET SCALE = T/2^NS ...
!
      DO I = 1,M
         WSP(I) = 0.0D0
      ENDDO
      DO J = 1,M
         DO I = 1,M
            WSP(I) = WSP(I) + ABS( H(I,J) )
         ENDDO
      ENDDO
      HNORM = 0.0D0
      DO I = 1,M
         HNORM = MAX( HNORM,WSP(I) )
      ENDDO
      HNORM = ABS( T*HNORM )
!     IF ( HNORM.EQ.0.0D0 ) STOP 'ERROR - NULL H IN INPUT OF DGPADM.'
      NS = MAX( 0,INT(LOG(HNORM)/LOG(2.0D0))+2 )
      SCALE = T / DBLE(2**NS)
      SCALE2 = SCALE*SCALE
!
!---  COMPUTE PADE COEFFICIENTS ...
!
      I = IDEG+1
      J = 2*IDEG+1
      WSP(ICOEF) = 1.0D0
      DO K = 1,IDEG
         WSP(ICOEF+K) = (WSP(ICOEF+K-1)*DBLE( I-K ))/DBLE( K*(J-K) )
      ENDDO
!
!---  H2 = SCALE2*H*H ...
!
      CALL DGEMM( 'N','N',M,M,M,SCALE2,H,LDH,H,LDH,0.0D0,WSP(IH2),M )
!
!---  INITIALIZE P (NUMERATOR) AND Q (DENOMINATOR) ...
!
      CP = WSP(ICOEF+IDEG-1)
      CQ = WSP(ICOEF+IDEG)
      DO J = 1,M
         DO I = 1,M
            WSP(IP + (J-1)*M + I-1) = 0.0D0
            WSP(IQ + (J-1)*M + I-1) = 0.0D0
         ENDDO
         WSP(IP + (J-1)*(M+1)) = CP
         WSP(IQ + (J-1)*(M+1)) = CQ
      ENDDO
!
!---  APPLY HORNER RULE ...
!
      IODD = 1
      K = IDEG - 1
 100  CONTINUE
      IUSED = IODD*IQ + (1-IODD)*IP
      CALL DGEMM( 'N','N',M,M,M, 1.0D0,WSP(IUSED),M,
     .             WSP(IH2),M, 0.0D0,WSP(IFREE),M )
      DO J = 1,M
         WSP(IFREE+(J-1)*(M+1)) = WSP(IFREE+(J-1)*(M+1))+WSP(ICOEF+K-1)
      ENDDO
      IP = (1-IODD)*IFREE + IODD*IP
      IQ = IODD*IFREE + (1-IODD)*IQ
      IFREE = IUSED
      IODD = 1-IODD
      K = K-1
      IF ( K.GT.0 )  GOTO 100
!
!---  OBTAIN (+/-)(I + 2*(P\Q)) ...
!
      IF ( IODD .EQ. 1 ) THEN
         CALL DGEMM( 'N','N',M,M,M, SCALE,WSP(IQ),M,
     .                H,LDH, 0.0D0,WSP(IFREE),M )
         IQ = IFREE
      ELSE
         CALL DGEMM( 'N','N',M,M,M, SCALE,WSP(IP),M,
     .                H,LDH, 0.0D0,WSP(IFREE),M )
         IP = IFREE
      ENDIF
      CALL DAXPY( MM, -1.0D0,WSP(IP),1, WSP(IQ),1 )
      CALL DGESV( M,M, WSP(IQ),M, IPIV, WSP(IP),M, IFLAG )
!     IF ( IFLAG.NE.0 ) STOP 'PROBLEM IN DGESV (WITHIN DGPADM)'
      CALL DSCAL( MM, 2.0D0, WSP(IP), 1 )
      DO J = 1,M
         WSP(IP+(J-1)*(M+1)) = WSP(IP+(J-1)*(M+1)) + 1.0D0
      ENDDO
      IPUT = IP
      IF ( NS.EQ.0 .AND. IODD.EQ.1 ) THEN
         CALL DSCAL( MM, -1.0D0, WSP(IP), 1 )
         GOTO 200
      ENDIF
!
!---  SQUARING : EXP(T*H) = (EXP(T*H))^(2^NS) ...
!
      IODD = 1
      DO K = 1,NS
         IGET = IODD*IP + (1-IODD)*IQ
         IPUT = (1-IODD)*IP + IODD*IQ
         CALL DGEMM( 'N','N',M,M,M, 1.0D0,WSP(IGET),M, WSP(IGET),M,
     .                0.0D0,WSP(IPUT),M )
         IODD = 1-IODD
      ENDDO
 200  CONTINUE
      IEXPH = IPUT
      END
!-----------------------------------------------------------------------!
