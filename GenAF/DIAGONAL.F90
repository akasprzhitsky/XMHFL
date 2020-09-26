! MODULE OF MATRIX DIAGONALIZATION SUB-PROGRAMS 
! VER 2.0 NEW 02 2006

module diagonal
 IMPLICIT REAL*8(A-H,O,P,R-Z)
      
	
      
 contains
      
 subroutine DIAGonalXX(N,A,B,Rval,Rvector,pi)
  !use msimsl
  integer(4)::N
  real(8)::pi
  real(8),dimension(:)::Rval
  real(8),dimension(:,:)::A,B,Rvector
      
  !call dgvcsp(N,A,N,B,N,Rval,Rvector,N)
  !pi=dgpisp(N,N,A,N,B,N,Rval,Rvector,N)
  return 
 end subroutine DIAGonalXX


 SUBROUTINE DIAG(NL,RZ,SR,RD,RZER)
   IMPLICIT REAL*8(A-H,O,P-Z)
   integer::NL
   real(8),dimension(:)::RD
   real(8),dimension(:,:)::RZ,SR,RZER
   real(8),allocatable,dimension(:,:)::ZZ,S,Rvector     
   real(8),allocatable,dimension(:)::DL,Rvalue
   ! WE ALWAYS A FIVE UNDER ARRAYS
   allocate(ZZ(size(RZ,dim=1),size(RZ,dim=2)),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "ZZ" IS NOT SELECTED'
	  stop 
   endif
   allocate(S(size(SR,dim=1),size(SR,dim=2)),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "S" IS NOT SELECTED'
	  stop 
   endif
   allocate(DL(size(SR,dim=1)),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "DL" IS NOT SELECTED'
	  stop 
   endif
   allocate(Rvector(size(SR,dim=1),size(SR,dim=2)),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "Rvector" IS NOT SELECTED'
	  stop 
   endif
   allocate(Rvalue(size(SR,dim=1)),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "Rvalue" IS NOT SELECTED'
	  stop 
   endif


   DO I=1,NL
      DO J=1,NL
         ZZ(I,J)=RZ(I,J)
         S(I,J)=SR(I,J)      
      ENDDO
   ENDDO	   
   KL=1
   N=NL
   NE=NL
   NK=NL
   call RTQLR(1,N,NE,NK,KL,Rvalue,Rvector,ZZ,S,DL)
   DO I=1,NL
      RD(I)=Rvalue(I)    
      DO J=1,NL
         RZER(I,J)=Rvector(I,J)  
      ENDDO
   ENDDO	  

      
  ! REMOVING ARRAYS FROM MEMORY
  deallocate(ZZ,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Diagonal THE FILE "ZZ" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(S,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Diagonal THE FILE "S" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(DL,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Diagonal THE FILE "DL" IS NOT REMOVED FROM MEMORY'
 	 stop 
  endif
   deallocate(Rvector,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Diagonal THE FILE "Rvector" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(Rvalue,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Diagonal THE FILE "Rvalue" IS NOT REMOVED FROM MEMORY'
 	 stop 
  endif
     
  RETURN
END SUBROUTINE

      SUBROUTINE RTQLR(IFLAG,N,NE,NK,KL,D,RZX,ZZ,S,DL)
      IMPLICIT REAL*8(A-H,O,P-Z)
	real(8),dimension(:,:):: ZZ,S,RZX
      real(8),dimension(:)::D,DL
      real(8),allocatable,dimension(:)::E
      
	! WE ALWAYS A FIVE UNDER ARRAYS
	allocate(E(size(D)),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'Diagonal MEMORY ON THE FILE "E" IS NOT SELECTED'
	stop 
	endif

      IF(IFLAG.EQ.1) CALL REDUC(N,ZZ,S,DL)
      CALL TRED2(ZZ,D,E,N,1.0D-50)
      CALL TQL2(N,1.0D-50,D,E,ZZ,KL)
      IF(IFLAG.EQ.1) CALL REBAK(N,ZZ,S,DL)
      RZX=ZZ 


      ! REMOVING ARRAYS FROM MEMORY
      deallocate(E,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'Diagonal THE FILE "E" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      RETURN
      END SUBROUTINE

 
      


    SUBROUTINE REDUC(N,A,B,DL)
      IMPLICIT REAL*8(A-H,O,P-Z)
	  real(8),dimension(:,:):: A,B
      real(8),dimension(:)::DL
      
      DO 1 I=1,N
      DO 2 J=I,N
      X=B(I,J)
      K1=I-1
      IF(K1.LE.0) GOTO 5
      DO 3 IK=1,K1
      K=K1+1-IK
3     X=X-B(I,K)*B(J,K)
5     IF(I.NE.J) GOTO 4
      IF(X.LE.0.D0) GOTO 7
      Y=DSQRT(X)
      DL(I)=Y
      GOTO 2
4     B(J,I)=X/Y
2     CONTINUE
1     CONTINUE
      DO 10 I=1,N
      Y=DL(I)
      DO 11 J=I,N
      X=A(I,J)
      K1=I-1
      IF(K1.LE.0) GOTO 12
      DO 13 IK=1,K1
      K=K1+1-IK
13    X=X-B(I,K)*A(J,K)
12    A(J,I)=X/Y
11    CONTINUE
10    CONTINUE
      DO 30 J=1,N
      DO 31 I=J,N
      X=A(I,J)
      K1=I-1
      IF(K1.LT.J) GOTO 32
      DO 33 IK=J,K1
      K=K1-IK+J
33    X=X-A(K,J)*B(I,K)
32    K1=J-1
      IF(K1.LT.1) GOTO 35
      DO 36 IK=1,K1
      K=K1+1-IK
36    X=X-A(J,K)*B(I,K)
35    A(I,J)=X/DL(I)
31    CONTINUE
30    CONTINUE
!      DO 40 I=1,N
!      DO 40 J=1,N
!40    A(I,J)=A(J,I)

!      DO 41 I=1,N
!      WRITE(6,100)(A(I,J),J=1,N)
!41    CONTINUE
!      DO 42 I=1,N
!      WRITE(6,100)(B(I,J),J=1,N)
!42    CONTINUE
100   FORMAT(10(1X,F10.4)) 
      RETURN
7     WRITE(11,875)
875   FORMAT(//20X,'OëíAçéÇ - èêéñÖÑìêÄ "REDUC"')
      STOP
      END SUBROUTINE

      SUBROUTINE REBAK(N,Z,B,DL)
      IMPLICIT REAL*8(A-H,O,P-Z)
	real(8),dimension(:,:):: Z,B
      real(8),dimension(:)::DL
     
      DO 1 J=1,N
      DO 2 I1=1,N
      I=N+1-I1
      X=Z(I,J)
      K1=I+1
      IF(K1.GT.N) GOTO 3
      DO 4 K=K1,N
4     X=X-B(K,I)*Z(K,J)
3     Z(I,J)=X/DL(I)
2     CONTINUE
1     CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE TQL2(N,MACH,D,E,Z,KL)
      IMPLICIT REAL*8(A-H,O,P-Z)
      REAL*8 MACH
     	real(8),dimension(:,:):: Z
      real(8),dimension(:)::E,D
      DO 1 I=2,N
1     E(I-1)=E(I)
      E(N)=0.D0
      B=0.D0
      F=0.D0
      DO 2 L=1,N
      J=0
      H=MACH*(DABS(D(L))+DABS(E(L)))
      IF(B.LT.H) B=H
      DO 3 M=L,N
      IF(DABS(E(M)).LE.B) GOTO 4
3     CONTINUE
4     IF(M.EQ.L) GOTO 5
6     IF(J.EQ.30) GOTO 7
      J=J+1
      G=D(L)
      P=(D(L+1)-G)/(2.*E(L))
      R=DSQRT(P*P+1.D0)
      IF(P.LT.0.D0) GOTO 8
      D(L)=E(L)/(P+R)
      GOTO 9
8     D(L)=E(L)/(P-R)
9     H=G-D(L)
      L1=L+1
      IF(L1.GT.N) GOTO 100
      DO 10 I=L1,N
      D(I)=D(I)-H
10    CONTINUE
100   F=F+H
      P=D(M)
      C=1.D0
      S=0.D0
      M1=M-1
      IF(L.GT.M1) GOTO 110
      DO 11 I=L,M1
      I1=L-I+M1
      G=C*E(I1)
      H=C*P
      IF(DABS(P).LT.DABS(E(I1))) GOTO 12
      C=E(I1)/P
      R=DSQRT(C*C+1.D0)
      E(I1+1)=S*P*R
      S=C/R
      C=1.D0/R
      GOTO 13
12    C=P/E(I1)
      R=DSQRT(C*C+1.D0)
      E(I1+1)=S*E(I1)*R
      S=1.D0/R
      C=C/R
13    P=C*D(I1)-S*G
      D(I1+1)=H+S*(C*G+S*D(I1))
      DO 11 K=1,N
      H=Z(K,I1+1)
      Z(K,I1+1)=S*Z(K,I1)+C*H
      Z(K,I1)=C*Z(K,I1)-S*H
11    CONTINUE
110   E(L)=S*P
      D(L)=C*P
      IF(DABS(E(L)).GT.B) GOTO 6
5     D(L)=D(L)+F
2     CONTINUE
      IF(KL.EQ.0) GOTO 16
      DO 14 I=1,N
      K=I
      P=D(I)
      I1=I+1
      IF(I1.GT.N) GOTO 150
      DO 15 J=I1,N
      IF(D(J).GE.P) GOTO 15
      K=J
      P=D(J)
15    CONTINUE
150   IF(K.EQ.I) GOTO 14
      D(K)=D(I)
      D(I)=P
      DO 140 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
      Z(J,K)=P
140   CONTINUE
14    CONTINUE
16    RETURN
7     WRITE(11,677)
677   FORMAT(//20X,' éëíéçéÇ - èêéñÖÑìêÄ "TQL2"')
      STOP
      END SUBROUTINE

      SUBROUTINE TRED2(Z,D,E,N,TOL)
      IMPLICIT REAL*8(A-H,O-Z)
      real(8),dimension(:,:):: Z
      real(8),dimension(:)::E,D
      DO 1 I=2,N
      I1=N+2-I
      L=I1-2
      F=Z(I1,I1-1)
      G=0.D0
      IF(L.LT.1) GOTO 20
      DO 2 K=1,L
      G=G+Z(I1,K)*Z(I1,K)
2     CONTINUE
20    H=G+F*F
      IF(G.GT.TOL) GOTO 3
      E(I1)=F
      H=0.D0
      GOTO 15
3     L=L+1
      G=-DSQRT(H)
      IF(F) 4,5,5
4     G=-G
5     E(I1)=G
      H=H-F*G
      Z(I1,I1-1)=F-G
      F=0.D0
      IF(L.LT.1) GOTO 80
      DO 8 J=1,L
      Z(J,I1)=Z(I1,J)/H
      G=0.D0
      DO 9 K=1,J
9     G=G+Z(J,K)*Z(I1,K)
      J1=J+1
      IF(L.LT.J1) GOTO 100
      DO 10 K=J1,L
      G=G+Z(K,J)*Z(I1,K)
10    CONTINUE
100   E(J)=G/H
      F=F+G*Z(J,I1)
8     CONTINUE
80    HH=F/(H+H)
      IF(L.LT.1) GOTO 60
      DO 6 J=1,L
      F=Z(I1,J)
      G=E(J)-HH*F
      E(J)=G
      DO 11 K=1,J 
11    Z(J,K)=Z(J,K)-F*E(K)-G*Z(I1,K)
6     CONTINUE
60    CONTINUE
15    D(I1)=H
1     CONTINUE
      D(1)=0.D0
      E(1)=0.D0
      DO 70 I=1,N
      L=I-1
      IF(DABS(D(I)).LT.1.0D-10) GOTO 120
      IF(L.LT.1) GOTO 120
      DO 12 J=1,L
      G=0.D0
      DO 13 K=1,L
13    G=G+Z(I,K)*Z(K,J)
      DO 12 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
12    CONTINUE
120   CONTINUE
      D(I)=Z(I,I)
      Z(I,I)=1.D0
      IF(L.LT.1) GOTO 70
      DO 7 J=1,L
      Z(I,J)=0.D0
      Z(J,I)=0.D0
7     CONTINUE
70    CONTINUE
      RETURN
      END SUBROUTINE
      end module diagonal