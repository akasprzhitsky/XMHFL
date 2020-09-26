      
! THE MODULE OF CALCULATION 3J-SYMBOLS VER 3.0 11.2019
! VER 3.0 NEW  11,2019 √Œƒ

!! MODULE FOR CALCULATION OF 3J-SYMBOLS VERSION 1.0 


	
module mc3js
 implicit none
    
	
 contains

 !! The SUB-SUBPROGRAMME FORMATES THE MASSIVE FOR CALCULATION OF 3J-SYMBOLS
 subroutine C3JS_VAR(NFF,FF,F) 
  implicit none
  integer::NFF
  real(8)::FF
  real(8),dimension(:)::F
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::ISD
 !! ZANULYAYEM BEFORE PAYMENT
  F=0.D0
  FF=1.D0
  F(1)=0.1D0
  DO ISD=2,NFF
     F(ISD)=F(ISD-1)*FLOAT(ISD)*0.1D0 
  ENDDO

  return
 end subroutine C3JS_VAR



!! SUBPROGRAMME OF CORNER PARTS CALCULATION OF CULONIAN INTEGRAL IN MEDIUM FIELD
† !! L1-ORBITAL MOMENT OF THE FIRST FUNCTION
† !! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
† !! K-ORBITAL MOMENT OF SPHERICAL HARMONIO
† !! Q1, Q2-NUMBERS OF FILLING
† !! AUXILIARY PARAMETERS FOR CALCULATION OF FF, FF (XX)
 real(8) function C3JS_GKK1(L1,L2,K,Q1,Q2,FF,F)
  implicit none
  integer::L1,L2,K,Q1,Q2
  real(8)::FF
  real(8),dimension(:)::F
      C3JS_GKK1=C3JS_CKK(L1,L2,K,FF,F)**2*float(Q1*Q2)/float((4*L1+2)*(2*L2+1))
  return
 end function C3JS_GKK1
 
    



	 
      
!! SUBPROGRAMME OF CALCULATION OF A MATRIX ELEMENT OF SPHERICAL HARMONIC
† !! L1-ORBITAL MOMENT OF THE FIRST FUNCTION
† !! M1-ORBITAL FIRST FUNCTION
† !! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
† !! M2-ORBITAL MOMENT OF THE SECOND FUNCTION
† !! K-ORBITAL MOMENT OF SPHERICAL HARMONIO
† !! The projection of a spherical harmonic is determined by the equality q = M1-M2
† !! The other projections lead to the zonation of the matrix element
† !! AUXILIARY PARAMETERS FOR CALCULATION OF FF, F (50)
 real(8) function C3JS_Ckq(L1,M1,L2,M2,K,FF,F)
  implicit none
  integer::L1,M1,L2,M2,K
  real(8)::FF
  real(8),dimension(:)::F
  !!!!!!!!!!!!!!!!!!!!!!!
  integer::II
  real(8)::A1,A2   
 
  II= L1+(L1+K+L2)/2.D0+1.D-1
  ! –¿—◊≈“ 3J-—»Ã¬ŒÀ¿,  Œ›‘‘»÷»≈Õ“¿ ¬»√Õ≈–¿ 
  A1=C3JS_W3JA(L1,K,L2,-M1,M1-M2,M2,FF,F)
  ! –¿—◊≈“ œ–»¬≈ƒ≈ÕÕŒ√Œ Ã¿“–»◊ÕŒ√Œ ›À≈Ã≈Õ“¿ Œœ≈–¿“Œ–¿ —‘≈–»◊≈— Œ… √¿–ÃŒÕ» » 
  A2=C3JS_CKK(L1,L2,K,FF,F)
  ! –≈«≈À‹“¿“
  C3JS_Ckq=(-1.D0)**(L1-M1+II)*A1*A2
  return
 end function C3JS_Ckq


      
  
!! PROGRAM OF CALCULATION OF 3J-SYMBOL
!! I1, I2, I3, I4, I5, I6-PARAMETERS OF THE SYMBOL
!! AUXILIARY PARAMETERS FOR CALCULATION OF FF, FF (50)
real(8) function C3JS_W3JA(I1,I2,I3,I4,I5,I6,FF,F)
 implicit none
 integer::I1,I2,I3,I4,I5,I6
 real(8)::FF
 real(8),dimension(:)::F
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::MIN,MAX,I,J
 integer::III1,III2,III3,III4,III5,III6,III7,III8,III9,III10
 real(8)::A1,A2,A3,A4,A5,A6,SUM,D1,D2,RRFUN,RRDSQRT
 
 455 FORMAT(1X,'VVV',100(1X,I4)) 
 456 FORMAT(1X,'NNN',100(1X,I4)) 

 ! œ¿–¿Ã≈“–€ 
 A1=DBLE(I1) 
 A2=DBLE(I2)
 A3=DBLE(I3)
 A4=DBLE(I4)
 A5=DBLE(I5)
 A6=DBLE(I6)
      
	
 C3JS_W3JA=0.D0

 IF(IDINT(A4+A5+A6+.1D0).NE.0) return
 IF(DABS(A1+A2+A3-DINT(A1+A2+A3+.1D0)).GT..1D0) return

 D1=DMAX1(A2-A3-A4,A1-A3+A5,0.D0)
 D2=DMIN1(A1+A2-A3,A1-A4,A2+A5)
 MIN=IDINT(D1)
 MAX=IDINT(D2) 
 IF(MIN.GT.MAX) return
 
 III1=IDINT(A1+A4+.1D0)
 III2=IDINT(A1-A4+.1D0)
 III3=IDINT(A2+A5+.1D0)
 III4=IDINT(A2-A5+.1D0)
 III5=IDINT(A3+A6+.1D0)
 III6=IDINT(A3-A6+.1D0)
 III7=IDINT(A1+A2-A3+.1D0)
 III8=IDINT(A1-A2+A3+.1D0)
 III9=IDINT(A2+A3-A1+.1D0)
 III10=IDINT(A1+A2+A3+1.1D0)
 !WRITE(100,*) 
 !WRITE(100,455) III1,III2,III3,III4,III5
 !WRITE(100,455) III6,III7,III8,III9,III10

 RRDSQRT=DSQRT(C3JS_F21(III1,FF,F))*DSQRT(C3JS_F21(III2,FF,F))*DSQRT(C3JS_F21(III3,FF,F))*DSQRT(C3JS_F21(III4,FF,F))*DSQRT(C3JS_F21(III5,FF,F))*DSQRT(C3JS_F21(III6,FF,F))*DSQRT(C3JS_F21(III7,FF,F))*DSQRT(C3JS_F21(III8,FF,F))*DSQRT(C3JS_F21(III9,FF,F))*DSQRT(0.1D0)/DSQRT(C3JS_F21(III10,FF,F))
 
 SUM=0.D0
 DO J=MIN,MAX
    III1=IDINT(A3-A2+A4+.1D0+J)
    III2=IDINT(A3-A1-A5+.1D0+J) 
    III3=IDINT(A1+A2-A3+.1D0-J)
    III4=IDINT(A1-A4+.1D0-J)
    III5=IDINT(A2+A5+.1D0-J)
	!WRITE(100,456) J,III1,III2,III3,III4,III5
    RRFUN=C3JS_F21(J,FF,F)*C3JS_F21(III1,FF,F)*C3JS_F21(III2,FF,F)*C3JS_F21(III3,FF,F)*C3JS_F21(III4,FF,F)*C3JS_F21(III5,FF,F)
	SUM=SUM +(-1.D0)**J*RRDSQRT/RRFUN
 ENDDO


 I=IDINT(DABS(A1-A2-A6)+.1D0)

 ! –≈«”À‹“¿“ –¿—◊≈“¿ 
 C3JS_W3JA=(-1.D0)**I*SUM
 
 return
end function C3JS_W3JA



! AUXILIARY PROCEDURE 
real(8) function C3JS_F21(IZX,FF,F)
 implicit none
 integer::IZX
 real(8)::FF
 real(8),dimension(:)::F
    
 IF(IZX.EQ.0) THEN
   C3JS_F21=FF
  ELSE
   C3JS_F21=F(IZX)
 ENDIF 
 return
end function C3JS_F21
      


!! CALCULATION OF THE APPLIED MATRIX ELEMENT OF SPHERICAL HARMONIC
!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
!! L1-ORBITAL MOMENT OF THE FIRST HARMONIC
!! L2-ORBITAL MOMENT OF THE SECOND HARMONIC
!! K-ORBITAL MOMENT OF SPHERICAL HARMONIC
real(8) function C3JS_CKK(L1,L2,K,FF,F)
 implicit none
 integer::L1,L2,K
 real(8)::FF
 real(8),dimension(:)::F
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::M,N
 real(8)::A,B,C,DGH,FHF


 C3JS_CKK=0.D0
 M=L1+L2+K+2
 N=M/2
 IF((M-N-N).EQ.0) THEN 
    A=DBLE(L1)
    B=DBLE(L2)
    C=DBLE(K)
    B=C3JS_DELTA1(A,B,C,FF,F)
   ELSE
   return
 ENDIF

 IF(B.GT.0.D0) THEN 
   DGH=B*DSQRT(DBLE((L1+L1+1)*(L2+L2+1))*0.1D0)*C3JS_F31(N,FF,F)
   FHF=C3JS_F31(N-L1,FF,F)*C3JS_F31(N-L2,FF,F)*C3JS_F31(N-K,FF,F)
   C3JS_CKK=DGH/FHF
 ENDIF 

 return
end function C3JS_CKK

 real(8) function C3JS_DELTA1(A,B,C,FF,F)
  implicit none
  real(8)::A,B,C,FF
  real(8),dimension(:)::F
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::I,J,K,M
  real(8)::D,DD,RRFG

  C3JS_DELTA1=0.D0

  D=A+B+C+1.1D0
  I=IDINT(D-C-C)
  J=IDINT(D-B-B)
  K=IDINT(D-A-A)
      
  IF(I*J*K.GT.0) THEN 
      M=IDINT(D)
      DD=M+.1D0-D
     ELSE 
      return
  ENDIF

  IF(IDINT(100D0*DD).GE.0) THEN 
     RRFG=C3JS_F31(I,FF,F)*C3JS_F31(J,FF,F)*C3JS_F31(K,FF,F)
     C3JS_DELTA1=DSQRT(RRFG/C3JS_F31(M+1,FF,F))
  ENDIF

  return
 end function C3JS_DELTA1

!! AUXILIARY PROCEDURE
 real(8) function C3JS_F31(IZX,FF,F)
  implicit none
  integer::IZX
  real(8)::FF
  real(8),dimension(:)::F
    
  IF(IZX.EQ.1) THEN
     C3JS_F31=FF
    ELSE
     C3JS_F31=F(IZX-1)
  ENDIF 
  return
 end function C3JS_F31


end module mc3js
