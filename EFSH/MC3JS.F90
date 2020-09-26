      
! THE MODULE OF CALCULATION 3J-SYMBOLS VER 3.0 11.2005
! VER 3.0 NEW  11,2017 цнд
! ver 4.0      03,2019 ЦНДЮ 

 


	
module mc3js
 implicit none
    
	
 contains

! AUXILIARY PROGRAM FORMS AN ARRAY FOR CALCULATING 3J CHARACTERS
 subroutine C3JS_VAR(NFF,FF,F) 
  implicit none
  integer::NFF
  real(8)::FF
  real(8),dimension(:)::F
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::ISD
  ! гюмскъел оепед пюявернл
  F=0.D0
  FF=1.D0
  F(1)=0.1D0
  DO ISD=2,NFF
     F(ISD)=F(ISD-1)*FLOAT(ISD)*0.1D0 
  ENDDO

  return
 end subroutine C3JS_VAR
	 
      
 ! SUB-PROGRAM FOR CALCULATING THE MATRIX ELEMENT OF THE SPHERICAL HARMONY
  ! L1-ORBITAL MOMENT OF FIRST FUNCTION
  ! M1-PROJECTION OF ORBITAL FIRST FUNCTION
  ! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
  ! M2-PROJECTION OF THE ORBITAL MOMENT OF THE SECOND FUNCTION
  ! K-ORBITAL MOMENT OF SPHERICAL HARMONY
  ! THE PROJECTION OF SPHERICAL HARMONICS IS DEFINED BY THE EQUALITY q = M1-M2
  ! OTHER PROJECTIONS LEAD TO ZEROING OF THE MATRIX ELEMENT
  ! AUXILIARY PARAMETERS FOR CALCULATION OF FF, FF (50)
 real(8) function C3JS_Ckq(L1,M1,L2,M2,K,FF,F)
  implicit none
  integer::L1,M1,L2,M2,K
  real(8)::FF
  real(8),dimension(:)::F
  !!!!!!!!!!!!!!!!!!!!!!!
  integer::II
  real(8)::A1,A2   
 
  II= L1+(L1+K+L2)/2.D0+1.D-1
  ! CALCULATION OF THE 3J SYMBOL, WIGNER RATIO
  A1=C3JS_W3JA(L1,K,L2,-M1,M1-M2,M2,FF,F)
 ! CALCULATION OF THE REDUCED MATRIX ELEMENT OF THE SPHERICAL HARMONIC OPERATOR
  A2=C3JS_CKK(L1,L2,K,FF,F)
  ! RESULT
  C3JS_Ckq=(-1.D0)**(L1-M1+II)*A1*A2
  return
 end function C3JS_Ckq

  ! SUB-PROGRAM FOR CALCULATING THE MATRIX ELEMENT OF THE SPHERICAL HARMONY
  ! L1-ORBITAL MOMENT OF FIRST FUNCTION
  ! M1-PROJECTION OF ORBITAL FIRST FUNCTION
  ! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
  ! M2-PROJECTION OF THE ORBITAL MOMENT OF THE SECOND FUNCTION
  ! K-ORBITAL MOMENT OF SPHERICAL HARMONY
  ! THE PROJECTION OF SPHERICAL HARMONICS IS DEFINED BY THE EQUALITY q = M1-M2
  ! OTHER PROJECTIONS LEAD TO ZEROING OF THE MATRIX ELEMENT
 real(8) function C3JS_CkqGG(L1,M1,L2,M2,K)
  implicit none
  integer::L1,M1,L2,M2,K
  !!!!!!!!!!!!!!!!!!!!!!!
  integer::II
  real(8)::A1,A2   
 
  II= L1+(L1+K+L2)/2.D0+1.D-1
  ! пюявер 3J-яхлбнкю, йнщттхжхемрю бхцмепю 
  A1=C3JS_W3JAGG(L1,K,L2,-M1,M1-M2,M2)
  ! пюявер опхбедеммнцн люрпхвмнцн щкелемрю ноепюрнпю ятепхвеяйни цюплнмхйх 
  A2=C3JS_CKKGG(L1,L2,K)
  ! пегекэрюр
  C3JS_CkqGG=(-1.D0)**(L1-M1+II)*A1*A2
  return
 end function C3JS_CkqGG
      
  
! 3J-SYMBOL CALCULATION SUB-PROGRAM
! I1, I2, I3, I4, I5, I6-SYMBOL PARAMETERS
! AUXILIARY PARAMETERS FOR CALCULATION OF FF, FF (50)
real(8) function C3JS_W3JA(I1,I2,I3,I4,I5,I6,FF,F)
 implicit none
 integer::I1,I2,I3,I4,I5,I6
 real(8)::FF,RR,RRRR
 real(8),dimension(:)::F
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::MIN,MAX,I,J
 integer::III1,III2,III3,III4,III5,III6,III7,III8,III9,III10
 real(8)::A1,A2,A3,A4,A5,A6,SUM,D1,D2,RRFUN,RRDSQRT



 ! оюпюлерпш 
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
 RRDSQRT=DSQRT(C3JS_F21(III1,FF,F))*DSQRT(C3JS_F21(III2,FF,F))*DSQRT(C3JS_F21(III3,FF,F))*DSQRT(C3JS_F21(III4,FF,F))*DSQRT(C3JS_F21(III5,FF,F))*DSQRT(C3JS_F21(III6,FF,F))*DSQRT(C3JS_F21(III7,FF,F))*DSQRT(C3JS_F21(III8,FF,F))*DSQRT(C3JS_F21(III9,FF,F))*DSQRT(0.1D0)/DSQRT(C3JS_F21(III10,FF,F))
 !WRITE(6,*) 'WW3'
 !WRITE(6,*) III1,III2,III3,III4,III5,III6,III7,III8,III9,III10
 SUM=0.D0
 DO J=MIN,MAX
    III1=IDINT(A3-A2+A4+.1D0+J)
    III2=IDINT(A3-A1-A5+.1D0+J) 
    III3=IDINT(A1+A2-A3+.1D0-J)
    III4=IDINT(A1-A4+.1D0-J)
    III5=IDINT(A2+A5+.1D0-J)
	!WRITE(6,*) 'SIK',III1,III2,III3,III4,III5
	RRFUN=C3JS_F21(J,FF,F)*C3JS_F21(III1,FF,F)*C3JS_F21(III2,FF,F)*C3JS_F21(III3,FF,F)*C3JS_F21(III4,FF,F)*C3JS_F21(III5,FF,F)
	SUM=SUM +(-1.D0)**J*RRDSQRT/RRFUN
 ENDDO


 I=IDINT(DABS(A1-A2-A6)+.1D0)

 ! пегскэрюр пюяверю 
 C3JS_W3JA=(-1.D0)**I*SUM
 
 return
end function C3JS_W3JA


! 3J-SYMBOL CALCULATION SUB-PROGRAM
! BY THE METHOD OF SIMPLIFICATION
! I1, I2, I3, I4, I5, I6-SYMBOL PARAMETERS
! AUXILIARY PARAMETERS FOR CALCULATION OF FF, FF (50)
real(8) function C3JS_W3JAGG(I1,I2,I3,I4,I5,I6)
 implicit none
 integer::I1,I2,I3,I4,I5,I6
 real(8)::RR,RRRR
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::MIN,MAX,I,J,III,IMAXI,IMAXJ,IIMAX,IIX
 integer::III1,III2,III3,III4,III5,III6,III7,III8,III9,III10
 integer,dimension(10)::IndexI,IndexJ
 real(8)::A1,A2,A3,A4,A5,A6,SUM,D1,D2,RRFUN,RRDSQRT,RCOFF
 real(8),dimension(10)::RndexI,RndexJ



 ! оюпюлерпш 
 A1=DBLE(I1) 
 A2=DBLE(I2)
 A3=DBLE(I3)
 A4=DBLE(I4)
 A5=DBLE(I5)
 A6=DBLE(I6)
      
	
 C3JS_W3JAGG=0.D0

 IF(IDINT(A4+A5+A6+.1D0).NE.0) THEN
    return
 ENDIF
 IF(DABS(A1+A2+A3-DINT(A1+A2+A3+.1D0)).GT..1D0) THEN
    return
 ENDIF
 D1=DMAX1(A2-A3-A4,A1-A3+A5,0.D0)
 D2=DMIN1(A1+A2-A3,A1-A4,A2+A5)
 MIN=IDINT(D1)
 MAX=IDINT(D2)
  
 IF(MIN.GT.MAX) THEN 
    return
 ENDIF

 
 IndexI(1)=IDINT(A1+A4+.1D0)
 IndexI(2)=IDINT(A1-A4+.1D0)
 IndexI(3)=IDINT(A2+A5+.1D0)
 IndexI(4)=IDINT(A2-A5+.1D0)
 IndexI(5)=IDINT(A3+A6+.1D0)
 IndexI(6)=IDINT(A3-A6+.1D0)
 IndexI(7)=IDINT(A1+A2-A3+.1D0)
 IndexI(8)=IDINT(A1-A2+A3+.1D0)
 IndexI(9)=IDINT(A2+A3-A1+.1D0)
 IndexI(10)=IDINT(A1+A2+A3+1.1D0)
  

 IndexJ(1)=IDINT(A3-A2+A4+.1D0+MIN)
 IndexJ(2)=IDINT(A3-A1-A5+.1D0+MIN) 
 IndexJ(3)=IDINT(A1+A2-A3+.1D0-MIN)
 IndexJ(4)=IDINT(A1-A4+.1D0-MIN)
 IndexJ(5)=IDINT(A2+A5+.1D0-MIN)
 IndexJ(6)=MIN
 ! нопедекъел люйяхлюкэмне гмювемхе
 IMAXI=IndexI(1)
 DO III=1,10
    IF(IMAXI.LT.IndexI(III)) THEN
       IMAXI=IndexI(III)
	ENDIF
 ENDDO
 IMAXJ=IndexJ(1)
 DO III=1,6
    IF(IMAXJ.LT.IndexJ(III)) THEN
       IMAXJ=IndexJ(III)
	ENDIF
 ENDDO

 IF(IMAXI.GE.IMAXJ) THEN
     IIMAX=IMAXI 
    ELSE
     IIMAX=IMAXJ
 ENDIF
 
 RCOFF=DSQRT(0.1D0)
 ! нясыеярбкъел пюявер оепбнцн якюцюелнцн
 DO III=1,IIMAX
    ! жхйк он йнщттхжхемрюл бундъыхл б ярпсйрспс оепбнцн якюцюелнцн
	DO IIX=1,10
	   ! опнбепъел нркхвем нр мскъ хмдейя
	   IF(IndexI(IIX).NE.0) THEN
           RndexI(IIX)=FLOAT(IndexI(IIX))*0.1D0
		   ! слемэьюел хмдейя мю едемхжс
           IndexI(IIX)=IndexI(IIX)-1
          ELSE
		   RndexI(IIX)=1.D0 	   
	   ENDIF    
	ENDDO 
	DO IIX=1,6
	   ! опнбепъел нркхвем нр мскъ хмдейя
	   IF(IndexJ(IIX).NE.0) THEN
           RndexJ(IIX)=FLOAT(IndexJ(IIX))*0.1D0
		   ! слемэьюел хмдейя мю едемхжс
           IndexJ(IIX)=IndexJ(IIX)-1
          ELSE
		   RndexJ(IIX)=1.D0 	   
	   ENDIF    
	ENDDO 
    RCOFF=RCOFF*DSQRT(RndexI(1))*DSQRT(RndexI(2))*DSQRT(RndexI(3))*DSQRT(RndexI(4))*DSQRT(RndexI(5))*DSQRT(RndexI(6))*DSQRT(RndexI(7))*DSQRT(RndexI(8))*DSQRT(RndexI(9))/(DSQRT(RndexI(10))*RndexJ(1)*RndexJ(2)*RndexJ(3)*RndexJ(4)*RndexJ(5)*RndexJ(6))
 ENDDO
 !WRITE(150,*)
 !write(150,*) 'DDDIF',RCOFF
 
 SUM=(-1.D0)**MIN*RCOFF
 ! опнбепъел вхякн якюцюелшу  
 IF(MIN.NE.MAX) THEN
    DO J=MIN+1,MAX
       III1=IDINT(A3-A2+A4+.1D0+J)
       III2=IDINT(A3-A1-A5+.1D0+J) 
       III3=IDINT(A1+A2-A3+.1D0-J)+1
       III4=IDINT(A1-A4+.1D0-J)+1
       III5=IDINT(A2+A5+.1D0-J)+1
	   RCOFF=RCOFF*FLOAT(III3)*FLOAT(III4)*FLOAT(III5)/(FLOAT(J)*FLOAT(III1)*FLOAT(III2))
	   !write(150,*) 'DDDIFXX',RCOFF
	   SUM=SUM+(-1.D0)**J*RCOFF 
	ENDDO
 ENDIF

 I=IDINT(DABS(A1-A2-A6)+.1D0)

 ! пегскэрюр пюяверю 
 C3JS_W3JAGG=(-1.D0)**I*SUM
 
 return
end function C3JS_W3JAGG




! бяонлнцюрекэмюъ опнжедспю 
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

! CALCULATION OF THE REDUCED MATRIX ELEMENT OF THE SPHERICAL HARMONY
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL MOMENT OF FIRST HARMONIC
! L2-ORBITAL MOMENT OF THE SECOND HARMONIC
! K-ORBITAL MOMENT OF SPHERICAL HARMONIC
real(8) function C3JS_CKKGG(L1,L2,K)
 implicit none
 integer::L1,L2,K
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::M,N,I,J,KK,MM,III,IMAXI,IIX
 real(8)::A,B,C,RCOFF,D,DD
 integer,dimension(10)::IndexID
 real(8),dimension(10)::RndexID 
 34 FORMAT(1X,'CKKNEW',F18.10)

 C3JS_CKKGG=0.D0

 M=L1+L2+K+2
 N=M/2

 IF((M-N-N).EQ.0) THEN 
      A=DBLE(L1)
      B=DBLE(L2)
      C=DBLE(K)
      ! пюявер декэрю 
	  D=A+B+C+1.1D0
      I=IDINT(D-C-C)
      J=IDINT(D-B-B)
      KK=IDINT(D-A-A)
    
	  IF(I*J*KK.GT.0) THEN 
           MM=IDINT(D)
           DD=MM+.1D0-D
         ELSE 
           return
      ENDIF
      ! опнбепъел бшонкмемхъ сякнбхъ
	  IF(IDINT(100D0*DD).LT.0) THEN 
         return
	  ENDIF
      ! нясыеярбкъел пюявер опхбедеммнцн 3J-яхлбнкю
  	  IndexID(1)=I-1
	  IndexID(2)=J-1
	  IndexID(3)=KK-1
	  IndexID(4)=MM
	  IndexID(5)=N-1
	  IndexID(6)=N-L1-1
	  IndexID(7)=N-L2-1
      IndexID(8)=N-K-1

	  
	  ! нопедекъел люйяхлюкэмне гмювемхе
      IMAXI=IndexID(1)
      DO III=1,8
         IF(IMAXI.LT.IndexID(III)) THEN
            IMAXI=IndexID(III)
         ENDIF
      ENDDO
    	
      ! нясыеярбкъел пюявер 
      RCOFF=DSQRT(DBLE((L1+L1+1)*(L2+L2+1))*0.1D0)
	  DO III=1,IMAXI
         ! жхйк он йнщттхжхемрюл бундъыхл б ярпсйрспс оепбнцн якюцюелнцн
	     DO IIX=1,8
	        ! опнбепъел нркхвем нр мскъ хмдейя
	        IF(IndexID(IIX).NE.0) THEN
               RndexID(IIX)=FLOAT(IndexID(IIX))*0.1D0
		       ! слемэьюел хмдейя мю едемхжс
               IndexID(IIX)=IndexID(IIX)-1
              ELSE
		       RndexID(IIX)=1.D0 	   
	        ENDIF    
	     ENDDO 
         RCOFF=RCOFF*DSQRT(RndexID(1))*DSQRT(RndexID(2))*DSQRT(RndexID(3))*RndexID(5)/(DSQRT(RndexID(4))*RndexID(6)*RndexID(7)*RndexID(8))
	  ENDDO
   
      C3JS_CKKGG=RCOFF
	  !WRITE(150,34) C3JS_CKKGG
   ELSE
      return
 ENDIF

 
 return
end function C3JS_CKKGG
      


! CALCULATION OF THE REDUCED MATRIX ELEMENT OF THE SPHERICAL HARMONY
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL MOMENT OF FIRST HARMONIC
! L2-ORBITAL MOMENT OF THE SECOND HARMONIC
! K-ORBITAL MOMENT OF SPHERICAL HARMONIC
real(8) function C3JS_CKK(L1,L2,K,FF,F)
 implicit none
 integer::L1,L2,K
 real(8)::FF
 real(8),dimension(:)::F
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::M,N
 real(8)::A,B,C,DGH,FHF
 35 FORMAT(1X,'CKKOLD',F18.10)

 C3JS_CKK=0.D0
 M=L1+L2+K+2
 N=M/2
 IF((M-N-N).EQ.0) THEN 
    A=DBLE(L1)
    B=DBLE(L2)
    C=DBLE(K)
    B=C3JS_DELTA1(A,B,C,FF,F)
	!WRITE(150,*) 'BOLD',B
   ELSE
    return
 ENDIF

 IF(B.GT.0.D0) THEN 
   DGH=B*DSQRT(DBLE((L1+L1+1)*(L2+L2+1))*0.1D0)*C3JS_F31(N,FF,F)
   FHF=C3JS_F31(N-L1,FF,F)*C3JS_F31(N-L2,FF,F)*C3JS_F31(N-K,FF,F)
   C3JS_CKK=DGH/FHF
 ENDIF 
 !WRITE(150,35) C3JS_CKK
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
  !WRITE(150,*) 'IJK',I,J,K   
  IF(I*J*K.GT.0) THEN 
     M=IDINT(D)
     DD=M+.1D0-D
    ELSE 
     return
  ENDIF
  
 
  !WRITE(150,*) 'DD',IDINT(100D0*DD)
  IF(IDINT(100D0*DD).GE.0) THEN 
     RRFG=C3JS_F31(I,FF,F)*C3JS_F31(J,FF,F)*C3JS_F31(K,FF,F)
     !WRITE(150,*) 'RRFG',RRFG
     !WRITE(150,*) 'DADA',C3JS_F31(M+1,FF,F)
	 C3JS_DELTA1=DSQRT(RRFG/C3JS_F31(M+1,FF,F))
     !WRITE(150,*) 'CCCF',C3JS_DELTA1
  ENDIF

  return
end function C3JS_DELTA1

! бяонлнцюрекэмюъ опнжедспю 
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


  ! SUB-PROGRAM FOR CALCULATING THE MATRIX ELEMENT OF THE SPHERICAL HARMONY
  ! DIRECT INTEGRATION OF SPHERICAL HARMONICS IS CARRIED OUT
  ! L1-ORBITAL MOMENT OF FIRST FUNCTION
  ! M1-PROJECTION OF ORBITAL FIRST FUNCTION
  ! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
  ! M2-PROJECTION OF THE ORBITAL MOMENT OF THE SECOND FUNCTION
  ! K-ORBITAL MOMENT OF SPHERICAL HARMONY
  ! THE PROJECTION OF SPHERICAL HARMONICS IS DEFINED BY THE EQUALITY q = M1-M2
  ! OTHER PROJECTIONS LEAD TO ZEROING OF THE MATRIX ELEMENT
 real(8) function C3JS_CkqISH(L1,M1,L2,M2,K)
  implicit none
  integer::L1,M1,L2,M2,K
  !!!!!!!!!!!!!!!!!!!!!!!
  integer::II,INDEX
  
  
  C3JS_CkqISH=0.D0
  II= L1+(L1+K+L2)/2.D0+1.D-1
  ! опнбепъел нркхвхе нр мскъ йнщттхжхемрю бхцмепю  (3J-яхлбнкю)
  INDEX=C3JS_W3JAISH(L1,K,L2,-M1,M1-M2,M2)
  ! опнбепъел нркхвхе нр мскъ
  IF(INDEX.EQ.0) THEN
     return
  ENDIF
  ! опнбепъел нркхвем нр мскъ опхбедеммши люрпхвмши щкелемр ноепюрнпю ятепхвеяйни цюплнмхйх 
  INDEX=C3JS_CKKISH(L1,L2,K)
  ! опнбепъел нркхвхе нр мскъ
  IF(INDEX.EQ.0) THEN
     return
  ENDIF
  ! пегскэрюр
  C3JS_CkqISH=CalculationCkq(L1,M1,L2,M2,K)
  return
 end function C3JS_CkqISH




 ! SUB-PROGRAM FOR CALCULATING THE MATRIX ELEMENT OF THE SPHERICAL HARMONY
  ! DONE BY DIRECT INTEGRATION
  ! L1-ORBITAL MOMENT
  ! M1-PROJECTION OF ORBITAL MOMENT
  ! L2-ORBITAL MOMENT
  ! M2-PROJECTION OF THE ORBITAL MOMENT
  ! HARMONIC K-NUMBER
 real(8) function CalculationCkq(L1,M1,L2,M2,K)
  implicit none
  integer::L1,M1,L2,M2,K
  !!!!!!!!!!!!!!!!!!!!!!!
  integer::ML1,ML2,ML3,IIP,IMAX,ierr,IXX,IIY1,IIY2,IIY3
  real(8)::Rcoff,RME
  integer,dimension(3)::NPP
  real(8),allocatable,dimension(:)::TL1,TL2,TL3,XXA,YYA
  real(8),allocatable,dimension(:,:)::RcoffSin

  ML1=IABS(M1)
  ML2=IABS(M2)
  ML3=IABS(M1-M2)

  CalculationCkq=0.D0
  ! опнбепъел яннрберярбхе опнейжхх х нпахрюкэмнцн лнлемрю
  IF(ML1.GT.L1) THEN
     return
  ENDIF
  IF(ML2.GT.L2) THEN
     return
  ENDIF
  IF(ML3.GT.K) THEN
     return
  ENDIF
   
  Rcoff=1.D0
  IF(M1.LT.0) THEN
     Rcoff=Rcoff*(-1.D0)**ML1     
  ENDIF
  IF(M2.LT.0) THEN
     Rcoff=Rcoff*(-1.D0)**ML2    
  ENDIF
  IF((M1-M2).LT.0) THEN
     Rcoff=Rcoff*(-1.D0)**ML3   
  ENDIF

  ! нопедекъел люйяхлюкэмне вхякн рнвей
  NPP(1)=L1-ML1+1
  NPP(2)=L2-ML2+1
  NPP(3)=K-ML3+1
  
  IMAX=NPP(1)
  DO IIP=1,3
     IF(IMAX.LT.NPP(IIP)) THEN
        IMAX=NPP(IIP)
	 ENDIF
  ENDDO

  ! нясыеярбкъел онксвемхъ ярпсйрсп рщрю-тсмйжхи
  allocate(TL1(NPP(1)),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "TL1" IS NOT SELECTED'
     stop 
  endif
  allocate(TL2(NPP(2)),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "TL2" IS NOT SELECTED'
     stop 
  endif
  allocate(TL3(NPP(3)),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "TL3" IS NOT SELECTED'
     stop 
  endif
  allocate(XXA(IMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "XXA" IS NOT SELECTED'
     stop 
  endif
  allocate(YYA(IMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "YYA" IS NOT SELECTED'
     stop 
  endif  
  allocate(RcoffSin(2,ML1+ML2+ML3+2),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'MEMORY ON THE FILE "RcoffSin" IS NOT SELECTED'
     stop 
  endif  

  ! нопедекъел гмювемхе ятепхвеяйху цюплнмхй дкъ ху юопнйяхлюжхх
  call VALUE_TETA_FUNCTION(L1,ML1,NPP(1),XXA,YYA)
  ! нопедекъел йнщттхжхемрш юопнйяхлхпнбюммни тсмйжхх 
  call COEFFICIENT_TETA_FUNCTION_APRO(ML1,NPP(1),XXA,YYA,TL1)
  ! нопедекъел гмювемхе ятепхвеяйху цюплнмхй дкъ ху юопнйяхлюжхх
  call VALUE_TETA_FUNCTION(L2,ML2,NPP(2),XXA,YYA)
  ! нопедекъел йнщттхжхемрш юопнйяхлхпнбюммни тсмйжхх 
  call COEFFICIENT_TETA_FUNCTION_APRO(ML2,NPP(2),XXA,YYA,TL2) 
  ! нопедекъел гмювемхе ятепхвеяйху цюплнмхй дкъ ху юопнйяхлюжхх
  call VALUE_TETA_FUNCTION(K,ML3,NPP(3),XXA,YYA)
  ! нопедекъел йнщттхжхемрш юопнйяхлхпнбюммни тсмйжхх 
  call COEFFICIENT_TETA_FUNCTION_APRO(ML3,NPP(3),XXA,YYA,TL3)
  ! нопедекъел йнщттхжхемрш опедярюбкемхъ SIN б  ML1+ML2+ML3+1 яреоемх
  call COEFFICIENT_TRANSFORMATION_POLINOM_SINN(ML1+ML2+ML3+1,RcoffSin)

  RME=0.D0
  ! нясыеярбкъел пюявер хмрецпюкю
  ! жхйк он яреоемъл яхмсяю 
  DO IXX=0,ML1+ML2+ML3+1
     ! жхйк он якюцюелшл цюплнмхй
	 DO IIY1=0,L1-ML1
	    DO IIY2=0,L2-ML2
		   DO IIY3=0,K-ML3
		      RME=RME+TL1(IIY1+1)*TL2(IIY2+1)*TL3(IIY3+1)*(RcoffSin(1,IXX+1)*RintSIN(IXX,IIY1,IIY2,IIY3)+RcoffSin(2,IXX+1)*RintCOS(IXX,IIY1,IIY2,IIY3))
           ENDDO
		ENDDO
	 ENDDO
  ENDDO
  
  ! гюохяшбюел пюявхрюммне гмювемхе люрпхвмнцн щкелемрю ятепхвеяйни цюплнмхйх
  CalculationCkq=DSQRT(2.D0/FLOAT(2*K+1))*Rcoff*RME


  ! сдюкемхе люяяхбнб хг оълърх  
  deallocate(TL1,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "TL1" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(TL2,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "TL2" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(TL3,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "TL3" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(XXA,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "XXA" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(YYA,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "YYA" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  deallocate(RcoffSin,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CalculationCkq'
     write(*,*) 'THE FILE "RcoffSin" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  
  return
 end function CalculationCkq


 ! ондопнцпюллю пюяверю хмрецпюкю 
 ! нясыеярбкъеряъ пюявер хмрецпюкю я ондхмрецпюкэмни тсмйжхи
 ! SIN(IXX*X)*COS(IIY1*X)*COS(IIY2*X)*COS(IIY3*X)
 real(8) function RintSIN(IXX,IIY1,IIY2,IIY3)
  implicit none
  integer::IXX,IIY1,IIY2,IIY3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::INDEXSD,ISX11,ISX12,ISX21,ISX22,ISX31,ISX32,ISX41,ISX42
  
  
  IF(IXX.EQ.0) THEN
     RintSIN=0.D0
	 return   
  ENDIF 
  
  INDEXSD=0
  IF(IIY1.EQ.0) THEN
     IF(IIY2.EQ.0) THEN
	   IF(IIY3.EQ.0) THEN
	      INDEXSD=1
	   ENDIF 
	 ENDIF 
  ENDIF
 
  ! йняхмсяш пюбмш едхмхже
  IF(INDEXSD.EQ.1) THEN
     IF((1-(-1)**IXX).NE.0) THEN
	     RintSIN=2.D0/FLOAT(IXX)
       ELSE
	     RintSIN=0.D0
	 ENDIF
  ENDIF

  IF(INDEXSD.EQ.0) THEN
     ISX11=IXX+IABS(IIY1+IIY2-IIY3)
	 ISX21=IXX+IABS(-IIY1+IIY2+IIY3)
	 ISX31=IXX+IABS(IIY1-IIY2+IIY3)
	 ISX41=IXX+IABS(IIY1+IIY2+IIY3)
     ISX12=IXX-IABS(IIY1+IIY2-IIY3)
	 ISX22=IXX-IABS(-IIY1+IIY2+IIY3)
	 ISX32=IXX-IABS(IIY1-IIY2+IIY3)
	 ISX42=IXX-IABS(IIY1+IIY2+IIY3)
     
     RintSIN=(RINTRSIN(ISX11)+RINTRSIN(ISX21)+RINTRSIN(ISX31)+RINTRSIN(ISX41)+RINTRSIN(ISX12)+RINTRSIN(ISX22)+RINTRSIN(ISX32)+RINTRSIN(ISX42))/8.D0 
  ENDIF




  return
 end function

 real(8) function RINTRSIN(Q)
   implicit none
   integer::Q

   IF(Q.EQ.0) THEN
      RINTRSIN=0.D0
   ENDIF
   IF(Q.NE.0) THEN
      IF((1-(-1)**Q).NE.0) THEN
	     RINTRSIN=2.D0/FLOAT(Q)
        ELSE
	     RINTRSIN=0.D0
	  ENDIF
   ENDIF
   return
 end function


  ! INTEGRAL CALCULATION SUBPROGRAM
  ! INTEGRAL CALCULATION WITH SUB-INTEGRAL FUNCTION IS PERFORMED
  ! COS (IXX * X) * COS (IIY1 * X) * COS (IIY2 * X) * COS (IIY3 * X)
 real(8) function RintCOS(IXX,IIY1,IIY2,IIY3)
  implicit none
  integer::IXX,IIY1,IIY2,IIY3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::INDEXSD,ISX11,ISX12,ISX21,ISX22,ISX31,ISX32,ISX41,ISX42
  
  
   
  INDEXSD=0
  IF(IXX.EQ.0) THEN
     IF(IIY1.EQ.0) THEN
        IF(IIY2.EQ.0) THEN
	       IF(IIY3.EQ.0) THEN
	          INDEXSD=1
	       ENDIF 
	    ENDIF 
     ENDIF
  ENDIF
  ! йняхмсяш пюбмш едхмхже
  IF(INDEXSD.EQ.1) THEN
     RintCOS=3.14159265358979D0
  ENDIF

  IF(INDEXSD.EQ.0) THEN
     ISX11=IXX+IABS(IIY1+IIY2-IIY3)
	 ISX21=IXX+IABS(-IIY1+IIY2+IIY3)
	 ISX31=IXX+IABS(IIY1-IIY2+IIY3)
	 ISX41=IXX+IABS(IIY1+IIY2+IIY3)
     ISX12=IXX-IABS(IIY1+IIY2-IIY3)
	 ISX22=IXX-IABS(-IIY1+IIY2+IIY3)
	 ISX32=IXX-IABS(IIY1-IIY2+IIY3)
	 ISX42=IXX-IABS(IIY1+IIY2+IIY3)
     
     RintCOS=(RINTRCOS(ISX11)+RINTRCOS(ISX21)+RINTRCOS(ISX31)+RINTRCOS(ISX41)+RINTRCOS(ISX12)+RINTRCOS(ISX22)+RINTRCOS(ISX32)+RINTRCOS(ISX42))/8.D0 
  ENDIF




  return
 end function

 real(8) function RINTRCOS(Q)
   implicit none
   integer::Q

   IF(Q.EQ.0) THEN
       RINTRCOS=3.14159265358979D0
      ELSE
	   RINTRCOS=0.D0 
   ENDIF
   
   return
 end function




  ! SUBPROGRAM FOR OBTAINING THE POLYNOMY TRANSFORMATION COEFFICIENTS
  ! sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! N-MAXIMUM POLYNOMA DEGREE
  ! Rcoff (2, N + 1) COEFFICIENT ARRAY REPRESENTATION sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
  ! Rcoff (1, N + 1) COEFFICIENTS AT SIN (PX)
  ! Rcoff (2, N + 1) COEFFICIENTS AT COS (PX)
 subroutine COEFFICIENT_TRANSFORMATION_POLINOM_SINN(N,Rcoff)
	implicit none

    integer::N,IYU,JYU,ierr
    real(8),dimension(:,:)::Rcoff
	real(8),allocatable,dimension(:,:,:)::RSP

	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RSP(2,N+2,N+2),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_SINN'
	write(*,*) 'MEMORY ON THE FILE "RSP" IS NOT SELECTED'
	stop 
	endif

    ! гюмскъел оепед пюявернл
	Rcoff=0.D0
    RSP=0.D0
      
	! гюохяшбюел йнщттхжхемрш дкъ яксвюъ sin(x)^0
    ! SIN
	RSP(1,1,1)=1.D0
    ! COS
	RSP(2,1,1)=0.D0
	! гюохяшбюел йнщттхжхемрш дкъ яксвюъ sin(x)^1
    ! SIN
	RSP(1,2,1)=0.D0
    RSP(1,2,2)=1.D0
    ! COS
	RSP(2,2,1)=0.D0
    RSP(2,2,2)=0.D0


	
	! жхйк он яреоемх SIN(X)^IYU
      DO IYU=2,N
         ! жхйк он йнщттхжхемрюл мнбнцн онкхмнлю (жхйк он яреоемъл) SIN
	     RSP(1,IYU+1,2)=RSP(1,IYU,1)+RSP(2,IYU,1)
	     DO JYU=1,IYU
	        RSP(1,IYU+1,JYU+2)=RSP(1,IYU+1,JYU+2)+RSP(2,IYU,JYU+1)*0.5D0
            RSP(1,IYU+1,JYU)=RSP(1,IYU+1,JYU)-RSP(2,IYU,JYU+1)*0.5D0
	     ENDDO
	     ! жхйк он йнщттхжхемрюл мнбнцн онкхмнлю (жхйк он яреоемъл) COS
	     DO JYU=1,IYU
	        RSP(2,IYU+1,JYU+2)=RSP(2,IYU+1,JYU+2)-RSP(1,IYU,JYU+1)*0.5D0
            RSP(2,IYU+1,JYU)=RSP(2,IYU+1,JYU)+RSP(1,IYU,JYU+1)*0.5D0
	     ENDDO    
	  ENDDO

	! гюохяшбюел гмювемхъ дкъ онксвеммнцн онкхмнлю
	DO IYU=1,N+1 
       ! SIN-йнщттхжхемрш
	   Rcoff(1,IYU)=RSP(1,N+1,IYU)
	   ! COS-йнщттхжхемрш
       Rcoff(2,IYU)=RSP(2,N+1,IYU)
    ENDDO

	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(RSP,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_SIN'
      write(*,*) 'THE FILE "RSP" IS NOT REMOVED FROM MEMORY'
	stop 
	endif


      
	return
  end subroutine COEFFICIENT_TRANSFORMATION_POLINOM_SINN



  ! SUB-PROGRAM FOR OBTAINING THE VALUE OF THETA FUNCTION (RECURRENT METHOD)
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! L-ORBITAL MOMENT FUNCTION
  ! M-PROJECTION OF THE ORBITAL MOMENTUM (THE PROGRAM CONSIDERS THE CASE (M> = 0)
  ! N-NUMBER OF POINTS WHICH WE FIND VALUES (IN THE INTERVAL (0, PI))
  ! X (N) -MASSIF OF ARGUMENT VALUES
  ! Y (N) -ARRAY OF FUNCTION VALUES
 subroutine VALUE_TETA_FUNCTION(L,M,N,X,Y)
     implicit none
     integer::L,M,N
     real(8),dimension(:)::X,Y
     !!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IAA1,IAA2,NNNMP
     real(8)::RHR,Rcoff1,Rcoff2,Rcoff3,RSDXF,RSDXX

     ! гюмскъел оепед пюявернл
     X=0.D0
	 Y=0.D0

     IF(N.EQ.1) THEN
        NNNMP=N
	   ELSE
        NNNMP=N-1
	 ENDIF
     ! тнплхпсел люяяхб гмювемхи юпцслемрю
     RHR=3.14159265358979D0/float(N+1)         
     DO IAA1=2,N+1   
        X(IAA1-1)=RHR*float(IAA1-1)
	 ENDDO
	  
	 ! жхйк он рнвйюл
	 DO IAA1=1,N
        ! якюцюелне Teta M,M
	    RSDXF=DSQRT(FACTORIAL(2*IABS(M)))/FACTORIAL(IABS(M))
	    RSDXX=(DSIN(X(IAA1)*0.5D0)*DCOS(X(IAA1)*0.5D0))**IABS(M)
        RSDXX=RSDXX*RSDXF/DSQRT(2.D0)
	    Rcoff1=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+1))*RSDXX
        !  якюцюелне Teta M+1,M 
	    RSDXF=DSQRT(FACTORIAL(2*IABS(M)+1))/FACTORIAL(IABS(M))
	    RSDXX=(DSIN(X(IAA1)*0.5D0)*DCOS(X(IAA1)*0.5D0))**IABS(M)
        RSDXX=RSDXX*RSDXF*DCOS(X(IAA1))/DSQRT(2.D0)
	    Rcoff2=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+3))*RSDXX
	    IF(IABS(M).EQ.L) Rcoff3=Rcoff1
        IF(IABS(M).EQ.L-1) Rcoff3=Rcoff2
        ! жхйк он лнлемрюл
	    DO IAA2=IABS(M),L-2
           RSDXF=SQRT(float(2*IAA2+3)*float(2*IAA2+5))
           RSDXF=RSDXF/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
         
		   RSDXX=SQRT(float(2*IAA2+5))
           RSDXX=RSDXX*SQRT(float(IAA2+1-IABS(M))*float(IAA2+1+IABS(M)))
           RSDXX=RSDXX/SQRT(float(2*IAA2+1))
           RSDXX=RSDXX/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
           ! пюявер якедсчыецн щкелемрю
		   Rcoff3=RSDXF*DCOS(X(IAA1))*Rcoff2-RSDXX*Rcoff1
		   ! ондцнрнбйю й пюяверс якедсчыецн щкелемрю
           Rcoff1=Rcoff2
	       Rcoff2=Rcoff3
        ENDDO
        ! гюохяшбюел гмювемхе тсмйжхх
	    Y(IAA1)=Rcoff3
     ENDDO


   return  
 end subroutine VALUE_TETA_FUNCTION

   ! SUBPROGRAM FOR APPROXIMATING THETA FUNCTIONS BY FUNCTIONS OF THE FORM sin (x) ^ M * SUM (Qi * cos (i * x))
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! M-PROJECTION OF ORBITAL MOMENT
   ! N-NUMBER OF DOTS
   ! X (N) -MASSIF OF VALUES OF THE TETE FUNCTION ARGUMENT
   ! Y (N) -ARRAY OF VALUES OF THEETA FUNCTION
   ! Qcoff (N) -MASSIVE OF VALUES OF THE APPROXIMATION COEFFICIENTS
  subroutine COEFFICIENT_TETA_FUNCTION_APRO(M,N,X,Y,Qcoff) 
   implicit none
   integer::N,M,ierr 
   real(8),dimension(:)::X,Y,Qcoff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::JZCVX,IZCVX
   real(8)::RXCVB
   real(8),allocatable,dimension(:,:)::AXXZ

   !бшдекъел оюлърэ дкъ люяяхбнб
   allocate(AXXZ(N,N),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION_APRO'
   	  write(*,*) 'MEMORY ON THE FILE "AXXZ" IS NOT SELECTED'
	  stop 
   endif


   ! гюмскъел оепед пюявернл
   AXXZ=0.D0


   ! гюонкмъел люяяхб
   DO JZCVX=1,N
      RXCVB=X(JZCVX)
      DO IZCVX=0,N-1
         AXXZ(JZCVX,IZCVX+1)=DSIN(RXCVB)**IABS(M)*DCOS(FLOAT(IZCVX)*RXCVB)
      ENDDO
   ENDDO
   ! ондопнцпюллю мюунфдемхъ йнпмеи яхярелш кхмеимшу спюбмемхи
   ! лернд цюсяяю я опхлемемхел яуелш вюярхвмнцн бшанпю
   call SYSTEM_LINEAR_EQUATIONS(N,AXXZ,Y,Qcoff)


      
	
   ! сдюкемхе люяяхбнб хг оълърх   
   deallocate(AXXZ,stat=ierr)
   if(ierr/=0) then
  	  write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION_APRO'
      write(*,*) 'THE FILE "AXXZ" IS NOT REMOVED FROM MEMORY'
	  stop 
   endif
   return
  end subroutine COEFFICIENT_TETA_FUNCTION_APRO

! SUBPROGRAM FOR SOLVING A SYSTEM OF LINEAR EQUATIONS A * Xs = Y
! THE GAUSS METHOD USING A PARTIAL SELECTION SCHEME
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! N-NUMBER OF UNKNOWN
! A-MATRIX OF A SYSTEM OF LINEAR EQUATIONS
! Y-ARRAY FUNCTION VALUE
! Xs-ARRAY OF ROOT SYSTEM VALUES
 subroutine SYSTEM_LINEAR_EQUATIONS(N,A,Y,Xs) 
	implicit none
      integer::N,ierr 
	real(8),dimension(:)::Y,Xs
      real(8),dimension(:,:)::A
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer,allocatable,dimension(:)::ipvt
      real(8),allocatable,dimension(:,:)::fac
      
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(ipvt(N),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_SYSTEM_LINEAR_EQUATIONS'
	write(*,*) 'MEMORY ON THE FILE "ipvt" IS NOT SELECTED'
	stop 
	endif
	allocate(fac(N,N),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_SYSTEM_LINEAR_EQUATIONS'
	write(*,*) 'MEMORY ON THE FILE "fac" IS NOT SELECTED'
	stop 
	endif




      ! ОНДОПНЦПЮЛЛЮ LU-пюгкнфемхъ люрпхжш A  
	call DLufac2Z(N,A,fac,ipvt)
      ! пЕЬЮЕЛ КХМЕИМСЧ ЯХЯРЕЛС AXs=Y
	call useLU2Z(N,fac,ipvt,Y,Xs)  

	
	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(ipvt,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_SYSTEM_LINEAR_EQUATIONS'
      write(*,*) 'THE FILE "ipvt" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(fac,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_SYSTEM_LINEAR_EQUATIONS'
      write(*,*) 'THE FILE "fac" IS NOT REMOVED FROM MEMORY'
	stop 
	endif

	return
end subroutine SYSTEM_LINEAR_EQUATIONS


      ! бяонлнцюрекэмюъ ондпнцпюллю
      subroutine useLU2Z(N,fac,ipvt,b,x)
	implicit none
      integer::N,k,i 
	real(8)::s,hold
	integer,dimension(:)::ipvt
	real(8),dimension(:)::b,x
      real(8),dimension(:,:)::fac
	  
      do i=1,N
         if(ipvt(i)/=i) then
            hold=b(ipvt(i))
	      b(ipvt(i))=b(i)
	      b(i)=hold
	   endif
	enddo

	do k=2,N
        b(k)=b(k)+dot_product(fac(k,1:k-1),b(1:k-1))
	enddo

	x(N)=b(N)/fac(N,N)
	do k=N-1,1,-1
         s=sum(fac(k,k+1:n)*x(k+1:n))
         x(k)=(b(k)-s)/fac(k,k)
	enddo
      
	
	return
      end subroutine useLU2Z
  
  subroutine DLufac2Z(n,a,fac,ipvt)
	implicit none
	integer::ierr,n,i,k,j
	real(8)::rhold  
	integer,dimension(:)::ipvt
	real(8),dimension(:,:)::a,fac
      real(8),allocatable,dimension(:)::p
      ! бшдекъел оълърэ онд люяяхбш
	allocate(p(n),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma MEMORY ON THE FILE "p" IS NOT SELECTED'
	stop 
	endif  
	
	fac=a
	ipvt=(/(i,i=1,n)/) 
	do k=1,n-1
	ipvt(k:k)=k+maxloc(dabs(fac(k:n,k)))-1
	
	if(ipvt(k)/=k)then
	p(k:n)=fac(k,k:n) 
      fac(k,k:n)=fac(ipvt(k),k:n)
	fac(ipvt(k),k:n)=p(k:n)

	do j=1,k-1
      rhold=fac(ipvt(k),j)
	fac(ipvt(k),j)=fac(k,j)
	fac(k,j)=rhold
	enddo
      endif
      p(k+1:n)=fac(k+1:n,k)/fac(k,k)
	do i=k+1,n
	fac(i,k+1:n)=fac(i,k+1:n)-fac(k,k+1:n)*p(i)
	enddo
	fac(k+1:n,k)=-p(k+1:n)
	enddo

	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(p,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma THE FILE "p" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	return
	end subroutine DLufac2Z

	   ! ондопнцпюллю пюяверю тюйрнпхюкю
	! нохяюмхе оюпюлерпнб
      ! NOPI-жекне онкнфхрекэмне вхякн 
	real(8) function FACTORIAL(NOPI)
      implicit none

	integer:: NOPI,IBKL
      real(8):: PROSS
      
	IF(NOPI.LT.0) THEN
        WRITE(*,*) 'ERROR FACTORIAL (N<0)'
        READ(*,*)
	  STOP
	ENDIF

      IF(NOPI.EQ.0) THEN
      FACTORIAL=1.D0
      return
	ENDIF
      
	PROSS=1.D0
	do IBKL=1,NOPI
         PROSS=PROSS*float(IBKL)
	enddo

	FACTORIAL=PROSS

   
      return
      end function FACTORIAL
	   

! CALCULATION OF THE REDUCED MATRIX ELEMENT OF THE SPHERICAL HARMONY
! WE CARRY OUT A CHECK FOR DIFFERENCE FROM ZERO
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL MOMENT OF FIRST HARMONIC
! L2-ORBITAL MOMENT OF THE SECOND HARMONIC
! K-ORBITAL MOMENT OF SPHERICAL HARMONIC
integer function C3JS_CKKISH(L1,L2,K)
 implicit none
 integer::L1,L2,K
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::M,N,I,J,KK,MM,III,IMAXI,IIX
 real(8)::A,B,C,RCOFF,D,DD

 C3JS_CKKISH=0

 M=L1+L2+K+2
 N=M/2

 IF((M-N-N).EQ.0) THEN 
      A=DBLE(L1)
      B=DBLE(L2)
      C=DBLE(K)
      ! пюявер декэрю 
	  D=A+B+C+1.1D0
      I=IDINT(D-C-C)
      J=IDINT(D-B-B)
      KK=IDINT(D-A-A)
    
	  IF(I*J*KK.GT.0) THEN 
           MM=IDINT(D)
           DD=MM+.1D0-D
         ELSE 
           return
      ENDIF
      ! опнбепъел бшонкмемхъ сякнбхъ
	  IF(IDINT(100D0*DD).LT.0) THEN 
         return
	  ENDIF
      
	  ! 3J-яхлбнк нркхвем нр мскъ
      C3JS_CKKISH=1
   ELSE
      return
 ENDIF

 
 return
end function C3JS_CKKISH

! 3J-SYMBOL CALCULATION SUB-PROGRAM
! CHECKING DIFFERENT FROM ZERO OR NOT
! I1, I2, I3, I4, I5, I6-SYMBOL PARAMETERS
 integer function C3JS_W3JAISH(I1,I2,I3,I4,I5,I6)
 implicit none
 integer::I1,I2,I3,I4,I5,I6
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::MIN,MAX
 real(8)::A1,A2,A3,A4,A5,A6,D1,D2




 ! оюпюлерпш 
 A1=DBLE(I1) 
 A2=DBLE(I2)
 A3=DBLE(I3)
 A4=DBLE(I4)
 A5=DBLE(I5)
 A6=DBLE(I6)
      
	
 C3JS_W3JAISH=0

 IF(IDINT(A4+A5+A6+.1D0).NE.0) THEN
    return
 ENDIF
 IF(DABS(A1+A2+A3-DINT(A1+A2+A3+.1D0)).GT..1D0) THEN
    return
 ENDIF
 D1=DMAX1(A2-A3-A4,A1-A3+A5,0.D0)
 D2=DMIN1(A1+A2-A3,A1-A4,A2+A5)
 MIN=IDINT(D1)
 MAX=IDINT(D2)
  
 IF(MIN.GT.MAX) THEN 
    return
 ENDIF
 
 ! пегскэрюр пюяверю 
 C3JS_W3JAISH=1
 
 return
end function C3JS_W3JAISH


end module mc3js
