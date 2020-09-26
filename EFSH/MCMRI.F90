! THE MODULE OF CALCULATION OF MOLECULAR RADIAL INTEGRALS VER 1.0 11.2019
! VER 1.0 NEW  11,2019 цнд

! MODULE FOR CALCULATION OF MOLECULAR RADIAL INTEGRALS VERSION 1.0
	

module mcmri
 implicit none
	
      
 contains


! SUB-PROGRAM FORMS ARRAYS FOR DESCRIBING DIFFERENT GRIDS
 ! Nsetka-NET TYPE
 ! Nsetka = 1-ATOMIC GRID RO = ALFA * R + BETTA * LN (R)
 ! Nsetka = 2-MOLECULAR GRID RO = ALFA * R + BETTA * LN (R) + ATANG ((R-Ral) / GAMMA)
 ! Npoint-NUMBER OF POINTS
 ! H-STEP
 ! ALFA-REPLACEMENT PARAMETERS
 ! BETTA-REPLACEMENT PARAMETERS
 ! GAMMA REPLACEMENT PARAMETERS
 ! Ral-OFFSET RELATIVE TO THE BEGINNING OF THE COORDINATES OF THE "SYSTEM IN WHICH THE DECOMPOSITION OCCURS"
 ! R (Npoint) - ARRAY OF RADIUS VALUES
 ! RO1 (Npoint) - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF THE NEW VARIABLE IN R (OLD VARIABLE)
 ! RO2 (Npoint) - ARRAY OF VALUES OF THE SECOND DERIVATIVE OF THE NEW VARIABLE IN R (OLD VARIABLE)
 ! RO3 (Npoint) - ARRAY OF VALUES OF THE THIRD DERIVATIVE OF THE NEW VARIABLE IN R (OLD VARIABLE)
 subroutine CMRI_VAR(NtipSetky,N,H,ALFA,BET,GAMMA,Ral,R,RO1X,RO2X,RO3X)
        implicit none
        integer::I1,I,N,NtipSetky
	    real(8)::H,ALFA,BET,GAMMA,RO1,RO2,RT,ROT,A,Ral
	    real(8),dimension(:)::R,RO1X,RO2X,RO3X
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		real(8)::R2,R3,RR,RR2,RDELTA,RDELTA2,RDELTA3,XXDF
          
		   
       ! ATOMIC GRID
        IF(NtipSetky.EQ.1) THEN
           RO1=-14.D0
		   RO2=RO1+FLOAT(N)*H 
           RT=(RO2-BET*DLOG(RO2))/ALFA
          DO I1=1,N
             I=N-I1+1
             ROT=RO1+FLOAT(I)*H
  3          A=(ALFA*RT+BET*DLOG(RT)-ROT)/(ALFA*RT+BET)
             RT=RT*(1.D0-A)
             IF(DABS(A)-1.D-12) 2,2,3  
  2          R(I)=RT
             RO1X(I)=(ALFA*RT+BET)/RT
             RO2X(I)=-BET/RT**2
             RO3X(I)=2.D0*BET/RT**3
          ENDDO

        ENDIF
        
		! MOLECULAR GRID
        IF(NtipSetky.EQ.2) THEN
          RO1=-15.D0
          RO2=RO1+FLOAT(N)*H  
		  RT=(RO2-BET*DLOG(RO2)-DATAN2(RO2-Ral,GAMMA))/ALFA
          DO I1=1,N
             I=N-I1+1
             ROT=RO1+FLOAT(I)*H
   30        RR=RT-Ral
	         RR2=RR*RR
	         RDELTA=GAMMA**2+RR2
			 A=(ALFA*RT+BET*DLOG(RT)+DATAN2(RT-Ral,GAMMA)-ROT)*RDELTA/((ALFA*RT+BET)*RDELTA+GAMMA*RT)
             RT=RT*(1.D0-A)
             IF(DABS(A)-1.D-12) 20,20,30  
   20        R(I)=RT
             R2=R(I)*R(I)
	         R3=R(I)*R(I)*R(I)
	         RR=R(I)-Ral
	         RR2=RR*RR
	         RDELTA=GAMMA**2+RR2
             RDELTA2=RDELTA*RDELTA
	         RDELTA3=RDELTA*RDELTA*RDELTA
             XXDF=(ALFA*R(I)+BET)*RDELTA+GAMMA*R(I)
             RO1X(I)=XXDF/(RDELTA*R(I))
	         XXDF=BET*RDELTA2+2.D0*GAMMA*RR*R2
		     RO2X(I)=-XXDF/(RDELTA2*R2)
             XXDF=2.D0*BET*RDELTA3-2.D0*GAMMA*R3*RDELTA
             XXDF=XXDF+8.D0*GAMMA*RR2*R3
		     RO3X(I)=XXDF/(RDELTA3*R3)
          ENDDO

        ENDIF

    

      return 
      end subroutine CMRI_VAR






 !subroutine CMRI_VAR(Nsetka,Npoint,H,ALFA,BETTA,GAMMA,Ral,R,RO1,RO2,RO3)
 ! implicit none
 ! integer::Nsetka,Npoint
 ! real(8)::H,ALFA,BETTA,GAMMA,Ral
 ! real(8),dimension(:)::R,RO1,RO2,RO3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! integer::III,I1
 ! real(8)::ROO1,ROO2,RT,R2,R3,RR,RR2,RDELTA,RDELTA2,RDELTA3
 ! real(8)::XXDF,ROT,AA
      
  ! гюмскъел оепед пюявернл
 ! R=0.D0
 ! RO1=0.D0
 ! RO2=0.D0
 ! RO3=0.D0



      ! юрнлмюъ яерйю
   !   IF(Nsetka.EQ.1) THEN
   !      ROO1=-14.D0
   !      ROO2=ROO1+FLOAT(Npoint)*H
  !       RT=(ROO2-BETTA*DLOG(ROO2))/ALFA

   !      DO III=Npoint,1,-1
   !         ROT=ROO1+FLOAT(III)*H
!	      AA=(ALFA*RT+BETTA*DLOG(RT)-ROT)/(ALFA*RT+BETTA)
!	      DO WHILE(DABS(AA).GT.1.D-6)
 !            AA=(ALFA*RT+BETTA*DLOG(RT)-ROT)/(ALFA*RT+BETTA)
 !            RT=RT*(1.D0-AA)
 !           ENDDO
	   
	      ! тнплхпсел люяяхбш
!	      R(III)=RT
 !           RO1(III)=(ALFA*RT+BETTA)/RT
!            RO2(III)=-BETTA/RT**2
 !           RO3(III)=2.D0*BETTA/RT**3     
!	   ENDDO

!      ENDIF
      ! лнкейскъпмюъ яерйю
!	IF(Nsetka.EQ.2) THEN
!        ! цпюмхжш яерйх он RO
!	    ROO1=-15.D0
!        ROO2=ROO1+FLOAT(Npoint)*H
!        ! жхйк он рнвйюл RO
!	  DO I1=1,Npoint
           !  рнвйю RO
!           RT=ROO1+H*FLOAT(I1)
!	     IF(I1.EQ.1) THEN 
!	call CMRI_SEARCH_ROOT(ALFA,BETTA,GAMMA,Ral,0.D0,0.D0,R(I1),RT)
!             ELSE
!	call CMRI_SEARCH_ROOT(ALFA,BETTA,GAMMA,Ral,R(I1-1),0.D0,R(I1),RT)
!	     ENDIF 
           ! гюонкмъел люяяхбш
!	     R2=R(I1)*R(I1)
!	     R3=R(I1)*R(I1)*R(I1)
!	     RR=R(I1)-Ral
!	     RR2=RR*RR
!	     RDELTA=GAMMA**2+RR2
!           RDELTA2=RDELTA*RDELTA
!	     RDELTA3=RDELTA*RDELTA*RDELTA
!           XXDF=(ALFA*R(I1)+BETTA)*RDELTA+GAMMA*R(I1)
!           RO1(I1)=XXDF/(RDELTA*R(I1))
!	     XXDF=BETTA*RDELTA2+2.D0*GAMMA*RR*R2
!		 RO2(I1)=-XXDF/(RDELTA2*R2)
!           XXDF=2.D0*BETTA*RDELTA3-2.D0*GAMMA*R3*RDELTA
!           XXDF=XXDF+8.D0*GAMMA*RR2*R3
!		 RO3(I1)=XXDF/(RDELTA3*R3)
!        ENDDO!
!	ENDIF


!  return  
! end subroutine CMRI_VAR


	
! EQUATION ROOT SEARCH SUBPROGRAM
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! ALFA, BET, GAMMA, Ral-REPLACEMENT PARAMETERS
! XH-START OF INTERVAL
! XK-END OF INTERVAL
! FOR XK = 0, THE END OF THE INTERVAL IS DETERMINED BY THE SUBPROGRAM
! AT XH = 0-START OF THE INTERVAL DETERMINED BY THE SUB-PROGRAM
! YR-CONSTANT
      subroutine CMRI_SEARCH_ROOT(ALFA,BET,GAMMA,Ral,XHx,XKx,Xc,YR)
      implicit none
      real(8)::ALFA,BET,GAMMA,Ral,Xc,XHx,XKx,YR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8)::X1,X2,REPS,FFR1,FFR2,FFRc,Xz1,Xz2,XH,XK

      
	! рнвмнярэ пюяверю
      REPS=10.D0**(-9)

      XH=XHx
	XK=XKx


	! нопедекъел мхфмчч цпюмхжс еякх нмю ме гюдюмю
     	IF(XH.EQ.0.D0) THEN
         XH=10.D0**(-50)
      ENDIF 
      
	! нопедекъел бепумчч цпюмхжс еякх нмю ме сярюмнбкемю
      IF(XK.EQ.0.D0) THEN
	  FFR1=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,XH,YR)
        XK=1.1D0*XH
        FFR2=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,XK,YR)
	  IF(FFR1.GT.0.D0) THEN
	     DO WHILE(FFR2.GT.0.D0)
              XK=XK+0.1D0
			FFR2=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,XK,YR)
	     ENDDO
	   ELSE
           DO WHILE(FFR2.LT.0.D0)
              XK=XK+0.1D0
              FFR2=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,XK,YR)
	     ENDDO
	  ENDIF
      ENDIF 
      
      X2=XK
	X1=XH
	  ! мюундхл йнпемэ лернднл ахяейжхи 
        DO WHILE(DABS(X2-X1).GT.REPS)
	     ! мюундхл япедмчч рнвйс
	     Xc=(X2+X1)*0.5D0
           ! бшъямъел б йюйнл хг мнбшу хмрепбюкнб мюундхряъ йнпемэ  
           FFR1=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,X1,YR)
		 FFRc=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,Xc,YR)
		 FFR2=CMRI_ROFUN(ALFA,BET,GAMMA,Ral,X2,YR)
           IF(FFR1.GT.0.D0.AND.FFRc.LT.0.D0) THEN
              X2=Xc
	     ENDIF
           IF(FFR1.LT.0.D0.AND.FFRc.GT.0.D0) THEN
              X2=Xc
	     ENDIF
           IF(FFRc.GT.0.D0.AND.FFR2.LT.0.D0) THEN
              X1=Xc
	     ENDIF
           IF(FFRc.LT.0.D0.AND.FFR2.GT.0.D0) THEN
              X1=Xc
	     ENDIF
        ENDDO
      ! ОНКСВЕММШИ ЙНПЕМЭ    
	Xc=(X2+X1)*0.5D0
     
      return
      end subroutine CMRI_SEARCH_ROOT


      ! тсмйжхъ RO (гюлемю б лнкейскъпмни яерйх)
	
	! ALFA,BET,GAMMA-оюпюлерпш гюлемш
	! Ral-ялеыемхъ юрнлмни яхярелш йннпдхмюр хг мювюкю йннпдхмюр "мнбни яхярелш йнндхмюр"
	! X-гмювемхе юпцслемрю
	! YR-онярнъммюъ
      real(8) function CMRI_ROFUN(ALFA,BET,GAMMA,Ral,X,YR)
      implicit none
      real(8)::ALFA,BET,GAMMA,Ral,X,YR
	
	CMRI_ROFUN=0.D0
	IF(X.NE.0.D0) THEN
	CMRI_ROFUN=ALFA*X+BET*DLOG(X)+DATAN2(X-Ral,GAMMA)-YR
      ENDIF
   
      return
      end function CMRI_ROFUN
    

	
     ! SUB-PROGRAM FOR CALCULATING THE ENERGY OF INTERACTION WITH THE NUCLEUS AT THE BEGINNING OF THE COORDINATES
     ! DESCRIPTION OF SUBPROGRAM PARAMETERS
     ! RO0-LEFT INTEGRATION BOUNDARY VALUE
     ! H-STEP
     ! Npoint-NUMBER OF POINTS
     ! R (Npoint) - ARGUMENT VALUE ARRAY
     ! Ffun1 (Npoint) - F-TYPE FIRST WAVE FUNCTION
     ! Ffun2 (Npoint) -SECOND F-TYPE WAVE FUNCTION
     ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of a new variable)
  real(8) function CMRI_CENT(RO0,H,Npoint,R,Ffun1,Ffun2,RO1)
    use mcbi,only:CBI_Integral_First_Type 
    implicit none
    integer::Npoint
    real(8)::Z,RO0,H
    real(8),dimension(:)::R,Ffun1,Ffun2,RO1
    
      
    ! пюявер опнхгбедем б юрнлмшу едемхжюу ноепюрнп бхдю 1/R 
    CMRI_CENT=CBI_Integral_First_Type(1,RO0,0.D0,-1,Npoint,H,R,Ffun1,Ffun2,1.D0,RO1)

    return
  end function  CMRI_CENT


   ! SUB-PROGRAM FOR CALCULATING THE ENERGY OF INTERACTION WITH THE NUCLEUS
   ! AT DISTANCE R FROM THE BEGINNING OF COORDINATES
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! Kgar-HARMONIC NUMBER
   ! Ral-DISTANCE FROM THE BEGINNING OF CORDINATES, VALUE OF RADIUS
   ! Roral-DISTANCE IN NEW VARIABLE, VALUE OF NEW VARIABLE
   ! RO0-LEFT INTEGRATION BOUNDARY VALUE
   ! H-STEP
   ! Npoint-NUMBER OF POINTS
   ! R (Npoint) - ARGUMENT VALUE ARRAY
   ! Ffun1 (Npoint) - F-TYPE FIRST WAVE FUNCTION
   ! Ffun2 (Npoint) -SECOND F-TYPE WAVE FUNCTION
   ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of new variable)
  real(8) function CMRI_CENT_R(Kgar,Ral,Roral,RO0,H,Npoint,R,Ffun1,Ffun2,RO1)
   use mcbi,only:CBI_Integral_First_Type 
   implicit none
   integer::Kgar,Npoint
   real(8)::RO0,H,Ral,Roral
   real(8),dimension(:)::R,Ffun1,Ffun2,RO1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)::RINT1,RINT2
  
   ! бшвхякъел оепбши хмрецпюк
   RINT1=CBI_Integral_First_Type(3,RO0,Roral,Kgar,Npoint,H,R,Ffun1,Ffun2,1.D0/Ral**(Kgar+1),RO1)
   ! бшвхякъел брнпни хмрецпюк 
   RINT2=CBI_Integral_First_Type(2,RO0,Roral,-(Kgar+1),Npoint,H,R,Ffun1,Ffun2,Ral**Kgar,RO1) 
   ! пюявер опнхгбедем б юрнлмшу едемхжюу 
   CMRI_CENT_R=RINT1+RINT2
     
   return
  end function  CMRI_CENT_R




   ! CENTRIFUGAL MEMBER CALCULATION SUB-PROGRAM
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! RO0-LEFT INTEGRATION BOUNDARY VALUE
   ! H-STEP
   ! Npoint-NUMBER OF POINTS
   ! R (Npoint) - ARGUMENT VALUE ARRAY
   ! Ffun1 (Npoint) - F-TYPE FIRST WAVE FUNCTION
   ! Ffun2 (Npoint) -SECOND F-TYPE WAVE FUNCTION
   ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of a new variable)
  real(8) function CMRI_KEN(RO0,H,Npoint,R,Ffun1,Ffun2,RO1)
   use mcbi,only:CBI_Integral_First_Type 
   implicit none
   integer::Npoint
   real(8)::RO0,H
   real(8),dimension(:)::R,Ffun1,Ffun2,RO1
           
   ! пюявер опнхгбедем б юрнлмшу едемхжюу
   CMRI_KEN=CBI_Integral_First_Type(1,RO0,0.D0,-2,Npoint,H,R,Ffun1,Ffun2,1.D0,RO1)
  
   return
  end function CMRI_KEN 


  ! SUBPROGRAM FOR CALCULATING THE SECOND DERIVATIVE
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! RO0-LEFT INTEGRATION BOUNDARY VALUE
  ! H-STEP
  ! Npoint-NUMBER OF POINTS
  ! R (Npoint) - ARGUMENT VALUE ARRAY
  ! Ffun1 (Npoint) - F-TYPE FIRST WAVE FUNCTION
  ! Ffun2 (Npoint) -SECOND F-TYPE WAVE FUNCTION
  ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of new variable)
  ! RO2 (Npoint) -FUNCTIONS ARRAY (second derivative of the new variable)
  ! RO3 (Npoint) -FUNCTION CONVERSION ARRAY (third derivative of the new variable)
 real(8) function CMRI_NABLA(Npoint,RO0,H,R,Ffun1,Ffun2,RO1,RO2,RO3)
  use mcbi,only:CBI_Integral_Second_Type
  implicit none
  integer::K,Npoint
  real(8)::RO0,H
  real(8),dimension(:)::R,Ffun1,Ffun2,RO1,RO2,RO3
  
  ! пюявер опнхгбедем б юрнлмшу едемхжюу
  CMRI_NABLA=CBI_Integral_Second_Type(RO0,H,Npoint,R,Ffun1,Ffun2,RO1,RO2,RO3)
  return
 end function CMRI_NABLA

 ! ондопнцпюллю пюяверю йскнмнбяйнцн хмрецпюкю
 real(8) function CMRI_COULOMB_INTEGRAL(Kgar,H,Npoint,R,Ffun1,Ffun2,Ffun3,Ffun4,RO1)
  use mcbi,only:CBI_Integral_First_Type 
  implicit none
  integer::Kgar,Npoint
  real(8)::H
  real(8),dimension(:)::R,Ffun1,Ffun2,Ffun3,Ffun4,RO1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::ierr,IIOOPP
  real(8)::RINT1,RINT2
  real(8),allocatable,dimension(:)::RPOT
    	
  !бшдекъел оюлърэ дкъ люяяхбнб
  allocate(RPOT(Npoint),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CMRI_COULOMB_INTEGRAL'
  	 write(*,*) 'MEMORY ON THE FILE "RPOT" IS NOT SELECTED'
     stop 
  endif

  ! нясыеярбкъел пюявер онремжхюкю
  RPOT=0.D0
  DO IIOOPP=1,Npoint
     ! бшвхякъел оепбши хмрецпюк
     RINT1=CBI_Integral_First_Type(3,-14.D0,-14.D0+FLOAT(IIOOPP)*H,Kgar,Npoint,H,R,Ffun2,Ffun4,R(IIOOPP)**(-Kgar),RO1)
     ! бшвхякъел брнпни хмрецпюк 
     RINT2=CBI_Integral_First_Type(2,-14.D0,-14.D0+FLOAT(IIOOPP)*H,-(Kgar+1),Npoint,H,R,Ffun2,Ffun4,R(IIOOPP)**(Kgar+1),RO1) 
     ! пюявер опнхгбедем б юрнлмшу едемхжюу 
     RPOT(IIOOPP)=RINT1+RINT2
  ENDDO
  
  ! нясыеярбкъел пюявер йскнмнбяйнцн хмрецпюкю
  DO IIOOPP=1,Npoint
     RPOT(IIOOPP)=RPOT(IIOOPP)*Ffun3(IIOOPP)
  ENDDO
  ! пюявер йскнмнбяйнцн хмрецпюкю б юрнлмшу едемжюу
  CMRI_COULOMB_INTEGRAL=CBI_Integral_First_Type(1,-14.D0,0.D0,-1,Npoint,H,R,Ffun1,RPOT,1.D0,RO1)


	
  ! сдюкемхе люяяхбнб хг оълърх 
  deallocate(RPOT,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'CMRI_COULOMB_INTEGRAL'
     write(*,*) 'THE FILE "RPOT" IS NOT REMOVED FROM MEMORY'
	 stop 
  endif
  return
 end function  CMRI_COULOMB_INTEGRAL



   ! CAPACITY CALCULATION SUB-PROGRAM
   ! Ffun1 (RO) * Ffun2 (RO) * R (RO) ** K * RRcoff / (RO1 (RO)) ** 2
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! K-DEGREE OF THE MULTIPLIER OF THE UNITTEGRAL FUNCTION
   ! RO0-BOTTOM BORDER
   ! H-STEP
   ! Npoint-NUMBER OF POINTS
   ! R (Npoint) - ARGUMENT VALUE ARRAY
   ! Ffun1 (Npoint) -FIRST FUNCTION OF F-TYPE
   ! Ffun2 (Npoint) - F-TYPE SECOND FUNCTION
   ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of a new variable)
   ! POT (Npoint) -ARRAY OF POTENTIAL VALUES
  real(8) function CMRI_POT_COULOMB_INTEGRAL(K,RO0,Npoint,H,R,Ffun1,Ffun2,RO1,POT)
   use mcbi,only:CBI_Integral_First_Type 
   implicit none
   integer::K,Npoint,NTYPE
   real(8)::H,RO0,XRO,RRcoff,CBI_Integral_First_Type_Betta
   real(8),dimension(:)::R,Ffun1,Ffun2,RO1,POT
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::ICBI,ierr,NintAAS,Ny1,Ny2,Ny3,Nyel0
   real(8)::SUMDDS,AR,BR,CR
   real(8),allocatable,dimension(:)::ZINT,UINT,FINT1,FINT2

      
   !бшдекъел оюлърэ дкъ люяяхбнб
   allocate(FINT1(Npoint),stat=ierr)
   if(ierr/=0) then
     write(*,*) 'CBI_Integral_First_Type_Betta'
	 write(*,*) 'MEMORY ON THE FILE "FINT1" IS NOT SELECTED'
	 stop 
   endif  
   allocate(FINT2(Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'CBI_Integral_First_Type_Betta'
  	  write(*,*) 'MEMORY ON THE FILE "FINT2" IS NOT SELECTED'
	  stop 
   endif  
   allocate(ZINT(Npoint),stat=ierr)
   if(ierr/=0) then
     write(*,*) 'CBI_Integral_First_Type_Betta'
	 write(*,*) 'MEMORY ON THE FILE "ZINT" IS NOT SELECTED'
	 stop 
   endif  
   allocate(UINT(Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'CBI_Integral_First_Type_Betta'
  	  write(*,*) 'MEMORY ON THE FILE "UINT" IS NOT SELECTED'
	  stop 
   endif  
	
	
	CMRI_POT_COULOMB_INTEGRAL=0.D0
    ! щрюо 0. тнплхпсел люяяхбш гмювемхи ондхмрецпюкэмшу тсмйжхх
    DO ICBI=1,Npoint 
       FINT1(ICBI)=Ffun1(ICBI)*Ffun2(ICBI)*R(ICBI)**K/RO1(ICBI)**2
       FINT2(ICBI)=Ffun1(ICBI)*Ffun2(ICBI)/(R(ICBI)**(K+1)*RO1(ICBI)**2)
    ENDDO

 
	! гюмскъел оепед пюявернл
    ZINT=0.D0
	UINT=0.D0
      
	! щрюо 1. тнплхпсел люяяхб гмювемхи онремжхюкю     
	ZINT(1)=0.D0
    ZINT(2)=(FINT1(1)+FINT1(2))*H
    UINT(1)=CBI_Integral_First_Type(1,-14.D0,XRO,K,Npoint,H,R,Ffun1,Ffun2,RRcoff,RO1)
	DO ICBI=1,Npoint-2
	   ! гмювемхе тсмйжхх
       ZINT(ICBI+2)=ZINT(ICBI)+(FINT1(ICBI)+4.D0*FINT1(ICBI+1)+FINT1(ICBI+2))*H/3.D0
	ENDDO 

	! щрюо 2. нопедекъел рхо хмрецпхпнбюмхъ


      ! гюохяшбюел пегскэрюр пюяверю
        !CBI_Integral_First_Type=SUMDDS
     







	
	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(ZINT,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'CBI_Integral_First_Type'
      write(*,*) 'THE FILE "RINT" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(UINT,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'CBI_Integral_First_Type'
      write(*,*) 'THE FILE "XINT" IS NOT REMOVED FROM MEMORY'
	stop 
	endif

    return
   end function CMRI_POT_COULOMB_INTEGRAL



  ! SUBPROGRAM FOR CALCULATING THE ORTHOGONALITY INTEGRAL
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! RO0-LEFT INTEGRATION BOUNDARY VALUE
  ! H-STEP
  ! Npoint-NUMBER OF POINTS
  ! Ffun1 (Npoint) - F-TYPE FIRST WAVE FUNCTION
  ! Ffun2 (Npoint) -SECOND F-TYPE WAVE FUNCTION
  ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of new variable)
 real(8) function CMRI_INT_ORT(Npoint,RO0,H,Ffun1,Ffun2,RO1)
  use mcbi,only:CBI_Integral_First_Type 
  implicit none
  integer::Npoint
  real(8)::RO0,H
  real(8),dimension(:)::Ffun1,Ffun2,RO1
   
  ! пюявер опнхгбедем б юрнлмшу едемхжюу
  CMRI_INT_ORT=CBI_Integral_First_Type(1,RO0,0.D0,0,Npoint,H,RO1,Ffun1,Ffun2,1.D0,RO1)
  return
 end function CMRI_INT_ORT

end module mcmri
