! THE MODULE OF CALCULATION OF BASE INTEGRALS VER 1.0 11.2005
! VER 1.0 NEW  11,2019 цнд

! MODULE FOR CALCULATION OF BASE INTEGRALS VERSION 1.0  
	

module mcbi
 implicit none
 
 contains


       
   ! SUB-PROGRAM FOR CALCULATING THE INTEGRAL OF THE FIRST TYPE FROM THE FUNCTION
   ! Ffun1 (RO) * Ffun2 (RO) * R (RO) ** K * RRcoff / (RO1 (RO)) ** 2
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! K-DEGREE OF THE MULTIPLIER OF THE UNITTEGRAL FUNCTION
   ! NTYPE-PARAMETER OF INTEGRATION INTERVAL
   ! NTYPE = 1- (RO0, INFINITY)
   ! NTYPE = 2- (XRO, INFINITY)
   ! NTYPE = 3- (RO0, XRO)
   ! RO0-LOWER INTERVAL BOUNDARY
   ! XRO-BORDER OF INTEGRATION INTERVAL
   ! H-STEP
   ! Npoint-NUMBER OF POINTS
   ! R (Npoint) - ARGUMENT VALUE ARRAY
   ! Ffun1 (Npoint) -FIRST FUNCTION OF F-TYPE
   ! Ffun2 (Npoint) - F-TYPE SECOND FUNCTION
   ! RRcoff-RADIAL PARAMETER
   ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of a new variable)
  real(8) function CBI_Integral_First_Type(NTYPE,RO0,XRO,K,Npoint,H,R,Ffun1,Ffun2,RRcoff,RO1)
   implicit none
   integer::K,Npoint,NTYPE
   real(8)::H,RO0,XRO,RRcoff
   real(8),dimension(:)::R,Ffun1,Ffun2,RO1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::ICBI,ierr,NintAAS,Ny1,Ny2,Ny3,Nyel0
   real(8)::SUMDDS,AR,BR,CR
   real(8),allocatable,dimension(:)::XINT,RINT

      
   !бшдекъел оюлърэ дкъ люяяхбнб
   allocate(RINT(Npoint),stat=ierr)
   if(ierr/=0) then
     write(*,*) 'CBI_Integral_First_Type'
	 write(*,*) 'MEMORY ON THE FILE "RINT" IS NOT SELECTED'
	 stop 
   endif  
   allocate(XINT(Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'CBI_Integral_First_Type'
  	  write(*,*) 'MEMORY ON THE FILE "XINT" IS NOT SELECTED'
	  stop 
   endif  
	
	
	CBI_Integral_First_Type=0.D0
 
	! гюмскъел оепед пюявернл
    RINT=0.D0
	XINT=0.D0
      
	! щрюо 1. тнплхпсел люяяхб гмювемхи ондхмрецпюкэмни тсмйжхх х юпцслемрю     
	DO ICBI=1,Npoint
	   ! гмювемхе юпцслемрю
      XINT(ICBI)=RO0+FLOAT(ICBI)*H
	   ! гмювемхе тсмйжхх
      RINT(ICBI)=RRcoff*Ffun1(ICBI)*Ffun2(ICBI)*R(ICBI)**K/RO1(ICBI)**2
	ENDDO 

	! щрюо 2. нопедекъел рхо хмрецпхпнбюмхъ

	! оепбши рхо хмрецпхпнбюмхъ
      IF(NTYPE.EQ.1) THEN
        ! сярюмюбкхбюел вхякн рнвей(вермне хкх мевермне)
        IF(2*(Npoint/2).EQ.Npoint) THEN
            ! вермне вхякн рнвей
		  ! пюявер дкъ бяецн хмрепбюкю йпнле онякедмецн
		  SUMDDS=CBI_SIMPSON_INT(1,Npoint-1,H,RINT)
	      ! днаюбкъел хмрецпюк лефдс опедонякедмеи х онякедмеи рнвйни
		  ! мюундхл йнщттхжхемрш
	        Ny1=Npoint-2
			Ny2=Npoint-1
			Ny3=Npoint
	      call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		    ! хмрецпюкэмюъ ясллю
			SUMDDS=SUMDDS+AR*(XINT(Ny3)**3-XINT(Ny2)**3)/3.D0
            SUMDDS=SUMDDS+BR*(XINT(Ny3)**2-XINT(Ny2)**2)/2.D0
            SUMDDS=SUMDDS+CR*(XINT(Ny3)-XINT(Ny2))  
	    ELSE
	       ! мевермне вхякн рнвей
           SUMDDS=CBI_SIMPSON_INT(1,Npoint,H,RINT)
	  ENDIF
        ! гюохяшбюел пегскэрюр пюяверю
        CBI_Integral_First_Type=SUMDDS
      ENDIF


      ! брнпни рхо хмрецпхпнбюмхъ
      IF(NTYPE.EQ.2) THEN
        ! нопедекъел онкнфемхъ мхфмецн опедекю
	  Nyel0=0
	  DO ICBI=1,Npoint-1
	     IF(XRO.GE.XINT(ICBI).AND.XRO.LE.XINT(ICBI+1)) THEN
             Nyel0=ICBI
		   EXIT 
	     ENDIF
	  ENDDO 
	  ! опнбепъел мюидемю мхфмъъ цпюмхжю  хмрепбюкю хкх мер
	  IF(Nyel0.EQ.0) THEN
	     WRITE(*,*) 'NTYPE=2,THE BOTTOM BORDER IS NOT FOUND' 
	     READ(*,*) 
	     STOP
	  ENDIF 
         ! сярюмюбкхбюел вхякн рнвей б дюммнл хмрепбюке 
	   ! (сякнбхе опнбепъер бнглнфмнярэ хмрецпхпнбюрэ я онлныэч тнплскш яхлоянмю)
	   ! (оепбюъ вюярэ хмрепбюкю, йсдю оноюдюер мхфмъъ цпюмхжю хмрецхпсеряъ нрдекэмн б мювюке) 
	   IF((Npoint-Nyel0).GE.3) THEN
             ! хмрецпхпсел оепбши хмрепбюк (нр мхфмеи цпюмхжш дн оепбни рнвйх, хмрепбюк (XRO,XINT(Nyel0+1)) )
		   ! мнлепю рнвей
		   Ny1=Nyel0
		   Ny2=Nyel0+1
		   Ny3=Nyel0+2
           call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
             ! хмрецпюкэмюъ ясллю
		   SUMDDS=AR*(XINT(Ny2)**3-XRO**3)/3.D0
             SUMDDS=SUMDDS+BR*(XINT(Ny2)**2-XRO**2)/2.D0
             SUMDDS=SUMDDS+CR*(XINT(Ny2)-XRO) 
		   ! хмрецпхпсеряъ нярюкэмюъ накюярэ
		   ! сярюмюбкхбюел вхякн рнвей(вермне хкх мевермне)
             IF(2*((Npoint-Nyel0)/2).EQ.(Npoint-Nyel0)) THEN
                ! вермне вхякн рнвей
		      ! пюявер дкъ бяецн хмрепбюкю йпнле онякедмецн
		      SUMDDS=SUMDDS+CBI_SIMPSON_INT(Nyel0+1,Npoint-1,H,RINT)
	          ! днаюбкъел хмрецпюк лефдс опедонякедмеи х онякедмеи рнвйни
		      ! мюундхл йнщттхжхемрш
	          Ny1=Npoint-2
			  Ny2=Npoint-1
			  Ny3=Npoint
	     call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		      ! хмрецпюкэмюъ ясллю
			  SUMDDS=SUMDDS+AR*(XINT(Ny3)**3-XINT(Ny2)**3)/3.D0
                SUMDDS=SUMDDS+BR*(XINT(Ny3)**2-XINT(Ny2)**2)/2.D0
                SUMDDS=SUMDDS+CR*(XINT(Ny3)-XINT(Ny2))  
	        ELSE
	         ! мевермне вхякн рнвей
               SUMDDS=SUMDDS+CBI_SIMPSON_INT(Nyel0+1,Npoint,H,RINT)
	       ENDIF
           ELSE
            
		   ! мнлепю рнвей
		   Ny1=Npoint-2
		   Ny2=Npoint-1
		   Ny3=Npoint
	     call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		   ! хмрецпюкэмюъ ясллю
		   SUMDDS=AR*(XINT(Ny3)**3-XRO**3)/3.D0
             SUMDDS=SUMDDS+BR*(XINT(Ny3)**2-XRO**2)/2.D0
             SUMDDS=SUMDDS+CR*(XINT(Ny3)-XRO)  
 
         ENDIF

        ! гюохяшбюел пегскэрюр пюяверю
        CBI_Integral_First_Type=SUMDDS
	ENDIF


      ! рперхи рхо хмрецпхпнбюмхъ
      IF(NTYPE.EQ.3) THEN
        ! нопедекъел онкнфемхе бепумецн опедекю
	  Nyel0=0
	  DO ICBI=1,Npoint-1
	     IF(XRO.GE.XINT(ICBI).AND.XRO.LE.XINT(ICBI+1)) THEN
             Nyel0=ICBI
		   EXIT 
	     ENDIF
	  ENDDO 
	  ! опнбепъел мюидемю бепумъъ цпюмхжю  хмрепбюкю хкх мер
	  IF(Nyel0.EQ.0) THEN
	     WRITE(*,*) 'NTYPE=3,THE TOP BORDER IS NOT FOUND' 
	     READ(*,*) 
	     STOP
	  ENDIF 
	      
	
         ! сярюмюбкхбюел вхякн рнвей б дюммнл хмрепбюке 
	   ! (сякнбхе опнбепъер бнглнфмнярэ хмрецпхпнбюрэ я онлныэч тнплскш яхлоянмю)
	   ! (онякедмъъ вюярэ хмрепбюкю, йсдю оноюдюер бепумъъ цпюмхжю хмрецхпсеряъ нрдекэмн б йнмже) 
	   IF(Nyel0.GE.3) THEN
             ! хмрецпхпсеряъ оепбюъ накюярэ
		   ! сярюмюбкхбюел вхякн рнвей(вермне хкх мевермне)
             IF(2*(Nyel0/2).EQ.Nyel0) THEN
                ! вермне вхякн рнвей
		      ! пюявер дкъ бяецн хмрепбюкю йпнле онякедмецн
		      SUMDDS=CBI_SIMPSON_INT(1,Nyel0-1,H,RINT)
	          ! днаюбкъел хмрецпюк лефдс опедонякедмеи х онякедмеи рнвйни
		      ! мюундхл йнщттхжхемрш
	          Ny1=Nyel0-2
			  Ny2=Nyel0-1
			  Ny3=Nyel0
	     call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		      ! хмрецпюкэмюъ ясллю
			  SUMDDS=SUMDDS+AR*(XINT(Ny3)**3-XINT(Ny2)**3)/3.D0
                SUMDDS=SUMDDS+BR*(XINT(Ny3)**2-XINT(Ny2)**2)/2.D0
                SUMDDS=SUMDDS+CR*(XINT(Ny3)-XINT(Ny2))  
	        ELSE
	         ! мевермне вхякн рнвей
               SUMDDS=CBI_SIMPSON_INT(1,Nyel0,H,RINT)
	       ENDIF

	       ! хмрецпхпсел онякедмхи хмрепбюк (XINT(Nyel0),XRO)
		   ! мнлепю рнвей
		   IF(Nyel0.EQ.(Npoint-1)) THEN
		      Ny1=Nyel0-1
		      Ny2=Nyel0
		      Ny3=Nyel0+1
               ELSE
                Ny1=Nyel0
		      Ny2=Nyel0+1
		      Ny3=Nyel0+2
             ENDIF
           call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
             ! хмрецпюкэмюъ ясллю
		   SUMDDS=SUMDDS+AR*(XRO**3-XINT(Nyel0)**3)/3.D0
             SUMDDS=SUMDDS+BR*(XRO**2-XINT(Nyel0)**2)/2.D0
             SUMDDS=SUMDDS+CR*(XRO-XINT(Nyel0)) 
           ELSE
            
		   ! мнлепю рнвей
		   Ny1=1
		   Ny2=2
		   Ny3=3
	     call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		   ! хмрецпюкэмюъ ясллю
		   SUMDDS=AR*(XRO**3-XINT(1)**3)/3.D0
             SUMDDS=SUMDDS+BR*(XRO**2-XINT(1)**2)/2.D0
             SUMDDS=SUMDDS+CR*(XRO-XINT(1)) 
 
         ENDIF


      ! гюохяшбюел пегскэрюр пюяверю
        CBI_Integral_First_Type=SUMDDS
      ENDIF







	
	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(RINT,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'CBI_Integral_First_Type'
      write(*,*) 'THE FILE "RINT" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(XINT,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'CBI_Integral_First_Type'
      write(*,*) 'THE FILE "XINT" IS NOT REMOVED FROM MEMORY'
	stop 
	endif

      return
      end function CBI_Integral_First_Type



 






	! SUB-PROGRAM FOR CALCULATION OF THE SECOND TYPE INTEGRAL
    ! Ffun1 (R) d ^ 2Ffun2 (R) / (dR)
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! RO0-LOWER INTERVAL BOUNDARY
    ! H-STEP
    ! Npoint-NUMBER OF POINTS
    ! R (Npoint) - ARGUMENT VALUE ARRAY
    ! Ffun1 (Npoint) -FIRST FUNCTION OF F-TYPE
    ! Ffun2 (Npoint) - F-TYPE SECOND FUNCTION
    ! RO1 (Npoint) -FUNCTIONS ARRAY (first derivative of new variable)
    ! RO2 (Npoint) -FUNCTIONS ARRAY (second derivative of the new variable)
    ! RO3 (Npoint) -FUNCTION CONVERSION ARRAY (third derivative of the new variable)
	real(8) function CBI_Integral_Second_Type(RO0,H,Npoint,R,Ffun1,Ffun2,RO1,RO2,RO3)
      implicit none
      integer::Npoint
      real(8)::H,RO0
      real(8),dimension(:)::R,Ffun1,Ffun2,RO1,RO2,RO3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::ICBI,ierr,Ny1,Ny2,Ny3
      real(8)::ROIF,RFHH,RO21,RO31,SUMDDS,AR,BR,CR
	real(8),allocatable,dimension(:)::XINT,RINT

      
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RINT(Npoint),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'CBI_Integral_Second_Type'
	write(*,*) 'MEMORY ON THE FILE "RINT" IS NOT SELECTED'
	stop 
	endif  
      allocate(XINT(Npoint),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'CBI_Integral_Second_Type'
	write(*,*) 'MEMORY ON THE FILE "XINT" IS NOT SELECTED'
	stop 
	endif
	
	
	
	! гюмскъел оепед пюанрни
      RINT=0.D0
      XINT=0.D0

	! щрюо 1. тнплхпсел люяяхб гмювемхи ондхмрецпюкэмни тсмйжхх х юпцслемрю     
      DO ICBI=1,Npoint
	   ! гмювемхе юпцслемрю
         XINT(ICBI)=RO0+FLOAT(ICBI)*H
	   ! гмювемхе тсмйжхх
         RO21=(RO2(ICBI)/RO1(ICBI))**2
	   RO31=RO3(ICBI)/RO1(ICBI)
	   ROIF=0.5D0*(RO31-1.5D0*RO21)/RO1(ICBI)**2
	   RFHH=PRO_FUN2(ICBI,Npoint,H,Ffun2)-Ffun2(ICBI)*ROIF
         RINT(ICBI)=Ffun1(ICBI)*RFHH 
      ENDDO

      
	! нясыеярбкъел хмрецпхпнбюмхе

	! сярюмюбкхбюел вхякн рнвей(вермне хкх мевермне)
      IF(2*(Npoint/2).EQ.Npoint) THEN
          ! вермне вхякн рнвей
		! пюявер дкъ бяецн хмрепбюкю йпнле онякедмецн
		SUMDDS=CBI_SIMPSON_INT(1,Npoint-1,H,RINT)
	    ! днаюбкъел хмрецпюк лефдс опедонякедмеи х онякедмеи рнвйни
		! мюундхл йнщттхжхемрш
	    Ny1=Npoint-2
		Ny2=Npoint-1
		Ny3=Npoint
	    call CBI_COEFFICIENT_POLINOM(Ny1,Ny2,Ny3,XINT,RINT,AR,BR,CR)
		! хмрецпюкэмюъ ясллю
		SUMDDS=SUMDDS+AR*(XINT(Ny3)**3-XINT(Ny2)**3)/3.D0
          SUMDDS=SUMDDS+BR*(XINT(Ny3)**2-XINT(Ny2)**2)/2.D0
          SUMDDS=SUMDDS+CR*(XINT(Ny3)-XINT(Ny2))  
	   ELSE
	    ! мевермне вхякн рнвей
          SUMDDS=CBI_SIMPSON_INT(1,Npoint,H,RINT)
	ENDIF

      ! гюохяшбюел пегскэрюр пюяверю
      CBI_Integral_Second_Type=SUMDDS
    
   return
  end function CBI_Integral_Second_Type

      
! SUBPROGRAM FOR CALCULATING THE SECOND DERIVATIVE FUNCTION
! SECOND DERIVATIVE AT THREE POINTS (STEP CONSTANT)
! SECOND DERIVATIVE AT FIVE POINTS (STEP CONSTANT)
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! IAS-NUMBER OF THE POINT AT WHICH WE CALCULATE THE DERIVATIVE
! Npoint-NUMBER OF POINTS
! H-STEP
! Fun (Npoint) - ARRAY OF FUNCTION VALUES
  real(8)  function PRO_FUN2(IAS,Npoint,H,Fun)
   implicit none
   integer::IAS,Npoint
   real(8)::H
   real(8),dimension(:)::Fun
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer,parameter::NproPoint=9
   
   ! яксвюи рпеу рнвей
   IF(NproPoint.EQ.3) THEN
      IF(IAS.EQ.1) THEN
        PRO_FUN2=(Fun(1)-2.D0*Fun(2)+Fun(3))/H**2
 	  ENDIF
      IF(IAS.EQ.Npoint) THEN
         PRO_FUN2=(Fun(Npoint-2)-2.D0*Fun(Npoint-1)+Fun(Npoint))/H**2
	  ENDIF
      IF(IAS.NE.1.AND.IAS.NE.Npoint) THEN
        PRO_FUN2=(Fun(IAS-1)-2.D0*Fun(IAS)+Fun(IAS+1))/H**2
	  ENDIF
   ENDIF

   ! яксвюи оърх рнвей
   IF(NproPoint.EQ.5) THEN   
      IF(IAS.EQ.1) THEN
        PRO_FUN2=(70.D0*Fun(1)-208.D0*Fun(2)+228.D0*Fun(3)-112.D0*Fun(4)+22.D0*Fun(5))/(24.D0*H**2)
	  ENDIF
      IF(IAS.EQ.2) THEN
        PRO_FUN2=(22.D0*Fun(1)-40.D0*Fun(2)+12.D0*Fun(3)+8.D0*Fun(4)-2.D0*Fun(5))/(24.D0*H**2)
	  ENDIF
      IF(IAS.GT.2.AND.IAS.LT.(Npoint-1)) THEN
        PRO_FUN2=(-2.D0*Fun(IAS-2)+32.D0*Fun(IAS-1)-60.D0*Fun(IAS)+32.D0*Fun(IAS+1)-2.D0*Fun(IAS+2))/(24.D0*H**2)
	  ENDIF
	  IF(IAS.EQ.(Npoint-1)) THEN
        PRO_FUN2=(-2.D0*Fun(Npoint-4)+8.D0*Fun(Npoint-3)+12.D0*Fun(Npoint-2)-40.D0*Fun(Npoint-1)+22.D0*Fun(Npoint))/(24.D0*H**2)
	  ENDIF
      IF(IAS.EQ.Npoint) THEN
        PRO_FUN2=(22.D0*Fun(Npoint-4)-112.D0*Fun(Npoint-3)+228.D0*Fun(Npoint-2)-208.D0*Fun(Npoint-1)+70.D0*Fun(Npoint))/(24.D0*H**2)
	  ENDIF
   
   ENDIF
   
   ! яксвюи яелх рнвей
   IF(NproPoint.EQ.7) THEN   
     
	  IF(IAS.EQ.1) THEN
        PRO_FUN2=(203.D0*Fun(1)/45.D0-261.D0*Fun(2)/15.D0+117.D0*Fun(3)/4.D0-254.D0*Fun(4)/9.D0+33.D0*Fun(5)/2.D0-81.D0*Fun(6)/15.D0+137.D0*Fun(7)/180.D0)/H**2
	  ENDIF
      IF(IAS.EQ.2) THEN
        PRO_FUN2=(548.D0*Fun(1)/720.D0-98.D0*Fun(2)/120.D0-68.D0*Fun(3)/48.D0+94.D0*Fun(4)/36.D0-76.D0*Fun(5)/48.D0+62.D0*Fun(6)/120.D0-52.D0*Fun(7)/720.D0)/H**2
	  ENDIF
      IF(IAS.EQ.3) THEN
        PRO_FUN2=(-52.D0*Fun(1)/720.D0+152.D0*Fun(2)/120.D0-112.D0*Fun(3)/48.D0+40.D0*Fun(4)/36.D0+4.D0*Fun(5)/48.D0-8.D0*Fun(6)/120.D0+8.D0*Fun(7)/720.D0)/H**2
	  ENDIF
	  IF(IAS.GT.3.AND.IAS.LT.(Npoint-2)) THEN
        PRO_FUN2=(8.D0*(Fun(IAS-3)+Fun(IAS+3))/720.D0-18.D0*(Fun(IAS-2)+Fun(IAS+2))/120.D0+72.D0*(Fun(IAS-1)+Fun(IAS+1))/48.D0-98.D0*Fun(IAS)/36.D0)/H**2
	  ENDIF
      IF(IAS.EQ.(Npoint-2)) THEN
        PRO_FUN2=(8.D0*Fun(Npoint-6)/720.D0-8.D0*Fun(Npoint-5)/120.D0+4.D0*Fun(Npoint-4)/48.D0+40.D0*Fun(Npoint-3)/36.D0-112.D0*Fun(Npoint-2)/48.D0+152.D0*Fun(Npoint-1)/120.D0-52.D0*Fun(Npoint)/720.D0)/H**2
	  ENDIF
	  IF(IAS.EQ.(Npoint-1)) THEN
         PRO_FUN2=(-52.D0*Fun(Npoint-6)/720.D0+62.D0*Fun(Npoint-5)/120.D0-76.D0*Fun(Npoint-4)/48.D0+94.D0*Fun(Npoint-3)/36.D0-68.D0*Fun(Npoint-2)/48.D0-98.D0*Fun(Npoint-1)/120.D0+548.D0*Fun(Npoint)/720.D0)/H**2
	  ENDIF
      IF(IAS.EQ.Npoint) THEN
        PRO_FUN2=(548.D0*Fun(Npoint-6)/720.D0-648.D0*Fun(Npoint-5)/120.D0+792.D0*Fun(Npoint-4)/48.D0-1016.D0*Fun(Npoint-3)/36.D0+1404.D0*Fun(Npoint-2)/48.D0-2088.D0*Fun(Npoint-1)/120.D0+3248.D0*Fun(Npoint)/720.D0)/H**2
	  ENDIF
   
   ENDIF

   ! яксвюи дебърх рнвей
   IF(NproPoint.EQ.9) THEN   
     
	  IF(IAS.EQ.1) THEN
        PRO_FUN2=(236248.D0*Fun(1)/40320.D0-138528.D0*Fun(2)/5040.D0+89424.D0*Fun(3)/1440.D0-64096.D0*Fun(4)/720.D0+49752.D0*Fun(5)/576.D0-40608.D0*Fun(6)/720.D0+34288.D0*Fun(7)/1440.D0-29664.D0*Fun(8)/5040.D0+26136.D0*Fun(9)/40320.D0)/H**2
	  ENDIF
      IF(IAS.EQ.2) THEN
        PRO_FUN2=(26136.D0*Fun(1)/40320.D0+128.D0*Fun(2)/5040.D0-5976.D0*Fun(3)/1440.D0+5508.D0*Fun(4)/720.D0-4232.D0*Fun(5)/576.D0+3384.D0*Fun(6)/720.D0-2808.D0*Fun(7)/1440.D0+2396.D0*Fun(8)/5040.D0-2088.D0*Fun(9)/40320.D0)/H**2
	  ENDIF
      IF(IAS.EQ.3) THEN
        PRO_FUN2=(-2088.D0*Fun(1)/40320.D0+5616.D0*Fun(2)/5040.D0-2648.D0*Fun(3)/1440.D0+144.D0*Fun(4)/720.D0+648.D0*Fun(5)/576.D0-592.D0*Fun(6)/720.D0+504.D0*Fun(7)/1440.D0-432.D0*Fun(8)/5040.D0+376.D0*Fun(9)/40320.D0)/H**2
	  ENDIF
	  IF(IAS.EQ.4) THEN
        PRO_FUN2=(376.D0*Fun(1)/40320.D0-684.D0*Fun(2)/5040.D0+2088.D0*Fun(3)/1440.D0-1888.D0*Fun(4)/720.D0+792.D0*Fun(5)/576.D0-36.D0*Fun(6)/720.D0-56.D0*Fun(7)/1440.D0+72.D0*Fun(8)/5040.D0-72.D0*Fun(9)/40320.D0)/H**2
	  ENDIF
	  IF(IAS.GT.4.AND.IAS.LT.(Npoint-3)) THEN
        PRO_FUN2=(-72.D0*Fun(IAS-4)/40320.D0+128.D0*Fun(IAS-3)/5040.D0-288.D0*Fun(IAS-2)/1440.D0+1152.D0*Fun(IAS-1)/720.D0-1640.D0*Fun(IAS)/576.D0+1152.D0*Fun(IAS+1)/720.D0-288.D0*Fun(IAS+2)/1440.D0+128.D0*Fun(IAS+3)/5040.D0-72.D0*Fun(IAS+4)/40320.D0)/H**2
	  ENDIF
      IF(IAS.EQ.(Npoint-3)) THEN
        PRO_FUN2=(-72.D0*Fun(Npoint-8)/40320.D0+72.D0*Fun(Npoint-7)/5040.D0-56.D0*Fun(Npoint-6)/1440.D0-36.D0*Fun(Npoint-5)/720.D0+792.D0*Fun(Npoint-4)/576.D0-1888.D0*Fun(Npoint-3)/720.D0+2088.D0*Fun(Npoint-2)/1440.D0-684.D0*Fun(Npoint-1)/5040.D0+376.D0*Fun(Npoint)/40320.D0)/H**2
	  ENDIF
      IF(IAS.EQ.(Npoint-2)) THEN
        PRO_FUN2=(376.D0*Fun(Npoint-8)/40320.D0-432.D0*Fun(Npoint-7)/5040.D0+504.D0*Fun(Npoint-6)/1440.D0-592.D0*Fun(Npoint-5)/720.D0+648.D0*Fun(Npoint-4)/576.D0+144.D0*Fun(Npoint-3)/720.D0-2648.D0*Fun(Npoint-2)/1440.D0+5616.D0*Fun(Npoint-1)/5040.D0-2088.D0*Fun(Npoint)/40320.D0)/H**2
	  ENDIF
	  IF(IAS.EQ.(Npoint-1)) THEN
        PRO_FUN2=(-2088.D0*Fun(Npoint-8)/40320.D0+2396.D0*Fun(Npoint-7)/5040.D0-2808.D0*Fun(Npoint-6)/1440.D0+3384.D0*Fun(Npoint-5)/720.D0-4232.D0*Fun(Npoint-4)/576.D0+5508.D0*Fun(Npoint-3)/720.D0-5976.D0*Fun(Npoint-2)/1440.D0+128.D0*Fun(Npoint-1)/5040.D0+26136.D0*Fun(Npoint)/40320.D0)/H**2
	  ENDIF
      IF(IAS.EQ.Npoint) THEN
        PRO_FUN2=(26136.D0*Fun(Npoint-8)/40320.D0-29664.D0*Fun(Npoint-7)/5040.D0+34288.D0*Fun(Npoint-6)/1440.D0-40608.D0*Fun(Npoint-5)/720.D0+49752.D0*Fun(Npoint-4)/576.D0-64096.D0*Fun(Npoint-3)/720.D0+89424.D0*Fun(Npoint-2)/1440.D0-138528.D0*Fun(Npoint-1)/5040.D0+236248.D0*Fun(Npoint)/40320.D0)/H**2
	  ENDIF
   
   ENDIF

    



     	
   return
  end function

  









   ! INTEGRAL CALCULATION SUBPROGRAM
   ! USING THE GENERALIZED SIMPSON FORMULA
   ! THE NUMBER OF DOTS MUST BE ODD
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! NH-NUMBER OF THE FIRST POINT OF INTEGRATION
   ! NK-LAST POINT OF INTEGRATION
   ! H-STEP
   ! A (NN) -ARRAY OF FUNCTION VALUES
  real(8)  function CBI_SIMPSON_INT(NH,NK,H,A)
   implicit none
   integer::NH,NK
   real(8)::H
   real(8),dimension(:)::A
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIKOP
   integer,parameter::Ntypp=1
   real(8)::SUMASD
      
    ! тнплскю яхлянмю  
    IF(Ntypp.EQ.1) THEN
       SUMASD=A(NH)+A(NK)
       do IIKOP=1,NK-NH-1
	      ! опнбепъел вермнярэ хмдейяю
          IF(2*(IIKOP/2).EQ.IIKOP) THEN
             SUMASD=SUMASD+2.D0*A(NH+IIKOP)
            ELSE
	         SUMASD=SUMASD+4.D0*A(NH+IIKOP)
	  	  ENDIF    
       enddo
       CBI_SIMPSON_INT=SUMASD*H/3.D0
	ENDIF
    ! тнплскю опълнсцнкэмхйю
    IF(Ntypp.EQ.2) THEN    
       SUMASD=0.D0
	   do IIKOP=NH,NK
          SUMASD=SUMASD+A(IIKOP)
       enddo
       CBI_SIMPSON_INT=SUMASD*H
    ENDIF

   return
  end function CBI_SIMPSON_INT



   ! SUBPROGRAM FOR CALCULATING THE COEFFICIENTS OF THE POLYNOMA OF THE 2ND ORDER
   ! BY CRAMERA METHOD
   ! DESCRIPTION OF SUBPROGRAM PARAMETERS
   ! Np1-FIRST POINT
   ! Np2-SECOND POINT
   ! Np3-THIRD POINT
   ! R (NN) - ARGUMENT VALUE ARRAY
   ! FUN (NN) -ARRAY OF FUNCTION VALUES
   ! A, B, C-COEFFICIENTS OF POLYNOMA Ar ^ 2 + Br + C
	subroutine CBI_COEFFICIENT_POLINOM(Np1,Np2,Np3,R,FUN,A,B,C)
      implicit none
      
      integer::Np1,Np2,Np3
	real(8)::A,B,C
	real(8),dimension(:)::R,FUN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8)::X1R,X2R,X3R,Y1R,Y2R,Y3R,RDX1,RDX2,RDX3,RDX4
     


      ! мюундхл йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
      X1R=R(Np1)
	X2R=R(Np2)
	X3R=R(Np3)
	Y1R=FUN(Np1)
	Y2R=FUN(Np2)
	Y3R=FUN(Np3)
	  
	! нопедекхрекэ яхярелш
	RDX1=X1R**2*(X2R-X3R)+X2R**2*(X3R-X1R)+X3R**2*(X1R-X2R)
      ! нопедекхрекэ йнщттхжхемрю Acoff(1)
	RDX2=Y1R*(X2R-X3R)+Y2R*(X3R-X1R)+Y3R*(X1R-X2R)
      ! нопедекхрекэ йнщттхжхемрю Acoff(2)
	RDX3=Y1R*(X3R**2-X2R**2)+Y2R*(X1R**2-X3R**2)+Y3R*(X2R**2-X1R**2)
      ! нопедекхрекэ йнщттхжхемрю Acoff(3)
	RDX4=X1R**2*(X2R*Y3R-Y2R*X3R)-X1R*(X2R**2*Y3R-Y2R*X3R**2)+Y1R*(X2R**2*X3R-X3R**2*X2R)
	
      ! гюохяшбюел йнщттхжхемрш
      A=RDX2/RDX1
      B=RDX3/RDX1
      C=RDX4/RDX1
         
	return  
	end subroutine CBI_COEFFICIENT_POLINOM




	
	
	end module mcbi
