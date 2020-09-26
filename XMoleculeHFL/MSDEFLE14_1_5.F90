      
!THE MODULE OF THE SOLUTION OF SYSTEM OF THE DIFFERENTIAL EQUATIONS VER 2.0 01.2020
! VER 2.0  01,2020 цнд

!! MODULE FOR SOLVING THE SYSTEM OF DIFFERENTIAL EQUATIONS VERSION 2.0


	
module msde
 implicit none
    
	
 contains

!! SUBPROGRAMME OF THE SOLUTION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS (PROGON + NUMEROV METHOD)
═!! OBTAINING THE WAVE FUNCTION OF A MOLECULAR SHELL
═!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
═!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
═!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
═!! NNEL-MAIN QUANTUM NUMBER OF MOLECULAR ORBITALS
═!! Niteration-Iteration Number
═!! NMO-NUMBER OF MOLECULAR SHELL IN CONFIGURATION
═!! NumeroMin-FIRST HARMONY NUMBER
═!! NumeroMax-NUMBER OF LAST HARMONIC
═!! NLST-KEY INDICATED ON THE PROCESS OF MIXING THE PREVIOUS CONFIGURATION TO THE DATA
═!! IreshimSO-PARAMETER INDICATING THE TYPE OF HARMONIZATION OF MOLECULAR ORBITALS
═!! IreshimSO = 0-STANDARD MODE OF CODE
═!! IreshimSO = 1- "STRENGTHENING" HARMONIZATION MODE THE BINDING COEFFICIENTS ARE EXPRESSED
═!! Npoint-NUMBER OF POINTS
═!! H-STEP
═!! EPS-ACCURACY DETERMINATION OF SINGLE-ELECTRON ENERGY
═!! TETMAX-MAXIMUM DEVIATION OF PREVIOUS ITERATION FROM THIS DATA
═!! TETA (NMO) -MASSIVE DEVIATIONS OF EACH FUNCTION
═!! MAXCNV-NUMBER OF THE POINT APPROPRIATE TO MAXIMUM REJECTION
═!! E3 (NMO) -MASSIVE OF SINGLE-ELECTRON ENERGIES RECEIVED ON NITERATION-1 ITERATION
═!! E4 (NMO) -MASSIVE OF SINGLE-ELECTRONIC ENERGIES OBTAINED ON NITERATION-2 ITERATIONS
═!! E5 (NMO) -MASSIVE OF SINGLE-ELECTRON ENERGIES RECEIVED ON NITERATION-3 ITERATIONS
═!! DeltaEnegry (NMO, 2) -MASSIVE INDICATORY DEVIATION OF SINGLE-ELECTRON ENERGY RECEIVED BY SOLUTION OF A HOMOGENEOUS AND INHOMOGENEOUS EQUATION
═!! RnormMOXFII-NORMAL FACTOR FOR THIS MOLECULAR ORBITAL FOR THIS ITERATION
═!! RnormMOXF (NMO) -MASSIVE NORMALITY COEFFICIENTS FOR EACH MOLECULAR ORBITAL OBTAINED AT THE FIRST STEP OF CRITERIA
═!! DETERMINING RECEIPT OF RIGHT FUNCTIONS AT EACH STAGE OF HARMONIZATION
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! Qz (NumbreGarmMO (NMO), Npoint + 2) -MASSIVE of VALUES of the Q (ro) -function into the equation in the equation Y '' = P * Y + Q-NONLOCAL PART OF POTENTIAL
═!! RFunMO (NMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
═!! FunMOOld (NMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration for writing the forward iteration
═!! FunMONew (NMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of harmonics of the molecular configuration orbitals for writing this iteration not mixed with the previous one
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! ENN0-LOWER ENERGY BORDER
═!! NumreISzam (NumbreMO) -MASSIVE NUMBER OF SHELLS OF EQUIVALENT DATA
═!! ISCalculZZ (NumbreMO, NumbreMO) -MASIVE OF EQUIVALENT DATA SHEET NUMBERS
═!! IKLExchange-KEY TO THE TYPE OF CALCULATION
═!! IKLExchange = 0- CALCULATION WITHOUT EXCHANGE
═!! IKLExchange = 1- CALCULATION WITH EXCHANGE
═!! IRETURN-KEY APPROPRIATE TO THE RESULT OF THE WORK OF THE PROGRAM
═!! IRETURN = 0-REMOVAL SOLVED AND ROOF FOUND
═!! IRETURN = 1-EQUATION NOT SOLVED ROOTS NOT AVAILABLE
═!! RparamterCalcul () - MASSIVE OF INTERMEDIATE PARAMETERS OF CALCULATION
═!! NumeroMaxLimid-MAXIMUM NUMBER HARMONIC FOR THIS MOLECULAR ORBITAL
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! NumbreRegion (NumbreMO) -MASSIVE NUMBER OF REGIONS IN WHICH A ROOT REDUCTION IS CARRIED OUT
═!! EnergyZeroZ (NumbreMO, NNEL + 1,2) -MASSIVE OF THE ROOTS OF A HOMOGENEOUS EQUATION
═!! EnergyRegion (NumbreMO, NNEL + 1,2) -MASSIVE OF AREAS IN WHICH STEP DECREASES
═!! IndexRegionHMIN (NumbreMO, NNEL + 1) -MASSIVE INDEX OF INDICATORS INDEXING STEP FOR DETERMINING BORDERS FOR THIS MINIMUM
═!! EnergyRegionHMIN (NumbreMO, NNEL + 1,2) -MASSIVE OF INITIAL STEPS FOR DETERMINATION OF THE BOUNDARIES OF MINIMUMS
 subroutine SolutionSystemDifferentialEquations(IndexExit,NNEL,Niteration,NMO,NumeroMin,NumeroMax,NLST,IreshimSO,Npoint,H,EPS,TETMAX,TETA,MAXCNV,E3,E4,E5,DeltaEnegry,RnormMOXFII,RnormMOXF,Pz,Qz,RFunMO,FunMOOld,FunMONew,RO1X,ENN0,NumreISzam,ISCalculZZ,IKLExchange,IRETURN,RparamterCalcul,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO,NumbreRegion,EnergyZeroZ,EnergyRegion,IndexRegionHMIN,EnergyRegionHMIN)
   implicit none
   integer::IndexExit,NNEL,Niteration,NMO,NumeroMin,NumeroMax,NLST,Npoint,MAXCNV,IKLExchange,IRETURN,NumeroMaxLimid,IreshimSO
   real(8)::H,EPS,TETMAX,ENN0,RnormMOXFII
   integer,dimension(:)::NumreISzam,NumbreLigand,NumbreRegion
   integer,dimension(:,:)::ISCalculZZ,NLigands,NFunLigands,IndexRegionHMIN
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::TETA,E3,E4,E5,RnormMOXF,RO1X,RparamterCalcul
   real(8),dimension(:,:)::Qz,DeltaEnegry
   real(8),dimension(:,:,:)::Pz,RFunMO,FunMOOld,FunMONew,AlfaCoffMO,EnergyZeroZ,EnergyRegion,EnergyRegionHMIN 
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IndexM,IndexML,IIX,ierr,Kpoint,NEpsIter,IYXDS,IIK,IIO,IparametrDeterEnergy,INEXDSEXF,IZNAK,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IndexSUMZK
   integer::IIXXYYT,IICFHK,IcorenZZ,INlevo,INspravo,IIIOPL,IntParXX,IntParXXZ,IKLSDCalcul,IndexCorrectSET,IndexSolutionsSET,INDEXGRNYT,IndexExitRerern
   integer::IISSXDDZS,IparDeterEnergy,IparDeter,Indexfr,IIISDD,INDEXZZSFA,INDEXEXITRF,IndexIIUMN,IndexMinDifini
   real(8)::RcofXZSD,Rcof,Exf,Henergy,T,RnormMO,C,EZR,BBL,RRnormcoff,ExfZero,ErrorRP,ExfLevo,ExfSpravo,Rkriter,EXXX1,Emin,RHS,RX0,RX1,RX2,RnormX0,RnormX1,RnormX2
   real(8)::Rfun0XZ,Xfun0XZ,Rfun1XZ,Xfun1XZ,Rfun2XZ,Xfun2XZ,RfunMIN,RgarMIN,REmin,E1X,E2X,E3X,E4X,D1X,D2X,D3X,D4X,E1XX,E2XX,E3XX,E4XX,D1XX,D2XX,D3XX,D4XX,ExfTemp,RGHK
   real(8)::HenergySET,ExfSETT,EZZmax,EZZmin,EZZ0cert,EESXZ0,REleft,REright,EEDP,EminGranL,EminGranR
   real(8),dimension(3)::ExfSolutions 
   real(8),dimension(10)::DetM,ExfM,DetMA,ExfMA,Egran,RdeltaZnak,ExfSET,DetSET
   integer,allocatable,dimension(:)::IndexMin  
   integer,allocatable,dimension(:,:)::IXM,IndexGranyMin  
   real(8),allocatable,dimension(:)::AX,BX,CX,DX,UNout,UNin,BXL,EEMmin,DDMmin
   real(8),allocatable,dimension(:,:)::E,AZ,BZ,FZ,DZ,FX,Uout,Uin,VNout,VNin,FunMO,EEGranymin
   real(8),allocatable,dimension(:,:,:)::A,B,AA,BB,Vout,Vin
   
   ! гюмскъел оепед пюявернл 
   ExfLevo=0.D0
   ExfSpravo=0.D0 
   

   ! нопедекъел пюглеп люяяхбю
   IndexM=NumeroMax-NumeroMin+1
   IndexML=NumeroMaxLimid-NumeroMin+1
   
   !бшдекъел оюлърэ дкъ люяяхбнб 
   allocate(EEMmin(NNEL),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "EEMmin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(EEGranymin(NNEL,2),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "EEGranymin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(IndexGranyMin(NNEL,2),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "IndexGranyMin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(DDMmin(NNEL),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "DDMmin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(IndexMin(NNEL-1),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "IndexMin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(FunMO(size(RFunMO,dim=2),size(RFunMO,dim=3)),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "FunMO" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(AX(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "AX" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(BX(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "BX" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(BXL(IndexML),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "BXL" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(CX(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "CX" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(DX(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "DX" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(E(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "E" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(FZ(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "FZ" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(AZ(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "AZ" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(BZ(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "BZ" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(DZ(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "DZ" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif     
   allocate(IXM(IndexM,2),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "IXM" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   !!!!!!!!!!! люяяхбш дкъ пюяверю !!!!!!!!!!!!!!!
   allocate(A(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "A" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(B(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "B" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(FX(IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "FX" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(AA(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "AA" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(BB(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "BB" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   allocate(Vout(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "Vout" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(Vin(IndexM,IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "Vin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(VNout(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "VNout" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(VNin(IndexM,IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "VNin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(Uout(IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "Uout" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(Uin(IndexM,Npoint),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "Uin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(UNout(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "UNout" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif
   allocate(UNin(IndexM),stat=ierr)
   if(ierr/=0) then
      write(*,*)'SolutionSystemDifferentialEquations'
      write(*,*)'MEMORY ON THE FILE "UNin" IS NOT SELECTED'
      read(*,*)
	  stop 
   endif

   !write(6,*)
   !write(6,*) 'local criteriy',Niteration,NMO
 
   34678 FORMAT(2X,100(F15.10,1X))
   35678 FORMAT(2X,100(1X,F20.16))
   89787 FORMAT(5X,I3,5X,I2,1X,F10.5,1X,F10.5,1X,1X,F10.5,1X,F10.5,2X,F7.5)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! щрюо 0. ондцнрнбхрекэмши !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! тнплхпсел едхмхвмсч люрпхжс
   
   Rcof=H*H/12.D0
   E=0.D0
   !DO IIX=1,IndexM
   !   E(IIX,IIX)=1.D0 
   !ENDDO 
   FORALL (IIX=1:IndexM) E(IIX,IIX)=1.D0 
     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! щрюо 1. тнплхпсел люяяхбш свюбярбсчыхе б пюявере!! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! нясыеярбкъел пюявер кнйюкэмшу якюцюелшу
   DO IIX=1,Npoint
      call ReadWriteMassiv2X(IIX,NumeroMin,NumeroMax,Pz,FZ) 
      
	  AZ=E-Rcof*FZ
	  BZ=2.D0*E+10.D0*Rcof*FZ
      	       
	  call ReadWriteMassiv2(1,IIX,IndexM,A,AZ) 
	  call ReadWriteMassiv2(1,IIX,IndexM,B,BZ)  
   ENDDO
      
   ! опнбепъел йкчв дкъ сярюмнбкеме менаундхлнярх б пюявере мендмнпндмнцн якюцюелнцн (налеммнцн якюцюелнцн)
   IF(IKLExchange.EQ.1) THEN
      IntParXX=1
   ENDIF
   IF(IKLExchange.EQ.0) THEN
      IntParXX=0
   ENDIF

   IF(IntParXX.EQ.1) THEN
      DO IIX=2,Npoint-1
         call ReadWriteMassiv1X(2,IIX-1,NumeroMin,NumeroMax,Qz,AX)
		 call ReadWriteMassiv1X(2,IIX,NumeroMin,NumeroMax,Qz,BX) 
         call ReadWriteMassiv1X(2,IIX+1,NumeroMin,NumeroMax,Qz,CX) 
         DX=Rcof*(AX+10.D0*BX+CX)
         call ReadWriteMassiv1(1,IIX,IndexM,FX,DX) 
      ENDDO 
   ENDIF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! щрюо 2. мюундхл ндмнщкейрпнммсч щмепцхч!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! щрюо 2.1  мюундхл ндмнщкейрпнммсч щмепцхч пеьемхел ндмнпндмни яхярелш спюбмемхи !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! сярюмюбкхбюел рхо пюяверю 
   IF(IKLExchange.EQ.1) THEN
      IntParXX=1
   ENDIF
   IF(IKLExchange.EQ.0) THEN
      IntParXX=0
   ENDIF
   ! йкчв сйюгшбючыхи мю пюявер мнплш
   IISSXDDZS=0
   
   ! опнбепъел рхо пюяверю
   IF(Niteration.EQ.1.AND.IntParXX.EQ.0) THEN

        ! дкъ яксвюъ оепбни хрепюжхх х аег налемю дкъ бяеу пеьемхи ндмнпндмнцн спюбмемхи хяонкэгсеряъ ярюмдюпрмши лернд
	    ! ярюмдюпрмши лернд пюяверю
	    WRITE(17,*) 'Energy zero',Niteration,NMO
		call SearchEnergySolutionHomogeneousSystemStandard(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfZero,NumbreRegion,EnergyZeroZ,EnergyRegion) 
	
	    ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
	    Exf=ExfZero
	 ELSE
            
        ! тхйяхпсел лернд пюяверю
        IF(NNEL.LT.4) THEN
		  
           ! ярюмдюпрмши лернд пюяверю
		   WRITE(17,*) 'Energy zero',Niteration,NMO
           call SearchEnergySolutionHomogeneousSystemStandard(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfZero,NumbreRegion,EnergyZeroZ,EnergyRegion) 
            
		   ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
		   Exf=ExfZero 
        ENDIF

        IF(NNEL.GT.3) THEN
		   ! сяйнпеммши лернд пюяверю
		   WRITE(17,*) 'Energy zero',Niteration,NMO
           ! дкъ оепбшу рпеу хрепюжхи янупюмъел бяе лерндш 
		   IF(Niteration.LT.4) THEN
              
			   call SearchEnergySolutionHomogeneousSystemAccelerated(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,NumbreRegion,EnergyZeroZ,EnergyRegion) 
            
		       ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
               ExfZero=ExfSolutions(2)
               Exf=ExfSolutions(2) 
		      
			  ELSE
			   ! опнбепъел нясыеярбкъеряъ пеьемхе спюбмемхъ уюпрпх хкх уюпхпх-тнйю
		       ! пеьемхе спюбмемхъ уюпрпх
		       IF(IntParXX.EQ.0) THEN
                  call SearchEnergySolutionHomogeneousSystemAccelerated(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,NumbreRegion,EnergyZeroZ,EnergyRegion) 
            
		          ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
                  ExfZero=ExfSolutions(2)
                  Exf=ExfSolutions(2) 
               ENDIF
		       ! пеьемхе спюбмемхъ уюпрпх-тнйю
		       ! мнбши лернд пеьемхъ спюбмемхъ уюпрпх-тнйю нямнбюм мю нопедекемхъ лхмхлслнб лефдс йнрнпшлх мюундхряъ йнпемэ
		       ! дкъ гюосяйю опнжедспш менаундхлн нопедекхрэ оепбши йнпемэ
               IF(IntParXX.EQ.1) THEN
		          ! сйюгшбюел рн врн б дюммнл яксвюе мнплю ме пюявхршбюеряъ
				  IISSXDDZS=1
                  
				  call SearchEnergySolutionHomogeneousSystemAccelerated2(IndexExit,Niteration,NMO,1,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,NumbreRegion,EnergyZeroZ,EnergyRegion) 
            
		          ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
                  ExfZero=ExfSolutions(1)
                  Exf=ExfSolutions(1)  
		       ENDIF
           ENDIF
        ENDIF
   ENDIF


   ! сярюмюбкхбюел рхо пюяверю
   IF(IKLExchange.EQ.1) THEN
      IntParXX=1
   ENDIF
   IF(IKLExchange.EQ.0) THEN
      IntParXX=0
   ENDIF

   ! нопедекъел мнплхпнбнвмсч йнмярюмрс йпхрепхи нранпю пеьемхи
   ! опнбепъел рхо бшундю
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.1.AND.IntParXX.EQ.1) THEN 
         ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
         call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
    
	     ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
         call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
       
         ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
         call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
 
         ! мюундхл гмювемхъ цюплнмхй лнкейскъпмни нпахрюкх б рнвйе яьхбйх
         call CalculationValuesGarmonikMO(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Vin,FunMO,AZ,BZ,E,AX,BX)
    
         ! нясыеярбкъел бняярюмнбкемхе гмювемхи тсмйжхх б бяеу рнвйюу нр рнвйх яьхбйх х мнплхпсел онксвеммсч тсмйжхч дкъ яксвюъ ндмнпндмни яхярелш
         call FunctionRecovery(1,Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,BXL,CX,VNout,UNout,VNin,UNin,RO1X,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO) 
   
         ! гюохяшбюел мнплхпнбнвмши йнщттхжхемр (йпхрепхи нопедекемхъ опюбхкэмни ндмнщкейрпнммни щмепцххх)
         RnormMOXF(NMO)=FunMO(NumeroMin,2)
      ENDIF
   ENDIF
     



   ! нясыеярбкъел пюявер мнплхпнбнвмни йнмярюмрш дкъ дюммни лнкейскъпмни нпахрюкх дкъ дюммни хрепюжхх
   RnormMOXFII=0.D0
   IF(IISSXDDZS.EQ.0) THEN
      ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
    
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
 
      ! мюундхл гмювемхъ цюплнмхй лнкейскъпмни нпахрюкх б рнвйе яьхбйх
      call CalculationValuesGarmonikMO(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Vin,FunMO,AZ,BZ,E,AX,BX)
    
      ! нясыеярбкъел бняярюмнбкемхе гмювемхи тсмйжхх б бяеу рнвйюу нр рнвйх яьхбйх х мнплхпсел онксвеммсч тсмйжхч дкъ яксвюъ ндмнпндмни яхярелш
      call FunctionRecovery(1,Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,BXL,CX,VNout,UNout,VNin,UNin,RO1X,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO) 
   
      ! гюохяшбюел мнплхпнбнвмши йнщттхжхемр (йпхрепхи нопедекемхъ опюбхкэмни ндмнщкейрпнммни щмепцххх)
      RnormMOXFII=FunMO(NumeroMin,2)
   ENDIF
    
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! щрюо 2.2  мюундхл ндмнщкейрпнммсч щмепцхч пеьемхел мендмнпндмни яхярелш спюбмемхи !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! сярюмюбкхбюел рхо пюяверю
   IF(IKLExchange.EQ.1) THEN
      IntParXX=1
   ENDIF
   IF(IKLExchange.EQ.0) THEN
      IntParXX=0
   ENDIF


   ! мюундхл ндмнщкейрпнммсч щмепцхч б яксвюе пеьемхъ мендмнпндмнцн спюбмемхъ
   IF(IntParXX.EQ.1) THEN 

       IISSXDDZS=0
	   ! нясыеярбкъел бшанп лефдс лерндюлх нопедекемхъ ндмнщкейрпнммни щмепцхх
	   IF(NNEL.GT.3) THEN
	 	  IF(Niteration.LT.4) THEN
              ! дкъ оепбшу рпеу хрепюжхи янупюмъел бяе лерндш 
		 	  IISSXDDZS=1  
	         ELSE
              ! опхлемъел мнбши лернд дкъ онхяйю ндмнщкейрпнммни щмепцхх
			  IISSXDDZS=2 
	      ENDIF
	   ENDIF

       ! мнбши лернд нямнбюммши мю онхяйе лхмхлслнб лефдс йнрнпшлх мюундъряъ йнпмх 
       IF(IISSXDDZS.EQ.2) THEN
          
		  ! щрюо 1. нопедекъел оепбши лхмхлсл дкъ сярюмнбкемхъ хмрепбюкю щмепцхи б йнрнпнл асдср нопедекемш йнпмх
		 
		  ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
          Henergy=0.1D0
          ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
          !DO IICFHK=1,NNEL
          !   ! хглемемхе ьюцю 
	      !   Henergy=Henergy/FLOAT(IICFHK+1)  
          !ENDDO
      
	      ! нопедекъел лхмхлсл якебю  
		  E1X=ExfZero*(1.D0-Henergy) 
          E2X=ExfZero*(1.D0-2.D0*Henergy)
		  E3X=ExfZero*(1.D0-3.D0*Henergy) 
		  ! опюбюъ цпюмхжю
		  EEDP=E3X 
	      ! нопедекъел якебю х яопюбю нр гмювемхъ ндмнщкейрпнммни щмепцхх онксвеммни опх пеьемхх ндмнпндмнцн спюбмемхъ
          WRITE(17,*) 'MIN',Niteration,NMO 
          ! 1. нопедекъел оепбне лхмхлюкэмне гмювемхе мнплш  
		  call MinValueNormEnergyFunctionDD(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,E3X,E2X,E1X,REmin,RfunMIN,REleft,IndexRegionHMIN(NMO,1),EnergyRegionHMIN(NMO,1,1))
		  
		  ! щрюо 2. нопедекъел опхакхфеммше гмювемхъ б хмрепбюке (0,REmin)
          ! оепбши лхмхлсл
		  EEMmin(1)=REmin
		  DDMmin(1)=RfunMIN
		  
          ! гюмскъел оепед пюявернл 
		  IndexGranyMin=0
		  EEGranymin=0.D0 
          
		  ! тхйяхпсел кебсч х опюбсч цпюмхжш оепбнцн лхмхлслю
          ! кебюъ цпюмхжю
		  IndexGranyMin(1,1)=1
          EEGranymin(1,1)=REleft
		  ! опюбюъ цпюмхжю
          IndexGranyMin(1,2)=1
          EEGranymin(1,2)=EEDP


		  ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
          Henergy=0.001D0
          ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
          DO IICFHK=1,NNEL
             ! хглемемхе ьюцю 
	         Henergy=Henergy/FLOAT(IICFHK+1)  
          ENDDO
          	 	
          
          ! сярюмюбкхбюел мювюкэмше гмювемхъ дкъ онкнфемхи лхмхлслнб
		  DO IICFHK=2,NNEL
		     ! мювюкэмне гмювемхе дкъ нопедекемхе якедсчыецн лхмхлслю
			 REmin=EEMmin(IICFHK-1)           
		     
			 ! оюпюлерп бшундю хг жхйкю 
             IparDeterEnergy=1 
             DO WHILE(IparDeterEnergy.EQ.1)
		        ! нопедекъел мювюкэмне опхакхфемхъ дкъ онкнфемхъ лхмхлслю
				REmin=REmin*DBLE(NNEL-1)/DBLE(NNEL)  
				WRITE(17,*) 'PremierShagMIN',Niteration,IICFHK
				! нопедекъел лхмхлюкэмне гмювемхе
   		        call MinValueNormEnergyFunctionD2(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,Henergy,REmin,Rkriter,RfunMIN,IndexGranyMin(IICFHK-1,1),EEGranymin(IICFHK-1,1),EEMmin(IICFHK-1),DDMmin(IICFHK-1),REleft,REright,IndexMinDifini,IndexRegionHMIN(NMO,IICFHK),EnergyRegionHMIN(NMO,IICFHK,1),EnergyRegionHMIN(NMO,IICFHK,2))
		        
				! сярюмюбкхбюел онксвем мнбши лхмхлсл хкх мер
				IF(IndexMinDifini.EQ.0) THEN
                     ! рпеанбюмхъ бшонкмемш
				     EEMmin(IICFHK)=Rkriter
				     DDMmin(IICFHK)=RfunMIN
				     ! сякнбхе бшундю нясыеярбкемн
				     IparDeterEnergy=2   
					 
					 ! гюохяшбюел опюбсч цпюмхжс лхмхлслю
                     IndexGranyMin(IICFHK,2)=1 
					 EEGranymin(IICFHK,2)=REright
					                    
					 ! гюохяшбюел кебсч цпюмхжс лхмхлслю
                     IndexGranyMin(IICFHK,1)=1
					 EEGranymin(IICFHK,1)=REleft

				   ELSE
				     ! сякнбхе ме бшонкмемн гмювхрэ щмепцхъ нрмняхряъ й опедедсыелс лхмхлслх
					 IF(REmin.LT.EEMmin(IICFHK-1)) THEN
					    ! опнбепъел гюохяшбюкюяэ кх сфе кебюъ цпюмхжю
						IF(IndexGranyMin(IICFHK-1,1).EQ.0) THEN
                            ! гюохяшбюел цпюмхжс б оепбши пюг
                            IndexGranyMin(IICFHK-1,1)=1 
							EEGranymin(IICFHK-1,1)=REmin
						   ELSE
						    ! опнбепъел дюммюъ щмепцхъ кефхр мхфе хкх мер
							IF(REmin.LT.EEGranymin(IICFHK-1,1)) THEN
							   EEGranymin(IICFHK-1,1)=REmin
							ENDIF	   
						ENDIF
					 ENDIF 
					 
				ENDIF
			 ENDDO           
		  ENDDO 
		  
		  WRITE(17,*)
		  WRITE(17,*) 'Min Value Energy Zero approximation',Niteration,NMO
          DO IICFHK=1,NNEL
             WRITE(17,*) EEMmin(IICFHK),DDMmin(IICFHK)
		  ENDDO
          WRITE(17,*)
		   

		  ! нопедекъел онкнфемхе лхмхлслнб
          IndexMin=0
          ! оюпюлерп бшундю хг жхйкю 
          IparDeter=1 
          DO WHILE(IparDeter.EQ.1)
		     
             ! жхйк он хмрепбюкюл лефдс лхмхлслюлх
			 DO IICFHK=1,NNEL-1
			    ! опнбепъел дюммши хмрепбюк ме яндепфхр лхмхлслнб хкх мер
				Indexfr=0
				IF(IndexMin(IICFHK).EQ.0) THEN
                   ! пюяялюрпхбюел накюярэ лефдс лхмхлслюлх
				   IF(IndexGranyMin(IICFHK,1).EQ.1) THEN
                        EZZmax=EEGranymin(IICFHK,1)
					  ELSE
                        EZZmax=EEMmin(IICFHK)  
				   ENDIF
                   IF(IndexGranyMin(IICFHK+1,2).EQ.1) THEN
                        EZZmin=EEGranymin(IICFHK+1,2)
					  ELSE
                        EZZmin=EEMmin(IICFHK+1)  
				   ENDIF

				   EZZ0cert=0.25D0*EZZmax+0.75D0*EZZmin 
				   EESXZ0=0.5D0*(EZZmax+EZZmin) 
                   ! нопедекъел хлееряъ б дюммни накюярх днонкмхрекэмше лхмхлслш хкх мер
                   IparDeterEnergy=1
				   DO WHILE(IparDeterEnergy.EQ.1)
				      ! нопедекъел йюйнлс лхмхлслс яннрберярбсер дюммюъ рнвйю
					  WRITE(17,*) 'DefiniMIN',Niteration,IICFHK
					  WRITE(17,*) 'EnergyGranyIntervala'
					  WRITE(17,*) EZZmax,EZZmin,EZZ0cert
                      call MinValueNormEnergyFunctionD3(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,Henergy,EZZ0cert,Rkriter,RfunMIN,IndexGranyMin(IICFHK+1,2),EEMmin(IICFHK+1),EEGranymin(IICFHK+1,2),DDMmin(IICFHK+1),IndexGranyMin(IICFHK,1),EEGranymin(IICFHK,1),EEMmin(IICFHK),DDMmin(IICFHK),IndexMinDifini,EminGranL,EminGranR,IndexRegionHMIN(NMO,IICFHK+1),EnergyRegionHMIN(NMO,IICFHK+1,1),EnergyRegionHMIN(NMO,IICFHK+1,2))
                      Indexfr=0
					  ! опнбепъел дюммюъ рнвйю оноюкю б накюярэ кебнцн хкх опюбнцн лхмхлслю
					  ! нопедекъел яннрберярбхе MAX  
				      IF(IndexMinDifini.EQ.2) THEN
                         ! дюммюъ рнвйю яннрберярбсер лхмхлслс б рнвйе EZZmax
                         Indexfr=1
                         EZZmax=EZZ0cert
						 ! йнппейрхпсел цпюмхжс накюярх дюммнцн лхмхлслю
						 IF(IndexGranyMin(IICFHK,1).EQ.1) THEN
						     ! цпюмхжю ясыеярбсер нясыеярбкъел ее йнппейрхпнбйс 
                             IF(EZZ0cert.LT.EEGranymin(IICFHK,1)) THEN
							    EEGranymin(IICFHK,1)=EZZ0cert
							 ENDIF                                    
							ELSE
                             ! цпюмхжю менопедекемю тхйяхпсел цпюмхжс
							 IndexGranyMin(IICFHK,1)=1
							 EEGranymin(IICFHK,1)=EZZ0cert 
						 ENDIF 
				      ENDIF
					  ! нопедекъел яннрберярбхе MIN 
                      IF(IndexMinDifini.EQ.1) THEN
                         ! дюммюъ рнвйю яннрберярбсер лхмхлслс б рнвйе EZZmin
                         Indexfr=2
                         EZZmin=EZZ0cert
                         
						 ! йнппейрхпсел цпюмхжс накюярх дюммнцн лхмхлслю
						 IF(IndexGranyMin(IICFHK+1,2).EQ.1) THEN
						     ! цпюмхжю ясыеярбсер нясыеярбкъел ее йнппейрхпнбйс 
                             IF(EZZ0cert.GT.EEGranymin(IICFHK+1,2)) THEN
							    EEGranymin(IICFHK+1,2)=EZZ0cert
							 ENDIF                                    
							ELSE
                             ! цпюмхжю менопедекемю тхйяхпсел цпюмхжс
							 IndexGranyMin(IICFHK+1,2)=1
							 EEGranymin(IICFHK+1,2)=EZZ0cert 
						 ENDIF  
				      
					  ENDIF
					  WRITE(17,*) 'Value Interval'
                      WRITE(17,*) Indexfr,DABS(EZZmax-EZZmin),DABS(EESXZ0*Henergy)
                      ! мнбюъ япедмъъ рнвйю
					  EZZ0cert=0.25D0*EZZmax+0.75D0*EZZmin 
				      ! опнбепъел мюидем мнбши лхмхлсл хкх мер
                      IF(IndexMinDifini.EQ.0) THEN
				         ! мюидем мнбши лхмхлсл
					     ! оепеярпюхбюел онпъднб лхмхлслнб
                         IF((IICFHK+1).NE.NNEL) THEN
						    DO IIISDD=NNEL,IICFHK+2,-1
                               EEMmin(IIISDD)=EEMmin(IIISDD-1)
						       DDMmin(IIISDD)=DDMmin(IIISDD-1) 
							   IndexGranyMin(IIISDD,1)=IndexGranyMin(IIISDD-1,1) 
							   IndexGranyMin(IIISDD,2)=IndexGranyMin(IIISDD-1,2)   
					        ENDDO
					     ENDIF
						 ! гюохяшбюел мнбши лхмхлсл
                         EEMmin(IICFHK+1)=Rkriter
                         DDMmin(IICFHK+1)=RfunMIN
						 ! тхйяхпсел цпюмхжш дкъ дюмнцн лхмхлслю
                         IndexGranyMin(IICFHK+1,1)=1 
						 IndexGranyMin(IICFHK+1,2)=1
						 EEGranymin(IICFHK+1,1)=EminGranL
						 EEGranymin(IICFHK+1,2)=EminGranR  
                         ! оняке оепеярпнхйх онхяй лхмхлслнб менаундхлн гюосяхрэ гюмнбн
                         Indexfr=4
						 EXIT
				      ENDIF 
					  ! опнбепъел пюгмхжю лефдс цпюмхжюлх лемэье лхмхлюкэмнцн гмювемхъ
					  IF(DABS(EZZmax-EZZmin).LT.4.D0*DABS(EESXZ0*Henergy)) THEN
					     ! б дюммнл хмрепбюке нрясрярбсчр дпсцхе лхмхлслш
						 IndexMin(IICFHK)=1
						 EXIT  
					  ENDIF 
				   ENDDO
                ENDIF

				! опнбепъел нясыеярбкем оепеялнрп онякеднбюрекэмнярх лхмхлслнб
				IF(Indexfr.EQ.4) THEN
                   EXIT
				ENDIF
			 ENDDO

			 ! опнбепъел мюидемш бяе хмрепбюкш бмсрпх йнрнпшу мер лхмхлслнб
             Indexfr=0
			 DO IICFHK=1,NNEL-1
                IF(IndexMin(IICFHK).EQ.1) THEN
                   Indexfr=Indexfr+1 
				ENDIF
             ENDDO
			 ! опнбепъел нопедекемш бяе хмрепбюкш хкх мер
			 IF(Indexfr.EQ.(NNEL-1)) THEN
                IparDeter=2 
			 ENDIF
		  ENDDO

		  ! тхйяхпсел мювюкэмши ьюц дкъ дюкэмеиьецн пюяверю 
		  ! свхршбюъ рнр тюйр, врн нопедекемш бяе цпюмхжш дкъ бяеу лхмхлслнб рн мюундхл пюгмхжс лефдс лхлхлслнл х цпюмхжеи
          DO IICFHK=1,NNEL
		     ! тхйяхпсел, врн ьюц дкъ дюммнцн лхмхлслю нопедекем 
             IndexRegionHMIN(NMO,IICFHK)=1
			 EnergyRegionHMIN(NMO,IICFHK,1)=DABS(1.D0-EEGranymin(IICFHK,1)/EEMmin(IICFHK)) 
             EnergyRegionHMIN(NMO,IICFHK,2)=DABS(1.D0-EEMmin(IICFHK)/EEGranymin(IICFHK,2)) 
          ENDDO 

		  WRITE(17,*)
		  WRITE(17,*) 'Min Value Energy',Niteration,NMO
          DO IICFHK=1,NNEL
             WRITE(17,*) EEMmin(IICFHK),DDMmin(IICFHK)
		  ENDDO
          WRITE(17,*)

		  ! щрюо 3. нопедекъел йнпмх 
          ! опнбепъел хлееряъ йнпмх хкх мер
		  ! опнбепъел мнплс яопюбю
		  IF(RnormMOXF(NMO).GT.DDMmin(NNEL-1)) THEN
		       WRITE(17,*) 'Defini grani Energy right',Niteration,NMO
		       ! нясыеярбкъел юмюкхг накюярх лефдс лхмхлслюлх х сярюмюбкхбюел леярн оепелемш гмюйю опнхгбндмни
			   ! гюмскъел оепед пюявернл
			   DetM=0.D0
			   DetMA=0.D0
			   ExfM=0.D0
			   ! гюохяшбюел оепбне гмювемхе
			   ExfM(7)=EEMmin(NNEL-1)
               DetM(7)=DDMmin(NNEL-1)
               ! нопедекъел мювюкэмши ьюц
               Henergy=0.02D0*DABS(EEMmin(NNEL)-EEMmin(NNEL-1))
               ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
			   E1X=EEMmin(NNEL-1)*(1-Henergy) 
			   IndexSUMZK=0	
			   ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
			   IndexIIUMN=0	         			  
			   ! оюпюлерп бшундю хг жхйкю 
               IparametrDeterEnergy=1 
			   ! оюпюлерп бшундю хг жхйкю
			   INDEXEXITRF=0
			   INDEXZZSFA=0
			   RcofXZSD=1.D0
               DO WHILE(IparametrDeterEnergy.EQ.1) 
			      ! хмдейя жхйкю
				  IndexSUMZK=IndexSUMZK+1
		          ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
				  call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			      WRITE(17,*) E1X,D1X,Henergy
				  ! гюохяшбюел пегскэрюрш пюяверю
				  ExfM(1)=ExfM(2)
                  DetM(1)=DetM(2)
				  ExfM(2)=ExfM(3)
                  DetM(2)=DetM(3)
				  ExfM(3)=ExfM(4)
                  DetM(3)=DetM(4)
				  ExfM(4)=ExfM(5)
                  DetM(4)=DetM(5)
                  ExfM(5)=ExfM(6)
                  DetM(5)=DetM(6)
                  ExfM(6)=ExfM(7)
                  DetM(6)=DetM(7)
				  ExfM(7)=E1X
                  DetM(7)=D1X
							 
				  IF(IndexSUMZK.EQ.1) THEN
                     ! пюявхршбюел опнхгбндмсч
                     DetMA(7)=DetM(7)*(1.D0-DetM(6)/DetM(7))/(ExfM(7)*(1.D0-ExfM(6)/ExfM(7)))
				  ENDIF
                  IF(IndexSUMZK.GT.1) THEN
                     ! пюявхршбюел опнхгбндмсч
					 DetMA(1)=DetMA(2)
					 DetMA(2)=DetMA(3)
                     DetMA(3)=DetMA(4)
					 DetMA(4)=DetMA(5)
					 DetMA(5)=DetMA(6)
                     DetMA(6)=DetMA(7)
				     DetMA(7)=DetM(7)*(1.D0-DetM(6)/DetM(7))/(ExfM(7)*(1.D0-ExfM(6)/ExfM(7)))

                     ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
				     call CorrectHRightD(IndexSUMZK,INDEXZZSFA,IndexIIUMN,Henergy,RcofXZSD,DetMA)
		
					 
					 ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
					 IF(INDEXZZSFA.EQ.1) THEN
                        ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
                        E1X=ExfM(5) 
                        ExfM(7)=ExfM(5) 
                        DetM(7)=DetM(5)
                        ExfM(6)=ExfM(4) 
                        DetM(6)=DetM(4)
						ExfM(5)=ExfM(3) 
                        DetM(5)=DetM(3)
						ExfM(4)=ExfM(2) 
                        DetM(4)=DetM(2)
						ExfM(3)=ExfM(1) 
                        DetM(3)=DetM(1)
                        DetMA(7)=DetMA(5)
                        DetMA(6)=DetMA(4)
						DetMA(5)=DetMA(3)
						DetMA(4)=DetMA(2)
						DetMA(3)=DetMA(1)
						! намскъел бшунд хг жхйкю
						INDEXEXITRF=0
						! опнбепъел вхякн хрепюжхи
			            IF(IndexSUMZK.LT.3) THEN
                           ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				           IndexSUMZK=0 
			            ENDIF
					 ENDIF
					  
				  ENDIF

				  ! мнбне опхакхфемхе дкъ щмепцхх                 
                  E1X=E1X*(1-Henergy)
				  ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
				  IF(INDEXZZSFA.EQ.0) THEN
                     IF(D1X.GT.RnormMOXF(NMO)) THEN 
                        INDEXEXITRF=INDEXEXITRF+1
					 ENDIF
					 ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
					 IF(INDEXEXITRF.EQ.5) THEN
					    ! сякнбхе бшундю бшонкмемш
                        IparametrDeterEnergy=2
					 ENDIF
				  ENDIF
			   ENDDO

			   ! 4. нопедекъел йнпемэ 
               ! нопедекъел мювюкэмше гмювемхъ опнжеяяю онхяйю
               RnormMO=RnormMOXF(NMO) 
               ! кебюъ цпюмхжю 
			   E1X=ExfM(7)
               E2X=ExfM(6)
			   E3X=ExfM(5)
               ! опюбюъ цпюмхжю
			   E4X=EEMmin(NNEL-1)
			   D4X=DDMmin(NNEL-1)
               ! жхйк он нопедекемхч йнпмъ
			   WRITE(17,*) 'Energy right',Niteration,NMO,RnormMO
               ! оюпюлерп бшундю хг жхйкю 
               IparametrDeterEnergy=1 
               DO WHILE(IparametrDeterEnergy.EQ.1)
			      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
			  	  call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E2X
			      call CalculationNormFunctionEnergyXF(E2X,D2X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3X
			      call CalculationNormFunctionEnergyXF(E3X,D3X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E4X
			      call CalculationNormFunctionEnergyXF(E4X,D4X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				  ! нясыеярбкъел опнбепйс сярюмнбкемн мскхбне гмювемхе хкх мер
				  ! мюундхл жемрпюкэмши йнпемэ
				  E3XX=((RnormMO-D2X)*E4X-(RnormMO-D4X)*E2X)/(D4X-D2X) 
                  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3XX
				  call CalculationNormFunctionEnergyXF(E3XX,D3XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			      ! опнбепъел бшонкмъеряъ сякнбхе бшундю
                  IF(DABS(D3XX-RnormMO).LT.EPS) THEN
                     ! сякнбхе бшундю бшонкмемн
				  	 IparametrDeterEnergy=2
					 Exf=E3XX 
				  ENDIF
                    
				  ! нясыеярбкъел пюявер 
				  ! гюохяшбюел дюммне опхакхфемхе 
				  E1XX=E1X
				  D1XX=D1X

				  E2XX=E2X
				  D2XX=D2X
					 
				  E3XX=E3X
				  D3XX=D3X

				  E4XX=E4X
				  D4XX=D4X
					 
				  ! нясыеярбкъел пюявер мнбнцн опхакхфемхъ
                  E2X=E2XX-(D2XX-RnormMO)*(E3XX-E1XX)/(D3XX-D1XX)
                  E3X=E2X*(1.D0+Henergy) 
                  E1X=E2X/(1.D0+Henergy)  
				  E4X=(E2XX*(D4XX-RnormMO)-E4XX*(D2XX-RnormMO))/(D4XX-D2XX)  
				  WRITE(17,*) E2X,E4X,D3XX
			   ENDDO 
					   

               ! йнпемэ меидем
		       INspravo=1
		       ExfSpravo=Exf
		       ! тхйяхпсел йнпемэ якебю йюй мскхбни 
               ! еякх мюидем мхфе кефюыхи анкее бшянйн кефюыхи ме мсфмн хяйюрэ
		       INlevo=0 
		       ExfLevo=0.D0
               ExfZero=0.D0  
			 ELSE 
                             
			  ! йнпемэ ме ясыеярбсер
              INspravo=0
		      ExfSpravo=0.D0
			  ExfZero=0.D0  
		  ENDIF 
          
		  ! опнбепъел мюидем йнпемэ яопюбю хкх мер
          IF(INspravo.EQ.0) THEN   
		     ! опнбепъел мнплс якебю
             IF(RnormMOXF(NMO).GT.DDMmin(NNEL)) THEN
                  WRITE(17,*) 'Defini grani Energy left',Niteration,NMO
		          ! нясыеярбкъел юмюкхг накюярх лефдс лхмхлслюлх х сярюмюбкхбюел леярн оепелемш гмюйю опнхгбндмни
			      ! гюмскъел оепед пюявернл
			      DetM=0.D0
			      DetMA=0.D0
			      ExfM=0.D0
			      ! гюохяшбюел оепбне гмювемхе
			      ExfM(7)=EEMmin(NNEL)
                  DetM(7)=DDMmin(NNEL)
                  ! нопедекъел мювюкэмши ьюц
                  Henergy=0.02D0*DABS(EEMmin(NNEL)-EEMmin(NNEL-1))
                  ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
			      E1X=EEMmin(NNEL)*(1+Henergy) 
			      IndexSUMZK=0	
			      ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
			      IndexIIUMN=0	         			  
			      ! оюпюлерп бшундю хг жхйкю 
                  IparametrDeterEnergy=1 
			      ! оюпюлерп бшундю хг жхйкю
			      INDEXEXITRF=0
			      INDEXZZSFA=0
			      RcofXZSD=1.D0
                  DO WHILE(IparametrDeterEnergy.EQ.1) 
			         ! хмдейя жхйкю
				     IndexSUMZK=IndexSUMZK+1
		             ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
				     call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			         WRITE(17,*) E1X,D1X,Henergy
				     ! гюохяшбюел пегскэрюрш пюяверю
				     ExfM(1)=ExfM(2)
                     DetM(1)=DetM(2)
				     ExfM(2)=ExfM(3)
                     DetM(2)=DetM(3)
				     ExfM(3)=ExfM(4)
                     DetM(3)=DetM(4)
				     ExfM(4)=ExfM(5)
                     DetM(4)=DetM(5)
                     ExfM(5)=ExfM(6)
                     DetM(5)=DetM(6)
                     ExfM(6)=ExfM(7)
                     DetM(6)=DetM(7)
				     ExfM(7)=E1X
                     DetM(7)=D1X
							 
				     IF(IndexSUMZK.EQ.1) THEN
                        ! пюявхршбюел опнхгбндмсч
                        DetMA(7)=DetM(7)*(1.D0-DetM(6)/DetM(7))/(ExfM(7)*(1.D0-ExfM(6)/ExfM(7)))
				     ENDIF
                     IF(IndexSUMZK.GT.1) THEN
                        ! пюявхршбюел опнхгбндмсч
					    DetMA(1)=DetMA(2)
					    DetMA(2)=DetMA(3)
                        DetMA(3)=DetMA(4)
					    DetMA(4)=DetMA(5)
					    DetMA(5)=DetMA(6)
                        DetMA(6)=DetMA(7)
				        DetMA(7)=DetM(7)*(1.D0-DetM(6)/DetM(7))/(ExfM(7)*(1.D0-ExfM(6)/ExfM(7)))

                        ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
			            call CorrectHLeftD(IndexSUMZK,INDEXZZSFA,IndexIIUMN,Henergy,RcofXZSD,DetMA)
			   
					 
					    ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
					    IF(INDEXZZSFA.EQ.1) THEN
                           ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
                           E1X=ExfM(5) 
                           ExfM(7)=ExfM(5) 
                           DetM(7)=DetM(5)
                           ExfM(6)=ExfM(4) 
                           DetM(6)=DetM(4)
						   ExfM(5)=ExfM(3) 
                           DetM(5)=DetM(3)
						   ExfM(4)=ExfM(2) 
                           DetM(4)=DetM(2)
						   ExfM(3)=ExfM(1) 
                           DetM(3)=DetM(1)
                           DetMA(7)=DetMA(5)
                           DetMA(6)=DetMA(4)
						   DetMA(5)=DetMA(3)
						   DetMA(4)=DetMA(2)
						   DetMA(3)=DetMA(1)
						   ! намскъел бшунд хг жхйкю
						   INDEXEXITRF=0
						   ! опнбепъел вхякн хрепюжхи
			               IF(IndexSUMZK.LT.3) THEN
                              ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				              IndexSUMZK=0 
			               ENDIF
					    ENDIF
					  
				     ENDIF

				     ! мнбне опхакхфемхе дкъ щмепцхх                 
                     E1X=E1X*(1+Henergy)
				     ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
				     IF(INDEXZZSFA.EQ.0) THEN
                        IF(D1X.GT.RnormMOXF(NMO)) THEN 
                           INDEXEXITRF=INDEXEXITRF+1
					    ENDIF
					    ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
					    IF(INDEXEXITRF.EQ.5) THEN
					       ! сякнбхе бшундю бшонкмемш
                           IparametrDeterEnergy=2
					    ENDIF
				     ENDIF
			      ENDDO

                  !4. нопедекъел йнпемэ
				  ! нопедекъел мювюкэмше гмювемхъ опнжеяяю онхяйю
                  RnormMO=RnormMOXF(NMO) 
                  ! опюбюъ цпюмхжю 
				  E1X=ExfM(7)
                  E2X=ExfM(6)
				  E3X=ExfM(5)
                  ! кебюъ цпюмхжю
				  E4X=EEMmin(NNEL)
				  D4X=DDMmin(NNEL)
                  ! жхйк он нопедекемхч йнпмъ
				  WRITE(17,*) 'Energy left',Niteration,NMO,RnormMO
                  ! оюпюлерп бшундю хг жхйкю 
                  IparametrDeterEnergy=1 
                  DO WHILE(IparametrDeterEnergy.EQ.1)
			         ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
					 call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			         ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E2X
			         call CalculationNormFunctionEnergyXF(E2X,D2X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					 ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3X
			         call CalculationNormFunctionEnergyXF(E3X,D3X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					 ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E4X
			         call CalculationNormFunctionEnergyXF(E4X,D4X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				     ! нясыеярбкъел опнбепйс сярюмнбкемн мскхбне гмювемхе хкх мер
					 ! мюундхл жемрпюкэмши йнпемэ
					 E3XX=((RnormMO-D2X)*E4X-(RnormMO-D4X)*E2X)/(D4X-D2X) 
                     ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3XX
					 call CalculationNormFunctionEnergyXF(E3XX,D3XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			         ! опнбепъел бшонкмъеряъ сякнбхе бшундю
                     IF(DABS(D3XX-RnormMO).LT.EPS) THEN
                        ! сякнбхе бшундю бшонкмемн
					    IparametrDeterEnergy=2
						Exf=E3XX 
					 ENDIF


				     ! нясыеярбкъел пюявер 
				     ! гюохяшбюел дюммне опхакхфемхе 
				     E1XX=E1X
					 D1XX=D1X

					 E2XX=E2X
					 D2XX=D2X
					  
					 E3XX=E3X
					 D3XX=D3X

					 E4XX=E4X
					 D4XX=D4X
					 
					 ! нясыеярбкъел пюявер мнбнцн опхакхфемхъ
                     E2X=E2XX-(D2XX-RnormMO)*(E3XX-E1XX)/(D3XX-D1XX)
                     E3X=E2X*(1.D0-Henergy) 
                     E1X=E2X/(1.D0-Henergy)  
					 E4X=(E2XX*(D4XX-RnormMO)-E4XX*(D2XX-RnormMO))/(D4XX-D2XX) 
					 WRITE(17,*) E2X,E4X,D3XX 
				  ENDDO 


                  ! йнпемэ мюидем
		          INlevo=1
		          ExfLevo=Exf
				  ExfZero=0.D0
		        ELSE
                 ! йнмемэ ме мюидем
		         INlevo=0 
		         ExfLevo=0.D0
				 ExfZero=0.D0
             ENDIF 
          ENDIF

	   ENDIF
      

	   ! сяйнпеммши лернд пюяверю
       IF(IISSXDDZS.EQ.1) THEN
		   
		    ! мскхбне опхакхфемхе
			ExfZero=Exf			 
           
			! мюундхл ндмнщкейрпнммше щмепцхх яннрберярбсчыхе дюммнлс гмювемхч мнплхпнбнвмни йнмярюмрш
	        ! якебю х яопюбю нр гмювемхъ ндмнщкейрпнммни щмепцхх онксвеммни опх пеьемхх ндмнпндмнцн спюбмемхъ
            WRITE(17,*) 'MIN',Niteration,NMO 
            ! нопедекъел гмювемхе яопюбю 
			! 1. нопедекъел лхмхлюкэмне гмювемхе мнплш  б хмрепбюке  (NNEL-1, NNEL)
			call MinValueNormEnergyFunction(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfSolutions(2),ExfSolutions(1),REmin,RfunMIN)
			
					
			! 2. япюбмхбюел онксвеммне лхмхлюкэмне гмювемхе мнплш я мнплни пеьемхъ дкъ сярюмнбкемхъ мюкхвхъ йнпмъ
			! нопедекъел ндмнщкейрпнммне гмювемхе щмепцхх б яннрберябхх я хмрецпюкэмшл йпхрепхел 
			IF(RnormMOXF(NMO).GT.RfunMIN) THEN
                
				! йнпемэ ясыеярбсер
				! нясыеярбкъел онхяй йнпмъ   
               
			      ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
                  Henergy=0.1D0
                  ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
                  DO IICFHK=1,NNEL
                     ! хглемемхе ьюцю 
	                 Henergy=Henergy/FLOAT(IICFHK+2)  
                  ENDDO

                  ! нопедекъел менаундхл пюявер хкх мер
                  IntParXXZ=2
	              IF(Niteration.EQ.1) THEN  
                     ! пюявер аег налемю
		             IF(IKLExchange.EQ.0) THEN
                        IntParXXZ=2
                     ENDIF
                     ! пюявер я налемнл
                     IF(IKLExchange.EQ.1) THEN
                        IntParXXZ=1
                     ENDIF
                  ENDIF

      
	              ! опнбепъел бшонкмемхъ йпхрепхъ дкъ дюммнцн ьюцю
	              ! ьюц днкфем ашрэ лемэье нрмньемхъ нрйкнмемхъ й мскхбнлс опхакхфемхч щмепцхх опх щрнл ашрэ рюйнбшл, 
	              ! врнаш лнфмн ашкн ядекюрэ рпх ьюцю
	              IF(Niteration.NE.IntParXXZ) THEN 
	                 Rkriter=DeltaEnegry(NMO,2)/(3.D0*ExfZero)
	                 DO WHILE(Henergy.GT.Rkriter)
		                Henergy=Henergy*0.97D0   
                     ENDDO
	              ENDIF
                  
                  ! нопедекъел мювюкэмше гмювемхъ опнжеяяю онхяйю
                  RnormMO=RnormMOXF(NMO) 
                  ! кебюъ цпюмхжю 
				  E1X=ExfZero*(1.D0+Henergy) 
                  E2X=ExfZero*(1.D0+2.D0*Henergy)
				  E3X=ExfZero*(1.D0+3.D0*Henergy) 
                  ! опюбюъ цпюмхжю
				  E4X=REmin 
				  D4X=RfunMIN   
                  ! жхйк он нопедекемхч йнпмъ
				  WRITE(17,*) 'Energy right',Niteration,NMO,RnormMO
                  ! оюпюлерп бшундю хг жхйкю 
                  IparametrDeterEnergy=1 
                  DO WHILE(IparametrDeterEnergy.EQ.1)
			         ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
					 call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			         ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E2X
			         call CalculationNormFunctionEnergyXF(E2X,D2X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					 ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3X
			         call CalculationNormFunctionEnergyXF(E3X,D3X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					 ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E4X
			         call CalculationNormFunctionEnergyXF(E4X,D4X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				     ! нясыеярбкъел опнбепйс сярюмнбкемн мскхбне гмювемхе хкх мер
					 ! мюундхл жемрпюкэмши йнпемэ
					 E3XX=((RnormMO-D2X)*E4X-(RnormMO-D4X)*E2X)/(D4X-D2X) 
                 	 ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3XX
					 call CalculationNormFunctionEnergyXF(E3XX,D3XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			         ! опнбепъел бшонкмъеряъ сякнбхе бшундю
                     IF(DABS(D3XX-RnormMO).LT.EPS) THEN
                        ! сякнбхе бшундю бшонкмемн
						IparametrDeterEnergy=2
						Exf=E3XX 
					 ENDIF
                    
				     ! нясыеярбкъел пюявер 
				     ! гюохяшбюел дюммне опхакхфемхе 
				     E1XX=E1X
					 D1XX=D1X

					 E2XX=E2X
					 D2XX=D2X
					 
					 E3XX=E3X
					 D3XX=D3X

					 E4XX=E4X
					 D4XX=D4X
					 
					 ! нясыеярбкъел пюявер мнбнцн опхакхфемхъ
                     E2X=E2XX-(D2XX-RnormMO)*(E3XX-E1XX)/(D3XX-D1XX)
                     E3X=E2X*(1.D0+Henergy) 
                     E1X=E2X/(1.D0+Henergy)  
					 E4X=(E2XX*(D4XX-RnormMO)-E4XX*(D2XX-RnormMO))/(D4XX-D2XX)  
					 WRITE(17,*) E2X,E4X,D3XX
				  ENDDO 
					   

                   ! йнпемэ меидем
		           INspravo=1
		           ExfSpravo=Exf
		           ! тхйяхпсел йнпемэ якебю йюй мскхбни 
                   ! еякх мюидем мхфе кефюыхи анкее бшянйн кефюыхи ме мсфмн хяйюрэ
		           INlevo=0 
		           ExfLevo=0.D0
               ELSE 
                             
			    ! йнпемэ ме ясыеярбсер
                INspravo=0
		        ExfSpravo=0.D0  
			ENDIF


			! 3. опнбепъел мюидем йнпемэ яопюбю хкх мер
            IF(INspravo.EQ.0) THEN
			   ! йнпемэ яопюбю ме мюидем нясыеярбкъел онхяй йнпмъ якебю
			   ! нясыеярбкъел онхяй онякедмеи ндмнщкейрпнммни щмепцхх дкъ NNEL+1
			   WRITE(17,*) 'Energy zero Left',Niteration,NMO
			   call SearchEnergySolutionHomogeneousSystemAcceleratedLeft2(NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,INDEXGRNYT) 
              
			   ! 1. нопедекъел лхмхлюкэмне гмювемхе мнплш  б хмрепбюке  (NNEL, NNEL+1)
			   WRITE(17,*) 
			   IF(INDEXGRNYT.EQ.0) THEN
			      ! хмрепбюк онкмнярэч нопедекем (нопедекемю цпюмхжю) 
			      call MinValueNormEnergyFunction(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfSolutions(3),ExfSolutions(2),REmin,RfunMIN)
			   ENDIF
		       IF(INDEXGRNYT.EQ.1) THEN
			      ! нопедекемю рнвйю б хмрепбюке нрмняхрекэмн йнрнпни мюундхряъ лхмхлсл 
                  call MinValueNormEnergyFunctionD(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfSolutions(3),ExfSolutions(2),REmin,RfunMIN)
			   ENDIF

               ! 2. япюбмхбюел онксвеммне лхмхлюкэмне гмювемхе мнплш я мнплни пеьемхъ дкъ сярюмнбкемхъ мюкхвхъ йнпмъ
			   IF(RnormMOXF(NMO).GT.RfunMIN) THEN
                    ! йнпемэ ясыеярбсер
				    ! нясыеярбкъел онхяй йнпмъ   		
			        ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
                    Henergy=0.1D0
                    ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
                    DO IICFHK=1,NNEL
                       ! хглемемхе ьюцю 
	                   Henergy=Henergy/FLOAT(IICFHK+1)  
                    ENDDO
      
	                ! нопедекъел менаундхл пюявер хкх мер
                    IntParXXZ=2
	                IF(Niteration.EQ.1) THEN  
                       ! пюявер аег налемю
		               IF(IKLExchange.EQ.0) THEN
                          IntParXXZ=2
                       ENDIF
                       ! пюявер я налемнл
                       IF(IKLExchange.EQ.1) THEN
                          IntParXXZ=1
                       ENDIF
                    ENDIF



	                ! опнбепъел бшонкмемхъ йпхрепхъ дкъ дюммнцн ьюцю
	                ! ьюц днкфем ашрэ лемэье нрмньемхъ нрйкнмемхъ й мскхбнлс опхакхфемхч щмепцхх опх щрнл ашрэ рюйнбшл, 
	                ! врнаш лнфмн ашкн ядекюрэ рпх ьюцю
                    IF(Niteration.NE.IntParXXZ) THEN 
	                   Rkriter=DeltaEnegry(NMO,1)/(3.D0*ExfZero)
	                   DO WHILE(Henergy.GT.Rkriter)
		                  Henergy=Henergy*0.97D0   
                       ENDDO
	                ENDIF

	                ! мюундхл ндмнщкейрпнммсч щмепцхч якебю  
	                ! йнщттхжхемр мнплхпнбйх дкъ лнкейскъпмни нпахрюкх йпхрепхи онксвемхъ опюбхкэмни ндмнщкейрпнммни щмепцхх  
                   
				    ! нопедекъел мювюкэмше гмювемхъ опнжеяяю онхяйю
                    RnormMO=RnormMOXF(NMO) 
                    ! опюбюъ цпюмхжю 
				    E1X=ExfZero*(1.D0-Henergy) 
                    E2X=ExfZero*(1.D0-2.D0*Henergy)
				    E3X=ExfZero*(1.D0-3.D0*Henergy) 
                    ! кебюъ цпюмхжю
				    E4X=REmin 
				    D4X=RfunMIN   
                    ! жхйк он нопедекемхч йнпмъ
                    ! оюпюлерп бшундю хг жхйкю 
                    IparametrDeterEnergy=1 
                    DO WHILE(IparametrDeterEnergy.EQ.1)
			           ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
					   call CalculationNormFunctionEnergyXF(E1X,D1X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			           ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E2X
			           call CalculationNormFunctionEnergyXF(E2X,D2X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					   ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3X
			           call CalculationNormFunctionEnergyXF(E3X,D3X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
					   ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E4X
			           call CalculationNormFunctionEnergyXF(E4X,D4X,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
				       ! нясыеярбкъел опнбепйс сярюмнбкемн мскхбне гмювемхе хкх мер
					   ! мюундхл жемрпюкэмши йнпемэ
					   E3XX=((RnormMO-D2X)*E4X-(RnormMO-D4X)*E2X)/(D4X-D2X) 
                       ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E3XX
					   call CalculationNormFunctionEnergyXF(E3XX,D3XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
			           ! опнбепъел бшонкмъеряъ сякнбхе бшундю
                       IF(DABS(D3XX-RnormMO).LT.EPS) THEN
                          ! сякнбхе бшундю бшонкмемн
						  IparametrDeterEnergy=2
						  Exf=E3XX 
					   ENDIF


				       ! нясыеярбкъел пюявер 
				       ! гюохяшбюел дюммне опхакхфемхе 
				       E1XX=E1X
					   D1XX=D1X

					   E2XX=E2X
					   D2XX=D2X
					  
					   E3XX=E3X
					   D3XX=D3X

					   E4XX=E4X
					   D4XX=D4X
					 
					   ! нясыеярбкъел пюявер мнбнцн опхакхфемхъ
                       E2X=E2XX-(D2XX-RnormMO)*(E3XX-E1XX)/(D3XX-D1XX)
                       E3X=E2X*(1.D0-Henergy) 
                       E1X=E2X/(1.D0-Henergy)  
					   E4X=(E2XX*(D4XX-RnormMO)-E4XX*(D2XX-RnormMO))/(D4XX-D2XX)  
				    ENDDO 


                    ! йнпемэ мюидем
		            INlevo=1
		            ExfLevo=Exf
				  
				  ELSE
                    ! йнмемэ ме мюидем
		            INlevo=0 
		            ExfLevo=0.D0
               ENDIF

			ENDIF

	   ENDIF
      


	   ! ярюмдюпрмши лернд пюяверю
	   IF(NNEL.LT.4) THEN 
	    
		 ! мюундхл ндмнщкейрпнммше щмепцхх яннрберярбсчыхе дюммнлс гмювемхч мнплхпнбнвмни йнмярюмрш
	     ! якебю х яопюбю нр гмювемхъ ндмнщкейрпнммни щмепцхх онксвеммни опх пеьемхх ндмнпндмнцн спюбмемхъ
	     ExfZero=Exf
      

         !WRITE(6,*) 'Energy One-Electron NE ODNORODNOE LEVOE',Niteration,NMO,Exf,Kpoint
         !WRITE(*,*) 'EEEX NE ODNORODNOE LEVOE',ExfLevo
         !READ(*,*) 
      
         ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
         Henergy=0.1D0
         ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
         DO IICFHK=1,NNEL
            ! хглемемхе ьюцю 
	        Henergy=Henergy/FLOAT(IICFHK+2)  
         ENDDO

         ! нопедекъел менаундхл пюявер хкх мер
         IntParXXZ=2
	     IF(Niteration.EQ.1) THEN  
            ! пюявер аег налемю
		    IF(IKLExchange.EQ.0) THEN
               IntParXXZ=2
            ENDIF
            ! пюявер я налемнл
            IF(IKLExchange.EQ.1) THEN
               IntParXXZ=1
            ENDIF
         ENDIF

      
	     ! опнбепъел бшонкмемхъ йпхрепхъ дкъ дюммнцн ьюцю
	     ! ьюц днкфем ашрэ лемэье нрмньемхъ нрйкнмемхъ й мскхбнлс опхакхфемхч щмепцхх опх щрнл ашрэ рюйнбшл, 
	     ! врнаш лнфмн ашкн ядекюрэ рпх ьюцю
	     IF(Niteration.NE.IntParXXZ) THEN 
	        Rkriter=DeltaEnegry(NMO,2)/(3.D0*ExfZero)
	        DO WHILE(Henergy.GT.Rkriter)
		       Henergy=Henergy*0.97D0   
            ENDDO
	     ENDIF

         ! мюундхл ндмнщкейрпнммсч щмепцхч яопюбю  
	     ! йнщттхжхемр мнплхпнбйх дкъ лнкейскъпмни нпахрюкх йпхрепхи онксвемхъ опюбхкэмни ндмнщкейрпнммни щмепцхх  
         RnormMO=RnormMOXF(NMO) 
         DetM=0.D0
         ExfM=0.D0
         ExfM(2)=ExfZero*(1.D0+Henergy) 
	     Exf=ExfZero*(1.D0+Henergy)   
      
         ! оюпюлерп бшундю хг жхйкю он пюяверс щмепцхх
         IparametrDeterEnergy=1
         NEpsIter=0
         INEXDSEXF=0
	     IcorenZZ=0
         ! оюпюлерп сйюгшбючыхи мю йнпемэ йнрнпши хыеряъ яопюбю
         IZNAK=2
         ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
         IndexVilky=0
         INDEXDF=0
         IIKK=0
         IndexCorrect=0
		 IndexExitRerern=0
         WRITE(17,*) 'Energy right',Niteration,NMO
	     !WRITE(*,*) 'EEEX NOW PRAVOE',Exf 
         !READ(*,*) 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !жхйк нопедекемхе ндмнщкейрпнммни щмепцхх лнкейскъпмни нпахрюкх яопюбю!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO WHILE(IparametrDeterEnergy.EQ.1)
            NEpsIter=NEpsIter+1
	        ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	        call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
	        !WRITE(*,*) 'Kpoint',Kpoint
	        !WRITE(*,*) '4 Exf',Exf

	        ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	        call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
            !WRITE(*,*) '5 Exf',Exf
            ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
            call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            !WRITE(*,*) '6 Exf',Exf  
            ! йнппейрхпсел гмювемхе щмепцхх Exf
	        !write(*,*) 'eeee',Exf 
    	    call CorrectionEnergy(IZNAK,NEpsIter,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DetM,ExfM,Exf,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,EPS,RdeltaZnak,IndexVilky,INDEXDF,IIKK,IndexCorrect,IcorenZZ,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO,IndexExitRerern)  
			!write(*,*) '3 exf',Exf 
         ENDDO

	
      
	     ! сярюмюбкхбюел мюидем кх яопюбю йнпемэ
	     IF(IcorenZZ.EQ.1) THEN
             ! йнпемэ мюидем
		     INspravo=1
		     ExfSpravo=Exf
		     ! тхйяхпсел йнпемэ якебю йюй мскхбни 
             ! еякх мюидем мхфе кефюыхи анкее бшянйн кефюыхи ме мсфмн хяйюрэ
		     INlevo=0 
		     ExfLevo=0.D0
            ELSE
		     ! йнмемэ ме мюидем
		     INspravo=0
		     ExfSpravo=0.D0 
	     ENDIF 
         
		 IF(Niteration.EQ.11.AND.NMO.EQ.3) THEN
            WRITE(6,*) 'Energy right',Niteration,NMO 
			WRITE(6,*) 'ExfSpravo',INspravo,ExfSpravo

		 ENDIF


	     ! онхяй йнпмъ якебю нясыеярбкъел б яксвюе нрясрярбхъ йнпмъ яопюбю
	     IF(INspravo.EQ.0) THEN
         
	        ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
            Henergy=0.1D0
            ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
            DO IICFHK=1,NNEL
               ! хглемемхе ьюцю 
	           Henergy=Henergy/FLOAT(IICFHK+1)  
            ENDDO
      
	        ! нопедекъел менаундхл пюявер хкх мер
            IntParXXZ=2
	        IF(Niteration.EQ.1) THEN  
               ! пюявер аег налемю
		       IF(IKLExchange.EQ.0) THEN
                  IntParXXZ=2
               ENDIF
               ! пюявер я налемнл
               IF(IKLExchange.EQ.1) THEN
                  IntParXXZ=1
               ENDIF
            ENDIF



	        ! опнбепъел бшонкмемхъ йпхрепхъ дкъ дюммнцн ьюцю
	        ! ьюц днкфем ашрэ лемэье нрмньемхъ нрйкнмемхъ й мскхбнлс опхакхфемхч щмепцхх опх щрнл ашрэ рюйнбшл, 
	        ! врнаш лнфмн ашкн ядекюрэ рпх ьюцю
            IF(Niteration.NE.IntParXXZ) THEN 
	           Rkriter=DeltaEnegry(NMO,1)/(3.D0*ExfZero)
	           DO WHILE(Henergy.GT.Rkriter)
		          Henergy=Henergy*0.97D0   
               ENDDO
	        ENDIF

	        ! мюундхл ндмнщкейрпнммсч щмепцхч якебю  
	        ! йнщттхжхемр мнплхпнбйх дкъ лнкейскъпмни нпахрюкх йпхрепхи онксвемхъ опюбхкэмни ндмнщкейрпнммни щмепцхх  
            RnormMO=RnormMOXF(NMO) 
            DetM=0.D0
            ExfM=0.D0
            ExfM(2)=ExfZero*(1.D0-Henergy) 
	        Exf=ExfZero*(1.D0-Henergy)   
      
            ! оюпюлерп бшундю хг жхйкю он пюяверс щмепцхх
            IparametrDeterEnergy=1
            NEpsIter=0
            INEXDSEXF=0
	        IcorenZZ=0
            ! оюпюлерп сйюгшбючыхи мю йнпемэ йнрнпши хыеряъ якебю
            IZNAK=1
            ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
            IndexVilky=0
            INDEXDF=0
            IIKK=0
            IndexCorrect=0
			IndexExitRerern=0
            WRITE(17,*) 'Energy left',Niteration,NMO
            !WRITE(*,*) 'EEEX NOW LEVOE',Exf 
            !READ(*,*) 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! жхйк нопедекемхе ндмнщкейрпнммни щмепцхх лнкейскъпмни нпахрюкх якебю!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO WHILE(IparametrDeterEnergy.EQ.1)
               NEpsIter=NEpsIter+1
	           ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	           call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
	           !WRITE(*,*) 'Kpoint',Kpoint
	           !WRITE(*,*) '4 Exf',Exf

	           ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	           call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
               !WRITE(*,*) '5 Exf',Exf
               ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
               call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
               !WRITE(*,*) '6 Exf',Exf  
               ! йнппейрхпсел гмювемхе щмепцхх Exf
               !write(*,*) 'eeee',Exf 
    	       call CorrectionEnergy(IZNAK,NEpsIter,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DetM,ExfM,Exf,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,EPS,RdeltaZnak,IndexVilky,INDEXDF,IIKK,IndexCorrect,IcorenZZ,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO,IndexExitRerern)  
               !write(*,*) '3 exf',Exf 
            ENDDO
      
	        ! сярюмюбкхбюел мюидем кх якебю йнпемэ
	        IF(IcorenZZ.EQ.1) THEN
                 ! йнпемэ меидем
		         INlevo=1
		         ExfLevo=Exf
               ELSE
		         ! йнмемэ ме мюидем
		         INlevo=0 
		         ExfLevo=0.D0
	        ENDIF 
         ENDIF
         
		
      ENDIF



      !WRITE(6,*) 'Energy One-Electron NE ODNORODNOE PRAVOE',Niteration,NMO,Exf,Kpoint
      !WRITE(*,*) 'EEEX NE ODNORODNOE PRAVOE',ExfSpravo
      !READ(*,*) 
              
	  ! опнбепъел мюкхвхе йнпмеи
      IF(INlevo.EQ.0.AND.INspravo.EQ.0) THEN
         ! йнпмх нрясрярбсчр
		 WRITE(7,*) ' ERROR. Solution not  Niteration=',Niteration,' NMO=',NMO 
	     WRITE(17,*) ' ERROR. Solution ont  Niteration=',Niteration,' NMO=',NMO  
	     
		 WRITE(*,*) ' ERROR. Solution not  Niteration=',Niteration,' NMO=',NMO  
		 !WRITE(*,*) 'Writing spectrum all solutions? (yes=1/no=2/exit=3):' 
		 !READ(*,*) IKLSDCalcul
		 ! тхйяхпсел йкчв
		 IKLSDCalcul=2 
		  
		 IF(IKLSDCalcul.EQ.3) THEN
            ! нярюмнбйю опнцпюллш
			STOP 
		 ENDIF
		 ! опнбепъел ядекюммши бшанп
         IF(IKLSDCalcul.EQ.1) THEN
		    WRITE(*,*) 'Writing Emin in Ry:'
			READ(*,*)  Emin  
		    ! опхбндхл гюбхяхлнярэ хмрецпюкю нпрнцнмюкэмнярх 
			! лнкейскъпмни нпахрюкх нр гмювемхи ндмнщкейрпнммни щмепцхх 
		    
			! оюпюлерп бшундю хг жхйкю он щмепцхх
            IparametrDeterEnergy=1
            ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
            Exf=ENN0  
			WRITE(18,*) 'Integral Ort(Energy)',Niteration,NMO
	        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !жхйк нопедекемхе ндмнщкейрпнммни щмепцхх лнкейскъпмни нпахрюкх!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO WHILE(IparametrDeterEnergy.EQ.1)
               ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	           call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
	 		   ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	           call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	           ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
               call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
               ! ондопнцпюллю гюохях гюбхяхлнярх хмрецпюкю нпрнцнмюкэмнярх нр ндмнщкейрпнммни щмепцхх
			   call SpectrumEnergy(Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,Exf,Emin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)  
			ENDDO 
    	 ENDIF

		 ! тхйяхпсел нрясрярбхе йнпмъ
		 IRETURN=1
		 ! хглемъел мнлеп хмрепюжхх
		 ! рюй йюй оепеундхл й янцкюянбюмхч б дпсцнл пефхле дкъ рни фе хрепюжхх 
         Niteration=Niteration-1
         return
	  ENDIF


	  IF(Niteration.EQ.15.AND.NMO.EQ.106) THEN 
         ! йнпмх нрясрярбсчр
		 !WRITE(7,*) ' ERROR. Solution not  Niteration=',Niteration,' NMO=',NMO 
	     !WRITE(17,*) ' ERROR. Solution ont  Niteration=',Niteration,' NMO=',NMO  
	     
		 !WRITE(*,*) ' ERROR. Solution not  Niteration=',Niteration,' NMO=',NMO  
		 !WRITE(*,*) 'Writing spectrum all solutions? (yes=1/no=2/exit=3):' 
		 !READ(*,*) IKLSDCalcul
		 ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
         Henergy=0.05D0
         ! нопедекъел ьюц дкъ дюммни ндмнщкейрпнммни тсмйжхх 
         DO IICFHK=1,NNEL
            ! хглемемхе ьюцю 
	        Henergy=Henergy/FLOAT(IICFHK+1)  
         ENDDO
		 ! тхйяхпсел йкчв
		 IKLSDCalcul=1 
		  
		 IF(IKLSDCalcul.EQ.3) THEN
            ! нярюмнбйю опнцпюллш
			STOP 
		 ENDIF
		 ! опнбепъел ядекюммши бшанп
         IF(IKLSDCalcul.EQ.1) THEN
		    Emin=0.01d0 
		    !WRITE(*,*) 'Writing Emin in Ry:'
			!READ(*,*)  Emin  
		    ! опхбндхл гюбхяхлнярэ хмрецпюкю нпрнцнмюкэмнярх 
			! лнкейскъпмни нпахрюкх нр гмювемхи ндмнщкейрпнммни щмепцхх 
		    
			! оюпюлерп бшундю хг жхйкю он щмепцхх
            IparametrDeterEnergy=1
            ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
            Exf=ENN0  
			WRITE(18,*) 'Integral Ort(Energy)',Niteration,NMO
	        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !жхйк нопедекемхе ндмнщкейрпнммни щмепцхх лнкейскъпмни нпахрюкх!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO WHILE(IparametrDeterEnergy.EQ.1)
               ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	           call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
	 		   ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	           call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	           ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
               call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
               ! ондопнцпюллю гюохях гюбхяхлнярх хмрецпюкю нпрнцнмюкэмнярх нр ндмнщкейрпнммни щмепцхх
			   call SpectrumEnergy(Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,Exf,Emin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)  
			ENDDO 
    	 ENDIF

		 !! тхйяхпсел нрясрярбхе йнпмъ
		 !IRETURN=1
		 !! хглемъел мнлеп хмрепюжхх
		 !! рюй йюй оепеундхл й янцкюянбюмхч б дпсцнл пефхле дкъ рни фе хрепюжхх 
         !Niteration=Niteration-1
         !return
	  ENDIF

      ! гюохяшбюел нрйкнмемхе ндмнщкейрпнммни щмепцхх онксвеммни пеьемхел ндмнпндмнцн х мендмнпндмнцн спюбмемхъ
      IF(INlevo.EQ.0) THEN
         ! кебши йнпемэ нрясрярбсер 
         DeltaEnegry(NMO,1)=DABS(ExfSpravo-ExfZero)
	     DeltaEnegry(NMO,2)=DABS(ExfSpravo-ExfZero) 
      ENDIF
	  IF(INspravo.EQ.0) THEN
         ! опюбши йнпемэ нрясрярбсер 
         DeltaEnegry(NMO,1)=DABS(ExfLevo-ExfZero)
	     DeltaEnegry(NMO,2)=DABS(ExfLevo-ExfZero) 
      ENDIF
      IF(INlevo.NE.0.AND.INspravo.NE.0) THEN
         DeltaEnegry(NMO,1)=DABS(ExfLevo-ExfZero)
	     DeltaEnegry(NMO,2)=DABS(ExfSpravo-ExfZero)    
      ENDIF

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! опнбндх юмюкхг онксвеммшу ндмнщкейрпнммшу щмепцхи !!!!!!!!!!
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(INlevo.EQ.0) THEN
         ! кебши йнпемэ нрясрярбсер 
         Exf=ExfSpravo 
      ENDIF
	  IF(INspravo.EQ.0) THEN
         ! опюбши йнпемэ нрясрярбсер 
         Exf=ExfLevo
      ENDIF
	  IF(INlevo.NE.0.AND.INspravo.NE.0) THEN
         Exf=ExfSpravo 
      ENDIF
   ENDIF


   ! пюявер ндмнщкейрпнммни щмепцхх гюбепьем
   EXXX1=Exf
   !WRITE(6,*) 'Energy One-Electron NE ODNORODNOE PRAVOE',Niteration,NMO,Exf,Kpoint
   !WRITE(*,*) 'EEEX FINAL',Exf 
   !READ(*,*) 

     
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! щрюо 3. мюундхл лнкейскъпмсч нпахрюкэ дкъ мюидемни ндмнщкейрпнммни щмепцхх уюпрпх-тнйю (оюпжхюкэмше цюплнмхйх)! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! нопедекъел рхо пюяверю
   IntParXXZ=2
   IF(Niteration.EQ.1) THEN  
        IF(IKLExchange.EQ.0) THEN
           IntParXXZ=1
        ENDIF
      ELSE
	    IF(IKLExchange.EQ.0) THEN
           IntParXXZ=1
        ENDIF 
   ENDIF 
     
   
   ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
   call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
   !WRITE(6,*) 'Kpoint=',Kpoint,'IndexM=',IndexM

   ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
   call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
       
   ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх
   call CalculationParameterProgonky(IntParXXZ,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
 
   ! мюундхл гмювемхъ цюплнмхй лнкейскъпмни нпахрюкх б рнвйе яьхбйх
   IF(IntParXXZ.EQ.1) THEN   
       call CalculationValuesGarmonikMO(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Vin,FunMO,AZ,BZ,E,AX,BX)
      ELSE
       call CalculationValuesGarmonikMOFull(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,FunMO,E,AX,BX,IXM) 
   ENDIF 
   
   ! нясыеярбкъел бняярюмнбкемхе гмювемхи тсмйжхх б бяеу рнвйюу нр рнвйх яьхбйх х мнплхпсел онксвеммсч тсмйжхч
   call FunctionRecovery(IntParXXZ,Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,BXL,CX,VNout,UNout,VNin,UNin,RO1X,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO) 
   ! гюохяшбюел мнплхпнбнвмши йнщттхжхемр (йпхрепхи нопедекемхъ опюбхкэмни ндмнщкейрпнммни щмепцххх)
   ! опнбепъел рхо бшундю
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.1.AND.IntParXXZ.EQ.1) THEN
         !WRITE(6,*) 'WRITING NORM MO=',NMO
	     RnormMOXF(NMO)=FunMO(NumeroMin,2)
      ENDIF
      IF(Niteration.EQ.3.AND.IntParXXZ.EQ.1) THEN
         !WRITE(6,*) 'WRITING NORM MO=',NMO
	     RnormMOXF(NMO)=FunMO(NumeroMin,2)
      ENDIF
   ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! щрюо 4. гюохяэ онксвеммни тсмйжхх дюммни хрепюжхх (опнжедспю опнбепйх янцкюянбюммнярх)!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   IF(Niteration.EQ.1) THEN
        ! оюпюлерп пюяянцкюянбюммнярх
		T=1.D0-RFunMO(NMO,NumeroMin,1)/Exf  
		TETA(NMO)=T
		TETMAX=TETMAX+DABS(T)
        ! гюохяшбюел ндмнщкейрпнммше щмепцхх онксвеммше дкъ дюкэмеиьецн сверю  
        E3(NMO)=RFunMO(NMO,NumeroMin,1)  
        E5(NMO)=Exf
		E4(NMO)=Exf

		! гюохяшбюел б люяяхб гмювемхе ндмнщкейрпнммни щмепцхх
		! ю рюйфе мнплхпнбнвмши йнщттхжхемр
        DO IIX=1,NumeroMaxLimid  !NumeroMax
		   RFunMO(NMO,IIX,1)=Exf
		   RFunMO(NMO,IIX,2)=FunMO(NumeroMin,2)
           FunMOOld(NMO,IIX,1)=Exf
           FunMOOld(NMO,IIX,2)=FunMO(NumeroMin,2)
        ENDDO
		! гюмскъел опефде цюплнмхйх мюундъыхеяъ мхфе NumeroMin
		DO IIO=1,Npoint 
		   DO IIX=1,NumeroMin
              RFunMO(NMO,IIX,IIO+2)=0.D0
		   ENDDO
        ENDDO
       
        ! гюохяшбюел онксвеммсч лнкейскъпмсч нпахрюкэ 
		DO IIO=1,Npoint 
		   DO IIX=NumeroMin,NumeroMaxLimid   !NumeroMax
		      RFunMO(NMO,IIX,IIO+2)=FunMO(IIX,IIO+2)
			  FunMOOld(NMO,IIX,IIO+2)=FunMO(IIX,IIO+2)
		   ENDDO
        ENDDO

		! гюохяшбюел лнкейскъпмше нпахрюкх щйбюбнкемрмше дюммни
        DO IIXXYYT=1,NumreISzam(NMO)
           ! оюпюлерп пюяянцкюянбюммнярх
		   TETA(ISCalculZZ(NMO,IIXXYYT+1))=T
		   TETMAX=TETMAX+DABS(T)
           ! гюохяшбюел ндмнщкейрпнммше щмепцхх онксвеммше дкъ дюкэмеиьецн сверю  
           E3(ISCalculZZ(NMO,IIXXYYT+1))=RFunMO(ISCalculZZ(NMO,IIXXYYT+1),NumeroMin,1)  
           E4(ISCalculZZ(NMO,IIXXYYT+1))=Exf
           ! гюохяшбюел б люяяхб гмювемхе ндмнщкейрпнммни щмепцхх
		   ! ю рюйфе мнплхпнбнвмши йнщттхжхемр
           DO IIX=1,NumeroMaxLimid  !NumeroMax
		      RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)=Exf
		      RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)=FunMO(NumeroMin,2)
			  FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)=Exf
              FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)=FunMO(NumeroMin,2)
           ENDDO
           ! гюмскъел опефде цюплнмхйх мюундъыхеяъ мхфе NumeroMin
		   DO IIO=1,Npoint 
		      DO IIX=1,NumeroMin
                 RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=0.D0
		      ENDDO
           ENDDO
           ! гюохяшбюел онксвеммсч лнкейскъпмсч нпахрюкэ 
		   DO IIO=1,Npoint 
		      DO IIX=NumeroMin,NumeroMaxLimid  !NumeroMax
		         RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=FunMO(IIX,IIO+2)
		         FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=FunMO(IIX,IIO+2) 
			  ENDDO
           ENDDO
        ENDDO

      ELSE

        ! оюпюлерп пюяянцкюянбюммнярх
	    T=1.D0-RFunMO(NMO,NumeroMin,1)/Exf  
	
	    ! опнбепъел ялеьхбюмхе нярясрярбсер хкх мер
	    IF(NLST.EQ.0) THEN
	       T=0.D0
	    ENDIF 
       	! ярюмдюпрмши пефхл янцкюянбюмхъ
        IF(IreshimSO.EQ.0) THEN
	       ! сярюмюбкхбюел йнщттхжхемр ялеьхбюмхъ
	       IF((Niteration-2).LT.0) THEN
               C=0.D0
              ELSE
               IF((-1)**Niteration.LE.0) THEN
		          C=AnalysCoffMixing(Niteration,E5(NMO),E4(NMO),Exf)
		         ELSE
			      EZR=Exf-E4(NMO)-RFunMO(NMO,NumeroMin,1)+E3(NMO)
                  IF(DABS(EZR).LT.1.D-8) THEN
			         EZR=1.D0
                  ENDIF
                  C=(Exf-E4(NMO))/EZR
				  IF(C.LT.0.D0.OR.C.GT.1.D0) THEN 
			         C=AnalysCoffMixing(Niteration,E5(NMO),E4(NMO),Exf)
			      ENDIF
		       ENDIF     
           ENDIF
	    ENDIF

        ! пефхл янцкюянбюмхъ я нцпюмхвемхел
        IF(IreshimSO.EQ.1) THEN
	       ! сярюмюбкхбюел йнщттхжхемр ялеьхбюмхъ
	       IF((Niteration-2).LT.0) THEN
               C=0.D0
              ELSE
               IF((-1)**Niteration.LE.0) THEN
		          C=AnalysCoffMixing(Niteration,E5(NMO),E4(NMO),Exf)
		         ELSE
			      EZR=Exf-E4(NMO)-RFunMO(NMO,NumeroMin,1)+E3(NMO)
                  IF(DABS(EZR).LT.1.D-8) THEN
			         EZR=1.D0
                  ENDIF
                  C=(Exf-E4(NMO))/EZR
				  IF(C.LT.0.D0.OR.C.GT.1.D0) THEN 
			         C=AnalysCoffMixing(Niteration,E5(NMO),E4(NMO),Exf)
			      ENDIF
		       ENDIF 
			   IF(C.LT.0.7D0) THEN
			      C=0.7D0
			   ENDIF    
           ENDIF
	    ENDIF

		
		! сяхкеммши пефхл янцкюянбюмхъ 
        IF(IreshimSO.EQ.2) THEN
	       ! сярюмюбкхбюел йнщттхжхемр ялеьхбюмхъ
	       IF((Niteration-2).LT.0) THEN
               C=0.D0
              ELSE
               IF((-1)**Niteration.LE.0) THEN
		           C=0.9D0
		          ELSE
				   C=AnalysCoffMixingStrong(Niteration,E5(NMO),E4(NMO),Exf) 
		       ENDIF     
           ENDIF
        ENDIF

		! сяхкеммши пефхл янцкюянбюмхъ 2
        IF(IreshimSO.EQ.3) THEN
	       ! сярюмюбкхбюел йнщттхжхемр ялеьхбюмхъ
	       IF((Niteration-2).LT.0) THEN
               C=0.D0
              ELSE
               IF((-1)**Niteration.LE.0) THEN
		           C=0.97D0
		          ELSE
				   C=AnalysCoffMixingStrong2(Niteration,E5(NMO),E4(NMO),Exf) 
		       ENDIF     
           ENDIF
        ENDIF
        
		! сяхкеммши пефхл янцкюянбюмхъ 3
		IF(IreshimSO.EQ.4) THEN
	       ! сярюмюбкхбюел йнщттхжхемр ялеьхбюмхъ
	       IF((Niteration-2).LT.0) THEN
               C=0.D0
              ELSE
               IF((-1)**Niteration.LE.0) THEN
		           C=0.98D0
		          ELSE
				   C=AnalysCoffMixingStrong3(Niteration,E5(NMO),E4(NMO),Exf) 
		       ENDIF     
           ENDIF
        ENDIF
        
	    E3(NMO)=RFunMO(NMO,NumeroMin,1)  
		E5(NMO)=E4(NMO)
        E4(NMO)=Exf
        IF(NLST.EQ.0) THEN
	       ! ялеьхбюмхъ мер 
	       C=0.D0
	    ENDIF  

      
	    ! гюохяэ ндмнщкейрпнммнцн гмювемхъ щмепцхх "мнбни" дюммни хрепюжхх
        ! гюохяшбюел б люяяхб гмювемхе ндмнщкейрпнммни щмепцхх
	    DO IIX=1,NumeroMaxLimid !NumeroMax
		   ! ГЮОХЯШБЮЕЛ ОПЕДШДСЫСЧ ХРЕПЮЖХЧ
		   FunMOOld(NMO,IIX,1)=RFunMO(NMO,IIX,1)
		   ! ГЮОХЯШБЮЕЛ МНПЛС ОПЕДЕДСЫЕИ ХРЕПЮЖХХ
           FunMOOld(NMO,IIX,2)=RFunMO(NMO,IIX,2)
		   ! ГЮОХЯШБЮЕЛ ДЮММСЧ ХРЕПЮЖХЧ
		   FunMONew(NMO,IIX,1)=Exf
		   RFunMO(NMO,IIX,1)=C*RFunMO(NMO,IIX,1)+(1.D0-C)*Exf
        ENDDO
        
	    ! гюохяшбюел онксвеммсч лнкейскъпмсч нпахрюкэ
	    DO IIO=1,Npoint  
	       DO IIX=NumeroMin,NumeroMaxLimid   !NumeroMax
	          BBL=RFunMO(NMO,IIX,IIO+2)-FunMO(IIX,IIO+2)
              IF(Niteration.EQ.3.AND.IndexExit.EQ.0) THEN
                 BBL=RFunMO(NMO,IIX,IIO+2)*DSQRT(FunMO(IIX,2)/RFunMO(NMO,IIX,2))-FunMO(IIX,IIO+2)
			  ENDIF
	  	      IF(DABS(BBL).GT.DABS(T)) THEN
			     MAXCNV=IIO
              ENDIF
			  IF(DABS(BBL).GT.DABS(T)) THEN
			     T=BBL
              ENDIF 
              ! ГЮОХЯШБЮЕЛ ОПЕДШДСЫСЧ ХРЕПЮЖХЧ
		      FunMOOld(NMO,IIX,IIO+2)=RFunMO(NMO,IIX,IIO+2)
		      ! ГЮОХЯШБЮЕЛ ДЮММСЧ ХРЕПЮЖХЧ
		      FunMONew(NMO,IIX,IIO+2)=FunMO(IIX,IIO+2)
			  IF(Niteration.EQ.3.AND.IndexExit.EQ.0) THEN
			       RFunMO(NMO,IIX,IIO+2)=C*RFunMO(NMO,IIX,IIO+2)*DSQRT(FunMO(IIX,2)/RFunMO(NMO,IIX,2))+(1.D0-C)*FunMO(IIX,IIO+2)
                 ELSE
				   RFunMO(NMO,IIX,IIO+2)=C*RFunMO(NMO,IIX,IIO+2)+(1.D0-C)*FunMO(IIX,IIO+2)  
		      ENDIF
		   ENDDO
        ENDDO

        ! нясыеярбкъел пюявер мнплш онксвеммни тсмйжхх
        RRnormcoff=CalculationNormalizationConstaResult(NMO,NumeroMin,NumeroMaxLimid,Npoint,H,RO1X,RFunMO) 
        
	    ! гюохяшбюел  мнплхпнбнвмши йнщттхжхемр
	    DO IIX=1,NumeroMaxLimid !NumeroMax
	       RFunMO(NMO,IIX,2)=RRnormcoff
        ENDDO
		
		! нясыеярбкъел пюявер мнплш дкъ тсмйжхи
        RRnormcoff=CalculationNormalizationConstaResult(NMO,NumeroMin,NumeroMaxLimid,Npoint,H,RO1X,FunMONew) 
        
        ! гюохяшбюел  мнплхпнбнвмши йнщттхжхемр
	    DO IIX=1,NumeroMaxLimid !NumeroMax
	       FunMONew(NMO,IIX,2)=RRnormcoff
        ENDDO
        


	    !WRITE(200,*) 'FUNN',Niteration,NMO,RFunMO(NMO,NumeroMin,1)
	    !DO IIO=1,Npoint
        !   WRITE(100*Niteration+NMO,35678) RR(IIO),(RFunMO(NMO,IIX,IIO+2)/DSQRT(RFunMO(NMO,IIX,2)*RO1X(IIO)),IIX=1,NumeroMax) !NumeroMin,NumeroMax)
        !ENDDO
        
	    TETA(NMO)=T
	    TETMAX=TETMAX+DABS(T)

		! гюохяшбюел лнкейскъпмше нпахрюкх щйбюбнкемрмше дюммни
        DO IIXXYYT=1,NumreISzam(NMO)
           ! оюпюлерп пюяянцкюянбюммнярх
		   TETA(ISCalculZZ(NMO,IIXXYYT+1))=T
		   TETMAX=TETMAX+DABS(T)
           ! гюохяшбюел ндмнщкейрпнммше щмепцхх онксвеммше дкъ дюкэмеиьецн сверю  
           E3(ISCalculZZ(NMO,IIXXYYT+1))=RFunMO(ISCalculZZ(NMO,IIXXYYT+1),NumeroMin,1)  
           E4(ISCalculZZ(NMO,IIXXYYT+1))=Exf
		   ! гюохяэ ндмнщкейрпнммнцн гмювемхъ щмепцхх "мнбни" дюммни хреппюжхх
           ! гюохяшбюел б люяяхб гмювемхе ндмнщкейрпнммни щмепцхх
		   DO IIX=1,NumeroMaxLimid !NumeroMax
              ! ГЮОХЯШБЮЕЛ ОПЕДШДСЫСЧ ХРЕПЮЖХЧ
		      FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)=RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)
		      ! ГЮОХЯШБЮЕЛ МНПЛС ОПЕДЕДСЫЕИ ХРЕПЮЖХХ
              FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)=RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)
		      ! ГЮОХЯШБЮЕЛ ДЮММСЧ ХРЕПЮЖХЧ
		      FunMONew(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)=Exf
			  RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)=C*RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,1)+(1.D0-C)*Exf
		   ENDDO
           ! гюохяшбюел онксвеммсч лнкейскъпмсч нпахрюкэ 
		   DO IIO=1,Npoint 
		      DO IIX=NumeroMin,NumeroMaxLimid !NumeroMax
			     ! ГЮОХЯШБЮЕЛ ОПЕДШДСЫСЧ ХРЕПЮЖХЧ
		         FunMOOld(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)  
		         ! ГЮОХЯШБЮЕЛ ДЮММСЧ ХРЕПЮЖХЧ
		         FunMONew(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=FunMO(IIX,IIO+2)
				 IF(Niteration.EQ.3.AND.IndexExit.EQ.0) THEN
				      RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=C*RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)*DSQRT(FunMO(IIX,2)/RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,2))+(1.D0-C)*FunMO(IIX,IIO+2)
		            ELSE
					  RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)=C*RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,IIO+2)+(1.D0-C)*FunMO(IIX,IIO+2)  
                 ENDIF
			  ENDDO
           ENDDO
	       ! гюохяшбюел  мнплхпнбнвмши йнщттхжхемр
		   DO IIX=1,NumeroMaxLimid !NumeroMax
              FunMONew(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)=FunMONew(NMO,IIX,2)
		      RFunMO(ISCalculZZ(NMO,IIXXYYT+1),IIX,2)=RFunMO(NMO,IIX,2)
           ENDDO
        ENDDO
        
   ENDIF

   ! бшдювю опнлефсрнвмшу пегскэрюрнб янцкюянбюмхъ
   !WRITE(7,89787) Niteration,NMO,Exf,ExfLevo,ExfZero,ExfSpravo,C 
   ! мюидеммюъ ндмнщкейрпнммюъ щмепцхъ 
   RparamterCalcul(1)=Exf
   ! йнпмх спюбмемхъ ндмнпндмни х мендмнпндмни яхярелш
   RparamterCalcul(2)=ExfLevo
   RparamterCalcul(3)=ExfZero
   RparamterCalcul(4)=ExfSpravo
   ! йнщттхжхемр ялеьхбюмхъ
   RparamterCalcul(5)=C 
   ! тхйяхпсел мюкхвхе йнпмъ
   IRETURN=0
   
   !WRITE(*,*) 'PARAMETR'
   !WRITE(*,*) Niteration,C,TETA(NMO),TETMAX
   !READ(*,*)

   !WRITE(*,*) 'DD'
   !READ(*,*)
   !STOP


   ! сдюкъел люяяхбш хг оюлърх 
   deallocate(IndexGranyMin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "IndexGranyMin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif  
   deallocate(EEGranymin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "EEGranymin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif  
   deallocate(IndexMin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "IndexMin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(DDMmin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "DDMmin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(EEMmin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "EEMmin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(IXM,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "IXM" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(AX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "AX" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(BX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "BX" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif
   deallocate(BXL,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "BXL" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif  
   deallocate(CX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "CX" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(DX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "DX" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(UNout,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "UNout" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(UNin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "UNin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(E,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "E" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(AZ,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "AZ" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(BZ,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "BZ" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(FZ,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "FZ" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif
   deallocate(DZ,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "DZ" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(FX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "FX" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(Uout,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "Uout" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(Uin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "Uin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(VNout,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "VNout" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(VNin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "VNin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(FunMO,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "FunMO" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(A,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "A" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif
   deallocate(B,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "B" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(AA,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "AA" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif 
   deallocate(BB,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "BB" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif
   deallocate(Vout,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "Vout" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif  
   deallocate(Vin,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'SolutionSystemDifferentialEquations'
      write(*,*) 'THE MASSIV "Vin" IS NOT REMOVED FROM MEMORY'
      read(*,*)
	  stop 
   endif  

   return
 end subroutine SolutionSystemDifferentialEquations

 
 
!! SUBPROGRAMME OF THE CALCULATION OF THE NORM OF A FUNCTION OF THE RELEVANT DATA ON THE ENERGY
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! Rcof-EQUITY RATIO H * H / 12
═!! H-STEP
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (IndexM, IndexM, Npoint) -MASSIVE A
═!! BB (IndexM, IndexM, Npoint) -MASSIVE B
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! Exf-SINGLE-ELECTRONIC ENERGY OF HARTRI-FOCA (MOLECULAR ORBITAL)
═!! RnormFunXF-NORMING CONSTANT OF THE HARTRI-FOCO FUNCTION MO
 subroutine CalculationNormFunctionEnergyXF(Exf,RnormFunXF,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
   implicit none
   integer::NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax
   real(8)::Exf,RnormFunXF,Rcof,H
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
    
   ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1X
   !  нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
   call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
	       
   ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
   call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
   ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
   call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                 
   ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
   call CalculationNormFunE(RnormFunXF,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
   return
 end subroutine CalculationNormFunctionEnergyXF 

 ! ондопнцпюллю йнппейрхпсчыюъ ьюц 
 subroutine CorrectHRight(INDEXZZSFA,IndexIIUMN,Henergy,RcoffChengH,ProDet)
   implicit none
   integer::INDEXZZSFA,IndexIIUMN
   real(8)::Henergy,RcoffChengH
   real(8),dimension(:)::ProDet
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)::RcoffChengHOld

   ! йкчв сйюгшбючыхи менаундхлнярэ бнгбпюыемхъ мю ндмс рнвйс мюгюд
   INDEXZZSFA=0 

       
   ! бнгбпюыюел б хяундмне янярнъмхе ьюц
   Henergy=Henergy*RcoffChengH*(10.D0)**IndexIIUMN 

  
     
   ! опнбепъел еярэ хглемемхе гмюйю опнхгбндмни хкх мер 
   IF(ProDet(7).LT.0.D0) THEN
        ! гмюй опнхгбндмни ме хглемхкяъ
      
	    ! пюявхршбюел мнбши йнщттхжхемр
        RcoffChengH=DABS(ProDet(7)/ProDet(6))   
     ELSE
	    ! гмюй опнхгбндмни хглемхкяъ 
		! тхйяхпсел хглемемхе гмюйю опнхгбндмни
        INDEXZZSFA=1
		! сбекхвхбюел хмдейя слемэьемхъ ьюцю
        IndexIIUMN=IndexIIUMN+1
		! пюявхршбюел мнбши йнщттхжхемр (я свернл бнгбпюыемхъ мю дбе рнвйх мюгюд)
        RcoffChengH=DABS(ProDet(5)/ProDet(4))  		 
   ENDIF

   ! опнбепъел йпхрепхи
   IF(RcoffChengH.GT.2.1D0) THEN
       RcoffChengH=RcoffChengH*2.D0
   ENDIF 
   ! нцпюмхвемхе он бекхвхме слемэьемхъ ьюцю
   ! мер пефхлю слемэьемхъ ьюцю б 10*IndexIIUMN пюг
   IF(RcoffChengH.GT.200.D0.AND.IndexIIUMN.EQ.0) THEN
      RcoffChengH=200.D0
   ENDIF 
   ! еярэ пефхл слемэьемхъ ьюцю б 10*IndexIIUMN пюг 
   IF(RcoffChengH.GT.400.D0.AND.IndexIIUMN.NE.0) THEN
      RcoffChengH=400.D0
   ENDIF 

   IF(RcoffChengH.LT.1.D0) THEN
      RcoffChengH=1.D0
   ENDIF

   ! хглемъел ьюц
   Henergy=Henergy/(RcoffChengH*(10.D0)**IndexIIUMN)  
   
  
   return
 end subroutine CorrectHRight


 ! ондопнцпюллю йнппейрхпсчыюъ ьюц 
 subroutine CorrectHLeft(INDEXZZSFA,IndexIIUMN,Henergy,RcoffChengH,ProDet)
   implicit none
   integer::INDEXZZSFA,IndexIIUMN
   real(8)::Henergy,RcoffChengH
   real(8),dimension(:)::ProDet
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)::RcoffChengHOld

   ! йкчв сйюгшбючыхи менаундхлнярэ бнгбпюыемхъ мю ндмс рнвйс мюгюд
   INDEXZZSFA=0 
   
   ! бнгбпюыюел б хяундмне янярнъмхе ьюц
   Henergy=Henergy*RcoffChengH*(10.D0)**IndexIIUMN 
     
   ! опнбепъел еярэ хглемемхе гмюйю опнхгбндмни хкх мер 
   IF(ProDet(7).GT.0.D0) THEN
        ! гмюй опнхгбндмни ме хглемхкяъ
	    ! пюявхршбюел мнбши йнщттхжхемр
        RcoffChengH=DABS(ProDet(7)/ProDet(6))   
     ELSE
	    ! гмюй опнхгбндмни хглемхкяъ 
		! тхйяхпсел хглемемхе гмюйю опнхгбндмни
        INDEXZZSFA=1
		! сбекхвхбюел хмдейя слемэьемхъ ьюцю
        IndexIIUMN=IndexIIUMN+1
		! пюявхршбюел мнбши йнщттхжхемр (я свернл бнгбпюыемхъ мю дбе рнвйх мюгюд)
        RcoffChengH=DABS(ProDet(5)/ProDet(4))  		 
   ENDIF

   ! опнбепъел йпхрепхи
   IF(RcoffChengH.GT.2.1D0) THEN
       RcoffChengH=RcoffChengH*2.D0
   ENDIF 
   ! нцпюмхвемхе он бекхвхме слемэьемхъ ьюцю
   ! мер пефхлю слемэьемхъ ьюцю б 10*IndexIIUMN пюг
   IF(RcoffChengH.GT.200.D0.AND.IndexIIUMN.EQ.0) THEN
      RcoffChengH=200.D0
   ENDIF 
   ! еярэ пефхл слемэьемхъ ьюцю б 10*IndexIIUMN пюг 
   IF(RcoffChengH.GT.400.D0.AND.IndexIIUMN.NE.0) THEN
      RcoffChengH=400.D0
   ENDIF 

   IF(RcoffChengH.LT.1.D0) THEN
      RcoffChengH=1.D0
   ENDIF

   ! хглемъел ьюц
   Henergy=Henergy/(RcoffChengH*(10.D0)**IndexIIUMN)  

   return
 end subroutine CorrectHLeft 


  ! ондопнцпюллю йнппейрхпсчыюъ ьюц 
 subroutine CorrectHRightD(IndexSUMZKX,INDEXZZSFA,IndexIIUMN,Henergy,RcoffChengH,ProDet)
   implicit none
   integer::IndexSUMZKX,INDEXZZSFA,IndexIIUMN
   real(8)::Henergy,RcoffChengH
   real(8),dimension(:)::ProDet
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)::RcoffChengHOld

   ! йкчв сйюгшбючыхи менаундхлнярэ бнгбпюыемхъ мю ндмс рнвйс мюгюд
   INDEXZZSFA=0 

       
   ! бнгбпюыюел б хяундмне янярнъмхе ьюц
   Henergy=Henergy*RcoffChengH*(10.D0)**IndexIIUMN 

  
     
   ! опнбепъел еярэ хглемемхе гмюйю опнхгбндмни хкх мер 
   IF(ProDet(7).LT.0.D0) THEN
        ! гмюй опнхгбндмни ме хглемхкяъ
	    ! пюявхршбюел мнбши йнщттхжхемр
        RcoffChengH=DABS(ProDet(7)/ProDet(6))   
     ELSE
	    ! гмюй опнхгбндмни хглемхкяъ 
		! тхйяхпсел хглемемхе гмюйю опнхгбндмни
        INDEXZZSFA=1
		! сбекхвхбюел хмдейя слемэьемхъ ьюцю
        IndexIIUMN=IndexIIUMN+1
		! опнбепъел мнлеп хрепюжхх 
		IF(IndexSUMZKX.GT.2) THEN
		    ! пюявхршбюел мнбши йнщттхжхемр (я свернл бнгбпюыемхъ мю дбе рнвйх мюгюд)
            RcoffChengH=DABS(ProDet(5)/ProDet(4))  		 
           ELSE
            ! оняйнкэйс вхякн хрепюжхи ме днярюрнвмн дкъ пюяверю йнщттхжхемрю тхйяхпсел 
			RcoffChengH=2.D0 
        ENDIF
   ENDIF

   ! опнбепъел йпхрепхи
   IF(RcoffChengH.GT.2.1D0) THEN
       RcoffChengH=RcoffChengH*2.D0
   ENDIF 
   ! нцпюмхвемхе он бекхвхме слемэьемхъ ьюцю
   ! мер пефхлю слемэьемхъ ьюцю б 10*IndexIIUMN пюг
   IF(RcoffChengH.GT.200.D0.AND.IndexIIUMN.EQ.0) THEN
      RcoffChengH=200.D0
   ENDIF 
   ! еярэ пефхл слемэьемхъ ьюцю б 10*IndexIIUMN пюг 
   IF(RcoffChengH.GT.400.D0.AND.IndexIIUMN.NE.0) THEN
      RcoffChengH=400.D0
   ENDIF 

   IF(RcoffChengH.LT.1.D0) THEN
      RcoffChengH=1.D0
   ENDIF

  
   ! хглемъел ьюц
   Henergy=Henergy/(RcoffChengH*(10.D0)**IndexIIUMN)  
   
 
   return
 end subroutine CorrectHRightD


 ! ондопнцпюллю йнппейрхпсчыюъ ьюц 
 subroutine CorrectHLeftD(IndexSUMZKX,INDEXZZSFA,IndexIIUMN,Henergy,RcoffChengH,ProDet)
   implicit none
   integer::IndexSUMZKX,INDEXZZSFA,IndexIIUMN
   real(8)::Henergy,RcoffChengH
   real(8),dimension(:)::ProDet
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)::RcoffChengHOld

   ! йкчв сйюгшбючыхи менаундхлнярэ бнгбпюыемхъ мю ндмс рнвйс мюгюд
   INDEXZZSFA=0 
   
   ! бнгбпюыюел б хяундмне янярнъмхе ьюц
   Henergy=Henergy*RcoffChengH*(10.D0)**IndexIIUMN 
     
   ! опнбепъел еярэ хглемемхе гмюйю опнхгбндмни хкх мер 
   IF(ProDet(7).GT.0.D0) THEN
        ! гмюй опнхгбндмни ме хглемхкяъ
	    ! пюявхршбюел мнбши йнщттхжхемр
        RcoffChengH=DABS(ProDet(7)/ProDet(6))   
     ELSE
	    ! гмюй опнхгбндмни хглемхкяъ 
		! тхйяхпсел хглемемхе гмюйю опнхгбндмни
        INDEXZZSFA=1
		! сбекхвхбюел хмдейя слемэьемхъ ьюцю
        IndexIIUMN=IndexIIUMN+1
		! опнбепъел мнлеп хрепюжхх 
		IF(IndexSUMZKX.GT.2) THEN
		    ! пюявхршбюел мнбши йнщттхжхемр (я свернл бнгбпюыемхъ мю дбе рнвйх мюгюд)
            RcoffChengH=DABS(ProDet(5)/ProDet(4))  		 
           ELSE
            ! оняйнкэйс вхякн хрепюжхи ме днярюрнвмн дкъ пюяверю йнщттхжхемрю тхйяхпсел 
			RcoffChengH=2.D0 
        ENDIF	 
   ENDIF

   ! опнбепъел йпхрепхи
   IF(RcoffChengH.GT.2.1D0) THEN
       RcoffChengH=RcoffChengH*2.D0
   ENDIF 
   ! нцпюмхвемхе он бекхвхме слемэьемхъ ьюцю
   ! мер пефхлю слемэьемхъ ьюцю б 10*IndexIIUMN пюг
   IF(RcoffChengH.GT.200.D0.AND.IndexIIUMN.EQ.0) THEN
      RcoffChengH=200.D0
   ENDIF 
   ! еярэ пефхл слемэьемхъ ьюцю б 10*IndexIIUMN пюг 
   IF(RcoffChengH.GT.400.D0.AND.IndexIIUMN.NE.0) THEN
      RcoffChengH=400.D0
   ENDIF 

   IF(RcoffChengH.LT.1.D0) THEN
      RcoffChengH=1.D0
   ENDIF
  
   ! хглемъел ьюц
   Henergy=Henergy/(RcoffChengH*(10.D0)**IndexIIUMN)  

 
   return
 end subroutine CorrectHLeftD 
 
 

 !! SUBPROGRAMME SEARCH FOR THE OWN ENERGY OF THE HOMOGENEOUS SYSTEM OF EQUATIONS BY "STANDARD" METHOD
═!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
═!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
═!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
═!! Niteration-NUMBER OF INTERACTION
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! NNEL- MAIN QUANTUM NUMBER
═!! Npoint-NUMBER OF POINTS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! Rcof-EQUITY RATIO H * H / 12
═!! EPS-ACCURACY
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! RO1X (Npoint) -MASSIVE FIRST DERIVATIVE NEW VARIABLE FOR THE OLD VARIABLE
═!! VNout, UNout, VNin, UNin, E, AZ, FZ, DZ, BZ, AX, BX-AUXILIARY MASSIVE
═!! ENN0-LOWER ENERGY BORDER
═!! ExfZero-SOLUTION OF THE HOMOGENEOUS SYSTEM OF EQUATIONS ONE-ELECTRONIC ENERGY
═!! NumbreRegion (NumbreMO) -MASSIVE NUMBER OF REGIONS IN WHICH A ROOT REDUCTION IS CARRIED OUT
═!! EnergyZeroZ (NumbreMO, NNEL + 1,2) -MASSIVE OF THE ROOTS OF A HOMOGENEOUS EQUATION
═!! EnergyRegion (NumbreMO, NNEL + 1,2) -MASSIVE OF AREAS IN WHICH STEP DECREASES
 subroutine SearchEnergySolutionHomogeneousSystemStandard(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfZero,NumbreRegion,EnergyZeroZ,EnergyRegion) 
   implicit none
   integer::IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM
   real(8)::Rcof,EPS,ENN0,ExfZero
   integer,dimension(:)::NumbreRegion
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,EnergyZeroZ,EnergyRegion 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IparametrDeterEnergy,NEpsIter,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexRegion,IndexSmena,IIFG,IKLKORX,IOPZ,IDFJK,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Henergy,Exf,Delta,RcoffInterval,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(10)::DetM,ExfM,ProDet 
    
   ! ярюмдюпрмши лернд пюяверю
   ! гюмскъел оепед пюявернл
   DetM=0.D0
   ExfM=0.D0
   ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
   Exf=ENN0  
   ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
   Henergy=0.1D0
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.1) THEN
         Henergy=Henergy/10.D0
      ENDIF
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=2
   ENDIF
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
      IF(Niteration.EQ.3) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=3
   ENDIF
   ! оюпюлерп бшундю хг жхйкю  
   IparametrDeterEnergy=1
   ! хмдейя жхйкю 
   NEpsIter=0
   ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
   IndexVilky=0
   INDEXDF=0
   IIKK=0
   IndexCorrect=0
   IndexRegion=1
   IndexSmena=0
   IKLKORX=0
   IOPZ=0
   IKLRXX=0
   INDEXAD=0
   INDEXRShag=0
   IIHD=0
   IIOXDPZ=0
   INDEXPeresc=0
   DeltaEnerPers=0.D0
   DeltaFunPers=0.D0
   RcoffChengH=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюяверю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
      ! хыел яннрберярбсчыхе йнпмх Exf
      call SearchSolution(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,EnergyZeroZ,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)  
      ! нясыеярбкъел йнмрпнкэ гю ьюцюл 
	  IF(Niteration.GT.IDFJK) THEN
	     call ControlH(NMO,IndexSmena,NumbreRegion(NMO),IndexRegion,EnergyRegion,Exf,Henergy)
      ENDIF 
   ENDDO

   ! гюохяшбюел вхякн йнпмеи пюяялнрпеммшу б унде онхяйю
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.1) THEN
         NumbreRegion(NMO)=IndexCorrect
      ENDIF
   ENDIF
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.EQ.2) THEN
         NumbreRegion(NMO)=IndexCorrect
      ENDIF
   ENDIF
   ! тнплхпсел накюярх цде ьюц асдер слемэьхм
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.GT.1) THEN
         ! жхйк он мюидеммшл йнпмъл
	     DO IIFG=1,NumbreRegion(NMO)
            Delta=DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
            ! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF

		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,IIFG,2)=EnergyZeroZ(NMO,IIFG,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,IIFG,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,IIFG,1)=EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta 
			ENDIF
			! опнбепъел врнаш накюярх йнпмеи ме оепеяхйюкхяэ мювхмюъ я брнпнцн йнпмъ
			IF(IIFG.GT.1) THEN
			   ! опнбепъел врнаш опюбюъ цпюмхжю дюммни накюярх ме ашкю анкэье кебни цпюмхжш опедшдсыеи накюярх
               IF(EnergyRegion(NMO,IIFG,2).GT.EnergyRegion(NMO,IIFG-1,1)) THEN
                  ! опхпюбмхбюел цпюмхжш опюбюъ цпюмхжю дюммни накюярх асдер янбоюдюрэ я кебни цпюмхжни опедшдсыеи накюярх
				  EnergyRegion(NMO,IIFG-1,1)=EnergyRegion(NMO,IIFG,2) 
			   ENDIF
			ENDIF
	     ENDDO
         !DO IIFG=1,NumbreRegion(NMO)
         !   Write(6,*) 'grani',NMO,EnergyRegion(NMO,IIFG,1),EnergyRegion(NMO,IIFG,2)
         !   Write(6,*) 'ZERO and Delta',NMO,EnergyZeroZ(NMO,IIFG,2),DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		 !ENDDO
	  ENDIF
   ENDIF
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.GT.2) THEN
         ! жхйк он мюидеммшл йнпмъл
	     DO IIFG=1,NumbreRegion(NMO)
            Delta=DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
			! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF

		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,IIFG,2)=EnergyZeroZ(NMO,IIFG,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,IIFG,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,IIFG,1)=EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta 
			ENDIF 
			! опнбепъел врнаш накюярх йнпмеи ме оепеяхйюкхяэ мювхмюъ я брнпнцн йнпмъ
			IF(IIFG.GT.1) THEN
			   ! опнбепъел врнаш опюбюъ цпюмхжю дюммни накюярх ме ашкю анкэье кебни цпюмхжш опедшдсыеи накюярх
               IF(EnergyRegion(NMO,IIFG,2).GT.EnergyRegion(NMO,IIFG-1,1)) THEN
                  ! опхпюбмхбюел цпюмхжш опюбюъ цпюмхжю дюммни накюярх асдер янбоюдюрэ я кебни цпюмхжни опедшдсыеи накюярх
				  EnergyRegion(NMO,IIFG-1,1)=EnergyRegion(NMO,IIFG,2) 
			   ENDIF
			ENDIF
	     ENDDO
      ENDIF
   ENDIF

   ! пеьемхе ндмнпндмни яхярелш спюбмемхи 
   ExfZero=Exf


   return
 end subroutine SearchEnergySolutionHomogeneousSystemStandard 


 !! SUBPROGRAMME SEARCH FOR THE OWN ENERGY OF THE HOMOGENEOUS SYSTEM OF EQUATIONS IN THE "ACCELERATED" METHOD
═!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
═!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
═!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
═!! Niteration-NUMBER OF INTERACTION
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! NNEL- MAIN QUANTUM NUMBER
═!! Npoint-NUMBER OF POINTS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! Rcof-EQUITY RATIO H * H / 12
═!! EPS-ACCURACY
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! RO1X (Npoint) -MASSIVE FIRST DERIVATIVE NEW VARIABLE FOR THE OLD VARIABLE
═!! VNout, UNout, VNin, UNin, E, AZ, FZ, DZ, BZ, AX, BX-AUXILIARY MASSIVE
═!! ENN0-LOWER ENERGY BORDER
═!! ExfSolutions (3) -SOLUTION OF A HOMOGENEOUS SYSTEM OF EQUATIONS SINGLE-ELECTRONIC ENERGY FOR NUMBER QUANTITY NUMERALS NNEL-1, NNEL, NNEL + 1
═!! NumbreRegion (NumbreMO) -MASSIVE NUMBER OF REGIONS IN WHICH A ROOT REDUCTION IS CARRIED OUT
═!! EnergyZeroZ (NumbreMO, NNEL + 1,2) -MASSIVE OF THE ROOTS OF A HOMOGENEOUS EQUATION
═!! EnergyRegion (NumbreMO, NNEL + 1,2) -MASSIVE OF AREAS IN WHICH STEP DECREASES
 subroutine SearchEnergySolutionHomogeneousSystemAccelerated(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,NumbreRegion,EnergyZeroZ,EnergyRegion) 
   implicit none
   integer::IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,IndexCorrectSET,IndexSolutionsSET
   real(8)::Rcof,EPS,ENN0,ExfZero,HenergySET,ExfSETT
   integer,dimension(:)::NumbreRegion 
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,ExfSolutions,ExfSET,DetSET
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,EnergyZeroZ,EnergyRegion 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IparametrDeterEnergy,NEpsIter,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexRegion,IndexSmena,IndexSolutions,IKLKORX,IOPZ,IKLRXX,IDFJK,IIFG,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Henergy,Exf,ExfTemp,Delta,RcoffInterval,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(10)::DetM,ExfM,DetMA,ExfMA,ProDet  
   
   
   ! сяйнпеммши лернд пюяверю
   ! гюмскъел оепед пюявернл
   DetM=0.D0
   ExfM=0.D0
   DetMA=0.D0
   ExfMA=0.D0
   ExfSolutions=0.D0
   ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
   Exf=ENN0  
   ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
   Henergy=0.1D0
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=2
   ENDIF
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
      IF(Niteration.EQ.3) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=3
   ENDIF

   ! оюпюлерп бшундю хг жхйкю  
   IparametrDeterEnergy=1
   ! хмдейя жхйкю 
   NEpsIter=0
   ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
   IndexVilky=0
   INDEXDF=0
   IIKK=0
   IndexCorrect=0
   IndexRegion=1
   IndexSmena=0
   IndexSolutions=0
   IKLKORX=0
   IOPZ=0
   IKLRXX=0
   INDEXAD=0
   INDEXRShag=0
   IIHD=0
   IIOXDPZ=0
   INDEXPeresc=0
   DeltaEnerPers=0.D0
   DeltaFunPers=0.D0
   RcoffChengH=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюяверю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
      ! хыел яннрберярбсчыхе йнпмх Exf
      call SearchSolutionAccelerated(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,IKLKORX,IOPZ,IKLRXX,EnergyZeroZ,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)  
      ! нясыеярбкъел йнмрпнкэ гю ьюцюл 
	  IF(Niteration.GT.IDFJK) THEN
	     call ControlH(NMO,IndexSmena,NumbreRegion(NMO),IndexRegion,EnergyRegion,Exf,Henergy)
      ENDIF 
   ENDDO

   ! гюохяшбюел вхякн йнпмеи пюяялнрпеммшу б унде онхяйю
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.EQ.2) THEN
         NumbreRegion(NMO)=IndexCorrect
      ENDIF
   ENDIF

   ! тнплхпсел накюярх цде ьюц асдер слемэьхм
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.GT.1) THEN
         ! жхйк он мюидеммшл йнпмъл
	     DO IIFG=1,NumbreRegion(NMO)
            Delta=DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
			! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF

		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,IIFG,2)=EnergyZeroZ(NMO,IIFG,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,IIFG,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,IIFG,1)=EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta 
			ENDIF
			! опнбепъел врнаш накюярх йнпмеи ме оепеяхйюкхяэ мювхмюъ я брнпнцн йнпмъ
			IF(IIFG.GT.1) THEN
			   ! опнбепъел врнаш опюбюъ цпюмхжю дюммни накюярх ме ашкю анкэье кебни цпюмхжш опедшдсыеи накюярх
               IF(EnergyRegion(NMO,IIFG,2).GT.EnergyRegion(NMO,IIFG-1,1)) THEN
                  ! опхпюбмхбюел цпюмхжш опюбюъ цпюмхжю дюммни накюярх асдер янбоюдюрэ я кебни цпюмхжни опедшдсыеи накюярх
				  EnergyRegion(NMO,IIFG-1,1)=EnergyRegion(NMO,IIFG,2) 
			   ENDIF
			ENDIF
	     ENDDO
         !DO IIFG=1,NumbreRegion(NMO)
         !   Write(6,*) 'grani',NMO,EnergyRegion(NMO,IIFG,1),EnergyRegion(NMO,IIFG,2)
		 !   Write(6,*) 'ZERO and Delta',NMO,EnergyZeroZ(NMO,IIFG,2),DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		 !ENDDO
	  ENDIF
   ENDIF
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.GT.2) THEN
         ! жхйк он мюидеммшл йнпмъл
	     DO IIFG=1,NumbreRegion(NMO)
            Delta=DABS(EnergyZeroZ(NMO,IIFG,2)-EnergyZeroZ(NMO,IIFG,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
			! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF
		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,IIFG,2)=EnergyZeroZ(NMO,IIFG,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,IIFG,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,IIFG,1)=EnergyZeroZ(NMO,IIFG,2)-RcoffInterval*Delta 
			ENDIF 
			! опнбепъел врнаш накюярх йнпмеи ме оепеяхйюкхяэ мювхмюъ я брнпнцн йнпмъ
			IF(IIFG.GT.1) THEN
			   ! опнбепъел врнаш опюбюъ цпюмхжю дюммни накюярх ме ашкю анкэье кебни цпюмхжш опедшдсыеи накюярх
               IF(EnergyRegion(NMO,IIFG,2).GT.EnergyRegion(NMO,IIFG-1,1)) THEN
                  ! опхпюбмхбюел цпюмхжш опюбюъ цпюмхжю дюммни накюярх асдер янбоюдюрэ я кебни цпюмхжни опедшдсыеи накюярх
				  EnergyRegion(NMO,IIFG-1,1)=EnergyRegion(NMO,IIFG,2) 
			   ENDIF
			ENDIF
	     ENDDO
      ENDIF
   ENDIF


          

   return
 end subroutine SearchEnergySolutionHomogeneousSystemAccelerated


 !! SUBPROGRAMME SEARCH FOR THE OWN ENERGY OF THE HOMOGENEOUS SYSTEM OF EQUATIONS IN THE "ACCELERATED" METHOD
═!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
═!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
═!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
═!! Niteration-NUMBER OF INTERACTION
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! NNEL- MAIN QUANTUM NUMBER
═!! Npoint-NUMBER OF POINTS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! Rcof-EQUITY RATIO H * H / 12
═!! EPS-ACCURACY
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! RO1X (Npoint) -MASSIVE FIRST DERIVATIVE NEW VARIABLE FOR THE OLD VARIABLE
═!! VNout, UNout, VNin, UNin, E, AZ, FZ, DZ, BZ, AX, BX-AUXILIARY MASSIVE
═!! ENN0-LOWER ENERGY BORDER
═!! ExfSolutions (3) -SOLUTION OF A HOMOGENEOUS SYSTEM OF EQUATIONS SINGLE-ELECTRONIC ENERGY FOR NUMBER QUANTITY NUMERALS NNEL-1, NNEL, NNEL + 1
═!! NumbreRegion (NumbreMO) -MASSIVE NUMBER OF REGIONS IN WHICH A ROOT REDUCTION IS CARRIED OUT
═!! EnergyZeroZ (NumbreMO, NNEL + 1,2) -MASSIVE OF THE ROOTS OF A HOMOGENEOUS EQUATION
═!! EnergyRegion (NumbreMO, NNEL + 1,2) -MASSIVE OF AREAS IN WHICH STEP DECREASES
 subroutine SearchEnergySolutionHomogeneousSystemAccelerated2(IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,NumbreRegion,EnergyZeroZ,EnergyRegion) 
   implicit none
   integer::IndexExit,Niteration,NMO,NNEL,Npoint,NumeroMin,Kpoint,IndexM,IndexCorrectSET,IndexSolutionsSET
   real(8)::Rcof,EPS,ENN0,ExfZero,HenergySET,ExfSETT
   integer,dimension(:)::NumbreRegion 
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,ExfSolutions,ExfSET,DetSET
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,EnergyZeroZ,EnergyRegion 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IparametrDeterEnergy,NEpsIter,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexRegion,IndexSmena,IndexSolutions,IKLKORX,IOPZ,IKLRXX,IDFJK,IIFG,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Henergy,Exf,ExfTemp,Delta,RcoffInterval,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(10)::DetM,ExfM,DetMA,ExfMA,ProDet  
   
   
   ! сяйнпеммши лернд пюяверю
   ! гюмскъел оепед пюявернл
   DetM=0.D0
   ExfM=0.D0
   DetMA=0.D0
   ExfMA=0.D0
   ExfSolutions=0.D0
   ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
   Exf=ENN0  
   ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
   Henergy=0.1D0
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=2
   ENDIF
   ! мю оепбшу дбсу хрепюжхъу слемэьюел ьюц б 10 пюг
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.EQ.2) THEN
         Henergy=Henergy/10.D0
      ENDIF
      IF(Niteration.EQ.3) THEN
         Henergy=Henergy/10.D0
      ENDIF
	  ! мнлеп хрепюжхх
	  IDFJK=3
   ENDIF

   ! оюпюлерп бшундю хг жхйкю  
   IparametrDeterEnergy=1
   ! хмдейя жхйкю 
   NEpsIter=0
   ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
   IndexVilky=0
   INDEXDF=0
   IIKK=0
   IndexCorrect=0
   IndexRegion=1
   IndexSmena=0
   IndexSolutions=0
   IKLKORX=0
   IOPZ=0
   IKLRXX=0
   INDEXAD=0
   INDEXRShag=0
   IIHD=0
   IIOXDPZ=0
   INDEXPeresc=0
   DeltaEnerPers=0.D0
   DeltaFunPers=0.D0
   RcoffChengH=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюяверю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
      ! хыел яннрберярбсчыхе йнпмх Exf
      call SearchSolutionAccelerated2(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,IKLKORX,IOPZ,IKLRXX,EnergyZeroZ,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)  
      ! нясыеярбкъел йнмрпнкэ гю ьюцюл 
	  IF(Niteration.GT.IDFJK) THEN
	     call ControlH(NMO,IndexSmena,NumbreRegion(NMO),IndexRegion,EnergyRegion,Exf,Henergy)
      ENDIF 
   ENDDO

  

   ! тнплхпсел накюярх цде ьюц асдер слемэьхм дкъ оепбнцн йнпмъ
   IF(IndexExit.EQ.0) THEN
      IF(Niteration.GT.1) THEN
       
            Delta=DABS(EnergyZeroZ(NMO,1,2)-EnergyZeroZ(NMO,1,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
			! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF

		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,1,2)=EnergyZeroZ(NMO,1,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,1,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,1,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,1,1)=EnergyZeroZ(NMO,1,2)-RcoffInterval*Delta 
			ENDIF
	  ENDIF
   ENDIF
   IF(IndexExit.EQ.1) THEN
      IF(Niteration.GT.2) THEN
            Delta=DABS(EnergyZeroZ(NMO,1,2)-EnergyZeroZ(NMO,1,1))
		    IF(Delta.LT.0.1D0) THEN
               Delta=0.1D0
		    ENDIF
			! нопедекъел йнщттхжхемр йнрнпши тнплхпсер хмрепбюк
            IF(Delta.LT.1.D0) THEN
               RcoffInterval=10.D0  
			ENDIF
            IF(Delta.GT.1.D0.AND.Delta.LT.5.D0) THEN
               RcoffInterval=5.D0 
			ENDIF
            IF(Delta.GT.10.D0.AND.Delta.LT.15.D0) THEN
               RcoffInterval=2.5D0 
			ENDIF
            IF(Delta.GT.15.D0) THEN
               RcoffInterval=1.D0 
			ENDIF
		    ! опюбюъ цпюмхжю 
		    EnergyRegion(NMO,1,2)=EnergyZeroZ(NMO,1,2)+RcoffInterval*Delta 
		    ! кебюъ цпюмхжю
			! опнбепъел кебсч цпюмхжс нмю ме лнфер гюирх б нрпхжюрекэмсч накюярэ
			IF((EnergyZeroZ(NMO,1,2)-RcoffInterval*Delta).LT.0.D0) THEN
			    EnergyRegion(NMO,1,1)=0.D0 
			   ELSE    
		        EnergyRegion(NMO,1,1)=EnergyZeroZ(NMO,1,2)-RcoffInterval*Delta 
			ENDIF 
      ENDIF
   ENDIF


          

   return
 end subroutine SearchEnergySolutionHomogeneousSystemAccelerated2


 !! SUBPROGRAMME SEARCH FOR THE OWN ENERGY OF THE HOMOGENEOUS SYSTEM OF EQUATIONS IN THE "ACCELERATED" METHOD FOR FINDING THE SOLUTION OF LEFT
═!! NNEL- MAIN QUANTUM NUMBER
═!! Npoint-NUMBER OF POINTS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! Rcof-EQUITY RATIO H * H / 12
═!! EPS-ACCURACY
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! RO1X (Npoint) -MASSIVE FIRST DERIVATIVE NEW VARIABLE FOR THE OLD VARIABLE
═!! VNout, UNout, VNin, UNin, E, AZ, FZ, DZ, BZ, AX, BX-AUXILIARY MASSIVE
═!! ENN0-LOWER ENERGY BORDER
═!! ExfSolutions (3) -SOLUTION OF A HOMOGENEOUS SYSTEM OF EQUATIONS SINGLE-ELECTRONIC ENERGY FOR NUMBER QUANTITY NUMERALS NNEL-1, NNEL, NNEL + 1
 subroutine SearchEnergySolutionHomogeneousSystemAcceleratedLeft(NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET) 
   implicit none
   integer::NNEL,Npoint,NumeroMin,Kpoint,IndexM,IndexCorrectSET,IndexSolutionsSET
   real(8)::Rcof,EPS,ENN0,ExfZero,HenergySET,ExfSETT
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,ExfSolutions,ExfSET,DetSET
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IparametrDeterEnergy,NEpsIter,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Henergy,Exf,ExfTemp,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(10)::DetM,ExfM,DetMA,ExfMA,ProDet   
   
   
   ! сяйнпеммши лернд пюяверю
   ! гюмскъел оепед пюявернл
   DetM=DetSET
   ExfM=ExfSET
   DetMA=0.D0
   ExfMA=0.D0
   
   ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
   Exf=ExfSETT
   ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
   ! свхршбюел рн тюйр, врн дкъ онхяйю дюммнцн йнпмъ слемэьюел ьюц б 10 пюг
   ! дюммши лернд хяонкэгнбюм бн бяеу ондопнцпюллю он онхяйс йнпмъ ндмнпндмнцн спюбмемхъ
   ! слемэьемхе ьюцю опнхгбндхряъ б накюярх бнглнфмнцн мюунфдемхъ йнпмъ
   ! оняйнкэйс дюммюъ опнжедспю гюосяйюеряъ ме мю йюфдни хрепюжхх
   ! ю рнкэйн опх сякнбхх нрясрярбхъ яопюбю рн накюярэч слемэьемхъ ьюцю явхрюел бяч накюярэ онхяйю    
   Henergy=HenergySET !*0.1D0

   ! оюпюлерп бшундю хг жхйкю  
   IparametrDeterEnergy=1
   ! хмдейя жхйкю 
   NEpsIter=3
   ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
   IndexVilky=0
   INDEXDF=0
   IIKK=0
   IndexCorrect=IndexCorrectSET
   IndexSolutions=IndexSolutionsSET
   IKLKORX=0
   IOPZ=0
   IKLRXX=0
   INDEXAD=0
   INDEXRShag=0
   IIHD=0
   IIOXDPZ=0
   INDEXPeresc=0
   DeltaEnerPers=0.D0
   DeltaFunPers=0.D0
   RcoffChengH=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюяверю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
      ! хыел яннрберярбсчыхе йнпмх Exf
      call SearchSolutionAcceleratedLeft(NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)  
   ENDDO
          

   return
 end subroutine SearchEnergySolutionHomogeneousSystemAcceleratedLeft

!! SUBPROGRAMME SEARCH FOR THE OWN ENERGY OF THE HOMOGENEOUS SYSTEM OF EQUATIONS IN THE "ACCELERATED" METHOD FOR FINDING THE SOLUTION OF LEFT
═!! NNEL- MAIN QUANTUM NUMBER
═!! Npoint-NUMBER OF POINTS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! IndexM-SIZE OF MASSIVE
═!! Rcof-EQUITY RATIO H * H / 12
═!! EPS-ACCURACY
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! RO1X (Npoint) -MASSIVE FIRST DERIVATIVE NEW VARIABLE FOR THE OLD VARIABLE
═!! VNout, UNout, VNin, UNin, E, AZ, FZ, DZ, BZ, AX, BX-AUXILIARY MASSIVE
═!! ENN0-LOWER ENERGY BORDER
═!! ExfSolutions (3) -SOLUTION OF A HOMOGENEOUS SYSTEM OF EQUATIONS SINGLE-ELECTRONIC ENERGY FOR NUMBER QUANTITY NUMERALS NNEL-1, NNEL, NNEL + 1
═!! INDEXGRNYT-TYPE OF DEFINED INTERVAL BORDER
 subroutine SearchEnergySolutionHomogeneousSystemAcceleratedLeft2(NNEL,Npoint,NumeroMin,Kpoint,IndexM,Rcof,EPS,Pz,RO1X,A,B,AA,BB,E,AZ,FZ,BZ,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DZ,AX,BX,ENN0,ExfSolutions,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,INDEXGRNYT) 
   implicit none
   integer::NNEL,Npoint,NumeroMin,Kpoint,IndexM,IndexCorrectSET,IndexSolutionsSET,INDEXGRNYT
   real(8)::Rcof,EPS,ENN0,ExfZero,HenergySET,ExfSETT
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,ExfSolutions,ExfSET,DetSET
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IparametrDeterEnergy,NEpsIter,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,INDEXSHAA,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Henergy,Exf,ExfTemp,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(10)::DetM,ExfM,DetMA,ExfMA,ProDet   
   
   
   ! сяйнпеммши лернд пюяверю  
   ! гюмскъел оепед пюявернл
   DetM=DetSET
   ExfM=ExfSET
   DetMA=0.D0
   ExfMA=0.D0
   
   ! мювюкэмне опхакхфемхе дкъ ндмнщкейрпнммни щмепцхх
   Exf=ExfSETT
   ! мювюкэмши ьюц он ндмнщкейрпнммни щмепцхх
   ! свхршбюел рн тюйр, врн дкъ онхяйю дюммнцн йнпмъ слемэьюел ьюц б 10 пюг
   ! дюммши лернд хяонкэгнбюм бн бяеу ондопнцпюллю он онхяйс йнпмъ ндмнпндмнцн спюбмемхъ
   ! слемэьемхе ьюцю опнхгбндхряъ б накюярх бнглнфмнцн мюунфдемхъ йнпмъ
   ! оняйнкэйс дюммюъ опнжедспю гюосяйюеряъ ме мю йюфдни хрепюжхх
   ! ю рнкэйн опх сякнбхх нрясрярбхъ яопюбю рн накюярэч слемэьемхъ ьюцю явхрюел бяч накюярэ онхяйю    
   Henergy=HenergySET !*0.1D0

   ! оюпюлерп бшундю хг жхйкю  
   IparametrDeterEnergy=1
   ! хмдейя жхйкю 
   NEpsIter=3
   ! оюпюлерп сйюгшбючыхи мю рн врн хмрепбюк цде мюундхряъ йнпемэ нопедекем
   IndexVilky=0
   INDEXDF=0
   IIKK=0
   IndexCorrect=IndexCorrectSET
   IndexSolutions=IndexSolutionsSET
   IKLKORX=0
   IOPZ=0
   IKLRXX=0
   INDEXAD=0
   INDEXRShag=0
   INDEXSHAA=0
   IIHD=0
   IIOXDPZ=0
   INDEXPeresc=0
   DeltaEnerPers=0.D0
   DeltaFunPers=0.D0
   RcoffChengH=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
      call PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
      ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюяверю оюпюлерпнб опнцнмйх
      call FormMassivAABB(Npoint,IndexM,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ ндмнпндмни яхярелш
      call CalculationParameterProgonky(1,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
      ! хыел яннрберярбсчыхе йнпмх Exf
      call SearchSolutionAcceleratedLeft2(NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,INDEXGRNYT,INDEXSHAA,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)  
   ENDDO
          

   return
 end subroutine SearchEnergySolutionHomogeneousSystemAcceleratedLeft2




!! SUBPROGRAMME DETERMINES THE MINIMUM VALUE OF NORM IN THE ENERGY INTERVAL
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! IndexM-SIZE OF MASSIVE
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! H-STEP
═!! Rcof-EQUITY RATIO H * H / 12
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! ExfMin, ExfMax-Energetic interval in which the minimum is searched
═!! REmin-Energy corresponding to a minimum
═!! RfunMIN-Value of the norm
 subroutine  MinValueNormEnergyFunction(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfMin,ExfMax,REmin,RfunMIN)
   implicit none
   integer::NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid
   real(8)::Rcof,ExfMin,ExfMax,H,RfunMIN,REmin
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NEpsIter,IparametrDeterEnergy 
   real(8)::RHS,RX0,RX1,RX2,Rfun0XZ,RnormX0,Xfun0XZ,Rfun1XZ,RnormX1,Xfun1XZ,Rfun2XZ,RnormX2,Xfun2XZ 
   
   ! мювюкэмши ьюц
   RHS=(ExfMax-ExfMin)/100.D0
   ! мювюкэмше рнвйх
   ! яепедхмю хмрепбюкю 
   RX0=(ExfMax+ExfMin)*0.5D0 
   ! анйнбше рнвйх
   RX1=RX0-RHS
   RX2=RX0+RHS
   NEpsIter=0
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX0
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX0,Pz,RO1X,Kpoint) 

	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX0,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                  
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX0,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
        
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX1,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX1,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	  
	  ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX1,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
            
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX2
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX2,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX2,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
             
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX2,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
	  
	  ! гюохяшбюел гмювемхе  
      IF(NEpsIter.EQ.1) THEN
           ! оепбюъ хрепюжхъ 
           Rfun0XZ=RnormX0
	       Xfun0XZ=RX0
           Rfun1XZ=RnormX1
           Xfun1XZ=RX1
		   Rfun2XZ=RnormX2
           Xfun2XZ=RX2 
		   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
           RX1=RX0-RHS
		   RX2=RX0+RHS
		 ELSE
           
		   WRITE(17,*) NEpsIter,RX0,RnormX0          
		   ! опнбепъел ндхмюйнбше гмювемхъ тсмйжхх лш онксвюел
		   ! якеднбюрекэмюъ оепбюъ опнхгбндмюъ пюбмю мскч	 
           IF(RnormX0.EQ.Rfun0XZ) THEN
              ! лхмхлюкэмне гмювемхе тсмйжхх
		      RfunMIN=RnormX0 
		      ! яннрберярбсчыхи юпцслемр
              REmin=RX0
			  ! сякнбхе бшундю бшонкмемн
              IparametrDeterEnergy=2 
           ENDIF 

		
		   ! опнбепъел онксвеммне мнбне опхакхфемхе лемэье опедедсыецн хкх мер
		   IF(RnormX0.LT.Rfun0XZ) THEN
                ! лхмхлюкэмне гмювемхе тсмйжхх
		        RfunMIN=RnormX0 
		        ! яннрберярбсчыхи юпцслемр
                REmin=RX0
              
			    ! опнбепъел сякнбхе бшундю (опнбепъел пюбемярбн мскч оепбни опнгбндмни)
			    IF(DABS((RfunMIN-Rfun0XZ)/(REmin-Xfun0XZ)).LT.10.D0**(-10)) THEN
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
			    ENDIF   
              
			    ! гюохяшбюел опедшдсыее гмювемхе
			    Rfun0XZ=RnormX0
			    Xfun0XZ=RX0
                Rfun1XZ=RnormX1
                Xfun1XZ=RX1
			    Rfun2XZ=RnormX2
                Xfun2XZ=RX2 

                ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! лхмхлюкэмне гмювемхе тсмйжхх
		             RfunMIN=RnormX0 
		             ! яннрберярбсчыхи юпцслемр
                     REmin=RX0
			         ! сякнбхе бшундю бшонкмемн
                     IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=RnormX0 
		           ! яннрберярбсчыхи юпцслемр
                   REmin=RX0
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

			
                ! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
			       RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
			
			  ELSE
				! гмювемхе тсмйжхх пюярер
				! лемъел ьюц
				RHS=-RHS*0.5D0
			

				! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

				! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
				   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
				
		   ENDIF
	  ENDIF
   ENDDO

 

   return
 end subroutine MinValueNormEnergyFunction

!! SUBPROGRAMME DETERMINES THE MINIMUM VALUE OF NORM IN THE ENERGY INTERVAL
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! IndexM-SIZE OF MASSIVE
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! H-STEP
═!! Rcof-EQUITY RATIO H * H / 12
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! ExfMin, ExfMax-Energetic interval in which the minimum is searched. The minimum is in the vicinity of ExfMin
═!! REmin-Energy corresponding to a minimum
═!! RfunMIN-Value of the norm
 subroutine  MinValueNormEnergyFunctionD(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfMin,ExfMax,REmin,RfunMIN)
   implicit none
   integer::NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid
   real(8)::Rcof,ExfMin,ExfMax,H,RfunMIN,REmin
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NEpsIter,IparametrDeterEnergy 
   real(8)::RHS,RX0,RX1,RX2,Rfun0XZ,RnormX0,Xfun0XZ,Rfun1XZ,RnormX1,Xfun1XZ,Rfun2XZ,RnormX2,Xfun2XZ 
   
   ! мювюкэмши ьюц
   RHS=(ExfMax-ExfMin)/100.D0
   ! мювюкэмше рнвйх
   RX0=ExfMin  
   ! анйнбше рнвйх
   RX1=RX0-RHS
   RX2=RX0+RHS
   NEpsIter=0
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX0
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX0,Pz,RO1X,Kpoint) 

	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX0,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                  
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX0,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
        
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX1,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX1,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	  
	  ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX1,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
            
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX2
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX2,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX2,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
             
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX2,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
	  
	  ! гюохяшбюел гмювемхе  
      IF(NEpsIter.EQ.1) THEN
           ! оепбюъ хрепюжхъ 
           Rfun0XZ=RnormX0
	       Xfun0XZ=RX0
           Rfun1XZ=RnormX1
           Xfun1XZ=RX1
		   Rfun2XZ=RnormX2
           Xfun2XZ=RX2 
		   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
           RX1=RX0-RHS
		   RX2=RX0+RHS
		 ELSE
           
		   WRITE(17,*) NEpsIter,RX0,RnormX0          
		   ! опнбепъел ндхмюйнбше гмювемхъ тсмйжхх лш онксвюел
		   ! якеднбюрекэмюъ оепбюъ опнхгбндмюъ пюбмю мскч	 
           IF(RnormX0.EQ.Rfun0XZ) THEN
              ! лхмхлюкэмне гмювемхе тсмйжхх
		      RfunMIN=RnormX0 
		      ! яннрберярбсчыхи юпцслемр
              REmin=RX0
			  ! сякнбхе бшундю бшонкмемн
              IparametrDeterEnergy=2 
           ENDIF 

		
		   ! опнбепъел онксвеммне мнбне опхакхфемхе лемэье опедедсыецн хкх мер
		   IF(RnormX0.LT.Rfun0XZ) THEN
                ! лхмхлюкэмне гмювемхе тсмйжхх
		        RfunMIN=RnormX0 
		        ! яннрберярбсчыхи юпцслемр
                REmin=RX0
              
			    ! опнбепъел сякнбхе бшундю (опнбепъел пюбемярбн мскч оепбни опнгбндмни)
			    IF(DABS((RfunMIN-Rfun0XZ)/(REmin-Xfun0XZ)).LT.10.D0**(-10)) THEN
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
			    ENDIF   
              
			    ! гюохяшбюел опедшдсыее гмювемхе
			    Rfun0XZ=RnormX0
			    Xfun0XZ=RX0
                Rfun1XZ=RnormX1
                Xfun1XZ=RX1
			    Rfun2XZ=RnormX2
                Xfun2XZ=RX2 

                ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! лхмхлюкэмне гмювемхе тсмйжхх
		             RfunMIN=RnormX0 
		             ! яннрберярбсчыхи юпцслемр
                     REmin=RX0
			         ! сякнбхе бшундю бшонкмемн
                     IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=RnormX0 
		           ! яннрберярбсчыхи юпцслемр
                   REmin=RX0
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

			
                ! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
			       RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
			
			  ELSE
				! гмювемхе тсмйжхх пюярер
				! лемъел ьюц
				RHS=-RHS*0.5D0
			

				! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

				! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
				   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
				
		   ENDIF
	  ENDIF
   ENDDO

 

   return
 end subroutine MinValueNormEnergyFunctionD


 !! SUBPROGRAMME DETERMINES THE MINIMUM VALUE OF NORM IN THE ENERGY INTERVAL
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! IndexM-SIZE OF MASSIVE
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! H-STEP
═!! Rcof-EQUITY RATIO H * H / 12
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! ExfMin, ExfMax-Energetic interval in which the minimum is searched. The minimum is in the vicinity of ExfMin
═!! REmin-Energy corresponding to a minimum
═!! RfunMIN-Value of the norm
═!! IndexX-INDEX INDICATES A BORDER OF A MINIMUM BORDER OR NOT
═!! IndexX = 0-NOT INSTALLED
═!! IndexX = 1-INSTALLED
═!! EminIndex-LEFT BORDER OF THE MINIMUM FIELD
═!! EmaxIndex-MINIMUM
═!! RFUNmaxIndex-The value of the norm in the minimum
═!! RELEFT-LEFT BORDER OF "NEW" MINIMUM
═!! REright-RIGHT BORDER OF "NEW" MINIMUM
═!! IndexMinDifini-INDEX INDICATOR DETERMINED BY NEW MINIMUM OR ORDERED BY OLD MINIMUM
═!! IndexMinDifini = 0-NEW MINIMUM DEFINED
═!! IndexMinDifini = 1-DEFINED AGAIN "OLD" MINIMUM
═!! IndexHDXX-INDEX INDICATOR DETERMINED THE INITIAL STEP TO DETERMINE BORDERS OF MINIMUM
═!! HXXZXL-STEP FOR LEFT ARE FROM MINIMUM
═!! HXXZXR-STEP FOR THE RIGHT OF THE RIGHT FROM THE MINIMUM
 subroutine  MinValueNormEnergyFunctionD2(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,HNN,ExfNN,REmin,RfunMIN,IndexX,EminIndex,EmaxIndex,RFUNmaxIndex,REleft,REright,IndexMinDifini,IndexHDXX,HXXZXL,HXXZXR)
   implicit none
   integer::NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,IndexX,IndexMinDifini,IndexHDXX
   real(8)::Rcof,HNN,ExfNN,H,RfunMIN,REmin,EminIndex,EmaxIndex,RFUNmaxIndex,REleft,REright,HXXZXL,HXXZXR
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NEpsIter,IparametrDeterEnergy,IndexRESX,IndexSUMZKX,IndexIIUMNX,INDEXEXITRFX,INDEXZZSFAX,IndexCorGRANF,IndexPovNorma,IndexNORMSDS 
   real(8)::RHS,RX0,RX1,RX2,Rfun0XZ,RnormX0,Xfun0XZ,Rfun1XZ,RnormX1,Xfun1XZ,Rfun2XZ,RnormX2,Xfun2XZ,EXDmin,EXDmax
   real(8)::RXDF,RfunLR,RcoffSD,HenergyX,E1XX,D1XX,RcofXZSDX,DetxfNN 
   real(8),dimension(10)::DetMX,ExfMX,DetMAX

   ! гюмскъел оепед пюявернл
   IndexMinDifini=0
   REleft=0.D0
   REright=0.D0
   
   ! мювюкэмше рнвйх
   EXDmin=ExfNN*(1.D0-HNN)
   EXDmax=ExfNN*(1.D0+HNN)

   ! мювюкэмши ьюц
   RHS=DABS(EXDmax-EXDmin)*0.5D0
   ! мювюкэмше рнвйх
   RX0=ExfNN 
   ! анйнбше рнвйх
   RX1=RX0-RHS
   RX2=RX0+RHS
   NEpsIter=0
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   DO WHILE(IparametrDeterEnergy.EQ.1)
      ! хмдейя рхою ялеыемхъ
	  IndexNORMSDS=0 
	  ! явервхй
      NEpsIter=NEpsIter+1
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX0
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX0,Pz,RO1X,Kpoint) 

	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX0,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                  
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX0,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
        
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX1,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX1,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	  
	  ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX1,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
            
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX2
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX2,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX2,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
             
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX2,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
	  
	  ! гюохяшбюел гмювемхе  
      IF(NEpsIter.EQ.1) THEN
           ! оепбюъ хрепюжхъ 
           ! сярюмюбкхбюел онбедемхе мнплш
		   IndexPovNorma=0
		   IF(RnormX1.GT.RnormX0.AND.RnormX0.GT.RnormX2) THEN
              IndexPovNorma=1
		   ENDIF
		   IF(RnormX2.GT.RnormX0.AND.RnormX0.GT.RnormX1) THEN
              IndexPovNorma=1
           ENDIF


           Rfun0XZ=RnormX0
	       Xfun0XZ=RX0
           Rfun1XZ=RnormX1
           Xfun1XZ=RX1
		   Rfun2XZ=RnormX2
           Xfun2XZ=RX2 
		   IF(IndexPovNorma.EQ.1) THEN
		       RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
              ELSE
			   RX0=Xfun0XZ-RHS*2.D0
		   ENDIF
		   RX1=RX0-RHS
		   RX2=RX0+RHS
		   ! гюохяшбюел гмювемхе мнплш дкъ дюкэмеиьецн пюяверю
		   DetxfNN=RnormX0 
		 ELSE
           
		   ! опнбепъел оноюдюер дюммне опхакхфемхе он щмепцхх б накюярэ сфе мюидеммнцн лхмхлслю 
           IF(IndexX.EQ.1) THEN
              ! кебюъ цпюмхжю накюярх нопедекемю
			  IF(RX0.GT.EminIndex.AND.RX0.LT.EmaxIndex) THEN
                 ! щмепцхъ оноюкю б накюярэ опедшдсыецн лхмхлслю
                 ! лхмхлюкэмне гмювемхе тсмйжхх
		         RfunMIN=RFUNmaxIndex
		         ! яннрберярбсчыхи юпцслемр
                 REmin=EmaxIndex
				 ! оюпюлерп сйюгшбючыхи, врн щмепцхъ оноюкю б накюярэ нопедекеммнцн лхмхлслю
                 IndexMinDifini=1
				 ! сякнбхе бшундю бшонкмемн
                 IparametrDeterEnergy=2 
				 WRITE(17,*) NEpsIter,REmin,RfunMIN
			  ENDIF
		   ENDIF
           
		   ! опнбепъел ме нясыеярбкем бшунд 
		   IF(IparametrDeterEnergy.NE.2) THEN
		      WRITE(17,*) NEpsIter,RX0,RnormX0          
		      ! опнбепъел ндхмюйнбше гмювемхъ тсмйжхх лш онксвюел
		      ! якеднбюрекэмюъ оепбюъ опнхгбндмюъ пюбмю мскч	 
              IF(RnormX0.EQ.Rfun0XZ) THEN
			     ! опнбепъел мнлеп хрепюжхх
				 IF(NEpsIter.EQ.2) THEN
                     IndexNORMSDS=1
					ELSE 
					 ! лхмхлюкэмне гмювемхе тсмйжхх
		             RfunMIN=RnormX0 
		             ! яннрберярбсчыхи юпцслемр
                     REmin=RX0
			         ! сякнбхе бшундю бшонкмемн
                     IparametrDeterEnergy=2 
                 ENDIF
			  ENDIF 

		
		      ! опнбепъел онксвеммне мнбне опхакхфемхе лемэье опедедсыецн хкх мер
		      IF(RnormX0.LT.Rfun0XZ) THEN
                  ! лхмхлюкэмне гмювемхе тсмйжхх
		          RfunMIN=RnormX0 
		          ! яннрберярбсчыхи юпцслемр
                  REmin=RX0
              
			      ! опнбепъел сякнбхе бшундю (опнбепъел пюбемярбн мскч оепбни опнгбндмни)
			      IF(DABS((RfunMIN-Rfun0XZ)/(REmin-Xfun0XZ)).LT.10.D0**(-10)) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
			         ENDIF 
				  ENDIF   
              
			      ! гюохяшбюел опедшдсыее гмювемхе
			      Rfun0XZ=RnormX0
			      Xfun0XZ=RX0
                  Rfun1XZ=RnormX1
                  Xfun1XZ=RX1
			      Rfun2XZ=RnormX2
                  Xfun2XZ=RX2 

                  ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                  IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
				         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=RnormX0 
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=RX0
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF
				  ENDIF 
				  !опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				  IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
					     ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=RnormX0 
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=RX0
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF
				  ENDIF

			
                  ! мнбне опхакхфемхе
				  IF(IparametrDeterEnergy.EQ.1) THEN
				     ! опнбепъел хмдейя рхою ялеыемхъ
					 IF(IndexNORMSDS.EQ.0) THEN
			             RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                        ELSE
                         RX0=Xfun0XZ-RHS*2.D0
					 ENDIF
					 RX1=RX0-RHS
			         RX2=RX0+RHS
			      ENDIF
			
			    ELSE
				  ! гмювемхе тсмйжхх пюярер
				  ! лемъел ьюц
				  RHS=-RHS*0.5D0
			

				  ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                  IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
				         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=Rfun0XZ
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=Xfun0XZ
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF
				  ENDIF 
				  !опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				  IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
                         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=Rfun0XZ
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=Xfun0XZ
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF 
				  ENDIF

				  ! мнбне опхакхфемхе
				  IF(IparametrDeterEnergy.EQ.1) THEN
				     ! опнбепъел хмдейя рхою ялеыемхъ
					 IF(IndexNORMSDS.EQ.0) THEN
				         RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                        ELSE
                         RX0=Xfun0XZ-RHS*4.D0
					 ENDIF 
					 RX1=RX0-RHS
			         RX2=RX0+RHS
			      ENDIF
			  ENDIF
           ENDIF
	  ENDIF
   ENDDO

   ! опнбепъел нопедекем мнбши лхмхлсл хкх нопедекем хгбеярмши лхмхлсл
   IF(IndexMinDifini.EQ.0) THEN
	  ! сярюмюбкхбюел йюйсч цпюмхжс мсфмн нопедекхрэ
	  IF(REmin.GT.ExfNN) THEN
           ! мсфмн нопедекхрэ опюбсч цпюмхжс
		   ! кебюъ цпюмхжю
           REleft=ExfNN 
		   ! оюпюлерп ядбхцю 
		   RcoffSD=1.D0
		 ELSE
           ! мсфмн нопедекхрэ кебсч цпюмхжс
		   ! опюбюъ цпюмхжю
           REright=ExfNN
		   ! оюпюлерп ядбхцю
           RcoffSD=-1.D0
	  ENDIF

	  ! сярюмюбкхбюел йпхрепхи дкъ нопедекемхъ кебни цпюмхжш лхмхлслю
      IF(0.2D0*RfunMIN.LT.0.5D0) THEN
           RXDF=1.2D0*RfunMIN
         ELSE
 	       RXDF=RfunMIN+0.2D0
      ENDIF

	  ! нопедекъел рс цпюмхжс йнрнпсч менопедекхкх
      ! нопедекъел кебсч цпюмхжс лхмхлслю
      WRITE(17,*)
	  IF(RcoffSD.LT.0.D0) THEN
         WRITE(17,*) 'Defeni grany left min'
      ENDIF
	  IF(RcoffSD.GT.0.D0) THEN
         WRITE(17,*) 'Defeni grany right min' 
	  ENDIF

      ! гюмскъел оепед пюявернл
      DetMX=0.D0
      DetMAX=0.D0
      ExfMX=0.D0
      ! нопедекъел мювюкэмши ьюц
	  IF(IndexHDXX.EQ.0) THEN
          HenergyX=0.02D0*DABS(ExfNN-REmin)
		  ! тхйяхпсел нцпюмхвемхе дкъ ьюцю
		  IF(HenergyX.GT.0.1D0) THEN
             HenergyX=0.1D0  
		  ENDIF
         ELSE
          IF(RcoffSD.LT.0.D0) THEN
              HenergyX=0.1D0*HXXZXL
			 ELSE
			  HenergyX=0.1D0*HXXZXR
		  ENDIF
	  ENDIF
	  ! гюохяшбюел оепбне гмювемхе
      ExfMX(7)=REmin
      DetMX(7)=RfunMIN
	  ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
      E1XX=REmin*(1+RcoffSD*HenergyX) 
      IndexSUMZKX=0	
      ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
      IndexIIUMNX=0	         			  
      ! оюпюлерп бшундю хг жхйкю 
      IparametrDeterEnergy=1 
      ! оюпюлерп бшундю хг жхйкю
      INDEXEXITRFX=0
      INDEXZZSFAX=0
      RcofXZSDX=1.D0
      DO WHILE(IparametrDeterEnergy.EQ.1) 
         ! хмдейя жхйкю
	     IndexSUMZKX=IndexSUMZKX+1
	     ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1XX
	     call CalculationNormFunctionEnergyXF(E1XX,D1XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
	     WRITE(17,*) E1XX,D1XX,HenergyX
	     ! гюохяшбюел пегскэрюрш пюяверю
	     ExfMX(1)=ExfMX(2)
         DetMX(1)=DetMX(2)
	     ExfMX(2)=ExfMX(3)
         DetMX(2)=DetMX(3)
	     ExfMX(3)=ExfMX(4)
         DetMX(3)=DetMX(4)
	     ExfMX(4)=ExfMX(5)
         DetMX(4)=DetMX(5)
         ExfMX(5)=ExfMX(6)
         DetMX(5)=DetMX(6)
         ExfMX(6)=ExfMX(7)
         DetMX(6)=DetMX(7)
	     ExfMX(7)=E1XX
         DetMX(7)=D1XX
							 
	     IF(IndexSUMZKX.EQ.1) THEN
            ! пюявхршбюел опнхгбндмсч
            DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))
         ENDIF
         IF(IndexSUMZKX.GT.1) THEN
            ! пюявхршбюел опнхгбндмсч
		    DetMAX(1)=DetMAX(2)
		    DetMAX(2)=DetMAX(3)
            DetMAX(3)=DetMAX(4)
		    DetMAX(4)=DetMAX(5)
		    DetMAX(5)=DetMAX(6)
            DetMAX(6)=DetMAX(7)
		    DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))

            ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
			IF(RcoffSD.LT.0.D0) THEN
		       call CorrectHRightD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			ENDIF
			IF(RcoffSD.GT.0.D0) THEN
			   call CorrectHLeftD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			ENDIF    		 
		    ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
		    IF(INDEXZZSFAX.EQ.1) THEN
               ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
               E1XX=ExfMX(5) 
               ExfMX(7)=ExfMX(5) 
               DetMX(7)=DetMX(5)
               ExfMX(6)=ExfMX(4) 
               DetMX(6)=DetMX(4)
			   ExfMX(5)=ExfMX(3) 
               DetMX(5)=DetMX(3)
			   ExfMX(4)=ExfMX(2) 
               DetMX(4)=DetMX(2)
			   ExfMX(3)=ExfMX(1) 
               DetMX(3)=DetMX(1)
               DetMAX(7)=DetMAX(5)
               DetMAX(6)=DetMAX(4)
			   DetMAX(5)=DetMAX(3)
			   DetMAX(4)=DetMAX(2)
			   DetMAX(3)=DetMAX(1)
			   ! намскъел бшунд хг жхйкю
			   INDEXEXITRFX=0
			   ! опнбепъел вхякн хрепюжхи
			   IF(IndexSUMZKX.LT.3) THEN
                  ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				  IndexSUMZKX=0 
			   ENDIF
	        ENDIF
	     ENDIF

	     ! мнбне опхакхфемхе дкъ щмепцхх                 
         E1XX=E1XX*(1+RcoffSD*HenergyX)
	     ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
	     IF(INDEXZZSFAX.EQ.0) THEN
            IF(D1XX.GT.RXDF) THEN 
               INDEXEXITRFX=INDEXEXITRFX+1
		    ENDIF
		    ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
		    IF(INDEXEXITRFX.EQ.5) THEN 
			   IF(RcoffSD.LT.0.D0) THEN 
                   ! кебюъ цпюмхжю нопедекемю
			       REleft=ExfMX(7)
			       RfunLR=DetMX(7)
			   ENDIF
			   IF(RcoffSD.GT.0.D0) THEN 
                   ! опюбюъ цпюмхжю нопедекемю
			       REright=ExfMX(7)
			       RfunLR=DetMX(7)
			   ENDIF    
		       ! сякнбхе бшундю бшонкмемш
               IparametrDeterEnergy=2
		    ENDIF
	     ENDIF
      ENDDO
      
	  ! опнбепъел цпюмхжс йнрнпюъ хгмювюкэмн ашкю нопедекемю
      ! опнбепъел дпсцюъ (кебюъ хкх опюбюъ) цпюмхжю яннрберярбсер сякнбхч хкх мер
      ! сярюмюбкхбюел йюйсч цпюмхжс мсфмн нопедекхрэ
	  ! йпхрепхи йнпейрхпнбйх цпюмхжш
	  IndexCorGRANF=0
	  IF(REmin.GT.ExfNN) THEN
	       ! опнбепъел бшонкмъеряъ сякнбхе дкъ кебни цпюмхжш
		   IF(DetxfNN.LT.RXDF) THEN
              ! мнплю б кебни цпюмхже лемэье йпхрепхъ менаундхлн йнппейрхпнбйю кебни цпюмхжш
			  ! оюпюлерп ядбхцю 
		      RcoffSD=-1.D0
			  IndexCorGRANF=1   
		   ENDIF		   
		 ELSE
           ! опнбепъел бшонкмъеряъ сякнбхе дкъ опюбни цпюмхжш
		   IF(DetxfNN.LT.RXDF) THEN
              ! мнплю б опюбни цпюмхже лемэье йпхрепхъ менаундхлн йнппейрхпнбйю опюбни цпюмхжш
			  ! оюпюлерп ядбхцю 
		      RcoffSD=1.D0
			  IndexCorGRANF=1   
		   ENDIF	
	  ENDIF

	  ! опнбепъел менаундхлн йнпейрхпнбйю цпюмхжш хкх мер
	  IF(IndexCorGRANF.EQ.1) THEN
	     
		 ! нопедекъел рс цпюмхжс йнрнпсч менопедекхкх
         ! нопедекъел кебсч цпюмхжс лхмхлслю
         WRITE(17,*)
	     IF(RcoffSD.LT.0.D0) THEN
            WRITE(17,*) 'Defeni grany left min'
         ENDIF
	     IF(RcoffSD.GT.0.D0) THEN
            WRITE(17,*) 'Defeni grany right min' 
	     ENDIF

         ! гюмскъел оепед пюявернл
         DetMX=0.D0
         DetMAX=0.D0
         ExfMX=0.D0
         ! гюохяшбюел оепбне гмювемхе
         ExfMX(7)=REmin
         DetMX(7)=RfunMIN
         ! нопедекъел мювюкэмши ьюц
         IF(IndexHDXX.EQ.0) THEN
             HenergyX=0.02D0*DABS(ExfNN-REmin)
			 ! тхйяхпсел нцпюмхвемхе дкъ ьюцю
		     IF(HenergyX.GT.0.1D0) THEN
                HenergyX=0.1D0  
		     ENDIF
            ELSE
             IF(RcoffSD.LT.0.D0) THEN
                 HenergyX=0.1D0*HXXZXL
			    ELSE
			     HenergyX=0.1D0*HXXZXR
		     ENDIF
	     ENDIF 
         ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
         E1XX=REmin*(1+RcoffSD*HenergyX) 
         IndexSUMZKX=0	
         ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
         IndexIIUMNX=0	         			  
         ! оюпюлерп бшундю хг жхйкю 
         IparametrDeterEnergy=1 
         ! оюпюлерп бшундю хг жхйкю
         INDEXEXITRFX=0
         INDEXZZSFAX=0
         RcofXZSDX=1.D0
         DO WHILE(IparametrDeterEnergy.EQ.1) 
            ! хмдейя жхйкю
	        IndexSUMZKX=IndexSUMZKX+1
	        ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1XX
	        call CalculationNormFunctionEnergyXF(E1XX,D1XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
	        WRITE(17,*) E1XX,D1XX,HenergyX
	        ! гюохяшбюел пегскэрюрш пюяверю
	        ExfMX(1)=ExfMX(2)
            DetMX(1)=DetMX(2)
	        ExfMX(2)=ExfMX(3)
            DetMX(2)=DetMX(3)
	        ExfMX(3)=ExfMX(4)
            DetMX(3)=DetMX(4)
	        ExfMX(4)=ExfMX(5)
            DetMX(4)=DetMX(5)
            ExfMX(5)=ExfMX(6)
            DetMX(5)=DetMX(6)
            ExfMX(6)=ExfMX(7)
            DetMX(6)=DetMX(7)
	        ExfMX(7)=E1XX
            DetMX(7)=D1XX
							 
	        IF(IndexSUMZKX.EQ.1) THEN
               ! пюявхршбюел опнхгбндмсч
               DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))
            ENDIF
            IF(IndexSUMZKX.GT.1) THEN
               ! пюявхршбюел опнхгбндмсч
		       DetMAX(1)=DetMAX(2)
		       DetMAX(2)=DetMAX(3)
               DetMAX(3)=DetMAX(4)
		       DetMAX(4)=DetMAX(5)
		       DetMAX(5)=DetMAX(6)
               DetMAX(6)=DetMAX(7)
		       DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))

            
			   ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
			   IF(RcoffSD.LT.0.D0) THEN
		          call CorrectHRightD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			   ENDIF
			   IF(RcoffSD.GT.0.D0) THEN
			      call CorrectHLeftD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			   ENDIF    	  		 
		       ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
		       IF(INDEXZZSFAX.EQ.1) THEN
                  ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
                  E1XX=ExfMX(5) 
                  ExfMX(7)=ExfMX(5) 
                  DetMX(7)=DetMX(5)
                  ExfMX(6)=ExfMX(4) 
                  DetMX(6)=DetMX(4)
			      ExfMX(5)=ExfMX(3) 
                  DetMX(5)=DetMX(3)
			      ExfMX(4)=ExfMX(2) 
                  DetMX(4)=DetMX(2)
			      ExfMX(3)=ExfMX(1) 
                  DetMX(3)=DetMX(1)
                  DetMAX(7)=DetMAX(5)
                  DetMAX(6)=DetMAX(4)
			      DetMAX(5)=DetMAX(3)
			      DetMAX(4)=DetMAX(2)
			      DetMAX(3)=DetMAX(1)
			      ! намскъел бшунд хг жхйкю
			      INDEXEXITRFX=0
				  ! опнбепъел вхякн хрепюжхи
			      IF(IndexSUMZKX.LT.3) THEN
                     ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				     IndexSUMZKX=0 
			      ENDIF
	           ENDIF
	        ENDIF

	        ! мнбне опхакхфемхе дкъ щмепцхх                 
            E1XX=E1XX*(1+RcoffSD*HenergyX)
	        ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
	        IF(INDEXZZSFAX.EQ.0) THEN
               IF(D1XX.GT.RXDF) THEN 
                  INDEXEXITRFX=INDEXEXITRFX+1
		       ENDIF
		       ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
		       IF(INDEXEXITRFX.EQ.5) THEN 
			      IF(RcoffSD.LT.0.D0) THEN 
                     ! кебюъ цпюмхжю нопедекемю
			         REleft=ExfMX(7)
			         RfunLR=DetMX(7)
			      ENDIF
			      IF(RcoffSD.GT.0.D0) THEN 
                     ! опюбюъ цпюмхжю нопедекемю
			         REright=ExfMX(7)
			         RfunLR=DetMX(7)
			      ENDIF    
		          ! сякнбхе бшундю бшонкмемш
                  IparametrDeterEnergy=2
		       ENDIF
	        ENDIF
         ENDDO
	  ENDIF  


   ENDIF


 

   return
 end subroutine MinValueNormEnergyFunctionD2


 
 !! SUBPROGRAMME DETERMINES THE MINIMUM VALUE OF NORM IN THE ENERGY INTERVAL
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! IndexM-SIZE OF MASSIVE
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! H-STEP
═!! Rcof-EQUITY RATIO H * H / 12
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! ExfMin, ExfMax-Energetic interval in which the minimum is searched. The minimum is in the vicinity of ExfMin
═!! REmin-Energy corresponding to a minimum
═!! RfunMIN-Value of the norm
═!! IndexXL-INDEX INDICATES A BORDER OF A MINIMUM BORDER OR NOT (LEFT MINIMUM)
═!! IndexXL = 0-NOT INSTALLED (LEFT MINIMUM)
═!! IndexXL = 1-INSTALLED (LEFT MINIMUM)
═!! EminIndexL-MINIMUM
═!! EmaxIndexL-BORDER OF A MINIMUM AREA
═!! RFUNmaxIndexL-MEANING OF THE NORM IN MINIMUM
═!! IndexXR-INDEX INDICATES A BORDER OF A MINIMUM BORDER OR NOT (RIGHT MINIMUM)
═!! IndexXR = 0-NOT INSTALLED (RIGHT MINIMUM)
═!! IndexXR = 1-INSTALLED (RIGHT MINIMUM)
═!! EminIndexR-BORDER OF A MINIMUM AREA
═!! EmaxIndexR-MINIMUM
═!! RFUNmaxIndexR-VALUE OF NORM IN MINIMUM
═!! IndexMinDifini-INDEX INDICATOR DEFINES A NEW MINIMUM OR FOUND AN RIGHT OR LEFT MINIMUM
═!! IndexMinDifini = 0-NEW MINIMUM DEFINED
═!! IndexMinDifini = 1-DEFINED AGAINST LEFT MINIMUM
═!! IndexMinDifini = 2-DETERMINED AGAIN RIGHT MINIMUM
═!! EminGranL-LEFT BORDER OF A NEW MINIMUM
═!! EminGranR-RIGHT BORDER OF A NEW MINIMUM
═!! IndexHDXX-INDEX INDICATOR DETERMINED THE INITIAL STEP TO DETERMINE BORDERS OF MINIMUM
═!! HXXZXL-STEP FOR LEFT ARE FROM MINIMUM
═!! HXXZXR-STEP FOR THE RIGHT OF THE RIGHT FROM THE MINIMUM
 subroutine MinValueNormEnergyFunctionD3(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,HNN,ExfNN,REmin,RfunMIN,IndexXL,EminIndexL,EmaxIndexL,RFUNmaxIndexL,IndexXR,EminIndexR,EmaxIndexR,RFUNmaxIndexR,IndexMinDifini,EminGranL,EminGranR,IndexHDXX,HXXZXL,HXXZXR)
   implicit none
   integer::NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,IndexXL,IndexXR,IndexMinDifini,IndexHDXX
   real(8)::Rcof,HNN,ExfNN,H,RfunMIN,REmin,EminIndexL,EmaxIndexL,RFUNmaxIndexL,EminIndexR,EmaxIndexR,RFUNmaxIndexR,EminGranL,EminGranR,HXXZXL,HXXZXR
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NEpsIter,IparametrDeterEnergy,IndexRESX,IndexSUMZKX,IndexIIUMNX,INDEXEXITRFX,INDEXZZSFAX,IndexCorGRANF,IndexPovNorma,IndexNORMSDS 
   real(8)::RHS,RX0,RX1,RX2,Rfun0XZ,RnormX0,Xfun0XZ,Rfun1XZ,RnormX1,Xfun1XZ,Rfun2XZ,RnormX2,Xfun2XZ,EXDmin,EXDmax,RcoffSD,RXDF,RfunLR,HenergyX,E1XX,D1XX,RcofXZSDX,DetxfNN   
   real(8),dimension(10)::DetMX,ExfMX,DetMAX

   ! гюмскъел дкъ мювюкю пюяверю
   IndexMinDifini=0
   EminGranL=0.D0
   EminGranR=0.D0
   
   ! мювюкэмше рнвйх
   EXDmin=ExfNN*(1.D0-HNN)
   EXDmax=ExfNN*(1.D0+HNN)

   ! мювюкэмши ьюц
   RHS=DABS(EXDmax-EXDmin)*0.5D0
   ! мювюкэмше рнвйх
   RX0=ExfNN 
   ! анйнбше рнвйх
   RX1=RX0-RHS
   RX2=RX0+RHS
   NEpsIter=0
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   DO WHILE(IparametrDeterEnergy.EQ.1)
      ! хмдейя рхою ялеыемхъ
	  IndexNORMSDS=0 
	  ! явервхй
      NEpsIter=NEpsIter+1
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX0
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX0,Pz,RO1X,Kpoint) 

	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX0,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                  
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX0,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
        
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX1,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX1,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	  
	  ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX1,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
            
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX2
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX2,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX2,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
             
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX2,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
	  
	  ! гюохяшбюел гмювемхе  
      IF(NEpsIter.EQ.1) THEN
           ! оепбюъ хрепюжхъ
		   WRITE(17,*) NEpsIter,RX0,RnormX0 
		   ! оепбюъ хрепюжхъ 
           ! сярюмюбкхбюел онбедемхе мнплш
		   IndexPovNorma=0
		   IF(RnormX1.GT.RnormX0.AND.RnormX0.GT.RnormX2) THEN
              IndexPovNorma=1
		   ENDIF
		   IF(RnormX2.GT.RnormX0.AND.RnormX0.GT.RnormX1) THEN
              IndexPovNorma=1
           ENDIF
           Rfun0XZ=RnormX0
	       Xfun0XZ=RX0
           Rfun1XZ=RnormX1
           Xfun1XZ=RX1
		   Rfun2XZ=RnormX2
           Xfun2XZ=RX2 
		   IF(IndexPovNorma.EQ.1) THEN
		       RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
              ELSE
			   RX0=Xfun0XZ-RHS*2.D0
		   ENDIF
           RX1=RX0-RHS
		   RX2=RX0+RHS
		   ! гюохяшбюел гмювемхе мнплш дкъ дюкэмеиьецн пюяверю
		   DetxfNN=RnormX0 
		 ELSE
           
		   ! опнбепъел оноюдюер дюммне опхакхфемхе он щмепцхх б накюярэ сфе мюидеммнцн лхмхлслю дкъ кебнцн лхмхлслю 
           IF(IndexXL.EQ.1) THEN
              ! цпюмхжю накюярх нопедекемю
			  IF(RX0.GT.EminIndexL.AND.RX0.LT.EmaxIndexL) THEN
                 ! щмепцхъ оноюкю б накюярэ опедшдсыецн лхмхлслю
                 ! лхмхлюкэмне гмювемхе тсмйжхх
		         RfunMIN=RFUNmaxIndexL
		         ! яннрберярбсчыхи юпцслемр
                 REmin=EminIndexL
				 ! сякнбхе бшундю бшонкмемн
                 IparametrDeterEnergy=2 
				 WRITE(17,*) NEpsIter,REmin,RfunMIN
				 ! тхйяхпсел, врн намюпсфем сфе хгбеярмши кебши лхмхлсл
                 IndexMinDifini=1
			  ENDIF
		   ENDIF
		   
		   ! опнбепъел оноюдюер дюммне опхакхфемхе он щмепцхх б накюярэ сфе мюидеммнцн лхмхлслю дкъ опюбнцн лхмхлслю 
           IF(IndexXR.EQ.1) THEN
              ! цпюмхжю накюярх нопедекемю
			  IF(RX0.GT.EminIndexR.AND.RX0.LT.EmaxIndexR) THEN
                 ! щмепцхъ оноюкю б накюярэ опедшдсыецн лхмхлслю
                 ! лхмхлюкэмне гмювемхе тсмйжхх
		         RfunMIN=RFUNmaxIndexR
		         ! яннрберярбсчыхи юпцслемр
                 REmin=EmaxIndexR
				 ! сякнбхе бшундю бшонкмемн
                 IparametrDeterEnergy=2 
				 WRITE(17,*) NEpsIter,REmin,RfunMIN
				 ! тхйяхпсел, врн намюпсфем сфе хгбеярмши опюбши лхмхлсл
                 IndexMinDifini=2
			  ENDIF
		   ENDIF
           
		   ! опнбепъел ме нясыеярбкем бшунд 
		   IF(IparametrDeterEnergy.NE.2) THEN
		      WRITE(17,*) NEpsIter,RX0,RnormX0          
		      ! опнбепъел ндхмюйнбше гмювемхъ тсмйжхх лш онксвюел
		      ! якеднбюрекэмюъ оепбюъ опнхгбндмюъ пюбмю мскч	 
              IF(RnormX0.EQ.Rfun0XZ) THEN
			     ! опнбепъел мнлеп хрепюжхх
				 IF(NEpsIter.EQ.2) THEN
                     IndexNORMSDS=1
					ELSE 
					 ! лхмхлюкэмне гмювемхе тсмйжхх
		             RfunMIN=RnormX0 
		             ! яннрберярбсчыхи юпцслемр
                     REmin=RX0
			         ! сякнбхе бшундю бшонкмемн
                     IparametrDeterEnergy=2 
                 ENDIF
              ENDIF 

		
		      ! опнбепъел онксвеммне мнбне опхакхфемхе лемэье опедедсыецн хкх мер
		      IF(RnormX0.LT.Rfun0XZ) THEN
                  ! лхмхлюкэмне гмювемхе тсмйжхх
		          RfunMIN=RnormX0 
		          ! яннрберярбсчыхи юпцслемр
                  REmin=RX0
              
			      ! опнбепъел сякнбхе бшундю (опнбепъел пюбемярбн мскч оепбни опнгбндмни)
			      IF(DABS((RfunMIN-Rfun0XZ)/(REmin-Xfun0XZ)).LT.10.D0**(-10)) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
			         ENDIF 
			      ENDIF   
              
			      ! гюохяшбюел опедшдсыее гмювемхе
			      Rfun0XZ=RnormX0
			      Xfun0XZ=RX0
                  Rfun1XZ=RnormX1
                  Xfun1XZ=RX1
			      Rfun2XZ=RnormX2
                  Xfun2XZ=RX2 

                  ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                  IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
				         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=RnormX0 
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=RX0
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF 
				  ENDIF 
				  !опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				  IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
				         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=RnormX0 
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=RX0
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF
				  ENDIF

			
                  ! мнбне опхакхфемхе
				  IF(IparametrDeterEnergy.EQ.1) THEN
				     ! опнбепъел хмдейя рхою ялеыемхъ
					 IF(IndexNORMSDS.EQ.0) THEN
			             RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                        ELSE
                         RX0=Xfun0XZ-RHS*2.D0
					 ENDIF
                     RX1=RX0-RHS
			         RX2=RX0+RHS
			      ENDIF
			
			    ELSE
				  ! гмювемхе тсмйжхх пюярер
				  ! лемъел ьюц
				  RHS=-RHS*0.5D0
			

				  ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                  IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
				         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=Rfun0XZ
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=Xfun0XZ
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF
				  ENDIF 
				  !опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				  IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
				     ! опнбепъел мнлеп хрепюжхх
				     IF(NEpsIter.EQ.2) THEN
                         IndexNORMSDS=1
					    ELSE 
                         ! лхмхлюкэмне гмювемхе тсмйжхх
		                 RfunMIN=Rfun0XZ
		                 ! яннрберярбсчыхи юпцслемр
                         REmin=Xfun0XZ
			             ! сякнбхе бшундю бшонкмемн
                         IparametrDeterEnergy=2 
				     ENDIF 
				  ENDIF

				  ! мнбне опхакхфемхе
				  IF(IparametrDeterEnergy.EQ.1) THEN
				     ! опнбепъел хмдейя рхою ялеыемхъ
					 IF(IndexNORMSDS.EQ.0) THEN
				         RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                        ELSE
                         RX0=Xfun0XZ-RHS*4.D0
					 ENDIF 
                     RX1=RX0-RHS
			         RX2=RX0+RHS
			      ENDIF
			  ENDIF
           ENDIF
	  ENDIF
   ENDDO

   ! опнбепъел мюидем мнбши лхмхлсл хкх намюпсфем ярюпши
   IF(IndexMinDifini.EQ.0) THEN
      ! сярюмюбкхбюел йюйсч цпюмхжс мсфмн нопедекхрэ
	  IF(REmin.GT.ExfNN) THEN
           ! мсфмн нопедекхрэ опюбсч цпюмхжс
		   ! кебюъ цпюмхжю
           EminGranL=ExfNN
		   ! оюпюлерп ядбхцю 
		   RcoffSD=1.D0
		 ELSE
           ! мсфмн нопедекхрэ кебсч цпюмхжс
		   ! опюбюъ цпюмхжю
           EminGranR=ExfNN
		   ! оюпюлерп ядбхцю
           RcoffSD=-1.D0
	  ENDIF

	  ! нопедекъел рс цпюмхжс йнрнпсч менопедекхкх
      ! нопедекъел кебсч цпюмхжс лхмхлслю
      WRITE(17,*)
	  IF(RcoffSD.LT.0.D0) THEN
         WRITE(17,*) 'Defeni grany left min'
      ENDIF
	  IF(RcoffSD.GT.0.D0) THEN
         WRITE(17,*) 'Defeni grany right min' 
	  ENDIF

	  
	  ! сярюмюбкхбюел йпхрепхи дкъ нопедекемхъ кебни цпюмхжш лхмхлслю
      IF(0.2D0*RfunMIN.LT.0.5D0) THEN
           RXDF=1.2D0*RfunMIN
         ELSE
 	       RXDF=RfunMIN+0.2D0
      ENDIF

	 
      ! гюмскъел оепед пюявернл
      DetMX=0.D0
      DetMAX=0.D0
      ExfMX=0.D0
      ! гюохяшбюел оепбне гмювемхе
      ExfMX(7)=REmin
      DetMX(7)=RfunMIN
      ! нопедекъел мювюкэмши ьюц
	  IF(IndexHDXX.EQ.0) THEN
          HenergyX=0.02D0*DABS(ExfNN-REmin)
		  ! тхйяхпсел нцпюмхвемхе дкъ ьюцю
		  IF(HenergyX.GT.0.1D0) THEN
             HenergyX=0.1D0  
		  ENDIF
         ELSE
          IF(RcoffSD.LT.0.D0) THEN
              HenergyX=0.1D0*HXXZXL
		     ELSE
		      HenergyX=0.1D0*HXXZXR
		  ENDIF
	  ENDIF 
      ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
      E1XX=REmin*(1+RcoffSD*HenergyX) 
      IndexSUMZKX=0	
      ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
      IndexIIUMNX=0	         			  
      ! оюпюлерп бшундю хг жхйкю 
      IparametrDeterEnergy=1 
      ! оюпюлерп бшундю хг жхйкю
      INDEXEXITRFX=0
      INDEXZZSFAX=0
      RcofXZSDX=1.D0
      DO WHILE(IparametrDeterEnergy.EQ.1) 
         ! хмдейя жхйкю
	     IndexSUMZKX=IndexSUMZKX+1
	     ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1XX
	     call CalculationNormFunctionEnergyXF(E1XX,D1XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
	     WRITE(17,*) E1XX,D1XX,HenergyX
	     ! гюохяшбюел пегскэрюрш пюяверю
	     ExfMX(1)=ExfMX(2)
         DetMX(1)=DetMX(2)
	     ExfMX(2)=ExfMX(3)
         DetMX(2)=DetMX(3)
	     ExfMX(3)=ExfMX(4)
         DetMX(3)=DetMX(4)
	     ExfMX(4)=ExfMX(5)
         DetMX(4)=DetMX(5)
         ExfMX(5)=ExfMX(6)
         DetMX(5)=DetMX(6)
         ExfMX(6)=ExfMX(7)
         DetMX(6)=DetMX(7)
	     ExfMX(7)=E1XX
         DetMX(7)=D1XX
							 
	     IF(IndexSUMZKX.EQ.1) THEN
            ! пюявхршбюел опнхгбндмсч
            DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))
         ENDIF
         IF(IndexSUMZKX.GT.1) THEN
            ! пюявхршбюел опнхгбндмсч
		    DetMAX(1)=DetMAX(2)
		    DetMAX(2)=DetMAX(3)
            DetMAX(3)=DetMAX(4)
		    DetMAX(4)=DetMAX(5)
		    DetMAX(5)=DetMAX(6)
            DetMAX(6)=DetMAX(7)
		    DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))

           	! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
			IF(RcoffSD.LT.0.D0) THEN
		       call CorrectHRightD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			ENDIF
			IF(RcoffSD.GT.0.D0) THEN
			   call CorrectHLeftD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			ENDIF   
		    ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
		    IF(INDEXZZSFAX.EQ.1) THEN
               ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
               E1XX=ExfMX(5) 
               ExfMX(7)=ExfMX(5) 
               DetMX(7)=DetMX(5)
               ExfMX(6)=ExfMX(4) 
               DetMX(6)=DetMX(4)
			   ExfMX(5)=ExfMX(3) 
               DetMX(5)=DetMX(3)
			   ExfMX(4)=ExfMX(2) 
               DetMX(4)=DetMX(2)
			   ExfMX(3)=ExfMX(1) 
               DetMX(3)=DetMX(1)
               DetMAX(7)=DetMAX(5)
               DetMAX(6)=DetMAX(4)
			   DetMAX(5)=DetMAX(3)
			   DetMAX(4)=DetMAX(2)
			   DetMAX(3)=DetMAX(1)
			   ! намскъел бшунд хг жхйкю
			   INDEXEXITRFX=0
			   ! опнбепъел вхякн хрепюжхи
			   IF(IndexSUMZKX.LT.3) THEN
                  ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				  IndexSUMZKX=0 
			   ENDIF
	        ENDIF
	     ENDIF

	     ! мнбне опхакхфемхе дкъ щмепцхх                 
         E1XX=E1XX*(1+RcoffSD*HenergyX)
	     ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
	     IF(INDEXZZSFAX.EQ.0) THEN
            IF(D1XX.GT.RXDF) THEN 
               INDEXEXITRFX=INDEXEXITRFX+1
		    ENDIF
		    ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
		    IF(INDEXEXITRFX.EQ.5) THEN              
			   IF(RcoffSD.LT.0.D0) THEN 
                  ! кебюъ цпюмхжю нопедекемю
			      EminGranL=ExfMX(7)
			      RfunLR=DetMX(7)
			   ENDIF
			   IF(RcoffSD.GT.0.D0) THEN 
                  ! опюбюъ цпюмхжю нопедекемю
			      EminGranR=ExfMX(7)
			      RfunLR=DetMX(7)
			   ENDIF    
		       ! сякнбхе бшундю бшонкмемш
               IparametrDeterEnergy=2
		    ENDIF
	     ENDIF
      ENDDO
      
      
	  ! опнбепъел цпюмхжс йнрнпюъ хгмювюкэмн ашкю нопедекемю
      ! опнбепъел дпсцюъ (кебюъ хкх опюбюъ) цпюмхжю яннрберярбсер сякнбхч хкх мер
      ! сярюмюбкхбюел йюйсч цпюмхжс мсфмн нопедекхрэ
	  ! йпхрепхи йнпейрхпнбйх цпюмхжш
	  IndexCorGRANF=0
	  IF(REmin.GT.ExfNN) THEN
	       ! опнбепъел бшонкмъеряъ сякнбхе дкъ кебни цпюмхжш
		   IF(DetxfNN.LT.RXDF) THEN
              ! мнплю б кебни цпюмхже лемэье йпхрепхъ менаундхлн йнппейрхпнбйю кебни цпюмхжш
			  ! оюпюлерп ядбхцю 
		      RcoffSD=-1.D0
			  IndexCorGRANF=1   
		   ENDIF		   
		 ELSE
           ! опнбепъел бшонкмъеряъ сякнбхе дкъ опюбни цпюмхжш
		   IF(DetxfNN.LT.RXDF) THEN
              ! мнплю б опюбни цпюмхже лемэье йпхрепхъ менаундхлн йнппейрхпнбйю опюбни цпюмхжш
			  ! оюпюлерп ядбхцю 
		      RcoffSD=1.D0
			  IndexCorGRANF=1   
		   ENDIF	
	  ENDIF

	  ! опнбепъел менаундхлн йнпейрхпнбйю цпюмхжш хкх мер
	  IF(IndexCorGRANF.EQ.1) THEN
	     
		 ! нопедекъел рс цпюмхжс йнрнпсч менопедекхкх
         ! нопедекъел кебсч цпюмхжс лхмхлслю
         WRITE(17,*)
	     IF(RcoffSD.LT.0.D0) THEN
            WRITE(17,*) 'Defeni grany left min'
         ENDIF
	     IF(RcoffSD.GT.0.D0) THEN
            WRITE(17,*) 'Defeni grany right min' 
	     ENDIF

         ! гюмскъел оепед пюявернл
         DetMX=0.D0
         DetMAX=0.D0
         ExfMX=0.D0
         ! гюохяшбюел оепбне гмювемхе
         ExfMX(7)=REmin
         DetMX(7)=RfunMIN
         ! нопедекъел мювюкэмши ьюц
		 IF(IndexHDXX.EQ.0) THEN
             HenergyX=0.02D0*DABS(ExfNN-REmin)
			 ! тхйяхпсел нцпюмхвемхе дкъ ьюцю
		     IF(HenergyX.GT.0.1D0) THEN
                HenergyX=0.1D0  
		     ENDIF
            ELSE
             IF(RcoffSD.LT.0.D0) THEN
                 HenergyX=0.1D0*HXXZXL
		        ELSE
		         HenergyX=0.1D0*HXXZXR
		     ENDIF
	     ENDIF 
         ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
         E1XX=REmin*(1+RcoffSD*HenergyX) 
         IndexSUMZKX=0	
         ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
         IndexIIUMNX=0	         			  
         ! оюпюлерп бшундю хг жхйкю 
         IparametrDeterEnergy=1 
         ! оюпюлерп бшундю хг жхйкю
         INDEXEXITRFX=0
         INDEXZZSFAX=0
         RcofXZSDX=1.D0
         DO WHILE(IparametrDeterEnergy.EQ.1) 
            ! хмдейя жхйкю
	        IndexSUMZKX=IndexSUMZKX+1
	        ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1XX
	        call CalculationNormFunctionEnergyXF(E1XX,D1XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
	        WRITE(17,*) E1XX,D1XX,HenergyX
	        ! гюохяшбюел пегскэрюрш пюяверю
	        ExfMX(1)=ExfMX(2)
            DetMX(1)=DetMX(2)
	        ExfMX(2)=ExfMX(3)
            DetMX(2)=DetMX(3)
	        ExfMX(3)=ExfMX(4)
            DetMX(3)=DetMX(4)
	        ExfMX(4)=ExfMX(5)
            DetMX(4)=DetMX(5)
            ExfMX(5)=ExfMX(6)
            DetMX(5)=DetMX(6)
            ExfMX(6)=ExfMX(7)
            DetMX(6)=DetMX(7)
	        ExfMX(7)=E1XX
            DetMX(7)=D1XX
							 
	        IF(IndexSUMZKX.EQ.1) THEN
               ! пюявхршбюел опнхгбндмсч
               DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))
            ENDIF
            IF(IndexSUMZKX.GT.1) THEN
               ! пюявхршбюел опнхгбндмсч
		       DetMAX(1)=DetMAX(2)
		       DetMAX(2)=DetMAX(3)
               DetMAX(3)=DetMAX(4)
		       DetMAX(4)=DetMAX(5)
		       DetMAX(5)=DetMAX(6)
               DetMAX(6)=DetMAX(7)
		       DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))

               
			   ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
			   IF(RcoffSD.LT.0.D0) THEN
		          call CorrectHRightD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			   ENDIF
			   IF(RcoffSD.GT.0.D0) THEN
			      call CorrectHLeftD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
			   ENDIF      		 
		       ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
		       IF(INDEXZZSFAX.EQ.1) THEN
                  ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
                  E1XX=ExfMX(5) 
                  ExfMX(7)=ExfMX(5) 
                  DetMX(7)=DetMX(5)
                  ExfMX(6)=ExfMX(4) 
                  DetMX(6)=DetMX(4)
			      ExfMX(5)=ExfMX(3) 
                  DetMX(5)=DetMX(3)
			      ExfMX(4)=ExfMX(2) 
                  DetMX(4)=DetMX(2)
			      ExfMX(3)=ExfMX(1) 
                  DetMX(3)=DetMX(1)
                  DetMAX(7)=DetMAX(5)
                  DetMAX(6)=DetMAX(4)
			      DetMAX(5)=DetMAX(3)
			      DetMAX(4)=DetMAX(2)
			      DetMAX(3)=DetMAX(1)
			      ! намскъел бшунд хг жхйкю
			      INDEXEXITRFX=0
				  ! опнбепъел вхякн хрепюжхи
			      IF(IndexSUMZKX.LT.3) THEN
                     ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
				     IndexSUMZKX=0 
			      ENDIF
	           ENDIF
	        ENDIF

	        ! мнбне опхакхфемхе дкъ щмепцхх                 
            E1XX=E1XX*(1+RcoffSD*HenergyX)
	        ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
	        IF(INDEXZZSFAX.EQ.0) THEN
               IF(D1XX.GT.RXDF) THEN 
                  INDEXEXITRFX=INDEXEXITRFX+1
		       ENDIF
		       ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
		       IF(INDEXEXITRFX.EQ.5) THEN 
			       IF(RcoffSD.LT.0.D0) THEN 
                      ! кебюъ цпюмхжю нопедекемю
			          EminGranL=ExfMX(7)
			          RfunLR=DetMX(7)
			       ENDIF
			       IF(RcoffSD.GT.0.D0) THEN 
                      ! опюбюъ цпюмхжю нопедекемю
			          EminGranR=ExfMX(7)
			          RfunLR=DetMX(7)
			       ENDIF    
		           ! сякнбхе бшундю бшонкмемш
                   IparametrDeterEnergy=2
		       ENDIF
	        ENDIF
         ENDDO
	  ENDIF  


   ENDIF 

   return
 end subroutine MinValueNormEnergyFunctionD3




!! SUBPROGRAMME DETERMINES THE MINIMUM VALUE OF NORM IN THE ENERGY INTERVAL
═!! NMO-NUMBER OF MOLECULAR ORBITALS
═!! Npoint-NUMBER OF POINTS
═!! IndexM-SIZE OF MASSIVE
═!! Kpoint-NUMBER OF POINT OF SOLVING SOLUTIONS
═!! NumeroMin-MINIMUM HARMONIC NUMBER
═!! H-STEP
═!! Rcof-EQUITY RATIO H * H / 12
═!! Pz (NumbreGarmMO (NMO), NumbreGarmMO (NMO), Npoint) - VALUE MAP of the P (ro) -particle of the log function P entering the equation Y '' = P * Y + Q-LOCAL POTENTIAL
═!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
═!! A (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! B (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING FUNCTION PARAMETERS
═!! AA (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! BB (Isize, Isize, Npoint) -MASSIVE FOR CALCULATING THE PARAMETERS OF THE PROLONG WITH ADDED ESSENTIAL SINGLE-ELECTRONIC ENERGY
═!! FX (IndexM, Npoint) -MASSIVE F
═!! Vout (IndexM, IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Uout (IndexM, Npoint) -MASSIVE OF FREQUENCY ON THE DIRECTION OF "OUTSIDE"
═!! Vin (IndexM, IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! Uin (IndexM, Npoint) -MASSIVE OF FIRING ON THE DIRECTION "INSIDE"
═!! VNout, UNout, VNin, UNIN, AZ, BZ, FZ, DZ, E, AX, BX, CX, IXM-AUXILIARY MASSIVE FOR CALCULATION
═!! FunMO-IN THIS CASE AUXILIARY MASSIVE
═!! NumeroMaxLimid-MAXIMUM NUMBER OF HARMONIC MOLECULAR ORBITALS
═!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
═!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
═!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
═!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
═!! ExfMin, ExfMax-Energetic interval in which the minimum is searched. The minimum is in the vicinity of ExfMin
═!! REmin-Energy corresponding to a minimum
═!! RfunMIN-Value of the norm
═!! REleft-Left Low Limit
═!! IndexHDXX-INDEX INDICATOR DETERMINED THE INITIAL STEP TO DETERMINE BORDERS OF MINIMUM
═!! HXXZXL-STEP FOR LEFT ARE FROM MINIMUM
 subroutine MinValueNormEnergyFunctionDD(NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,Rcof,H,Uout,Uin,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,AlfaCoffMO,RfunLigands,RO1X,FunMO,Pz,IXM,UNout,UNin,Vout,Vin,AX,BX,CX,E,AZ,FZ,BZ,FX,VNout,VNin,DZ,A,B,AA,BB,ExfMin,Exf0,ExfMax,REmin,RfunMIN,REleft,IndexHDXX,HXXZXL)
   implicit none
   integer::NMO,Npoint,IndexM,NumeroMin,NumeroMax,Kpoint,NumeroMaxLimid,IndexHDXX
   real(8)::Rcof,ExfMin,Exf0,ExfMax,H,RfunMIN,REmin,REleft,HXXZXL
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::RO1X,UNout,UNin,AX,BX,CX
   real(8),dimension(:,:)::E,AZ,FZ,BZ,FX,Uout,Uin,VNout,VNin,DZ,FunMO
   real(8),dimension(:,:,:)::Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NEpsIter,IparametrDeterEnergy,IndexRESX,IndexSUMZKX,IndexIIUMNX,INDEXEXITRFX,INDEXZZSFAX 
   real(8)::RHS,RX0,RX1,RX2,Rfun0XZ,RnormX0,Xfun0XZ,Rfun1XZ,RnormX1,Xfun1XZ,Rfun2XZ,RnormX2,Xfun2XZ,RcofXZSDX 
   real(8)::RfunLeft,RXDF,HenergyX,E1XX,D1XX  
   real(8),dimension(10)::DetMX,ExfMX,DetMAX
 
   ! мювюкэмши ьюц
   RHS=DABS(ExfMin-Exf0)
   ! мювюкэмше рнвйх
   RX0=Exf0 
   ! анйнбше рнвйх
   RX1=ExfMin
   RX2=ExfMax
   NEpsIter=0
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   DO WHILE(IparametrDeterEnergy.EQ.1)
      NEpsIter=NEpsIter+1
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX0
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX0,Pz,RO1X,Kpoint) 

	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX0,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
                  
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX0,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
        
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX1
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX1,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX1,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	  
	  ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
            
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX1,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
            
      ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх RX2
	  ! нопедекъел мнлеп рнвйх яьхбйх пеьемхи OUT AND IN
	  call PointIntersection(Npoint,NumeroMin,RX2,Pz,RO1X,Kpoint)
	  
	  ! ондопнцпюллю тнплхпсер люяяхбш дкъ пюявербю оюпюлерпнб опнцнмйх
	  call FormMassivAABB(Npoint,IndexM,Rcof,RX2,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
	 
      ! нясыеярбкъел пюявер оюпюлерпнб опнцнмйх дкъ мендмнпндмни яхярелш
      call CalculationParameterProgonky(2,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
             
      ! нопедекъел гмювемхе мнплш яннрберярбсчыее дюммнлс гмювемхч щмепцхх RX0
      call CalculationNormFunE(RnormX2,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)    
     
	  
	  ! гюохяшбюел гмювемхе  
      IF(NEpsIter.EQ.1) THEN
           ! оепбюъ хрепюжхъ 
           Rfun0XZ=RnormX0
	       Xfun0XZ=RX0
           Rfun1XZ=RnormX1
           Xfun1XZ=RX1
		   Rfun2XZ=RnormX2
           Xfun2XZ=RX2 
		   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
           RX1=RX0-RHS
		   RX2=RX0+RHS
		 ELSE
           
		   WRITE(17,*) NEpsIter,RX0,RnormX0          
		   ! опнбепъел ндхмюйнбше гмювемхъ тсмйжхх лш онксвюел
		   ! якеднбюрекэмюъ оепбюъ опнхгбндмюъ пюбмю мскч	 
           IF(RnormX0.EQ.Rfun0XZ) THEN
              ! лхмхлюкэмне гмювемхе тсмйжхх
		      RfunMIN=RnormX0 
		      ! яннрберярбсчыхи юпцслемр
              REmin=RX0
			  ! сякнбхе бшундю бшонкмемн
              IparametrDeterEnergy=2
           ENDIF 

		
		   ! опнбепъел онксвеммне мнбне опхакхфемхе лемэье опедедсыецн хкх мер
		   IF(RnormX0.LT.Rfun0XZ) THEN
                ! лхмхлюкэмне гмювемхе тсмйжхх
		        RfunMIN=RnormX0 
		        ! яннрберярбсчыхи юпцслемр
                REmin=RX0

			              
			    ! опнбепъел сякнбхе бшундю (опнбепъел пюбемярбн мскч оепбни опнгбндмни)
			    IF(DABS((RfunMIN-Rfun0XZ)/(REmin-Xfun0XZ)).LT.10.D0**(-10)) THEN
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
			    ENDIF   
              
			    ! гюохяшбюел опедшдсыее гмювемхе
			    Rfun0XZ=RnormX0
			    Xfun0XZ=RX0
                Rfun1XZ=RnormX1
                Xfun1XZ=RX1
			    Rfun2XZ=RnormX2
                Xfun2XZ=RX2 

                ! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				     ! лхмхлюкэмне гмювемхе тсмйжхх
		             RfunMIN=RnormX0 
		             ! яннрберярбсчыхи юпцслемр
                     REmin=RX0
			         ! сякнбхе бшундю бшонкмемн
                     IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=RnormX0 
		           ! яннрберярбсчыхи юпцслемр
                   REmin=RX0
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

			
                ! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
			       RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
			
			  ELSE
				! гмювемхе тсмйжхх пюярер
				! лемъел ьюц
				RHS=-RHS*0.5D0
			

				! опнбепъел пюбмю мскч оепбюъ опнхгбндмюъ хкх мер 
                IF(Rfun2XZ.EQ.Rfun1XZ) THEN
				   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF 
				!опнбепъел пюбмю мскч брнпюъ опнхгбндмюъ  опх щрнл пюявер нярюмюбкхбюеряъ
				IF(DABS((Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)/RHS**2).LT.10.D0**(-10)) THEN
                   ! лхмхлюкэмне гмювемхе тсмйжхх
		           RfunMIN=Rfun0XZ
		           ! яннрберярбсчыхи юпцслемр
                   REmin=Xfun0XZ
			       ! сякнбхе бшундю бшонкмемн
                   IparametrDeterEnergy=2 
				ENDIF

				! мнбне опхакхфемхе
				IF(IparametrDeterEnergy.EQ.1) THEN
				   RX0=Xfun0XZ-RHS*(Rfun2XZ-Rfun1XZ)*0.5D0/(Rfun2XZ-2.D0*Rfun0XZ+Rfun1XZ)
                   RX1=RX0-RHS
			       RX2=RX0+RHS
			    ENDIF
				
		   ENDIF
	  ENDIF
   ENDDO

   ! нопедекъел кебсч цпюмхжс лхмхлслю
   WRITE(17,*)
   WRITE(17,*) 'Defeni grany left min'
 
   ! сярюмюбкхбюел йпхрепхи дкъ нопедекемхъ кебни цпюмхжш лхмхлслю
   IF(0.2D0*RfunMIN.LT.0.5D0) THEN
        RXDF=1.2D0*RfunMIN
      ELSE
 	    RXDF=RfunMIN+0.2D0
   ENDIF
    
   ! нясыеярбкъел нопедекемхе кебни цпюмхжш хмрепбюкю 
   ! гюмскъел оепед пюявернл
   DetMX=0.D0
   DetMAX=0.D0
   ExfMX=0.D0
   ! гюохяшбюел оепбне гмювемхе
   ExfMX(7)=REmin
   DetMX(7)=RfunMIN
   ! нопедекъел мювюкэмши ьюц
   IF(IndexHDXX.EQ.0) THEN
       HenergyX=0.02D0*DABS(Exf0-REmin)
	   ! тхйяхпсел нцпюмхвемхе дкъ ьюцю
	   IF(HenergyX.GT.0.1D0) THEN
          HenergyX=0.1D0  
	   ENDIF
      ELSE
	   HenergyX=0.1D0*HXXZXL
   ENDIF
   ! нопедекъел мювюкэмне гмювемхе дкъ щмепцхх 
   E1XX=REmin*(1-HenergyX) 
   IndexSUMZKX=0	
   ! хмдейя сйюгшбюер яйнкэйн пюг нясыеярбкъкняэ слемэьемхе б 10 пюг ьюцю
   IndexIIUMNX=0	         			  
   ! оюпюлерп бшундю хг жхйкю 
   IparametrDeterEnergy=1 
   ! оюпюлерп бшундю хг жхйкю
   INDEXEXITRFX=0
   INDEXZZSFAX=0
   RcofXZSDX=1.D0
   DO WHILE(IparametrDeterEnergy.EQ.1) 
      ! хмдейя жхйкю
	  IndexSUMZKX=IndexSUMZKX+1
	  ! мюундхл гмювемхе мнплш дкъ гмювемхъ щмепцхх E1XX
	  call CalculationNormFunctionEnergyXF(E1XX,D1XX,NMO,NumeroMaxLimid,Npoint,IndexM,Kpoint,NumeroMin,NumeroMax,Rcof,H,NumbreLigand,IXM,NLigands,NFunLigands,NumbreFunctionLig,RO1X,UNout,UNin,AX,BX,CX,FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E,FunMO,Pz,A,B,AA,BB,Vout,Vin,AlfaCoffMO,RfunLigands)
	  WRITE(17,*) E1XX,D1XX,HenergyX
	  ! гюохяшбюел пегскэрюрш пюяверю
	  ExfMX(1)=ExfMX(2)
      DetMX(1)=DetMX(2)
	  ExfMX(2)=ExfMX(3)
      DetMX(2)=DetMX(3)
	  ExfMX(3)=ExfMX(4)
      DetMX(3)=DetMX(4)
	  ExfMX(4)=ExfMX(5)
      DetMX(4)=DetMX(5)
      ExfMX(5)=ExfMX(6)
      DetMX(5)=DetMX(6)
      ExfMX(6)=ExfMX(7)
      DetMX(6)=DetMX(7)
	  ExfMX(7)=E1XX
      DetMX(7)=D1XX
							 
	  IF(IndexSUMZKX.EQ.1) THEN
         ! пюявхршбюел опнхгбндмсч
         DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))
      ENDIF
      IF(IndexSUMZKX.GT.1) THEN
         ! пюявхршбюел опнхгбндмсч
		 DetMAX(1)=DetMAX(2)
		 DetMAX(2)=DetMAX(3)
         DetMAX(3)=DetMAX(4)
		 DetMAX(4)=DetMAX(5)
		 DetMAX(5)=DetMAX(6)
         DetMAX(6)=DetMAX(7)
		 DetMAX(7)=DetMX(7)*(1.D0-DetMX(6)/DetMX(7))/(ExfMX(7)*(1.D0-ExfMX(6)/ExfMX(7)))

         ! опнбндхл юмюкхг х йнппейрхпнбйс ьюцю
		 call CorrectHRightD(IndexSUMZKX,INDEXZZSFAX,IndexIIUMNX,HenergyX,RcofXZSDX,DetMAX)
		
					 
		 ! опнбепъел менаундхлн бепмсряъ мю дбе рнвйх мюгюд хкх мер
		 IF(INDEXZZSFAX.EQ.1) THEN
            ! менаундхлн бепмсряъ мю дбе рнвйх мюгюд
            E1XX=ExfMX(5) 
            ExfMX(7)=ExfMX(5) 
            DetMX(7)=DetMX(5)
            ExfMX(6)=ExfMX(4) 
            DetMX(6)=DetMX(4)
			ExfMX(5)=ExfMX(3) 
            DetMX(5)=DetMX(3)
			ExfMX(4)=ExfMX(2) 
            DetMX(4)=DetMX(2)
			ExfMX(3)=ExfMX(1) 
            DetMX(3)=DetMX(1)
            DetMAX(7)=DetMAX(5)
            DetMAX(6)=DetMAX(4)
			DetMAX(5)=DetMAX(3)
			DetMAX(4)=DetMAX(2)
			DetMAX(3)=DetMAX(1)
			! намскъел бшунд хг жхйкю
			INDEXEXITRFX=0
			! опнбепъел вхякн хрепюжхи
			IF(IndexSUMZKX.LT.3) THEN
               ! намскъел явервхй оняйнкэйс оняке бнгбпюыемхъ мю дбе рнвйх лш оноюдюел б мювюкэмсч рнвйс
			   IndexSUMZKX=0 
			ENDIF
	     ENDIF
	  ENDIF

	  ! мнбне опхакхфемхе дкъ щмепцхх                 
      E1XX=E1XX*(1-HenergyX)
	  ! опнбепъел бшонкмхкняэ сякнбхе хкх мер
	  IF(INDEXZZSFAX.EQ.0) THEN
         IF(D1XX.GT.RXDF) THEN 
            INDEXEXITRFX=INDEXEXITRFX+1
		 ENDIF
		 ! опнбепъел бшонкмемн сякнбхе бшундю хкх мер
		 IF(INDEXEXITRFX.EQ.5) THEN 
			! кебюъ цпюмхжю нопедекемю
			REleft=ExfMX(7)
            RfunLeft=DetMX(7)
		    ! сякнбхе бшундю бшонкмемш
            IparametrDeterEnergy=2
		 ENDIF
	  ENDIF
   ENDDO
    

 

   return
 end subroutine MinValueNormEnergyFunctionDD




 
 ! ондопнцпюллю йнмрпнкъ ьюцю
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! IndexSmena-хмдейя хглемемхъ ьюцю
 ! NumbreMaxRegion-люйяхлюкэмне вхякн хмрепбюкнб йнпмеи
 ! IndexRegion- мнлеп накюярх йнпмъ
 ! EnergyRegion(NMO,IndexRegion,2)-люяяхб хмрепбюкнб йнпмеи
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! Henergy-ьюц
 subroutine ControlH(NMO,IndexSmena,NumbreMaxRegion,IndexRegion,EnergyRegion,Exf,Henergy)
   implicit none
   integer::NMO,IndexSmena,NumbreMaxRegion,IndexRegion
   real(8)::Exf,Henergy
   real(8),dimension(:,:,:)::EnergyRegion
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IKLDF,IIIUU


   ! опнбепъел мер опхбшьемхъ он вхякс хмрепбюкнб йнпмеи
   IF(IndexRegion.GT.NumbreMaxRegion) THEN
      ! бшундхл хг-гю рнцн, врн бяе хмрепбюкш опнялнрпемш
	  return
   ENDIF
   
   IKLDF=0
   IF(IndexSmena.EQ.0) THEN
      IKLDF=1
   ENDIF
   IF(IndexSmena.EQ.1) THEN
      IKLDF=2
   ENDIF 
   
      
   ! опнбепъел щмепцхъ мюундхряъ бшье цпюмхжш хмрепбюкю
   IF(Exf.GT.EnergyRegion(NMO,IndexRegion,2)) THEN
      ! пефхл слемэьемхъ ьюцю нрйкчвхрэ
	  IndexSmena=0  
   ENDIF
   ! щмепцхъ мюундхряъ б хмрепбюке
   IF(Exf.LT.EnergyRegion(NMO,IndexRegion,2).AND.Exf.GT.EnergyRegion(NMO,IndexRegion,1)) THEN
      ! щмепцхъ оноюкю б хмрепбюк
      ! нясыеярбхрэ слемэьемхе ьюцю
      IndexSmena=1
   ENDIF

   ! щмепцхъ мюундхряъ мхфе цпюмхжш
   IIIUU=0
   IF(Exf.LT.EnergyRegion(NMO,IndexRegion,1)) THEN
      ! нрйкчвюел пефхл слемэьемхъ йнпмъ
      IndexSmena=0
	  ! накюярэ щмепцхх гюйнмвемю оепеундхл й якедсчыеи накюярх
      IndexRegion=IndexRegion+1
	  IIIUU=1
   ENDIF

   ! нясыеярбкъел слемэьемхе ьюцю
   IF(IndexSmena.EQ.1.AND.IKLDF.EQ.1) THEN
      Henergy=Henergy*0.1D0
   ENDIF

   ! бнгбпюыюел пефхл б хяундмне янярнъмхе
   IF(IndexSmena.EQ.0.AND.IKLDF.EQ.2) THEN
      Henergy=Henergy*10.D0
   ENDIF

   
   ! опнбепъел оепеькх й мнбни накюярх хкх мер
   IF(IIIUU.EQ.1) THEN
      ! опнбепъел врнаш вхякн хмрепбюкнб ме опебняундхкн люйяхлюкэмнцн вхякю гмювемхи 
      IF(IndexRegion.LE.NumbreMaxRegion) THEN
         ! опнбепъел оноюкх лш б мнбсч накюярэ хкх мер
	     IF(Exf.LT.EnergyRegion(NMO,IndexRegion,2).AND.Exf.GT.EnergyRegion(NMO,IndexRegion,1)) THEN
            ! щмепцхъ оноюкю б хмрепбюк
            ! нясыеярбхрэ слемэьемхе ьюцю
            IndexSmena=1
            Henergy=Henergy*0.1D0
	     ENDIF
      ENDIF
   ENDIF  
   
   return
 end subroutine ControlH

 ! ондопнцпюллю сярюмюбкхбюер йнщттхжхемр ялеьхбюмхъ
 ! N-мнлеп хрепюжхх
 ! E3-ндмнщкейрпнммюъ щмепцхъ мю N-2 хрепюжхх
 ! E2-ндмнщкейрпнммюъ щмепцхъ мю N-1 хрепюжхх
 ! E1-ндмнщкейрпнммюъ щмепцхъ мю N хрепюжхх
 real(8) function AnalysCoffMixing(N,E3,E2,E1)
   implicit none
   integer::N
   real(8)::E3,E2,E1
   !!!!!!!!!!!!!!!!!!!!!
   integer::NumbreStatus
   
   ! тхйяхпнбюммне гмювемхе 
   AnalysCoffMixing=0.8D0
    
   ! сярюмюбкхбюел пефхл янцкюянбюмхъ
   IF((E3-E2).GT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   IF((E3-E1).GT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   
   IF(NumbreStatus.EQ.1) THEN
      IF(N.LT.10) THEN
          AnalysCoffMixing=0.8D0
	     ELSE
	      AnalysCoffMixing=0.8D0 
	      IF(DABS(E3-E2).LT.0.1D0) THEN
             AnalysCoffMixing=0.7D0
	      ENDIF
          IF(DABS(E3-E2).LT.0.01D0) THEN
             AnalysCoffMixing=0.6D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.001D0) THEN
             AnalysCoffMixing=0.5D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.0001D0) THEN
             AnalysCoffMixing=0.45D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.00001D0) THEN
             AnalysCoffMixing=0.4D0
	      ENDIF
      ENDIF	  
   ENDIF

   IF(NumbreStatus.EQ.2) THEN
      AnalysCoffMixing=0.9D0
   ENDIF	   
   
   return
 end function AnalysCoffMixing

 ! ондопнцпюллю сярюмюбкхбюер йнщттхжхемр ялеьхбюмхъ (б сяхкеммнл пефхле)
 ! N-мнлеп хрепюжхх
 ! E3-ндмнщкейрпнммюъ щмепцхъ мю N-2 хрепюжхх
 ! E2-ндмнщкейрпнммюъ щмепцхъ мю N-1 хрепюжхх
 ! E1-ндмнщкейрпнммюъ щмепцхъ мю N хрепюжхх
 real(8) function AnalysCoffMixingStrong(N,E3,E2,E1)
   implicit none
   integer::N
   real(8)::E3,E2,E1
   !!!!!!!!!!!!!!!!!!!!!
   integer::NumbreStatus
   
   ! тхйяхпнбюммне гмювемхе 
   AnalysCoffMixingStrong=0.88D0
    
   ! сярюмюбкхбюел пефхл янцкюянбюмхъ
   IF((E3-E2).GT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   IF((E3-E1).GT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   
   IF(NumbreStatus.EQ.1) THEN
      IF(N.LT.10) THEN
          AnalysCoffMixingStrong=0.9D0
	     ELSE
	      AnalysCoffMixingStrong=0.9D0 
	      IF(DABS(E3-E2).LT.0.1D0) THEN
             AnalysCoffMixingStrong=0.86D0
	      ENDIF
          IF(DABS(E3-E2).LT.0.01D0) THEN
             AnalysCoffMixingStrong=0.83D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.001D0) THEN
             AnalysCoffMixingStrong=0.8D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.0001D0) THEN
             AnalysCoffMixingStrong=0.78D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.00001D0) THEN
             AnalysCoffMixingStrong=0.75D0
	      ENDIF
      ENDIF	  
   ENDIF

   IF(NumbreStatus.EQ.2) THEN
      AnalysCoffMixingStrong=0.95D0
   ENDIF	   
   
   return
 end function AnalysCoffMixingStrong 



 ! ондопнцпюллю сярюмюбкхбюер йнщттхжхемр ялеьхбюмхъ (б сяхкеммнл пефхле 2)
 ! N-мнлеп хрепюжхх
 ! E3-ндмнщкейрпнммюъ щмепцхъ мю N-2 хрепюжхх
 ! E2-ндмнщкейрпнммюъ щмепцхъ мю N-1 хрепюжхх
 ! E1-ндмнщкейрпнммюъ щмепцхъ мю N хрепюжхх
 real(8) function AnalysCoffMixingStrong2(N,E3,E2,E1)
   implicit none
   integer::N
   real(8)::E3,E2,E1
   !!!!!!!!!!!!!!!!!!!!!
   integer::NumbreStatus
   
   ! тхйяхпнбюммне гмювемхе 
   AnalysCoffMixingStrong2=0.93D0
    
   ! сярюмюбкхбюел пефхл янцкюянбюмхъ
   IF((E3-E2).GT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   IF((E3-E1).GT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   
   IF(NumbreStatus.EQ.1) THEN
      IF(N.LT.10) THEN
          AnalysCoffMixingStrong2=0.95D0
	     ELSE
	      AnalysCoffMixingStrong2=0.96D0 
	      IF(DABS(E3-E2).LT.0.1D0) THEN
             AnalysCoffMixingStrong2=0.98D0
	      ENDIF
          IF(DABS(E3-E2).LT.0.01D0) THEN
             AnalysCoffMixingStrong2=0.91D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.001D0) THEN
             AnalysCoffMixingStrong2=0.9D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.0001D0) THEN
             AnalysCoffMixingStrong2=0.89D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.00001D0) THEN
             AnalysCoffMixingStrong2=0.88D0
	      ENDIF
      ENDIF	  
   ENDIF

   IF(NumbreStatus.EQ.2) THEN
      AnalysCoffMixingStrong2=0.99D0
   ENDIF	   
   
   return
 end function AnalysCoffMixingStrong2 

 ! ондопнцпюллю сярюмюбкхбюер йнщттхжхемр ялеьхбюмхъ (б сяхкеммнл пефхле 3)
 ! N-мнлеп хрепюжхх
 ! E3-ндмнщкейрпнммюъ щмепцхъ мю N-2 хрепюжхх
 ! E2-ндмнщкейрпнммюъ щмепцхъ мю N-1 хрепюжхх
 ! E1-ндмнщкейрпнммюъ щмепцхъ мю N хрепюжхх
 real(8) function AnalysCoffMixingStrong3(N,E3,E2,E1)
   implicit none
   integer::N
   real(8)::E3,E2,E1
   !!!!!!!!!!!!!!!!!!!!!
   integer::NumbreStatus
   
   ! тхйяхпнбюммне гмювемхе 
   AnalysCoffMixingStrong3=0.97D0
    
   ! сярюмюбкхбюел пефхл янцкюянбюмхъ
   IF((E3-E2).GT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=1
   ENDIF
   IF((E3-E2).LT.0.D0.AND.(E2-E1).GT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   IF((E3-E1).GT.0.D0.AND.(E2-E1).LT.0.D0) THEN
      NumbreStatus=2
   ENDIF 
   
   IF(NumbreStatus.EQ.1) THEN
      IF(N.LT.10) THEN
          AnalysCoffMixingStrong3=0.98D0
	     ELSE
	      AnalysCoffMixingStrong3=0.975D0 
	      IF(DABS(E3-E2).LT.0.1D0) THEN
             AnalysCoffMixingStrong3=0.99D0
	      ENDIF
          IF(DABS(E3-E2).LT.0.01D0) THEN
             AnalysCoffMixingStrong3=0.96D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.001D0) THEN
             AnalysCoffMixingStrong3=0.955D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.0001D0) THEN
             AnalysCoffMixingStrong3=0.94D0
	      ENDIF  
          IF(DABS(E3-E2).LT.0.00001D0) THEN
             AnalysCoffMixingStrong3=0.92D0
	      ENDIF
      ENDIF	  
   ENDIF

   IF(NumbreStatus.EQ.2) THEN
      AnalysCoffMixingStrong3=0.99D0
   ENDIF	   
   
   return
 end function AnalysCoffMixingStrong3 

 ! ондопнцпюллю тнплхпсел люяяхбш свюбярбсчыхе б пюявере
 ! Npoint-вхякн рнвей 
 ! Isize-пюглеп люяяхбю
 ! Rcof-йнщттхжхемр пюбмши H*H/12
 ! A(Isize,Isize,Npoint)-люяяхб дкъ пюяверю оюпюлерпнб опнцнмйх
 ! B(Isize,Isize,Npoint)-люяяхб дкъ пюяверю оюпюлерпнб опнцнмйх
 ! AA(Isize,Isize,Npoint)-люяяхб дкъ пюяверю оюпюлерпнб опнцнмйх я днаюбкеммшл якюцюелшл ндмнщкейрпнммни щмепцхх 
 ! BB(Isize,Isize,Npoint)-люяяхб дкъ пюяверю оюпюлерпнб опнцнмйх я днаюбкеммшл якюцюелшл ндмнщкейрпнммни щмепцхх 
 ! RO1X(Npoint)-люяяхб оепбни опнхгбндмни мнбни оепелемни он ярюпни оепелеммни
 ! E,AZ,FZ,BZ-бяонлнцюрекэмше люяяхбш
 subroutine FormMassivAABB(Npoint,Isize,Rcof,Exf,A,B,AA,BB,RO1X,E,AZ,FZ,BZ)
   implicit none
   integer::Npoint,Isize
   real(8)::Rcof,Exf
   real(8),dimension(:)::RO1X
   real(8),dimension(:,:)::E,AZ,FZ,BZ
   real(8),dimension(:,:,:)::A,B,AA,BB
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIX,IYXDS
   
   ! тнплхпсел люяяхбш я ндмнщкейрпнммни щмепцхеи
   E=0.D0
   DO IIX=1,Npoint
      ! оепеявхршбюел гмювемхе я ндмнщкейрпнммни щмепцхеи         
  	  !DO IYXDS=1,Isize
      !    E(IYXDS,IYXDS)=Exf/RO1X(IIX)**2
      !ENDDO
	  FORALL (IYXDS=1:Isize) E(IYXDS,IYXDS)=Exf/RO1X(IIX)**2  
	    
	  ! тнплхпсел онкмши люяяхб A
      call ReadWriteMassiv2(2,IIX,Isize,A,FZ) 
      AZ=FZ-Rcof*E
	  call ReadWriteMassiv2(1,IIX,Isize,AA,AZ) 
	  ! тнплхпсел онкмши люяяхб B
	  call ReadWriteMassiv2(2,IIX,Isize,B,FZ) 
	  BZ=FZ+10.D0*Rcof*E
      call ReadWriteMassiv2(1,IIX,Isize,BB,BZ)  
   ENDDO

   return
 end subroutine FormMassivAABB

 ! ондопнцпюллю бняярюмнбкемхъ гмювемхи лнкейскъпмни нпахрюкх бн бяеу рнвйюу хг рнвйх яьхбйх
 ! NKLZ-йкчв сйюгбючыхи рхо пюяверю
 ! NKLZ=1-пюявер оюпюлерпнб дкъ яксвюъ ндмнпндмни яхярелш спюбмемхи
 ! NKLZ=2-пюявер оюпюлерпнб дкъ яксвюъ мендмнпндмни яхярелш спюбмемхи
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Npoint-вхякн рнвей
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! H-ьюц
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)-люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)-люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! FunMO(NGarmonik,Npoint+2)-люяяхб гмювемхи пюдхюкэмни вюярх тсмйжхх лнкейскъпмни нпахрюкх
 ! AX,BX,BXL,CX,VNout,UNout,VNin,UNin-бяонлнцюрекэмше люяяхбш
 ! RO1X(Npoint)-бяонлнцюрекэмши люяяхб  оепбюъ опнхгбндмюъ мнбни оепелеммни он ярюпни 
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMaxLimid-люйяхлюкэмне вхякн цюплнмхй лнкейскъпмни нпахрюкх
 ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
 ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
 ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
 subroutine FunctionRecovery(NKLZ,Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,BXL,CX,VNout,UNout,VNin,UNin,RO1X,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)
   implicit none
   integer::Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,NKLZ,NMO,NumeroMaxLimid
   real(8)::H
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::AX,BX,BXL,CX,UNout,UNin,RO1X
   real(8),dimension(:,:)::Uout,Uin,FunMO,VNout,VNin
   real(8),dimension(:,:,:)::Vout,Vin,AlfaCoffMO 
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IXXZ,IILV,IIDF,IOFG
   real(8)::RnormConstas
   
   10000 FORMAT(2X,100(1X,F20.16))

   IF(NKLZ.EQ.1) THEN
        ! онксвюел пеьемхе дкъ ндмнпндмни яхярелш  
        
		! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (1,Kpoint-1)
		DO IXXZ=Kpoint-1,1,-1 
           ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ+1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	       call ReadWriteMassiv1X(2,IXXZ+3,NumeroMin,NumeroMax,FunMO,AX)
	       !WRITE(6,*) 'POINTGGGGGG',IXXZ+1
	       !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
           ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	       call ReadWriteMassiv2(2,IXXZ,IndexM,Vout,VNout)
	       !WRITE(6,*) 'POINT',IXXZ
	       !DO IIDF=1,IndexM
           !   WRITE(6,10000) (VNout(IIDF,IOFG),IOFG=1,IndexM)
	       !ENDDO
	       !WRITE(6,*)
	       ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	       BX=MATMUL(VNout,AX)
	       ! гюохяшбюел пегскэрюр 
           call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
        ENDDO
   
        ! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (Kpoint+1,Npoint)
        DO IXXZ=Kpoint+1,Npoint 
           ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ-1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	       call ReadWriteMassiv1X(2,IXXZ+1,NumeroMin,NumeroMax,FunMO,AX)
	       !WRITE(6,*) 'POINTGGGGGG',IXXZ-1
	       !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
           ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	       call ReadWriteMassiv2(2,IXXZ,IndexM,Vin,VNin)
	       ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	       BX=MATMUL(VNin,AX)
	       ! гюохяшбюел пегскэрюр 
           call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
        ENDDO
   	  
   ENDIF	  

   IF(NKLZ.EQ.2) THEN

        ! онксвюел пеьемхе дкъ мендмнпндмни яхярелш
		! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (1,Kpoint-1)
		DO IXXZ=Kpoint-1,1,-1 
           ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ+1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	       call ReadWriteMassiv1X(2,IXXZ+3,NumeroMin,NumeroMax,FunMO,AX)
	       !WRITE(6,*) 'POINTGGGGGG',IXXZ+1
	       !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
           ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	       call ReadWriteMassiv2(2,IXXZ,IndexM,Vout,VNout)
	       call ReadWriteMassiv1(2,IXXZ,IndexM,Uout,UNout)  
           !WRITE(6,*) 'POINT',IXXZ
	       !DO IIDF=1,IndexM
           !   WRITE(6,10000) (VNout(IIDF,IOFG),IOFG=1,IndexM)
	       !ENDDO
	       !WRITE(6,*)
	       ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	       BX=UNout+MATMUL(VNout,AX)
	       ! гюохяшбюел пегскэрюр 
           call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
        ENDDO
   
        ! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (Kpoint+1,Npoint)
        DO IXXZ=Kpoint+1,Npoint 
           ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ-1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	       call ReadWriteMassiv1X(2,IXXZ+1,NumeroMin,NumeroMax,FunMO,AX)
	       !WRITE(6,*) 'POINTGGGGGG',IXXZ-1
	       !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
           ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	       call ReadWriteMassiv2(2,IXXZ,IndexM,Vin,VNin)
	       call ReadWriteMassiv1(2,IXXZ,IndexM,Uin,UNin)  
           ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	       BX=UNin+MATMUL(VNin,AX)
	       ! гюохяшбюел пегскэрюр 
           call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
        ENDDO
   ENDIF

   ! тнплхпсел ярпсйрспс бшяьху цюплнмхй лнкейскъпмни нпахрюкх
		
   ! нопедекъел йнщттхжхемрш юкэтю
   call CalculationCoffAlfaMO(NMO,Npoint,H,NumeroMax,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,FunMO,RfunLigands,AlfaCoffMO,RO1X) 
        
   ! тнплхпсел ярпсйрспс бшяьху цюплнмхй лнкейскъпмшу нпахрюкеи я свернл онксвеммшу йнщттхжхемрнб юкэтю
   call FormStructureMO(NMO,Npoint,H,NumeroMax,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,FunMO,RfunLigands,AlfaCoffMO) 



   
   ! мнплхпсел онксвеммсч лнкейскъпмсч нпахрюкэ
   RnormConstas=CalculationNormalizationConsta(NumeroMin,NumeroMaxLimid,Npoint,H,RO1X,FunMO) 

   !WRITE(*,*) 'NORMA',RnormConstas
   !READ(*,*)

   ! гюохяшбюел онксвеммсч мнплхпнбнвмсч йнмярюмрс 
   BXL=RnormConstas
   call ReadWriteMassiv1X(1,2,NumeroMin,NumeroMaxLimid,FunMO,BXL)

   !WRITE(*,*) 'NORMA',BX
   !READ(*,*)

   
   !DO IXXZ=1,Npoint
   !   WRITE(6,10000) (FunMO(IILV,IXXZ+2)/DSQRT(RnormConstas*RO1X(IXXZ)),IILV=NumeroMin,NumeroMax)
   !ENDDO

   
   return
 end subroutine FunctionRecovery


 ! ондпнцпюллю пюяверю мнплхпнбнвмнцн йнщттхжхемрю онксвеммни тсмйжхх
 ! NMO- мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! Npoint-вхякн рнвей
 ! H-ьюц
 ! RO1X(Npoint)-бяонлнцюрекэмши люяяхб  оепбюъ опнхгбндмюъ мнбни оепелеммни он ярюпни 
 ! FunMO(NMO,NGarmonik,Npoint+2)-люяяхб гмювемхи пюдхюкэмни вюярх тсмйжхх лнкейскъпмни нпахрюкх
 real(8) function CalculationNormalizationConstaResult(NMO,NumeroMin,NumeroMax,Npoint,H,RO1X,FunMO) 
   implicit none
   integer::NumeroMin,NumeroMax,Npoint,NMO
   real(8)::H
   real(8),dimension(:)::RO1X
   real(8),dimension(:,:,:)::FunMO
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::INNC,IDFZ
   real(8)::Snorm,SING,SING2,SING1,R0,RN1,RN

  
   Snorm=0.D0
   DO INNC=NumeroMin,NumeroMax
      SING=0.D0 
	  !DO IDFZ=1,Npoint2
      !   SING=SING+(FunMO(NMO,INNC,IDFZ+2)/RO1X(IDFZ))**2
	  !ENDDO
      !Snorm=Snorm+sum((FunMO(NMO,INNC,2+1:Npoint)/RO1X(1:Npoint))**2)
      
	  ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
      IF(Npoint/2..EQ.Npoint/2) THEN
	        ! пюявер он тнплске яхлоянмю
            ! вермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,Npoint-2,2
               SING2=SING2+(FunMO(NMO,INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,Npoint-3,2
               SING1=SING1+(FunMO(NMO,INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
		    ! пегскэрхпсчыюъ ясллю 
		    R0=(FunMO(NMO,INNC,3)/RO1X(1))**2
		    RN1=(FunMO(NMO,INNC,Npoint+1)/RO1X(Npoint-1))**2
		    RN=(FunMO(NMO,INNC,Npoint+2)/RO1X(Npoint))**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
         ELSE
            ! пювер он тнплске яхлоянмю
			! мевермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,Npoint-1,2
               SING2=SING2+(FunMO(NMO,INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,Npoint-2,2
               SING1=SING1+(FunMO(NMO,INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            ! пегскэрхпсчыюъ ясллю 
		    R0=(FunMO(NMO,INNC,3)/RO1X(1))**2
		    RN=(FunMO(NMO,INNC,Npoint+2)/RO1X(Npoint))**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
      ENDIF
      ! ясллхпсел хмрецпюкш
	  Snorm=Snorm+SING
   ENDDO

   CalculationNormalizationConstaResult=Snorm*H
  

   return
 end function CalculationNormalizationConstaResult




 
 
 ! ондпнцпюллю пюяверю мнплхпнбнвмнцн йнщттхжхемрю
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! Npoint-вхякн рнвей
 ! H-ьюц
 ! RO1X(Npoint)-бяонлнцюрекэмши люяяхб  оепбюъ опнхгбндмюъ мнбни оепелеммни он ярюпни 
 ! FunMO(NGarmonik,Npoint+2)-люяяхб гмювемхи пюдхюкэмни вюярх тсмйжхх лнкейскъпмни нпахрюкх
 real(8) function CalculationNormalizationConsta(NumeroMin,NumeroMax,Npoint,H,RO1X,FunMO) 
   implicit none
   integer::NumeroMin,NumeroMax,Npoint
   real(8)::H
   real(8),dimension(:)::RO1X
   real(8),dimension(:,:)::FunMO
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::INNC,IDFZ
   real(8)::Snorm,SING,SING2,SING1,R0,RN1,RN

   Snorm=0.D0
   DO INNC=NumeroMin,NumeroMax
      SING=0.D0 
	  !DO IDFZ=1,Npoint
      !   SING=SING+(FunMO(INNC,IDFZ+2)/RO1X(IDFZ))**2
	  !ENDDO
	  !Snorm=Snorm+sum((FunMO(INNC,2+1:Npoint)/RO1X(1:Npoint))**2)
      ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
      IF(Npoint/2..EQ.Npoint/2) THEN
	        ! пюявер он тнплске яхлоянмю
            ! вермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,Npoint-2,2
               SING2=SING2+(FunMO(INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,Npoint-3,2
               SING1=SING1+(FunMO(INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
		    ! пегскэрхпсчыюъ ясллю 
		    R0=(FunMO(INNC,3)/RO1X(1))**2
		    RN1=(FunMO(INNC,Npoint+1)/RO1X(Npoint-1))**2
		    RN=(FunMO(INNC,Npoint+2)/RO1X(Npoint))**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
         ELSE
            ! пювер он тнплске яхлоянмю
			! мевермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,Npoint-1,2
               SING2=SING2+(FunMO(INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,Npoint-2,2
               SING1=SING1+(FunMO(INNC,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            ! пегскэрхпсчыюъ ясллю 
		    R0=(FunMO(INNC,3)/RO1X(1))**2
		    RN=(FunMO(INNC,Npoint+2)/RO1X(Npoint))**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
      ENDIF
      ! ясллхпсел хмрецпюкш
	  Snorm=Snorm+SING
   ENDDO

   CalculationNormalizationConsta=Snorm*H

   return
 end function CalculationNormalizationConsta




  ! ондопнцпюллю нопедекъчыюъ йнщттхжхемр юкэтю бундъыхи б ярпсйрспс лнкейскъпмни нпахрюкх
  ! NMO-мнлеп лнкейскъпмни нпахрюкх
  ! Npoint-вхякн рнвей
  ! H-ьюц 
  ! NumbreGarmMO-вхякн цюплнмхй лнкейскъпмни нпахрюкх дн Lk бйкчвхрекэмн
  ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
  ! RFunMO(NumbreGarmMO,Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
  ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
  ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
  ! RO1X(Npoint)-бяонлнцюрекэмши люяяхб  оепбюъ опнхгбндмюъ мнбни оепелеммни он ярюпни 
  subroutine CalculationCoffAlfaMO(NMO,Npoint,H,NumbreGarmMO,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFunMO,RfunLigands,AlfaCoffMO,RO1X) 
    implicit none
    integer::NMO,Npoint,NumbreGarmMO
	real(8)::H
    integer,dimension(:)::NumbreLigand
    integer,dimension(:,:)::NLigands,NFunLigands
    integer,dimension(:,:,:)::NumbreFunctionLig
	real(8),dimension(:)::RO1X 
    real(8),dimension(:,:)::RFunMO
	real(8),dimension(:,:,:)::AlfaCoffMO
    real(8),dimension(:,:,:,:)::RfunLigands
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::NumbreFunMOZ,IIXX,IIXY,NumbreFunLigZX,NumbreLigandVV,NumbreGMOAS
	integer::NumbreLigandVV1,NumbreLigandVV2,NumbreLigandVV3,NumbreFunLigZX1,NumbreFunLigZX2,NumbreFunLigZX3,NumbreGMOAS1,NumbreGMOAS2,NumbreGMOAS3
    integer::ISUMINDEX,IXZZ1,IXZY1
	real(8)::Snorm,RintMO,RintFL,RintMO1,RintMO2,RintMO3,RintFL11,RintFL12,RintFL13,RintFL21,RintFL22,RintFL23,RintFL31,RintFL32,RintFL33
    real(8)::Rdelta,Rdelta1,Rdelta2,Rdelta3

	! нопедекъел вхякн тсмйжхи кхцюмднб нохяшбючыхе бшяьхе цюплнмхйх лнкейскъпмни нпахрюкх
    NumbreFunMOZ=0
	DO IIXX=1,NumbreLigand(NMO)
       NumbreFunMOZ=NumbreFunMOZ+NFunLigands(NMO,IIXX)
	ENDDO

	! опнбепъел вхякн сярюмнбкеммшу тсмйжхи
	IF(NumbreFunMOZ.GT.3) THEN
       WRITE(6,*) 'ERROR. Numbre function of ligands >3 NumbreFunMOZ= ',NumbreFunMOZ 
	   WRITE(7,*) 'ERROR. Numbre function of ligands >3 NumbreFunMOZ= ',NumbreFunMOZ  
	   WRITE(*,*) 'ERROR. Numbre function of ligands >3 NumbreFunMOZ= ',NumbreFunMOZ 
	   READ(*,*)
	   STOP
	ENDIF
    
	! яксвюи ндмни тсмйжхх
	IF(NumbreFunMOZ.EQ.1) THEN
       ! мнлеп цюплнмхй он йнрнпни нясыеярбъел нопедекемхе йнщттхжхемрнб
	   NumbreGMOAS=NumbreGarmMO
	   ! мнлеп кхцюмдю
	   NumbreLigandVV=NLigands(NMO,1)
       ! мнлеп тсмйжхх кхцюмдю
	   NumbreFunLigZX=NumbreFunctionLig(NMO,1,1)
	   RintMO=sum((RFunMO(NumbreGMOAS,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS,2) 
	   RintFL=sum(RFunMO(NumbreGMOAS,2+1:Npoint)*RfunLigands(NumbreLigandVV,NumbreFunLigZX,NumbreGMOAS,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS,2)*RfunLigands(NumbreLigandVV,NumbreFunLigZX,NumbreGMOAS,2)) 
	   ! пюявхршбюел йнщттхжхемр юкэтю
	   AlfaCoffMO(NMO,NumbreLigandVV,NumbreFunLigZX)=RintMO/RintFL
	ENDIF
    
	! яксвюи дбсу тсмйжхи
    IF(NumbreFunMOZ.EQ.2) THEN
       ! мнлеп цюплнмхй он йнрнпни нясыеярбъел нопедекемхе йнщттхжхемрнб
	   NumbreGMOAS1=NumbreGarmMO-1
	   NumbreGMOAS2=NumbreGarmMO
       ! опнбепъел вхякн кхцюмднб
	   IF(NumbreLigand(NMO).EQ.2) THEN
            ! мнлеп кхцюмднб
			NumbreLigandVV1=NLigands(NMO,1)
            NumbreLigandVV2=NLigands(NMO,2)
            ! мнлеп тсмйжхи
            NumbreFunLigZX1=NumbreFunctionLig(NMO,1,1)
            NumbreFunLigZX2=NumbreFunctionLig(NMO,2,1)
          ELSE
		    ! мнлеп кхцюмднб
			NumbreLigandVV1=NLigands(NMO,1)
            NumbreLigandVV2=NLigands(NMO,1)
            ! мнлеп тсмйжхи
			NumbreFunLigZX1=NumbreFunctionLig(NMO,1,1)
            NumbreFunLigZX2=NumbreFunctionLig(NMO,1,2)
       ENDIF
       RintMO1=sum((RFunMO(NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS1,2)        
       RintMO2=sum((RFunMO(NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS2,2)
       RintFL11=sum(RFunMO(NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS1,2)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2)) 
	   RintFL12=sum(RFunMO(NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS1,2)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2))
	   RintFL21=sum(RFunMO(NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS2,2)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2)) 
	   RintFL22=sum(RFunMO(NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS2,2)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2))
	   ! пюявхршбюел йнщттхжхемрш юкэтю
	   AlfaCoffMO(NMO,NumbreLigandVV2,NumbreFunLigZX2)=(RintFL21*RintMO1-RintFL11*RintMO2)/(RintFL12*RintFL21-RintFL11*RintFL22)  
       AlfaCoffMO(NMO,NumbreLigandVV1,NumbreFunLigZX1)=(RintMO1-AlfaCoffMO(NMO,NumbreLigandVV2,NumbreFunLigZX2)*RintFL12)/RintFL11
	  
	ENDIF

    ! яксвюи рпеу тсмйжхи
    IF(NumbreFunMOZ.EQ.3) THEN
	   ! мнлепю цюплнмхй он йнрнпни нясыеярбъел нопедекемхе йнщттхжхемрнб
	   NumbreGMOAS1=NumbreGarmMO-2
	   NumbreGMOAS2=NumbreGarmMO-1
	   NumbreGMOAS3=NumbreGarmMO
       
	   ISUMINDEX=0
	   ! жхйк он кхцюмдюл
	   DO IXZZ1=1,NumbreLigand(NMO)
	      ! жхйк он тсмйжхъл кхцюмдю
		  DO IXZY1=1,NFunLigands(NMO,IXZZ1)
		      ! мнбне якюцюелне
			  ISUMINDEX=ISUMINDEX+1
			  IF(ISUMINDEX.EQ.1) THEN
                 ! мнлеп кхцюмдю
				 NumbreLigandVV1=NLigands(NMO,IXZZ1) 
                 ! мнлеп тсмйжхх кхцюмдю
				 NumbreFunLigZX1=NumbreFunctionLig(NMO,IXZZ1,IXZY1)
			  ENDIF		 
              IF(ISUMINDEX.EQ.2) THEN
                 ! мнлеп кхцюмдю
				 NumbreLigandVV2=NLigands(NMO,IXZZ1) 
                 ! мнлеп тсмйжхх кхцюмдю
				 NumbreFunLigZX2=NumbreFunctionLig(NMO,IXZZ1,IXZY1)
			  ENDIF		 
              IF(ISUMINDEX.EQ.3) THEN
                 ! мнлеп кхцюмдю
				 NumbreLigandVV3=NLigands(NMO,IXZZ1) 
                 ! мнлеп тсмйжхх кхцюмдю
				 NumbreFunLigZX3=NumbreFunctionLig(NMO,IXZZ1,IXZY1)
			  ENDIF		 
		  ENDDO
       ENDDO
       ! ПЮЯВЕР ОПНБЕДЕМ ОН ЛЕРНДС ЙПЮЛЕПЯЮ (ПЕЬЕМХЪ ЯХЯРЕЛШ КХМЕИМШУ СПЮБМЕМХИ)
	   RintMO1=sum((RFunMO(NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintMO2=sum((RFunMO(NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintMO3=sum((RFunMO(NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintFL11=sum(RFunMO(NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H  
	   RintFL12=sum(RFunMO(NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL13=sum(RFunMO(NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL21=sum(RFunMO(NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL22=sum(RFunMO(NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL23=sum(RFunMO(NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL31=sum(RFunMO(NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL32=sum(RFunMO(NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL33=sum(RFunMO(NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   Rdelta=RintFL11*(RintFL22*RintFL33-RintFL23*RintFL32)-RintFL12*(RintFL21*RintFL33-RintFL23*RintFL31)+RintFL13*(RintFL21*RintFL32-RintFL22*RintFL31)
	   Rdelta1=RintMO1*(RintFL22*RintFL33-RintFL23*RintFL32)-RintFL12*(RintMO2*RintFL33-RintFL23*RintMO3)+RintFL13*(RintMO2*RintFL32-RintFL22*RintMO3)
	   Rdelta2=RintFL11*(RintMO2*RintFL33-RintFL23*RintMO3)-RintMO1*(RintFL21*RintFL33-RintFL23*RintFL31)+RintFL13*(RintFL21*RintMO3-RintMO2*RintFL31)
	   Rdelta3=RintFL11*(RintFL22*RintMO3-RintMO2*RintFL32)-RintFL12*(RintFL21*RintMO3-RintMO2*RintFL31)+RintMO1*(RintFL21*RintFL32-RintFL22*RintFL31)
	   ! пюявхршбюел йнщттхжхемрш юкэтю
	   AlfaCoffMO(NMO,NumbreLigandVV1,NumbreFunLigZX1)=Rdelta1/Rdelta
	   AlfaCoffMO(NMO,NumbreLigandVV2,NumbreFunLigZX2)=Rdelta2/Rdelta
       AlfaCoffMO(NMO,NumbreLigandVV3,NumbreFunLigZX3)=Rdelta3/Rdelta
	ENDIF 



    return
  end subroutine CalculationCoffAlfaMO


  ! ондопнцпюллю нясыеярбкъер тнплхпнбюмхе ярпсйрспш бяьху цюплнмхй лнкейскъпмни нпахрюкх
  ! NMO-мнлеп лнкейскъпмни нпахрюкх
  ! Npoint-вхякн рнвей
  ! H-ьюц 
  ! NumbreGarmMO-вхякн цюплнмхй лнкейскъпмни нпахрюкх дн Lk бйкчвхрекэмн
  ! NumbreGarmMOLimid-вхякн цюплнмхй лнкейскъпмни нпахрюкх оняке Lk дн люйяхлюкэмни цюплнмхйх 
  ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
  ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
  ! RFunMO(NumbreGarmMO,Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
  ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
  ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
  subroutine FormStructureMO(NMO,Npoint,H,NumbreGarmMO,NumbreGarmMOLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFunMO,RfunLigands,AlfaCoffMO) 
    implicit none
    integer::NMO,Npoint,NumbreGarmMO,NumbreGarmMOLimid
	real(8)::H
    integer,dimension(:)::NumbreLigand
    integer,dimension(:,:)::NLigands,NFunLigands
    integer,dimension(:,:,:)::NumbreFunctionLig
	real(8),dimension(:,:)::RFunMO
	real(8),dimension(:,:,:)::AlfaCoffMO
    real(8),dimension(:,:,:,:)::RfunLigands
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIXDGAR,IX1,IX2,IX3,NumbreLigandX,NumbreFunLigX,ierr
	real(8)::Rcoff
    real(8),allocatable,dimension(:)::RfunXD
    
	! бшдекъел оюлърэ онд люяяхб
	allocate(RfunXD(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*)'FormStructureMO'
      write(*,*)'MEMORY ON THE FILE "RfunXD" IS NOT SELECTED'
      read(*,*)
	  stop 
    endif


    
	! жхйк он цюплнмхйюл
	DO IIXDGAR=NumbreGarmMO+1,NumbreGarmMOLimid
       RfunXD=0.D0
	   ! жхйк он вхякс кхцюмднб
	   DO IX1=1,NumbreLigand(NMO)
	      ! мнлеп кхцюмдю 
	      NumbreLigandX=NLigands(NMO,IX1)
          ! жхйк он тсмйжхъл
		  DO IX2=1,NFunLigands(NMO,IX1)
             ! мнлеп тсмйжхх
			 NumbreFunLigX=NumbreFunctionLig(NMO,IX1,IX2)
             Rcoff=AlfaCoffMO(NMO,NumbreLigandX,NumbreFunLigX)
			 DO IX3=1,Npoint
			    RfunXD(IX3)=RfunXD(IX3)+Rcoff*RfunLigands(NumbreLigandX,NumbreFunLigX,IIXDGAR,2+IX3)
	         ENDDO
		  ENDDO
	   ENDDO
       ! гюохяшбюел ятнплхпнбюммсч цюплнмхйс
       RFunMO(IIXDGAR,2+1:Npoint)=RfunXD(1:Npoint)
	ENDDO

	
	! сдюкъел люяяхбш хг оюлърх 
    deallocate(RfunXD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'FormStructureMO'
       write(*,*) 'THE MASSIV "RfunXD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 

    return
  end subroutine FormStructureMO


 ! ондопнцпюллю пюяверю гмювемхи бнкмнбшу тсмйжхи б рнвйе яьхбйх (яксвюи хрепюжхх оняке оепбни)
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! FunMO(IndexM,Npoint+2)-пюдхюкэмюъ вюярэ бнкмнбни тсмйжхх лнкейскъпмни нпахрюкх 
 ! VNout,UNout,VNin,UNin,E,AX,BX,IXM-бяонлнцюрекэмше люяяхбш
 subroutine CalculationValuesGarmonikMOFull(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,FunMO,E,AX,BX,IXM)
   implicit none
   integer::Kpoint,IndexM,NumeroMin,NumeroMax
   integer,dimension(:,:)::IXM
   real(8),dimension(:)::UNout,UNin,AX,BX
   real(8),dimension(:,:)::Uout,Uin,VNout,VNin,E,FunMO
   real(8),dimension(:,:,:)::Vout,Vin
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIXD,IIOF,IIDD

   69685 FORMAT(2X,100(F20.16,1X))

   ! явхршбюел оюпюлерпш опнцнмйх дкъ тнплхпнбюмхъ яхярелш спюбмемхи
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,VNout)
   call ReadWriteMassiv1(2,Kpoint,IndexM,Uout,UNout)   
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,VNin) 
   call ReadWriteMassiv1(2,Kpoint+1,IndexM,Uin,UNin)  

   
   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   !DO IIDD=1,IndexM
   !   E(IIDD,IIDD)=1.D0
   !ENDDO
   FORALL (IIDD=1:IndexM) E(IIDD,IIDD)=1.D0 
   
   ! тнплхпсел люрпхжс яхярелш  
   E=E-MATMUL(VNout,VNin)
   
   ! тнплхпсел ярнкаеж ябнандмшу вкемнб
   AX=MATMUL(VNout,UNin)+UNout
  
   !WRITE(6,*) 'E'
   !DO IIXD=1,IndexM
   !   WRITE(6,69685) (E(IIXD,IIOF),IIOF=1,IndexM)
   !ENDDO
   
   !WRITE(6,*) 'AX'
   !WRITE(6,69685) (AX(IIOF),IIOF=1,IndexM) 



   ! пеьюел онксвеммсч яхярелс спюбмемхи
   call SolutionSistemLinearEquationFull(E,AX,IXM,BX) 
   
 

   !WRITE(6,*) 'POINT SHIVKA',Kpoint
   !WRITE(6,69685) (BX(IIOF),IIOF=1,IndexM)
   
   ! гюохяшбюел онксвеммше гмювемхъ лнкейскъпмни нпахрюкх б рнвйе яьхбйх  Kpoint+2=>мювюкн люяяхбю б рнвйе 3 
   call ReadWriteMassiv1X(1,Kpoint+2,NumeroMin,NumeroMax,FunMO,BX)
    
   return
 end subroutine CalculationValuesGarmonikMOFull





 ! ондопнцпюллю пюяверю гмювемхи бнкмнбшу тсмйжхи б рнвйе яьхбйх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! FunMO(IndexM,Npoint+2)-пюдхюкэмюъ вюярэ бнкмнбни тсмйжхх лнкейскъпмни нпахрюкх 
 ! AZ,BZ,AX,BX-бяонлнцюрекэмше люяяхбш
 subroutine CalculationValuesGarmonikMO(Kpoint,IndexM,NumeroMin,NumeroMax,Vout,Vin,FunMO,AZ,BZ,E,AX,BX)
   implicit none
   integer::Kpoint,IndexM,NumeroMin,NumeroMax
   real(8),dimension(:)::AX,BX
   real(8),dimension(:,:)::FunMO,AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIXD,IIOF

   69685 FORMAT(2X,100(F20.16,1X))

   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ)
   
   !WRITE(6,*) 'AZ'
   !DO IIXD=1,IndexM
   !  WRITE(6,69685) (AZ(IIXD,IIOF),IIOF=1,IndexM)
   !ENDDO
   !WRITE(6,*) 'BZ'
   !DO IIXD=1,IndexM
   !   WRITE(6,69685) (BZ(IIXD,IIOF),IIOF=1,IndexM)
   !ENDDO
  
   E=0.D0
   !DO IIXD=1,IndexM 
   !   E(IIXD,IIXD)=1.D0
   !ENDDO
   FORALL (IIXD=1:IndexM) E(IIXD,IIXD)=1.D0

   ! люрпхжю яхярелш спюбмемхи 
   E=E-MATMUL(AZ,BZ)
   !WRITE(6,*) 'E'
   !DO IIXD=1,IndexM
   !   WRITE(6,69685) (E(IIXD,IIOF),IIOF=1,IndexM)
   !ENDDO

   
   !READ(*,*)
   ! оепеярпюхбюел люрпхжс пюяонкнцюъ ярпнвйх б напюрмнл онпъдйе дкъ рнцн, врнаш опхпюбмърэ й едхмхже гмювемхе лхмхлюкэмни цюплнмхйх
   ! щрн ябъгюмн я реумхйни онксвемхъ пеьемхъ яхярелш спюбмемхи (напюрмши унд) хяонкэгсчр онякедмхи йнпемэ нм днкфем ашрэ хгбеярем
   ! оепеярюбкъел б напюрмнл онпъдйе ярнкажш 
   call InterchangeStolbsovMatrix(E)
   ! оепеярюбкъел ярпнвйх 
   call InterchangeStrokMatrix(E)
   ! опнбепъел дхюцнмюкэмше щкелемрш люрпхжш - нмх днкфмш сашбюрэ я сбекхвемхе онпъдйнбнцн мнлепю
   !call ControlMatrix(E)


   !WRITE(6,*) 'DAA-1'
   !DO IIXD=1,IndexM
   !   WRITE(6,69685) (E(IIXD,IIOF),IIOF=1,IndexM)
   !ENDDO
   

   ! мюундхл пеьемхе яхярелш спюбмемхи (онксвюел гмювемхе лнкейскъпмни нпахрюкх б рнвйе яьхбйх)
   ! нясыеярбкъеряъ оепеунд дкъ онксвеммнцн бейрнпю й нашвмнлс онпъдйс ярпнй (йнпмеи)
   call SolutionSistemLinearEquation(E,AX,BX)

   !WRITE(6,*) 'POINT SHIVKA',Kpoint
   !WRITE(6,69685) (AX(IIOF),IIOF=1,IndexM)
   !WRITE(6,*) 'Kpoint',Kpoint
  
   ! гюохяшбюел онксвеммше гмювемхъ ! Kpoint+2=>мювюкн люяяхбю б рнвйе 3 
   call ReadWriteMassiv1X(1,Kpoint+2,NumeroMin,NumeroMax,FunMO,AX)
   
   return
 end subroutine CalculationValuesGarmonikMO



 ! опнбепйю дхюцнмюкэмшу люрпхвмшу щкелемрнб (нмх днкфмш сашбюрэ я пнярнл мнлепю хмдейяю)
 subroutine  ControlMatrix(E)
   implicit none 
   real(8),dimension(:,:)::E
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NXX,IIXX

   NXX=size(E,dim=1)

   DO IIXX=1,NXX-1
      IF(DABS(E(IIXX,IIXX)).LT.DABS(E(IIXX+1,IIXX+1))) THEN
	     WRITE(6,*) 'First iteration'
		 WRITE(6,*) 'Irregular order matrix elements'
		 EXIT
	  ENDIF 
   ENDDO 

   return
 end subroutine ControlMatrix

 ! ондопнцпюллю пеьемхъ яхярелш спюбмемхи 
 ! E-йбюдпюрмюъ люрпхжю
 ! AX,BX-ндмнлепмши люяяхб (бейрнп пеьемхи)
 subroutine SolutionSistemLinearEquation(E,AX,P)
   implicit none 
   real(8),dimension(:)::AX,P
   real(8),dimension(:,:)::E
   !!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NXX,ierr,IIXX,IIYY,IKL
   real(8)::S
  
  
   NXX=size(E,dim=1)
   
   
   ! бшонкмъел опълни унд
   DO IIXX=1,NXX-1 
      P(IIXX+1:NXX)=E(IIXX+1:NXX,IIXX)/E(IIXX,IIXX)
	  DO IIYY=IIXX+1,NXX
	     DO IKL=1,NXX
            E(IIYY,IKL)=E(IIYY,IKL)-E(IIXX,IKL)*P(IIYY)
		 ENDDO
	  ENDDO
   ENDDO
   
   
       
    
 

   ! бшонкмъел напюрмши унд 
   AX(NXX)=1.D0 ! тхйяхпсел пеьемхе яхярелш опхпюбмхбюъ лхмхлюкэмни цюплнмхйе б рнвйе яьхйх едемхжс
   DO IIXX=NXX-1,1,-1
      S=sum(E(IIXX,IIXX+1:NXX)*AX(IIXX+1:NXX))
	  AX(IIXX)=-S/E(IIXX,IIXX)   
   ENDDO 

   ! нясыеярбкъел бняCрюмнбкемхе опюбхкэмни онякеднбюрекэмнярх
   IIYY=NXX
   DO IIXX=1,NXX 
      P(IIYY)=AX(IIXX)
      IIYY=IIYY-1
   ENDDO
   
   ! гюохяшбюел пегскэрюр
   AX=P

   
  
   
   return
 end subroutine SolutionSistemLinearEquation




 ! ондопнцпюллю пеьемхъ яхярелш кхмеимшу спюбмемхи Ax=B лерндюл цюсяяю я бшанпнл бедсыецн щкелемрю он бяеи люрпхже 
 ! A(N,N)-йбюдпюрмюъ люрпхжю яхярелш
 ! B(N)-ярнкаеж ябнандмшу вкемнб
 ! IS(N,2)-люрпхжю тюйрнпхгюжхх
 ! X(N)-онксвеммне пеьемхе яхярелш 
 subroutine SolutionSistemLinearEquationFull(A,B,IS,X)
   implicit none
   integer,dimension(:,:)::IS
   real(8),dimension(:)::B,X
   real(8),dimensioN(:,:)::A
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::S,N,P,F,I,J,K,L,M,BP,BI,I1,BH
   real(8)::BO,BK,BX,CC,CD,PR,Z,W,R  

         
   ! пюглепмнярэ яхярелш
   N=size(A,dim=1)

   IF(N.EQ.1) THEN 
      ! опнбепъел нркхвхе нр мскъ
      IF(A(1,1).NE.0.D0) THEN
         X(1)=B(1)/A(1,1) 
      ENDIF
	  
	  return
   ENDIF   

	
   ! тнйрюпхгюжхъ люрпхжш яхярелш
   F=0
   1  F=F+1
      BK=0.0D0
      DO J=F,N
         BP=F
         BO=0.D0
         DO I=F,N
            BX=DABS(A(I,J))
            IF(BO.GE.BX) THEN 
  	           CYCLE 
            ENDIF
		    BO=BX
            BP=I
         ENDDO 
      
	     IF(BK.GT.BO) THEN 
	        CYCLE 
         ENDIF   
	     BK=BO
         BH=BP
         BI=J
      ENDDO 
      
      DO I=1,N,1
         CC=A(I,BI)
         A(I,BI)=A(I,F)
         A(I,F)=CC
      ENDDO 
   
      DO J=1,N,1
         CD=A(BH,J)
         A(BH,J)=A(F,J)
         A(F,J)=CD
      ENDDO

      IS(F,1)=BH
      IS(F,2)=BI
   
      PR=1.D0/A(F,F)
      L=F+1
      DO I=L,N,1
         A(I,F)=-(PR*A(I,F))
      ENDDO
     
      DO J=L,N,1
         DO I=L,N,1
            A(I,J)=A(I,J)+A(F,J)*A(I,F)
         ENDDO
      ENDDO
   IF(L-N.LT.0) THEN 
      GOTO 1
   ENDIF 

   DO I=1,N
      X(I)=B(I)
   ENDDO
   
   M=N-1
   DO I=1,M
      I1=IS(I,1)
      Z=X(I1)
      X(I1)=X(I)
      X(I)=Z
   ENDDO
   
   DO J=1,M
      W=X(J)
      K=J+1
      DO I=K,N
         X(I)=X(I)+A(I,J)*W
      ENDDO
   ENDDO
   
   X(N)=X(N)/A(N,N)
   K=N
   DO I=1,M
      R=0.0D0
      L=K-1
      DO J=K,N
         R=R+A(L,J)*X(J)
      ENDDO
      K=N-I
      X(K)=(X(K)-R)/A(K,K)
   ENDDO
   
   DO J=1,M
      I=N-J
      I1=IS(I,2)
      Z=X(I)
      X(I)=X(I1)
      X(I1)=Z
   ENDDO 
   
   return 
 end subroutine SolutionSistemLinearEquationFull







 
 ! ондопнцпюллю оепеярюмнбйх ярнкажнб люрпхжш б напюрмнл онпъдйе
 ! A-йбюдпюрмюъ люрпхжю
 subroutine InterchangeStolbsovMatrix(A)
   implicit none
   real(8),dimension(:,:)::A
   !!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NNX,III1,III2,III3,ierr
   real(8),allocatable,dimension(:,:)::VC
   
   NNX=size(A,dim=1)
   
   ! бшдекъел оълърэ онд люяяхбш
   allocate(VC(NNX,NNX),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'InterchangeStrokMatrix'
      write(*,*) 'MEMORY ON THE FILE "VC" IS NOT SELECTED'
  	  read(*,*)
	  stop 
   endif  

   
   III3=NNX
   DO III1=1,NNX
      DO III2=1,NNX
         VC(III2,III3)=A(III2,III1) 
	  ENDDO
	  III3=III3-1
   ENDDO
   
   ! гюохяшбюел онксвеммсч люрпхжс (я напюрмшл онпъдйнл ярнкажнб)
   A=VC
   
   ! сдюкемхе люяяхбнб хг оълърх
   deallocate(VC,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'InterchangeStrokMatrix'
      write(*,*) 'THE FILE "VC" IS NOT REMOVED FROM MEMORY'
	  read(*,*) 
	  stop 
   endif

   return
 end subroutine InterchangeStolbsovMatrix


 ! ондопнцпюллю оепеярюмнбйх ярпнй люрпхжш б напюрмнл онпъдйе
 ! A-йбюдпюрмюъ люрпхжю
 subroutine InterchangeStrokMatrix(A)
   implicit none
   real(8),dimension(:,:)::A
   !!!!!!!!!!!!!!!!!!!!!!!!!
   integer::NNX,III1,III2,III3,ierr
   real(8),allocatable,dimension(:,:)::VC
   
   NNX=size(A,dim=1)
   
   ! бшдекъел оълърэ онд люяяхбш
   allocate(VC(NNX,NNX),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'InterchangeStrokMatrix'
      write(*,*) 'MEMORY ON THE FILE "VC" IS NOT SELECTED'
  	  read(*,*)
	  stop 
   endif  

   
   III3=NNX
   DO III1=1,NNX
      DO III2=1,NNX
         VC(III3,III2)=A(III1,III2) 
	  ENDDO
	  III3=III3-1
   ENDDO
   
   ! гюохяшбюел онксвеммсч люрпхжс (я напюрмшл онпъдйнл ярпнй)
   A=VC
   
   ! сдюкемхе люяяхбнб хг оълърх
   deallocate(VC,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'InterchangeStrokMatrix'
      write(*,*) 'THE FILE "VC" IS NOT REMOVED FROM MEMORY'
	  read(*,*) 
	  stop 
   endif

   return
 end subroutine InterchangeStrokMatrix


 ! ондопнцпюллю йнппейрхпнбйх ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх 
 ! дкъ яксвюъ пеьемхъ мендмнпндмни яхярелш спюбмемхи
 ! IZNAK-Cреоемэ сйюгшбючыюъ рхо йнпмъ йнрнпши хыеряъ кебши IZNAK=1 хкх опюбши IZNAK=2  
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! Npoint-вхякн рнвей
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи опнлефсрнвмни мнплхпнбйх
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! Exf-ндмнщкейрпнммюъ тсмйжъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM-бяонлнцюрекэмше люяяхбш
 ! FunMO-б дюммнл яксвюе бяонлнцюрекэмши люяяхб
 ! RO1X(Npoint)-люяяхб опнхгбндмни мнбни оепелеммни он ярюпни
 ! H-ьюц
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ
 ! RnormMO-мнплхпнбйю лнкейскъпмни нпахрюкх онксвеммни мю опедшдсыеи хрепюжхх (йпхрепхи опюбхкэмнярх онксвемхъ ндмнщкейрпнммни щмепцхх)  
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! RdeltaZnak(4)-люяяхб нопедекъчыхи гмюйх опхйнппейрхпнбйх ндмнщкейрпнммни щмепцхх
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя йнппейрхпнбйх хяонкэгсеряъ б пефхле йнппейрхпнбйх гмювемхи тсмйжхх мю цпюмхжюу хмрепбюкю  
 ! IcorenZZ-оюпюлерп сйюгшбючыхи мюидем йнпемэ хкх мер
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMaxLimid-люйяхлюкэмне вхякн цюплнмхй лнкейскъпмни нпахрюкх
 ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
 ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
 ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
 ! IndexExitRerern-хмдейя сйюгшбючыхи вепег яйнкэйн онбрнпемхи бшунд асдер нясыеярбкем 
 subroutine CorrectionEnergy(IZNAK,NEpsIter,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,DetM,ExfM,Exf,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,EPS,RdeltaZnak,IndexVilky,INDEXDF,IIKK,IndexCorrect,IcorenZZ,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO,IndexExitRerern)  
   implicit none
   integer::NEpsIter,Kpoint,IndexM,NumeroMin,NumeroMax,IparametrDeterEnergy,Npoint,IZNAK,IndexVilky,INDEXDF,IIKK,IndexCorrect,IcorenZZ,NMO,NumeroMaxLimid,IndexExitRerern
   real(8)::Exf,Henergy,RnormMO,EPS,H
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::UNout,UNin,DetM,ExfM,AX,BX,CX,RO1X,RdeltaZnak
   real(8),dimension(:,:)::Uout,Uin,VNout,VNin,AZ,BZ,E,FunMO
   real(8),dimension(:,:,:)::Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,ISUM,klkl,IndexRepete
   real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,DDD,DMtemp

   ! нясыеярбкъел пюявер мнплхпнбнвмнцн йнщттхжхемрю (дкъ йнппейрхпнбйх ндмнщкейрпнммни щмепцхх)
   ! явхршбюел оюпюлерпш опнцнмйх дкъ тнплхпнбюмхъ яхярелш спюбмемхи
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,VNout)
   call ReadWriteMassiv1(2,Kpoint,IndexM,Uout,UNout)   
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,VNin) 
   call ReadWriteMassiv1(2,Kpoint+1,IndexM,Uin,UNin)  
   
   80685 FORMAT(2X,100(I4,1X))
   79685 FORMAT(2X,100(F20.12,1X))
   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   !DO IIDD=1,IndexM
   !   E(IIDD,IIDD)=1.D0
   !ENDDO
   FORALL (IIDD=1:IndexM) E(IIDD,IIDD)=1.D0 
   
   ! тнплхпсел люрпхжс яхярелш  
   E=E-MATMUL(VNout,VNin)
   !WRITE(6,*) 'EEEEEX'
   !DO IIDD=1,IndexM
   !   WRITE(6,79685) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! тнплхпсел ярнкаеж ябнандмшу вкемнб
   AX=MATMUL(VNout,UNin)+UNout
   
   !WRITE(6,*) 'AAAAX'
   !WRITE(6,79685) (AX(IBBB),IBBB=1,IndexM)
   
  
   ! пеьюел онксвеммсч яхярелс спюбмемхи
   call SolutionSistemLinearEquationFull(E,AX,IXM,BX) 
    
  
   !WRITE(*,*) 'BBBBX'
   !WRITE(*,79685) (BX(IBBB),IBBB=1,IndexM)
   !READ(*,*)  
   !WRITE(400,*) Exf,BX(1)
   
   ! гюохяшбюел онксвеммше гмювемхъ лнкейскъпмни нпахрюкх б рнвйе яьхбйх
   call ReadWriteMassiv1X(1,Kpoint+2,NumeroMin,NumeroMax,FunMO,BX)

   ! бняярюмюбкхбюел гмювемхъ оюпжхюкэмшу цюплнмхй (лнкейскъпмни нпахрюкх) бн бяеу рнвйюу дкъ онксвемхъ мнплхпнбнвмни йнмярюмрш 
   call FunctionRecoveryCorrection(Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,CX,VNout,UNout,VNin,UNin,RO1X,RnormConstas,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)
   
     
   WRITE(17,*) Exf,RnormConstas,RnormMO
   !WRITE(*,*) Exf,RnormConstas,RnormMO
   !WRITE(*,*) 'DDD,klkl'
   !READ(*,*) DDD,klkl
   !IF(klkl.EQ.-1) THEN
   !   Exf=ExfTemp
   !  WRITE(*,*) 'CVCV',Exf
   !ENDIF
   !ExfTemp=Exf
   
   
  
   
   
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(2)=RnormConstas
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(3)=RnormConstas
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(4)=DetM(1)
            DetM(1)=DetM(2)
            DetM(2)=DetM(3)
            DetM(3)=RnormConstas
        ENDIF
	  ELSE

        IF(INDEXDF.EQ.0) THEN 
            ! йнппейрхпсел цпюмхжш хмрепбюкю
			! хмдейя сйюгшбючыхи, врн йнппейрхпнбйю нясыеярбкемю 
			IndexCorrect=1
            DetM(1)=RnormConstas
		   ELSE
		    IF(INDEXDF.EQ.1) THEN
                DetM(4)=DetM(1)
                DetM(1)=DetM(2)
                DetM(2)=DetM(3)
                DetM(3)=RnormConstas
		      ELSE
		       ! йнпемэ оноюк б оепбши хрепбюк
               IF(IIKK.EQ.1) THEN
                  DetM(2)=DetM(3)
                  DetM(3)=RnormConstas
		       ENDIF
               ! йнпемэ оноюк б брнпни хрепбюк
               IF(IIKK.EQ.2) THEN
                  DetM(1)=DetM(3)
                  DetM(3)=RnormConstas
		       ENDIF
            ENDIF
        ENDIF  
   ENDIF
   
   
   !WRITE(*,*) 'ZERO ITERASTION2'
   !write(*,80685)  NEpsIter,IndexVilky,INDEXDF,IIKK
   !write(*,79685)  DetM(1),DetM(3),DetM(2)
   !write(*,79685)  RnormMO,DABS(DetM(3)-RnormMO),DABS(DetM(2)-RnormMO),RcoffGran
   !write(*,79685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,79685)  Exf
   !read(*,*) 
    		    
      
   ! хмрепбюк йнпмъ 
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл юмюкхг онксвеммнцн гмювемхъ
      IF(NEpsIter.EQ.1) THEN
           IF((DetM(2)-RnormMO).GT.0.D0) THEN
	           Exf=Exf*(1.D0+(-1.D0)**IZNAK*Henergy)
               ExfM(3)=Exf
	       ENDIF
           IF((DetM(2)-RnormMO).LT.0.D0) THEN
	       	   Exf=Exf*(1.D0-(-1.D0)**IZNAK*Henergy)
               ExfM(3)=Exf
	       ENDIF
           IF((DetM(2)-RnormMO).EQ.0.D0) THEN
	           ! сякнбхе бшундю бшонкмемн
	           IparametrDeterEnergy=2 
			   ! йнпемэ мюидем
	           IcorenZZ=1 
	       ENDIF
         ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(2)-RnormMO)*(DetM(3)-RnormMO).LT.0.D0) THEN
		         !WRITE(*,*) 'VILKA'
				 !READ(*,*)
				 ! опнбепъел сцкнбни йнщттхжхемр дкъ рнцн, врн бшъямхрэ опюбхкэмюъ накюярэ хкх мер
				 ! нясыеярбкъел пюявер сцкнбнцн йнщттхжхемрю
                 Kcoff=(DetM(3)-DetM(2))/(ExfM(3)-ExfM(2))
				 IF(IZNAK.EQ.1) THEN
				    ! опнбепъел сцкнбни йнщттхжхемр дкъ кебнцн йнпмъ
                    IF(Kcoff.LT.0.D0) THEN
                       ! намюпсфеммши хмрепбюк йнпмъ мюундхряъ б накюярх дпсцни щмепцхх
					   ! йнпемэ ме мюидем
                       ! сякнбхе бшундю бшонкмемн
	                   IparametrDeterEnergy=2 
			           ! йнпемэ ме мюидем
	                   IcorenZZ=0 
					   !WRITE(*,*) 'EXIT1'
					   !READ(*,*) 
					   return
					ENDIF
				 ENDIF  
                 IF(IZNAK.EQ.2) THEN
				    ! опнбепъел сцкнбни йнщттхжхемр дкъ опюбнцн йнпмъ
                    IF(Kcoff.GT.0.D0) THEN
                       ! намюпсфеммши хмрепбюк йнпмъ мюундхряъ б накюярх дпсцни щмепцхх
					   ! йнпемэ ме мюидем
                       ! сякнбхе бшундю бшонкмемн
	                   IparametrDeterEnergy=2 
			           ! йнпемэ ме мюидем
	                   IcorenZZ=0 
					   !WRITE(*,*) 'EXIT2'
					   !READ(*,*) 
					   return
					ENDIF
				 ENDIF  
		         ! сярюмнбкем хмрепбюк йнпмъ оепеундх б пефхл онхяйю йнпмъ 
	             IndexVilky=1                    
			   ELSE

			     ! нясыеярбкъел пюявер сцкнбнцн йнщттхжхемрю
                 Kcoff=(DetM(3)-DetM(2))/(ExfM(3)-ExfM(2))

                 IF(IZNAK.EQ.1) THEN
				    ! опнбепъел сцкнбни йнщттхжхемр дкъ кебнцн йнпмъ
                    IF(Kcoff.LT.0.D0) THEN
                       ! онхяй бшьек гю цпюмхжш дюммни накюярх щмепцхи
					   ! йнпемэ ме мюидем
                       ! сякнбхе бшундю бшонкмемн
	                   IparametrDeterEnergy=2 
			           ! йнпемэ ме мюидем
	                   IcorenZZ=0 
					   return
					ENDIF
				 ENDIF  
                 IF(IZNAK.EQ.2) THEN
				    ! опнбепъел сцкнбни йнщттхжхемр дкъ опюбнцн йнпмъ
                    IF(Kcoff.GT.0.D0) THEN
                      ! онхяй бшьек гю цпюмхжш дюммни накюярх щмепцхи
					   ! йнпемэ ме мюидем
                       ! сякнбхе бшундю бшонкмемн
	                   IparametrDeterEnergy=2 
			           ! йнпемэ ме мюидем
	                   IcorenZZ=0 
					   return
					ENDIF
				 ENDIF  

                 ! сярюмюбкхбюел б йюйни накюярх мюундъряъ рнвйх  
		         IF(DetM(3).GT.RnormMO) THEN
                    ! сярюмюбкхбюел гмюй сцкнбнцн йнщттхжхемрю
					IF(Kcoff.GT.0.D0) THEN
                         Exf=Exf*(1.D0-Henergy)
                         ExfM(1)=ExfM(2)
                         ExfM(2)=ExfM(3)
			             ExfM(3)=Exf
                       ELSE
					     Exf=Exf*(1.D0+Henergy)
						 ExfM(1)=ExfM(2)
                         ExfM(2)=ExfM(3)
			             ExfM(3)=Exf
					ENDIF
                 ENDIF
				 IF(DetM(3).LT.RnormMO) THEN
				    ! сярюмюбкхбюел гмюй сцкнбнцн йнщттхжхемрю
					IF(Kcoff.GT.0.D0) THEN
                         Exf=Exf*(1.D0+Henergy)
						 ExfM(1)=ExfM(2)
                         ExfM(2)=ExfM(3)
			             ExfM(3)=Exf
                       ELSE
					     Exf=Exf*(1.D0-Henergy)
                         ExfM(1)=ExfM(2)
                         ExfM(2)=ExfM(3)
			             ExfM(3)=Exf
					ENDIF
				 ENDIF
				 IF(DetM(3).EQ.RnormMO) THEN
				    ! сякнбхе бшундю бшонкмемн
	                IparametrDeterEnergy=2 
				 ENDIF 
           ENDIF      
	  ENDIF
   ENDIF
   
   ! хмрепбюк дкъ йнпмъ нопедекем
   IF(IndexVilky.EQ.1) THEN
        ! тхйяхпсел хмдейя онбрнпемхи
        IndexRepete=0


        ! сярюмюбкхбюел б йюйнл хмрепбюке мюундхряъ йнпемэ
	    IF(IndexCorrect.EQ.1) THEN 
           ! сярюмюбкхбюел хмрепбюк б йнрнпнл мюундхряъ йнпемэ
		   IF((DetM(2)-RnormMO)*(DetM(1)-RnormMO).LT.0.D0) THEN
		       ExfM(3)=ExfM(1)
		  	   DetM(3)=DetM(1)   
		   ENDIF
           IF((DetM(1)-RnormMO)*(DetM(3)-RnormMO).LT.0.D0) THEN
		       ExfM(2)=ExfM(1)
			   DetM(2)=DetM(1)   
		   ENDIF
		
	       ! йнппейрхпнбйю опнбедемю ноепюжхъ гюйнмвемю
		   IndexCorrect=0 
        ENDIF
	    

        !  опнбепъел менаундхлю хкх мер йнппейрхпнбйю цпюмхж хмрепбюкю он Y (бднкэ нях тсмйжхх)
        ! ббндхл йнщттхжхемр сйюгшбючыхи мю нрйкнмемхе гмювемхи тсмйжхх мю цпюмхжюу хмрепбюкю    
	    RcoffGran=DABS((DetM(2)-RnormMO)/(DetM(3)-RnormMO))
		! еякх йнппейрхпнбйю опнькю х мювхкяъ онхяй йнпмъ б онксвеммнл хмрепбюке рн йнппейрхпнбйю ме днкфмю опнхяундхрэ
		IF(INDEXDF.NE.0) THEN
		   RcoffGran=1.D0
		ENDIF 
	    ! йнппейрхпнбйю опнхяундхр еякх цпюмхжш нркхвючяъ б 5 пюг
	    IF(RcoffGran.GT.5.D0.OR.RcoffGran.LT.0.2D0) THEN
             ! сярюмюбкхбюел пюяонкнфемхе цпюмхж 
		     IF(ExfM(2).LT.ExfM(3)) THEN
			     Exf=ExfM(2)+0.5D0*(ExfM(3)-ExfM(2)) 	
				 ExfM(1)=Exf	  
                ELSE
				 Exf=ExfM(3)+0.5D0*(ExfM(2)-ExfM(3)) 	
				 ExfM(1)=Exf	  
  	         ENDIF
           
          ELSE
		      !WRITE(*,*) 'SSSGT'
			  !READ(*,*) 
		      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
              INDEXDF=INDEXDF+1
              ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	          IF(INDEXDF.EQ.1) THEN
                   Exf=((RnormMO-DetM(2))*ExfM(3)-(RnormMO-DetM(3))*ExfM(2))/(DetM(3)-DetM(2))
		           ExfM(1)=ExfM(2)
                   ExfM(2)=ExfM(3)
		           ExfM(3)=Exf
		         ELSE
                   ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		           ! оепбши хмрепбюк
		           IF((DetM(1)-RnormMO)*(DetM(3)-RnormMO).LT.0.D0) THEN
                      IIKK=1
			          Exf=((RnormMO-DetM(1))*ExfM(3)-(RnormMO-DetM(3))*ExfM(1))/(DetM(3)-DetM(1))
		              ! опнбепъел еCрэ нркхвхе лефдс мнбшл опхакхфемхе х ярюпшл опхакхфемхел
					  IF(DABS(Exf-ExfM(3)).EQ.0.D0) THEN
					     Exf=Exf*0.99D0  
						 ! тхйяхпсел оноюдюмхе яхярелш б пефхл онбрнпемхъ 
					     IndexExitRerern=IndexExitRerern+1
						 ! пефхл онбрнпемхъ бйкчвем
						 IndexRepete=1
					  ENDIF
                      ! еякх оноюдемхе яхярелш б дюммши пефхл нясыеярбкемн 20 пюг
					  ! пеьемхе явхрюеряъ мюидеммшл
					  IF(IndexExitRerern.GT.20) THEN
                         Exf=Exf/0.99D0
						 ! пефхл онбрнпемхи нрйкчвем
                         IndexRepete=0
						 ! тхйяхпсел врн онксвемн пеьемхе  
						 DetM(3)=RnormMO  
					  ENDIF

					  ExfM(2)=ExfM(3)
                      ExfM(3)=Exf
		           ENDIF
		           ! брнпни хмрепбюк
                   IF((DetM(3)-RnormMO)*(DetM(2)-RnormMO).LT.0.D0) THEN
		              IIKK=2
                      Exf=((RnormMO-DetM(3))*ExfM(2)-(RnormMO-DetM(2))*ExfM(3))/(DetM(2)-DetM(3))
		              ! опнбепъел еCрэ нркхвхе лефдс мнбшл опхакхфемхе х ярюпшл опхакхфемхел
					  IF(DABS(Exf-ExfM(3)).EQ.0.D0) THEN
					     Exf=Exf*0.99D0 
						 ! тхйяхпсел онондюмхе яхярелш б пефхл онбрнпемхъ 
					     IndexExitRerern=IndexExitRerern+1
						 ! пефхл онбрнпемхъ бйкчвем
						 IndexRepete=1
					  ENDIF
                      
					  ! еякх оноюдемхе яхярелш б дюммши пефхл нясыеярбкемн 20 пюг
					  ! пеьемхе явхрюеряъ мюидеммшл
					  IF(IndexExitRerern.GT.20) THEN
                         Exf=Exf/0.99D0
						 ! тхйяхпсел врн онксвемн пеьемхе  
						 DetM(3)=RnormMO  
						 ! пефхл онбрнпемхи нрйкчвем
                         IndexRepete=0
					  ENDIF

					  ExfM(1)=ExfM(3)
                      ExfM(3)=Exf
		           ENDIF
	          ENDIF
        ENDIF
		
       
		! опнбепъел бшонкмхкняэ кх сякнбхе
        ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
        IF(DABS(DetM(3)-RnormMO).LT.EPS) THEN
		   ! опнбепъел бкчвем пефхл онбрнпемхи хкх мер
		   IF(IndexRepete.EQ.1) THEN
              ! пефхл бйкчвем. нрйкчвюел пефхл онбрнпемхи
              Exf=Exf/0.99D0
		   ENDIF
           !       WRITE(*,*) 'SSDD EXIT'
  	       ! сякнбхе бшундю бшонкмемн
           IparametrDeterEnergy=2
	       ! йнпемэ мюидем
	       IcorenZZ=1 
        ENDIF		  		   
   ENDIF

   
  
     !WRITE(*,*) 'ITERASTION2'
     !write(*,80685)  NEpsIter,IndexVilky,INDEXDF,IIKK
     !write(*,79685)  DetM(1),DetM(3),DetM(2)
     !write(*,79685)  RnormMO,DABS(DetM(3)-RnormMO),DABS(DetM(2)-RnormMO),RcoffGran
     !write(*,79685)  ExfM(1),ExfM(3),ExfM(2)
     !write(*,79685)  Exf
     !read(*,*) 
   return
 end subroutine CorrectionEnergy




 ! ондопнцпюллю нясыеярбкъел пюявер мнплш яннрберярбсчыеи дюммни щмепцхх
 ! дкъ яксвюъ пеьемхъ мендмнпндмни яхярелш спюбмемхи
 ! RnormConstas-онксвеммюъ мнплю тсмйжхх  
 ! Kpoint-рнвйю яьхбйх
 ! Npoint-вхякн рнвей
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи опнлефсрнвмни мнплхпнбйх
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! Exf-ндмнщкейрпнммюъ тсмйжъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM-бяонлнцюрекэмше люяяхбш
 ! FunMO-б дюммнл яксвюе бяонлнцюрекэмши люяяхб
 ! RO1X(Npoint)-люяяхб опнхгбндмни мнбни оепелеммни он ярюпни
 ! H-ьюц
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ
 ! RnormMO-мнплхпнбйю лнкейскъпмни нпахрюкх онксвеммни мю опедшдсыеи хрепюжхх (йпхрепхи опюбхкэмнярх онксвемхъ ндмнщкейрпнммни щмепцхх)  
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! RdeltaZnak(4)-люяяхб нопедекъчыхи гмюйх опхйнппейрхпнбйх ндмнщкейрпнммни щмепцхх
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя йнппейрхпнбйх хяонкэгсеряъ б пефхле йнппейрхпнбйх гмювемхи тсмйжхх мю цпюмхжюу хмрепбюкю  
 ! IcorenZZ-оюпюлерп сйюгшбючыхи мюидем йнпемэ хкх мер
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMaxLimid-люйяхлюкэмне вхякн цюплнмхй лнкейскъпмни нпахрюкх
 ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
 ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
 ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
 subroutine CalculationNormFunE(RnormConstas,Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)  
   implicit none
   integer::Kpoint,IndexM,NumeroMin,NumeroMax,Npoint,NMO,NumeroMaxLimid
   real(8)::RnormConstas,H
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::UNout,UNin,AX,BX,CX,RO1X
   real(8),dimension(:,:)::Uout,Uin,VNout,VNin,AZ,BZ,E,FunMO
   real(8),dimension(:,:,:)::Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,ISUM,klkl
   

   ! нясыеярбкъел пюявер мнплхпнбнвмнцн йнщттхжхемрю (дкъ йнппейрхпнбйх ндмнщкейрпнммни щмепцхх)
   ! явхршбюел оюпюлерпш опнцнмйх дкъ тнплхпнбюмхъ яхярелш спюбмемхи
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,VNout)
   call ReadWriteMassiv1(2,Kpoint,IndexM,Uout,UNout)   
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,VNin) 
   call ReadWriteMassiv1(2,Kpoint+1,IndexM,Uin,UNin)  
   
   80685 FORMAT(2X,100(I4,1X))
   79685 FORMAT(2X,100(F20.12,1X))
   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   !DO IIDD=1,IndexM
   !   E(IIDD,IIDD)=1.D0
   !ENDDO
   FORALL (IIDD=1:IndexM) E(IIDD,IIDD)=1.D0 
   
   ! тнплхпсел люрпхжс яхярелш  
   E=E-MATMUL(VNout,VNin)
   !WRITE(6,*) 'EEEEEX'
   !DO IIDD=1,IndexM
   !   WRITE(6,79685) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! тнплхпсел ярнкаеж ябнандмшу вкемнб
   AX=MATMUL(VNout,UNin)+UNout
   
   !WRITE(6,*) 'AAAAX'
   !WRITE(6,79685) (AX(IBBB),IBBB=1,IndexM)
   
  
   ! пеьюел онксвеммсч яхярелс спюбмемхи
   call SolutionSistemLinearEquationFull(E,AX,IXM,BX) 
    
  
   !WRITE(*,*) 'BBBBX'
   !WRITE(*,79685) (BX(IBBB),IBBB=1,IndexM)
   !READ(*,*)  
   !WRITE(400,*) Exf,BX(1)
   
   ! гюохяшбюел онксвеммше гмювемхъ лнкейскъпмни нпахрюкх б рнвйе яьхбйх
   call ReadWriteMassiv1X(1,Kpoint+2,NumeroMin,NumeroMax,FunMO,BX)

   ! бняярюмюбкхбюел гмювемхъ оюпжхюкэмшу цюплнмхй (лнкейскъпмни нпахрюкх) бн бяеу рнвйюу дкъ онксвемхъ мнплхпнбнвмни йнмярюмрш 
   call FunctionRecoveryCorrection(Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,CX,VNout,UNout,VNin,UNin,RO1X,RnormConstas,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)
   
   
   
    return
 end subroutine CalculationNormFunE





 

 ! ондопнцпюллю онксвюер гюбхяхлнярэ хмрецпюкю нпрнцнмюкэмнярх нр ндмнщкейрпнммни щмепцхх
 ! б хмрепбюке Exf>Emin  
 ! дкъ яксвюъ пеьемхъ мендмнпндмни яхярелш спюбмемхи
 ! Kpoint-рнвйю яьхбйх
 ! Npoint-вхякн рнвей
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Exf-ндмнщкейрпнммюъ тсмйжъ
 ! Emin-цпюмхжю хмрепбюкю
 ! VNout,UNout,VNin,UNin,AZ,BZ,E,AX,BX,CX,IXM-бяонлнцюрекэмше люяяхбш
 ! FunMO-б дюммнл яксвюе бяонлнцюрекэмши люяяхб
 ! RO1X(Npoint)-люяяхб опнхгбндмни мнбни оепелеммни он ярюпни
 ! H-ьюц
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ
 ! RnormMO-мнплхпнбйю лнкейскъпмни нпахрюкх онксвеммни мю опедшдсыеи хрепюжхх (йпхрепхи опюбхкэмнярх онксвемхъ ндмнщкейрпнммни щмепцхх)  
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMaxLimid-люйяхлюкэмне вхякн цюплнмхй лнкейскъпмни нпахрюкх
 ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
 ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
 ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
 subroutine SpectrumEnergy(Kpoint,Npoint,IndexM,NumeroMin,NumeroMax,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,Exf,Emin,AZ,BZ,E,AX,BX,CX,IXM,FunMO,RO1X,H,Henergy,RnormMO,IparametrDeterEnergy,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)  
   implicit none
   integer::Niteration,NEpsIter,Kpoint,IndexM,NumeroMin,NumeroMax,IparametrDeterEnergy,Npoint,IZNAK,IndexVilky,INDEXDF,IIKK,IndexCorrect,IcorenZZ,NMO,NumeroMaxLimid
   real(8)::Exf,Emin,Henergy,RnormMO,EPS,H
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::IXM,NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::UNout,UNin,AX,BX,CX,RO1X
   real(8),dimension(:,:)::Uout,Uin,VNout,VNin,AZ,BZ,E,FunMO
   real(8),dimension(:,:,:)::Vout,Vin,AlfaCoffMO
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,ISUM,klkl
   real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,DDD,ExfTemp,DMtemp

   ! нясыеярбкъел пюявер мнплхпнбнвмнцн йнщттхжхемрю (дкъ йнппейрхпнбйх ндмнщкейрпнммни щмепцхх)
   ! явхршбюел оюпюлерпш опнцнмйх дкъ тнплхпнбюмхъ яхярелш спюбмемхи
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,VNout)
   call ReadWriteMassiv1(2,Kpoint,IndexM,Uout,UNout)   
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,VNin) 
   call ReadWriteMassiv1(2,Kpoint+1,IndexM,Uin,UNin)  
   
   80685 FORMAT(2X,100(I4,1X))
   79685 FORMAT(2X,100(F20.12,1X))
   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш  
   E=E-MATMUL(VNout,VNin)
   !WRITE(6,*) 'EEEEEX'
   !DO IIDD=1,IndexM
   !   WRITE(6,79685) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! тнплхпсел ярнкаеж ябнандмшу вкемнб
   AX=MATMUL(VNout,UNin)+UNout
   
   !WRITE(6,*) 'AAAAX'
   !WRITE(6,79685) (AX(IBBB),IBBB=1,IndexM)
   
  
   ! пеьюел онксвеммсч яхярелс спюбмемхи
   call SolutionSistemLinearEquationFull(E,AX,IXM,BX) 
    
  
   !WRITE(*,*) 'BBBBX'
   !WRITE(*,79685) (BX(IBBB),IBBB=1,IndexM)
   !READ(*,*)  
   !WRITE(400,*) Exf,BX(1)
   
   ! гюохяшбюел онксвеммше гмювемхъ лнкейскъпмни нпахрюкх б рнвйе яьхбйх
   call ReadWriteMassiv1X(1,Kpoint+2,NumeroMin,NumeroMax,FunMO,BX)

   ! бняярюмюбкхбюел гмювемхъ оюпжхюкэмшу цюплнмхй (лнкейскъпмни нпахрюкх) бн бяеу рнвйюу дкъ онксвемхъ мнплхпнбнвмни йнмярюмрш 
   call FunctionRecoveryCorrection(Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,CX,VNout,UNout,VNin,UNin,RO1X,RnormConstas,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)
   
   WRITE(*,*) Exf,RnormConstas,RnormMO
   WRITE(18,*) Exf,RnormConstas,RnormMO
   Exf=Exf*(1.D0-Henergy) 
   
   ! сякнбхе бшундю хг ондопнцпюллш сйюгшбюел мхфмхч цпюмхжс   
   IF(Exf.LT.Emin) THEN
      ! сякнбхе бшундю бшонкмемн
      IparametrDeterEnergy=2
   ENDIF
   
   return
 end subroutine SpectrumEnergy




 ! ондопнцпюллю бняярюмнбкемхъ гмювемхи лнкейскъпмни нпахрюкх бн бяеу рнвйюу хг рнвйх яьхбйх 
 ! бонлнцюрекэмюъ ондопнцпюллю дкъ йнппейрхпнбйх ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Npoint-вхякн рнвей
 ! NumeroMin-мнлеп оепбни цюплнмхйх
 ! NumeroMax-мнлеп онякедмеи цюплнмхйх
 ! H-ьюц
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Uout(IndexM,Npoint)-люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! Uin(IndexM,Npoint)-люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! FunMO(NGarmonik,Npoint+2)-люяяхб гмювемхи пюдхюкэмни вюярх тсмйжхх лнкейскъпмни нпахрюкх
 ! AX,BX,CX,VNout,UNout,VNin,UNin-бяонлнцюрекэмше люяяхбш
 ! RO1X(Npoint)-бяонлнцюрекэмши люяяхб  оепбюъ опнхгбндмюъ мнбни оепелеммни он ярюпни 
 ! RnormConstas-мнплхпнбнвмюъ йнмярюмрю
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NumeroMaxLimid-люйяхлюкэмне вхякн цюплнмхй лнкейскъпмни нпахрюкх
 ! NumbreLigand(NumbreMO)-люяяхб вхякю кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб мнлепнб кхцюмднб тсмйжхх йнрнпшу нохяшбючр бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NFunLigands(NumbreMO,NumbreLigand(NumbreMO))-люяяхб вхякю тсмйжхи йюфднцн кхцюмдю нохяшбючыху бшяьхе цюплнмхйх дюммни лнкейскъпмни нпахрюкх
 ! NumbreFunctionLig(NumbreMO,NumbreLigand(NumbreMO),NFunLigands(NumbreMO,NumbreLigand(NumbreMO)))-люяяхб мнлепнб тсмйжхи кхцюмднб йнрнпше яннрберярбсчр дюммни лнкейскъпмни нпахрюкх
 ! RfunLigands(Nnuclei,NumbreFunctionLig,NumbreGarmLigand,Npoint+2)-люяяхб гмювемхи тсмйжхи кхцюмднб
 ! AlfaCoffMO(NumbreMO,NumbreLigand,NumbreFunLigand)-люяяхб йнщттхжхемрнб я йнрнпшлх тсмйжхи кхцюмдю бундър б ярпсйрспс лнкейскъпмни нпахрюкх 
 subroutine FunctionRecoveryCorrection(Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,H,Vout,Uout,Vin,Uin,FunMO,AX,BX,CX,VNout,UNout,VNin,UNin,RO1X,RnormConstas,NMO,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO)
   implicit none
   integer::Kpoint,IndexM,Npoint,NumeroMin,NumeroMax,NMO,NumeroMaxLimid
   real(8)::H,RnormConstas
   integer,dimension(:)::NumbreLigand
   integer,dimension(:,:)::NLigands,NFunLigands
   integer,dimension(:,:,:)::NumbreFunctionLig
   real(8),dimension(:)::AX,BX,CX,UNout,UNin,RO1X
   real(8),dimension(:,:)::Uout,Uin,FunMO,VNout,VNin
   real(8),dimension(:,:,:)::Vout,Vin,AlfaCoffMO 
   real(8),dimension(:,:,:,:)::RfunLigands
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IXXZ,IILV,IIDF,IOFG
   
   
   10000 FORMAT(2X,100(1X,F20.16))
   ! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (1,Kpoint-1)
   DO IXXZ=Kpoint-1,1,-1 
      ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ+1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	  call ReadWriteMassiv1X(2,IXXZ+3,NumeroMin,NumeroMax,FunMO,AX)
	  !WRITE(6,*) 'POINTGGGGGG',IXXZ+1
	  !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
      ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	  call ReadWriteMassiv2(2,IXXZ,IndexM,Vout,VNout)
	  call ReadWriteMassiv1(2,IXXZ,IndexM,Uout,UNout)  
      !WRITE(6,*) 'POINT',IXXZ
	  !DO IIDF=1,IndexM
      !   WRITE(6,10000) (VNout(IIDF,IOFG),IOFG=1,IndexM)
	  !ENDDO
	  !WRITE(6,*)
	  ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	  BX=UNout+MATMUL(VNout,AX)
	  ! гюохяшбюел пегскэрюр 
      call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
   ENDDO
   
   ! бняярюмюбкхбюел гмювемхе пюдхюкэмни вюярх лнкейскъпмни нпахрюкх б рнвйюу (Kpoint+1,Npoint)
   DO IXXZ=Kpoint+1,Npoint 
      ! явхршбюел гмювемхе тсмйжхх б рнвйе  IXXZ-1 +(2)-мювюкн мслепюжхх я 3 рнвйе б люяяхбе
	  call ReadWriteMassiv1X(2,IXXZ+1,NumeroMin,NumeroMax,FunMO,AX)
	  !WRITE(6,*) 'POINTGGGGGG',IXXZ-1
	  !WRITE(6,10000) (AX(IOFG),IOFG=1,IndexM)
      ! явхршбюел гмювемхе оюпюлерпнб опнцнмйх
	  call ReadWriteMassiv2(2,IXXZ,IndexM,Vin,VNin)
	  call ReadWriteMassiv1(2,IXXZ,IndexM,Uin,UNin)  
      ! нясыеярбкъел пюявер гмювемхи лнкейскъпмни нпахрюкх (оюпжхюкэмшу цюплнмхй) б рнвйе IXXZ 
	  BX=UNin+MATMUL(VNin,AX)
	  ! гюохяшбюел пегскэрюр 
      call ReadWriteMassiv1X(1,IXXZ+2,NumeroMin,NumeroMax,FunMO,BX)
   ENDDO   

   ! тнплхпсел ярпсйрспс бшяьху цюплнмхй лнкейскъпмни нпахрюкх
		
   ! нопедекъел йнщттхжхемрш юкэтю
   call CalculationCoffAlfaMO(NMO,Npoint,H,NumeroMax,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,FunMO,RfunLigands,AlfaCoffMO,RO1X) 
        
   ! тнплхпсел ярпсйрспс бшяьху цюплнмхй лнкейскъпмшу нпахрюкеи я свернл онксвеммшу йнщттхжхемрнб юкэтю
   call FormStructureMO(NMO,Npoint,H,NumeroMax,NumeroMaxLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,FunMO,RfunLigands,AlfaCoffMO) 


   
   ! мнплхпсел онксвеммсч лнкейскъпмсч нпахрюкэ
   RnormConstas=CalculationNormalizationConsta(NumeroMin,NumeroMaxLimid,Npoint,H,RO1X,FunMO) 


   !DO IXXZ=1,Npoint
   !   WRITE(6,10000) RR(IXXZ),(FunMO(IILV,IXXZ+2)/DSQRT(RnormConstas*RO1X(IXXZ)),IILV=NumeroMin,NumeroMax)
   !ENDDO


   !WRITE(*,*) 'NORMA',RnormConstas
   !READ(*,*)
      
   return
 end subroutine FunctionRecoveryCorrection

 ! ондопнцпюллю нясыеярбкъчыюъ онхяй хмрепбюкю мсфмнцн йнпмъ ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх дкъ яксвюъ оепбни хрепюжхх
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NNEL- цкюбмне йбюмрнбне вхякн
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи дереплхмюмрю
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! AZ,BZ,E-бяонлнцюрекэмше люяяхбш
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ (оняке нйнмвюмхъ опнцпюллш дюер мскхбне опхакхфемхе дкъ йнпмъ)
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! IndexVilky=2-пефхл онхяйю йнпмъ 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя сйюгшбючыхи вхякн йнпмеи мюидеммшу мю бяел хмрепбюке  
 ! EnergyZeroZ(NumbreMO,NNEL+1,2)-люяяхб йнпмеи ндмнпндмнцн спюбмемхъ 
 ! IKLKORX,IOPZ,IKLRXX-хмдейяш йнппейрхпнбйх дкъ сярюмюбкемхъ йнпемэ хкх пюгпшб
 ! INDEXAD-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу б яксвюе йнцдю хмрепбюк менопедекем 
 ! RcoffChengH-йнщттхжхемр хглемемхъ ьюцю
 ! ProDet(2)-люяяхб опнхгбндмшу
 ! INDEXRShag-пефхл хглемемхъ ьюцю
 ! IIHD-хмдейя слемэьемхъ ьюцю
 ! IIOXDPZ-йкчв сйюгшбючыхи врн лш оноюкх б пефхл слемэьемхъ ьюцю х йнмрпнкхпсел бнглнфмнярэ оепеяйнйю
 ! INDEXPeresc-хмдейя сйюгшбючыхи вюярх хрепюжхи б пефхле бнглнфмнцн оепеяйнйю
 ! DeltaEnerPers-хмрепбюк щмепцхи б накюярх оепеяйнйю
 ! DeltaFunPers-хмрепбюк тсмйжхи б накюярх оепеяйнйю 
 subroutine SearchSolution(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,EnergyZeroZ,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::NMO,NNEL,NEpsIter,Kpoint,IndexM,IparametrDeterEnergy,IndexVilky,INDEXDF,IIKK,IndexCorrect,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Exf,Henergy,EPS,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::DetM,ExfM,ProDet
   real(8),dimension(:,:)::AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin,EnergyZeroZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,IIII,klkl,IIIFKL,IKLCOMBEC
   real(8)::DeltaEnergy,DDD,Kcoff,RcoffGran,ExfTemp,DMtemp,Rdet
   !real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,ExfTemp,DMtemp,DDD
   ! нопедекъел дереплхмюмр яхярелш
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ) 
   10685 FORMAT(2X,100(I4,1X))
   19685 FORMAT(2X,100(F20.12,1X))

   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   !WRITE(100,*) 'AZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (AZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO 
   !WRITE(100,*) 'BZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (BZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   E=E-MATMUL(AZ,BZ)
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! гюохяшбюел пегскэрюр пюяверю
   ! нопедекъел дереплхмюмр яхярелш
   Rdet=Determenant(E)

   WRITE(17,*) Exf,Rdet,RcoffChengH    !Kpoint
   !WRITE(*,*) 'ENERGY',Exf,Rdet,Kpoint
   !WRITE(*,*) 'DDD,klkl'
   !READ(*,*) DDD,klkl
   !IF(klkl.EQ.-1) THEN
   !   Exf=ExfTemp
   !	  WRITE(*,*) 'CVCV',Exf
   !  ENDIF
   !ExfTemp=Exf
   !DDD=0.0001D0
   !Exf=Exf*(1.D0-DDD)
   !write(*,*) '0 Exf',Exf
   !DeltaEnergy=1.D0
   !IndexVilky=10
  
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(4)=Rdet
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(5)=Rdet
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(1)=DetM(2)
			DetM(2)=DetM(3)
            DetM(3)=DetM(4)
            DetM(4)=DetM(5)
            DetM(5)=Rdet
        ENDIF
		! хмдейя тхйяхпсчыхи вхякн хрепюжхи б яксвюе йнцдю хмрепбюк менопедекем
		INDEXAD=INDEXAD+1
		! б мювюке бнгбпюыюеляъ й опедшдсыелс ьюцс сахпюел йнщттхжхемр онксвеммши б опедедсыел хмрепбюке 
		IF(INDEXAD.EQ.4) THEN
		   Henergy=Henergy*RcoffChengH
        ENDIF 
		! мювхмюел йнппейрхпнбюрэ ьюц 
		IF(INDEXAD.EQ.5) THEN
		   ! бнгбпюыюел йнщттхжхемр б хяундмне янярнъмхе
		   RcoffChengH=1.D0
           ProDet(1)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
	    IF(INDEXAD.EQ.6) THEN
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
		IF(INDEXAD.GT.6) THEN
           ProDet(1)=ProDet(2)
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		   ! нясыеярбкъел пюявер йнщттхжхемрю х хглемемхе ьюцю
		   call ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
		   
		   ! опнбепъел оепеунд б пефхл слемэьемхъ ьюцю
		   IF(IKLCOMBEC.EQ.1) THEN
              ! бнгбпюыюеляъ мю ндмс рнвйс мюгюд
              Exf=ExfM(4)
			  ExfM(5)=ExfM(4)
              ExfM(4)=ExfM(3)
              ExfM(3)=ExfM(2)
              ExfM(2)=ExfM(1)
              DetM(5)=DetM(4)
              DetM(4)=DetM(3)
              DetM(3)=DetM(2)
              DetM(2)=DetM(1)   
		   ENDIF
		ENDIF
   ENDIF

   IF(IndexVilky.EQ.1) THEN	
      IF(INDEXDF.EQ.1) THEN 
         DetM(1)=DetM(2)
	     DetM(2)=DetM(3)
         DetM(3)=DetM(4)
         DetM(4)=DetM(5)
         DetM(5)=Rdet
      ENDIF 
	  IF(INDEXDF.EQ.2) THEN 
         ! опнбепъел йюйни хг пефхлнб бйкчвем
		 IF(IKLKORX.EQ.1) THEN
            ! йнппейрхпнбйю яопюбю
            ! опнбепъел врнаш гмювемхе яопюбю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(3).GT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(3).LT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
         IF(IKLKORX.EQ.2) THEN
            ! йнппейрхпнбйю якебю
            ! опнбепъел врнаш гмювемхе якебю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(4).GT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(4).LT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
		
      ENDIF   
   ENDIF 

   
   IF(IndexVilky.EQ.2) THEN
      IF(INDEXDF.EQ.1) THEN
           DetM(4)=DetM(1)
           DetM(1)=DetM(2)
           DetM(2)=DetM(3)
           DetM(3)=Rdet
	     ELSE
	       ! йнпемэ оноюк б оепбши хрепбюк
           IF(IIKK.EQ.1) THEN
              DetM(2)=DetM(3)
              DetM(3)=Rdet
	       ENDIF
           ! йнпемэ оноюк б брнпни хрепбюк
           IF(IIKK.EQ.2) THEN
              DetM(1)=DetM(3)
              DetM(3)=Rdet
 	       ENDIF
      ENDIF
   ENDIF 
   		    
      
   ! опнбндхл онхяй йнпмеи мю бяел щмепцхрхвеяйнл хмрепбюке
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл оепеанп гмювемхи щмепцхх
      IF(NEpsIter.EQ.1) THEN
	        ExfM(4)=Exf
            Exf=Exf*(1.D0-Henergy)
            ExfM(5)=Exf  
          ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(4)*DetM(5)).LT.0.D0) THEN
	             IndexVilky=1  
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи бме хмрепбюкю дкъ мнбнцн яксвюъ бме хмрепбюкю  
				 INDEXAD=0
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи я слемэьеммшл ьюцнл
				 INDEXRShag=0  
				 ! йнмрпнкэ яйювйю нрйкчвюел
				 IIOXDPZ=0                 
			   ELSE
			     Exf=Exf*(1.D0-Henergy)
                 ExfM(1)=ExfM(2)
                 ExfM(2)=ExfM(3)
				 ExfM(3)=ExfM(4)
                 ExfM(4)=ExfM(5)
				 ExfM(5)=Exf  
           ENDIF      
	  ENDIF
   ENDIF
   
   ! сярюмнбкемю накюярэ оепелемш гмюйю 
   IF(IndexVilky.EQ.1) THEN
        ! опнбепъел оепелемю гмюйю ябъгюмю я бнгмхймнбемел пюгпшбю б дюммни накюярх хкх я мюкхвхел йнпмъ
        ! хмдейя сйюгшбюер, врн б дюммнл пефхле  опнцпюллю пюанрюер дкъ дюммнцн хмрепбюкю б оепбши пюг
		INDEXDF=INDEXDF+1
        IF(INDEXDF.EQ.1) THEN
             Exf=Exf*(1.D0-Henergy/FLOAT(IndexCorrect+2))
             ExfM(1)=ExfM(2)
             ExfM(2)=ExfM(3)
		     ExfM(3)=ExfM(4)
             ExfM(4)=ExfM(5)
		     ExfM(5)=Exf 
			 ! тхйяхпсел IOPZ дкъ йнппейрхпнбйх йнпмеи
			 IOPZ=4 
			 ! тхйяхпсел вхякн рнвей йнрнпше ашкх свремш опх ясфемхх накюярх 
			 IKLRXX=0 
           ELSE   
             IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(2)).GT.DABS(DetM(3))) THEN
                  ! намюпсфем хмрепбюк йнпмъ 
                  IndexCorrect=IndexCorrect+1
                  ! хглемемхе ьюцю 
				  Henergy=Henergy/FLOAT(IndexCorrect+1)
				  ! тхйяхпсел бнглнфмши йнпемэ б дюммни накюярх
                  EnergyZeroZ(NMO,IndexCorrect,1)=EnergyZeroZ(NMO,IndexCorrect,2)
				  EnergyZeroZ(NMO,IndexCorrect,2)=ExfM(3)  
                  
				  ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			      ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			      ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			      IndexVilky=0 
                  INDEXDF=0 
                  Exf=ExfM(5)
				ELSE
				  ! сярюмюбкхбюел опхвхмю б пюгпшбе хкх мер
				  ! бшъямъел ясыеярбсер йнпемэ хкх хлеер леярн пюгпшб
				  ! опх щрнл, сярюмюбкхбюел рюйне онбедемхе ябъгюмн я йнкеаюмхълх б накюярх оепелемш гмюйю
				  ! нясыеярбкъел ясфемхе накюярх оепелемш гмюйю дкъ бмеяемхъ ъямнярх б онбедемхх нопедекхрекъ б накюярх мскъ
                   
				  IIIFKL=0
				  ! сярюмюбкхбюел опхвхмс нрясрярбхъ йнпмъ
				  ! опнбепъел бшонкмемхе оепбнцн сякнбхъ
			      IF(DABS(DetM(5)).GT.DABS(DetM(4))) THEN
				     ! брнпне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (яопюбю)     
				     IKLKORX=1
                     Exf=ExfM(3)-DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
				  ! опнбепъел бшонкмемхе брнпнцн сякнбхъ
				  IF(DABS(DetM(2)).GT.DABS(DetM(3))) THEN
				     ! оепбне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (якебю)     
				     IKLKORX=2
                     Exf=ExfM(4)+DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
                  
				  ! опнбепъел яйнкэйн пюг онксвемю мнбюъ рнвйю рн еярэ яйнкэйн пюг лш слемэьюкх накюярэ оепелемш гмюйю  
                  IF(IKLRXX.EQ.3) THEN
				     ! вхякн рнвей пюбмн рпел, опх щрнл ме опнхгнькн хглемемхи б онбедемхх нопедекхрекъ б акхгх рнвйх оепелемш гмюйю 
					 ! рюйхл напюгнл хлеер леярн пюгпшб 
					 IIIFKL=0 
				  ENDIF 

				  ! опнбепъел хлеер леярн пюгпшб хкх мер 
				  IF(IIIFKL.EQ.0) THEN
			         ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			         ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			         ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			         IndexVilky=0 
                     INDEXDF=0 
					 Exf=ExfM(5) 
				  ENDIF
			 ENDIF
			
		ENDIF 
	    
		! опнбепъел бшонкмхкняэ кх сякнбхе
        ! мнлеп мюидемнцн йнпмъ он онпъдйс(хмрепбюкю б йнрнпнл мюундхряъ йнпемэ) яннрберярбсер цкюбмнлс йбюмрнбнлс вхякс
        IF(IndexCorrect.EQ.NNEL) THEN
           ! тхйяхпсел пефхл онхяйю йнпмъ б хмрепбюке
		   IndexVilky=2 
		   !WRITE(*,*) 'ZERO',IndexVilky
		   !WRITE(*,*) DetM(5),DetM(4),DetM(2),DetM(3)
		   !WRITE(*,*) ExfM(5),ExfM(4),ExfM(2),ExfM(3)
		   !READ(*,*) 
		   ! гюохяшбюел цпюмхжш хмрепбюкю 
           ExfM(2)=ExfM(4)
		   DetM(2)=DetM(4)
		   ExfM(3)=ExfM(3)
		   DetM(3)=DetM(3)
		   !WRITE(*,*) DetM(2),DetM(3)
		   !WRITE(*,*) ExfM(2),ExfM(3)
		   !READ(*,*)   
		   ! гюмскъел оепед онхяйнл йнпмъ
		   ExfM(1)=0.D0
		   DetM(1)=0.D0
		   ExfM(4)=0.D0
           DetM(4)=0.D0
		   ExfM(5)=0.D0
		   DetM(5)=0.D0
		   INDEXDF=0 
       ENDIF
		  		    
   ENDIF

   ! пефхл онхяйю йнпмъ б хмрепбюке
   IF(IndexVilky.EQ.2) THEN 
      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
      INDEXDF=INDEXDF+1
      ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	  IF(INDEXDF.EQ.1) THEN
           Exf=(DetM(3)*ExfM(2)-DetM(2)*ExfM(3))/(DetM(3)-DetM(2))
		   ExfM(1)=ExfM(2)
           ExfM(2)=ExfM(3)
		   ExfM(3)=Exf
		 ELSE
           ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		   ! оепбши хмрепбюк
		   IF((DetM(1)*DetM(3)).LT.0.D0) THEN
               IIKK=1
		       Exf=(DetM(3)*ExfM(1)-DetM(1)*ExfM(3))/(DetM(3)-DetM(1))
			   ExfM(2)=ExfM(3)
               ExfM(3)=Exf
		   ENDIF
		   ! брнпни хмрепбюк
           IF((DetM(3)*DetM(2)).LT.0.D0) THEN
		       IIKK=2
               Exf=(DetM(2)*ExfM(3)-DetM(3)*ExfM(2))/(DetM(2)-DetM(3))
		       ExfM(1)=ExfM(3)
               ExfM(3)=Exf
		   ENDIF
      ENDIF		  		   
   

      ! опнбепъел бшонкмхкняэ кх сякнбхе
      ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
      ! мювхмюел опнбепърэ сякнбхе б яксвюе намюпсфемхъ хмрепбюкю цде мюундхряъ йнпемэ
      IF(INDEXDF.GT.1) THEN
         IF(DABS(DetM(3)).LT.EPS) THEN
            ! сякнбхе бшундю бшонкмемн
   	        IparametrDeterEnergy=2 
			! гюохяшбюел йнпемэ мюидеммши
			EnergyZeroZ(NMO,IndexCorrect,2)=ExfM(3)  
         ENDIF
      ENDIF  
   ENDIF
    

  
   !WRITE(*,*) 'ITERASTION ZERO'
   !write(*,10685)  NEpsIter,IndexVilky,IndexCorrect,NNEL
   !write(*,19685)  DetM(1),DetM(3),DetM(2)
   !write(*,19685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,19685)  Exf
   !read(*,*)
   return
 end subroutine SearchSolution




 ! ондопнцпюллю нясыеярбкъчыюъ онхяй хмрепбюкю мсфмнцн йнпмъ ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх дкъ яксвюъ оепбни хрепюжхх
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NNEL- цкюбмне йбюмрнбне вхякн
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи дереплхмюмрю
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! DetMA-люяяхб гмювемхи дереплхмюмрю дкъ онхяйю йнпмъ
 ! ExfMA-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх дкъ онхяйю йнпмъ
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! AZ,BZ,E-бяонлнцюрекэмше люяяхбш
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ (оняке нйнмвюмхъ опнцпюллш дюер мскхбне опхакхфемхе дкъ йнпмъ)
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! IndexVilky=2-пефхл онхяйю йнпмъ 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя сйюгшбючыхи вхякн йнпмеи мюидеммшу мю бяел хмрепбюке  
 ! IndexSolutions-хмдейя вхякю пеьемхе, менаундхлн мюирх пеьемхъ дкъ NNEL-1, NNEL, NNEL+1, рн еярэ дкъ бшундю 
 ! менаундхлн, врнаш IndexSolutions=3 
 ! ExfSolutions(IndexSolutions)-щмепцхх мюидеммшу йнпмеи  NNEL-1, NNEL, NNEL+1
 ! ExfTemp-бяонлнцюрекэмюъ оепелеммюъ
 ! IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET-дюммше менаундхлш дкъ пюяверю б яксвюе мюунфдемхъ йнпмъ якебю дкъ мендмнпндмни яхярелш
 ! IKLKORX,IOPZ,IKLRXX-хмдейяш йнппейрхпнбйх дкъ сярюмюбкемхъ йнпемэ хкх пюгпшб
 ! EnergyZeroZ(NumbreMO,NNEL+1,2)-люяяхб йнпмеи ндмнпндмнцн спюбмемхъ 
 ! INDEXAD-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу б яксвюе йнцдю хмрепбюк менопедекем 
 ! RcoffChengH-йнщттхжхемр хглемемхъ ьюцю
 ! ProDet(2)-люяяхб опнхгбндмшу
 ! INDEXRShag-пефхл хглемемхъ ьюцю
 ! IIHD-хмдейя слемэьемхъ ьюцю
 ! IIOXDPZ-йкчв сйюгшбючыхи врн лш оноюкх б пефхл слемэьемхъ ьюцю х йнмрпнкхпсел бнглнфмнярэ оепеяйнйю
 ! INDEXPeresc-хмдейя сйюгшбючыхи вюярх хрепюжхи б пефхле бнглнфмнцн оепеяйнйю
 ! DeltaEnerPers-хмрепбюк щмепцхи б накюярх оепеяйнйю
 ! DeltaFunPers-хмрепбюк тсмйжхи б накюярх оепеяйнйю 
 subroutine SearchSolutionAccelerated(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,IKLKORX,IOPZ,IKLRXX,EnergyZeroZ,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::NMO,NNEL,NEpsIter,Kpoint,IndexM,IparametrDeterEnergy,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IndexCorrectSET,IndexSolutionsSET,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Exf,Henergy,EPS,ExfTemp,ExfSETT,HenergySET,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::DetM,ExfM,DetMA,ExfMA,ExfSolutions,ExfSET,DetSET,ProDet
   real(8),dimension(:,:)::AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin,EnergyZeroZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,IIII,klkl,IKLPERZnak,ISolution,IIIFKL,IKLCOMBEC
   real(8)::DeltaEnergy,DDD,Kcoff,RcoffGran,DMtemp,Rdet
   !real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,ExfTemp,DMtemp,DDD
   ! нопедекъел дереплхмюмр яхярелш
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ) 
   10685 FORMAT(2X,100(I4,1X))
   19685 FORMAT(2X,100(F20.12,1X))

   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   !WRITE(100,*) 'AZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (AZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO 
   !WRITE(100,*) 'BZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (BZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   E=E-MATMUL(AZ,BZ)
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! гюохяшбюел пегскэрюр пюяверю
   ! нопедекъел дереплхмюмр яхярелш
   Rdet=Determenant(E)

   WRITE(17,*) Exf,Rdet,RcoffChengH
   
  
  
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(4)=Rdet
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(5)=Rdet
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(1)=DetM(2)
			DetM(2)=DetM(3)
            DetM(3)=DetM(4)
            DetM(4)=DetM(5)
            DetM(5)=Rdet
        ENDIF
		
		! хмдейя тхйяхпсчыхи вхякн хрепюжхи б яксвюе йнцдю хмрепбюк менопедекем
		INDEXAD=INDEXAD+1
		! б мювюке бнгбпюыюеляъ й опедшдсыелс ьюцс сахпюел йнщттхжхемр онксвеммши б опедедсыел хмрепбюке 
		IF(INDEXAD.EQ.4) THEN
		   Henergy=Henergy*RcoffChengH
        ENDIF 
		! мювхмюел йнппейрхпнбюрэ ьюц 
		IF(INDEXAD.EQ.5) THEN
		   ! бнгбпюыюел йнщттхжхемр б хяундмне янярнъмхе
		   RcoffChengH=1.D0
           ProDet(1)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
	    IF(INDEXAD.EQ.6) THEN
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
		IF(INDEXAD.GT.6) THEN
           ProDet(1)=ProDet(2)
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		   ! нясыеярбкъел пюявер йнщттхжхемрю х хглемемхе ьюцю
		   call ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
		   ! опнбепъел оепеунд б пефхл слемэьемхъ ьюцю
		   IF(IKLCOMBEC.EQ.1) THEN
              ! бнгбпюыюеляъ мю ндмс рнвйс мюгюд
              Exf=ExfM(4)
			  ExfM(5)=ExfM(4)
              ExfM(4)=ExfM(3)
              ExfM(3)=ExfM(2)
              ExfM(2)=ExfM(1)
              DetM(5)=DetM(4)
              DetM(4)=DetM(3)
              DetM(3)=DetM(2)
              DetM(2)=DetM(1)   
		   ENDIF
		ENDIF
   ENDIF

   IF(IndexVilky.EQ.1) THEN	
      IF(INDEXDF.EQ.1) THEN 
         DetM(1)=DetM(2)
	     DetM(2)=DetM(3)
         DetM(3)=DetM(4)
         DetM(4)=DetM(5)
         DetM(5)=Rdet
      ENDIF 
	  IF(INDEXDF.EQ.2) THEN 
         ! опнбепъел йюйни хг пефхлнб бйкчвем
		 IF(IKLKORX.EQ.1) THEN
            ! йнппейрхпнбйю яопюбю
            ! опнбепъел врнаш гмювемхе яопюбю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(3).GT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(3).LT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
         IF(IKLKORX.EQ.2) THEN
            ! йнппейрхпнбйю якебю
            ! опнбепъел врнаш гмювемхе якебю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(4).GT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
               ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4  
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(4).LT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
		
      ENDIF    
   ENDIF 

   
   IF(IndexVilky.EQ.2) THEN
      IF(INDEXDF.EQ.1) THEN
           DetMA(4)=DetMA(1)
           DetMA(1)=DetMA(2)
           DetMA(2)=DetMA(3)
           DetMA(3)=Rdet
	     ELSE
	       ! йнпемэ оноюк б оепбши хрепбюк
           IF(IIKK.EQ.1) THEN
              DetMA(2)=DetMA(3)
              DetMA(3)=Rdet
	       ENDIF
           ! йнпемэ оноюк б брнпни хрепбюк
           IF(IIKK.EQ.2) THEN
              DetMA(1)=DetMA(3)
              DetMA(3)=Rdet
 	       ENDIF
      ENDIF
   ENDIF 
   		    
      
   ! опнбндхл онхяй йнпмеи мю бяел щмепцхрхвеяйнл хмрепбюке
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл оепеанп гмювемхи щмепцхх
      IF(NEpsIter.EQ.1) THEN
	        ExfM(4)=Exf
            Exf=Exf*(1.D0-Henergy)
            ExfM(5)=Exf  
          ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(4)*DetM(5)).LT.0.D0) THEN
	             IndexVilky=1 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи бме хмрепбюкю дкъ мнбнцн яксвюъ бме хмрепбюкю  
				 INDEXAD=0  
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи я слемэьеммшл ьюцнл
				 INDEXRShag=0 
				 ! нрйкчвюел йнмрпнкэ яйювйю
				 IIOXDPZ=0                          
			   ELSE
			     Exf=Exf*(1.D0-Henergy)
                 ExfM(1)=ExfM(2)
                 ExfM(2)=ExfM(3)
				 ExfM(3)=ExfM(4)
                 ExfM(4)=ExfM(5)
				 ExfM(5)=Exf  
           ENDIF      
	  ENDIF
   ENDIF
   
   ! сярюмнбкемю накюярэ оепелемш гмюйю 
   IF(IndexVilky.EQ.1) THEN
        ! опнбепъел оепелемю гмюйю ябъгюмю я бнгмхймнбемел пюгпшбю б дюммни накюярх хкх я мюкхвхел йнпмъ
        ! хмдейя сйюгшбюер, врн б дюммнл пефхле  опнцпюллю пюанрюер дкъ дюммнцн хмрепбюкю б оепбши пюг
		INDEXDF=INDEXDF+1
		! йкчв сйюгшбюер мю оепелемс гмюйю
		IKLPERZnak=0
        IF(INDEXDF.EQ.1) THEN
             Exf=Exf*(1.D0-Henergy/FLOAT(IndexCorrect+2))
             ExfM(1)=ExfM(2)
             ExfM(2)=ExfM(3)
		     ExfM(3)=ExfM(4)
             ExfM(4)=ExfM(5)
		     ExfM(5)=Exf  
			 ExfTemp=Exf 
			 ! тхйяхпсел IOPZ дкъ йнппейрхпнбйх йнпмеи
			 IOPZ=4  
			 ! тхйяхпсел вхякн рнвей йнрнпше ашкх свремш опх ясфемхх накюярх 
			 IKLRXX=0 
           ELSE   
             IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(2)).GT.DABS(DetM(3))) THEN
                 ! тхйяхпсел рнр тюйр, врн намюпсфем йнпемэ
				 IKLPERZnak=1 
				 ! намюпсфем хмрепбюк йнпмъ 
                 IndexCorrect=IndexCorrect+1
				 ! тхйяхпсел бнглнфмши йнпемэ б дюммни накюярх
                 EnergyZeroZ(NMO,IndexCorrect,1)=EnergyZeroZ(NMO,IndexCorrect,2)
				 EnergyZeroZ(NMO,IndexCorrect,2)=ExfM(3)  
                 ! хглемемхе ьюцю 
				 Henergy=Henergy/FLOAT(IndexCorrect+1)  
                 ! гюохяшбюел дюммше дкъ якедсчыецн йнпмъ дкъ мюунфдемхъ йнпмъ мендмнпндмнцн спюбмемхъ якебю 
                 ExfSET=ExfM
				 DetSET=DetM
				 HenergySET=Henergy
                 IndexCorrectSET=IndexCorrect
                 
				 ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			     ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			     ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			     IndexVilky=0 
                 INDEXDF=0
				 Exf=ExfM(5) 
				ELSE
				  ! сярюмюбкхбюел опхвхмю б пюгпшбе хкх мер
				  ! бшъямъел ясыеярбсер йнпемэ хкх хлеер леярн пюгпшб
				  ! опх щрнл, сярюмюбкхбюел рюйне онбедемхе ябъгюмн я йнкеаюмхълх б накюярх оепелемш гмюйю
				  ! нясыеярбкъел ясфемхе накюярх оепелемш гмюйю дкъ бмеяемхъ ъямнярх б онбедемхх нопедекхрекъ б накюярх мскъ
                     
				  IIIFKL=0
				  ! сярюмюбкхбюел опхвхмс нрясрярбхъ йнпмъ
				  ! опнбепъел бшонкмемхе оепбнцн сякнбхъ
			      IF(DABS(DetM(5)).GT.DABS(DetM(4))) THEN
				     ! брнпне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (яопюбю)     
				     IKLKORX=1
					 Exf=ExfM(3)-DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
				  ! опнбепъел бшонкмемхе брнпнцн сякнбхъ
				  IF(DABS(DetM(2)).GT.DABS(DetM(3))) THEN
				     ! оепбне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (якебю)     
				     IKLKORX=2
                     Exf=ExfM(4)+DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
                  
				  ! опнбепъел яйнкэйн пюг онксвемю мнбюъ рнвйю рн еярэ яйнкэйн пюг лш слемэьюкх накюярэ оепелемш гмюйю  
                  IF(IKLRXX.EQ.3) THEN
				     ! вхякн рнвей пюбмн рпел, опх щрнл ме опнхгнькн хглемемхи б онбедемхх нопедекхрекъ б акхгх рнвйх оепелемш гмюйю 
					 ! рюйхл напюгнл хлеер леярн пюгпшб 
					 IIIFKL=0 
				  ENDIF 
  

				  ! опнбепъел хлеер леярн пюгпшб хкх мер 
				  IF(IIIFKL.EQ.0) THEN
			         ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			         ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			         ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			         IndexVilky=0 
                     INDEXDF=0 
					 Exf=ExfM(5) 
				  ENDIF

             ENDIF
			
		ENDIF 
        
        ! опнбепъел сярюмнбкем тюйр мюкхвхъ йнпмъ хкх мер
		IF(IKLPERZnak.EQ.1) THEN
           ! опнбепъел бшонкмхкняэ кх сякнбхе
           ! мнлеп мюидемнцн йнпмъ он онпъдйс(хмрепбюкю б йнрнпнл мюундхряъ йнпемэ) яннрберярбсер цкюбмнлс йбюмрнбнлс вхякс
           ISolution=0
           IF(IndexCorrect.EQ.(NNEL-1)) THEN
              ISolution=1
		   ENDIF
           IF(IndexCorrect.EQ.NNEL) THEN
              ISolution=1
		   ENDIF
          		   
		   
		   IF(ISolution.EQ.1) THEN
              ! тхйяхпсел пефхл онхяйю йнпмъ б хмрепбюке
		      IndexVilky=2 
		      !WRITE(*,*) 'ZERO',IndexVilky
		      !WRITE(*,*) DetM(5),DetM(4),DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(5),ExfM(4),ExfM(2),ExfM(3)
		      !READ(*,*) 
		      ! гюохяшбюел цпюмхжш хмрепбюкю 
              ExfMA(2)=ExfM(4)
		      DetMA(2)=DetM(4)
		      ExfMA(3)=ExfM(3)
		      DetMA(3)=DetM(3)
		      !WRITE(*,*) DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(2),ExfM(3)
		      !READ(*,*)   
		      ! гюмскъел оепед онхяйнл йнпмъ
		      ExfMA(1)=0.D0
		      DetMA(1)=0.D0
		      ExfMA(4)=0.D0
              DetMA(4)=0.D0
		      ExfMA(5)=0.D0
		      DetMA(5)=0.D0
		      INDEXDF=0 
           ENDIF
        ENDIF 
	
		  		    
   ENDIF

   ! пефхл онхяйю йнпмъ б хмрепбюке
   IF(IndexVilky.EQ.2) THEN 
      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
      INDEXDF=INDEXDF+1
      ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	  IF(INDEXDF.EQ.1) THEN
           Exf=(DetMA(3)*ExfMA(2)-DetMA(2)*ExfMA(3))/(DetMA(3)-DetMA(2))
		   ExfMA(1)=ExfMA(2)
           ExfMA(2)=ExfMA(3)
		   ExfMA(3)=Exf
		 ELSE
           ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		   ! оепбши хмрепбюк
		   IF((DetMA(1)*DetMA(3)).LT.0.D0) THEN
               IIKK=1
		       Exf=(DetMA(3)*ExfMA(1)-DetMA(1)*ExfMA(3))/(DetMA(3)-DetMA(1))
			   ExfMA(2)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
		   ! брнпни хмрепбюк
           IF((DetMA(3)*DetMA(2)).LT.0.D0) THEN
		       IIKK=2
               Exf=(DetMA(2)*ExfMA(3)-DetMA(3)*ExfMA(2))/(DetMA(2)-DetMA(3))
		       ExfMA(1)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
      ENDIF		  		   
   

      ! опнбепъел мюидем йнпемэ хкх мер
	  ! опнбепъел бшонкмхкняэ кх сякнбхе
      ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
      ! мювхмюел опнбепърэ сякнбхе б яксвюе намюпсфемхъ хмрепбюкю цде мюундхряъ йнпемэ
      IF(INDEXDF.GT.1) THEN
         IF(DABS(DetMA(3)).LT.EPS) THEN
            ! йнпемэ мюидем 
			! тхйяхпсел мюидеммши йнпемэ
			IndexSolutions=IndexSolutions+1
			IndexSolutionsSET=IndexSolutions
			! гюохяшбюел мюидеммши йнпемэ
			ExfSolutions(IndexSolutions)=Exf
			! гюохяшбюел йнпемэ мюидеммши
			EnergyZeroZ(NMO,IndexCorrect,2)=ExfMA(3)  
			! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			IndexVilky=0 
            INDEXDF=0 
			Exf=ExfTemp 
			ExfSETT=ExfTemp 
			!write(6,*) 'kor',NNEL,IndexSolutions,ExfSolutions(IndexSolutions)
         ENDIF
      ENDIF 
	  
	  ! бшунд хг жхйкю он онхяйс йнпмеи опнхгнидер опх сякнбхх
	  IF(IndexSolutions.EQ.2) THEN
	     ! сякнбхе бшундю бшонкмемн
   	     IparametrDeterEnergy=2 
	  ENDIF 
	   
   ENDIF
    

  
   !WRITE(*,*) 'ITERASTION ZERO'
   !write(*,10685)  NEpsIter,IndexVilky,IndexCorrect,NNEL
   !write(*,19685)  DetM(1),DetM(3),DetM(2)
   !write(*,19685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,19685)  Exf
   !read(*,*)
   return
 end subroutine SearchSolutionAccelerated


 ! ондопнцпюллю нясыеярбкъчыюъ онхяй хмрепбюкю мсфмнцн йнпмъ ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх дкъ яксвюъ оепбни хрепюжхх
 ! NMO-мнлеп лнкейскъпмни нпахрюкх
 ! NNEL- цкюбмне йбюмрнбне вхякн
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи дереплхмюмрю
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! DetMA-люяяхб гмювемхи дереплхмюмрю дкъ онхяйю йнпмъ
 ! ExfMA-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх дкъ онхяйю йнпмъ
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! AZ,BZ,E-бяонлнцюрекэмше люяяхбш
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ (оняке нйнмвюмхъ опнцпюллш дюер мскхбне опхакхфемхе дкъ йнпмъ)
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! IndexVilky=2-пефхл онхяйю йнпмъ 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя сйюгшбючыхи вхякн йнпмеи мюидеммшу мю бяел хмрепбюке  
 ! IndexSolutions-хмдейя вхякю пеьемхе, менаундхлн мюирх пеьемхъ дкъ NNEL-1, NNEL, NNEL+1, рн еярэ дкъ бшундю 
 ! менаундхлн, врнаш IndexSolutions=3 
 ! ExfSolutions(IndexSolutions)-щмепцхх мюидеммшу йнпмеи  NNEL-1, NNEL, NNEL+1
 ! ExfTemp-бяонлнцюрекэмюъ оепелеммюъ
 ! IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET-дюммше менаундхлш дкъ пюяверю б яксвюе мюунфдемхъ йнпмъ якебю дкъ мендмнпндмни яхярелш
 ! IKLKORX,IOPZ,IKLRXX-хмдейяш йнппейрхпнбйх дкъ сярюмюбкемхъ йнпемэ хкх пюгпшб
 ! EnergyZeroZ(NumbreMO,NNEL+1,2)-люяяхб йнпмеи ндмнпндмнцн спюбмемхъ 
 ! INDEXAD-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу б яксвюе йнцдю хмрепбюк менопедекем 
 ! RcoffChengH-йнщттхжхемр хглемемхъ ьюцю
 ! ProDet(2)-люяяхб опнхгбндмшу
 ! INDEXRShag-пефхл хглемемхъ ьюцю
 ! IIHD-хмдейя слемэьемхъ ьюцю
 ! IIOXDPZ-йкчв сйюгшбючыхи врн лш оноюкх б пефхл слемэьемхъ ьюцю х йнмрпнкхпсел бнглнфмнярэ оепеяйнйю
 ! INDEXPeresc-хмдейя сйюгшбючыхи вюярх хрепюжхи б пефхле бнглнфмнцн оепеяйнйю
 ! DeltaEnerPers-хмрепбюк щмепцхи б накюярх оепеяйнйю
 ! DeltaFunPers-хмрепбюк тсмйжхи б накюярх оепеяйнйю 
 subroutine SearchSolutionAccelerated2(NMO,NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IndexCorrectSET,IndexSolutionsSET,HenergySET,ExfSETT,ExfSET,DetSET,IKLKORX,IOPZ,IKLRXX,EnergyZeroZ,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::NMO,NNEL,NEpsIter,Kpoint,IndexM,IparametrDeterEnergy,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IndexCorrectSET,IndexSolutionsSET,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Exf,Henergy,EPS,ExfTemp,ExfSETT,HenergySET,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::DetM,ExfM,DetMA,ExfMA,ExfSolutions,ExfSET,DetSET,ProDet
   real(8),dimension(:,:)::AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin,EnergyZeroZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,IIII,klkl,IKLPERZnak,ISolution,IIIFKL,IKLCOMBEC
   real(8)::DeltaEnergy,DDD,Kcoff,RcoffGran,DMtemp,Rdet
   !real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,ExfTemp,DMtemp,DDD
   ! нопедекъел дереплхмюмр яхярелш
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ) 
   10685 FORMAT(2X,100(I4,1X))
   19685 FORMAT(2X,100(F20.12,1X))

   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   !WRITE(100,*) 'AZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (AZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO 
   !WRITE(100,*) 'BZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (BZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   E=E-MATMUL(AZ,BZ)
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! гюохяшбюел пегскэрюр пюяверю
   ! нопедекъел дереплхмюмр яхярелш
   Rdet=Determenant(E)

   WRITE(17,*) Exf,Rdet,RcoffChengH
   
  
  
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(4)=Rdet
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(5)=Rdet
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(1)=DetM(2)
			DetM(2)=DetM(3)
            DetM(3)=DetM(4)
            DetM(4)=DetM(5)
            DetM(5)=Rdet
        ENDIF
		
		! хмдейя тхйяхпсчыхи вхякн хрепюжхи б яксвюе йнцдю хмрепбюк менопедекем
		INDEXAD=INDEXAD+1
		! б мювюке бнгбпюыюеляъ й опедшдсыелс ьюцс сахпюел йнщттхжхемр онксвеммши б опедедсыел хмрепбюке 
		IF(INDEXAD.EQ.4) THEN
		   Henergy=Henergy*RcoffChengH
        ENDIF 
		! мювхмюел йнппейрхпнбюрэ ьюц 
		IF(INDEXAD.EQ.5) THEN
		   ! бнгбпюыюел йнщттхжхемр б хяундмне янярнъмхе
		   RcoffChengH=1.D0
           ProDet(1)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
	    IF(INDEXAD.EQ.6) THEN
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
		IF(INDEXAD.GT.6) THEN
           ProDet(1)=ProDet(2)
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		   ! нясыеярбкъел пюявер йнщттхжхемрю х хглемемхе ьюцю
		   call ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
		   ! опнбепъел оепеунд б пефхл слемэьемхъ ьюцю
		   IF(IKLCOMBEC.EQ.1) THEN
              ! бнгбпюыюеляъ мю ндмс рнвйс мюгюд
              Exf=ExfM(4)
			  ExfM(5)=ExfM(4)
              ExfM(4)=ExfM(3)
              ExfM(3)=ExfM(2)
              ExfM(2)=ExfM(1)
              DetM(5)=DetM(4)
              DetM(4)=DetM(3)
              DetM(3)=DetM(2)
              DetM(2)=DetM(1)   
		   ENDIF
		ENDIF
   ENDIF

   IF(IndexVilky.EQ.1) THEN	
      IF(INDEXDF.EQ.1) THEN 
         DetM(1)=DetM(2)
	     DetM(2)=DetM(3)
         DetM(3)=DetM(4)
         DetM(4)=DetM(5)
         DetM(5)=Rdet
      ENDIF 
	  IF(INDEXDF.EQ.2) THEN 
         ! опнбепъел йюйни хг пефхлнб бйкчвем
		 IF(IKLKORX.EQ.1) THEN
            ! йнппейрхпнбйю яопюбю
            ! опнбепъел врнаш гмювемхе яопюбю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(3).GT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(3).LT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
         IF(IKLKORX.EQ.2) THEN
            ! йнппейрхпнбйю якебю
            ! опнбепъел врнаш гмювемхе якебю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(4).GT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
               ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4  
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(4).LT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
		
      ENDIF    
   ENDIF 

   
   IF(IndexVilky.EQ.2) THEN
      IF(INDEXDF.EQ.1) THEN
           DetMA(4)=DetMA(1)
           DetMA(1)=DetMA(2)
           DetMA(2)=DetMA(3)
           DetMA(3)=Rdet
	     ELSE
	       ! йнпемэ оноюк б оепбши хрепбюк
           IF(IIKK.EQ.1) THEN
              DetMA(2)=DetMA(3)
              DetMA(3)=Rdet
	       ENDIF
           ! йнпемэ оноюк б брнпни хрепбюк
           IF(IIKK.EQ.2) THEN
              DetMA(1)=DetMA(3)
              DetMA(3)=Rdet
 	       ENDIF
      ENDIF
   ENDIF 
   		    
      
   ! опнбндхл онхяй йнпмеи мю бяел щмепцхрхвеяйнл хмрепбюке
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл оепеанп гмювемхи щмепцхх
      IF(NEpsIter.EQ.1) THEN
	        ExfM(4)=Exf
            Exf=Exf*(1.D0-Henergy)
            ExfM(5)=Exf  
          ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(4)*DetM(5)).LT.0.D0) THEN
	             IndexVilky=1 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи бме хмрепбюкю дкъ мнбнцн яксвюъ бме хмрепбюкю  
				 INDEXAD=0  
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи я слемэьеммшл ьюцнл
				 INDEXRShag=0 
				 ! нрйкчвюел йнмрпнкэ яйювйю
				 IIOXDPZ=0                          
			   ELSE
			     Exf=Exf*(1.D0-Henergy)
                 ExfM(1)=ExfM(2)
                 ExfM(2)=ExfM(3)
				 ExfM(3)=ExfM(4)
                 ExfM(4)=ExfM(5)
				 ExfM(5)=Exf  
           ENDIF      
	  ENDIF
   ENDIF
   
   ! сярюмнбкемю накюярэ оепелемш гмюйю 
   IF(IndexVilky.EQ.1) THEN
        ! опнбепъел оепелемю гмюйю ябъгюмю я бнгмхймнбемел пюгпшбю б дюммни накюярх хкх я мюкхвхел йнпмъ
        ! хмдейя сйюгшбюер, врн б дюммнл пефхле  опнцпюллю пюанрюер дкъ дюммнцн хмрепбюкю б оепбши пюг
		INDEXDF=INDEXDF+1
		! йкчв сйюгшбюер мю оепелемс гмюйю
		IKLPERZnak=0
        IF(INDEXDF.EQ.1) THEN
             Exf=Exf*(1.D0-Henergy/FLOAT(IndexCorrect+2))
             ExfM(1)=ExfM(2)
             ExfM(2)=ExfM(3)
		     ExfM(3)=ExfM(4)
             ExfM(4)=ExfM(5)
		     ExfM(5)=Exf  
			 ExfTemp=Exf 
			 ! тхйяхпсел IOPZ дкъ йнппейрхпнбйх йнпмеи
			 IOPZ=4  
			 ! тхйяхпсел вхякн рнвей йнрнпше ашкх свремш опх ясфемхх накюярх 
			 IKLRXX=0 
           ELSE   
             IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(2)).GT.DABS(DetM(3))) THEN
                 ! тхйяхпсел рнр тюйр, врн намюпсфем йнпемэ
				 IKLPERZnak=1 
				 ! намюпсфем хмрепбюк йнпмъ 
                 IndexCorrect=IndexCorrect+1
				 ! тхйяхпсел бнглнфмши йнпемэ б дюммни накюярх
                 EnergyZeroZ(NMO,IndexCorrect,1)=EnergyZeroZ(NMO,IndexCorrect,2)
				 EnergyZeroZ(NMO,IndexCorrect,2)=ExfM(3)  
                 ! хглемемхе ьюцю 
				 Henergy=Henergy/FLOAT(IndexCorrect+1)  
                 ! гюохяшбюел дюммше дкъ якедсчыецн йнпмъ дкъ мюунфдемхъ йнпмъ мендмнпндмнцн спюбмемхъ якебю 
                 ExfSET=ExfM
				 DetSET=DetM
				 HenergySET=Henergy
                 IndexCorrectSET=IndexCorrect
                 
				 ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			     ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			     ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			     IndexVilky=0 
                 INDEXDF=0
				 Exf=ExfM(5) 
				ELSE
				  ! сярюмюбкхбюел опхвхмю б пюгпшбе хкх мер
				  ! бшъямъел ясыеярбсер йнпемэ хкх хлеер леярн пюгпшб
				  ! опх щрнл, сярюмюбкхбюел рюйне онбедемхе ябъгюмн я йнкеаюмхълх б накюярх оепелемш гмюйю
				  ! нясыеярбкъел ясфемхе накюярх оепелемш гмюйю дкъ бмеяемхъ ъямнярх б онбедемхх нопедекхрекъ б накюярх мскъ
                     
				  IIIFKL=0
				  ! сярюмюбкхбюел опхвхмс нрясрярбхъ йнпмъ
				  ! опнбепъел бшонкмемхе оепбнцн сякнбхъ
			      IF(DABS(DetM(5)).GT.DABS(DetM(4))) THEN
				     ! брнпне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (яопюбю)     
				     IKLKORX=1
					 Exf=ExfM(3)-DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
				  ! опнбепъел бшонкмемхе брнпнцн сякнбхъ
				  IF(DABS(DetM(2)).GT.DABS(DetM(3))) THEN
				     ! оепбне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (якебю)     
				     IKLKORX=2
                     Exf=ExfM(4)+DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
                  
				  ! опнбепъел яйнкэйн пюг онксвемю мнбюъ рнвйю рн еярэ яйнкэйн пюг лш слемэьюкх накюярэ оепелемш гмюйю  
                  IF(IKLRXX.EQ.3) THEN
				     ! вхякн рнвей пюбмн рпел, опх щрнл ме опнхгнькн хглемемхи б онбедемхх нопедекхрекъ б акхгх рнвйх оепелемш гмюйю 
					 ! рюйхл напюгнл хлеер леярн пюгпшб 
					 IIIFKL=0 
				  ENDIF 
  

				  ! опнбепъел хлеер леярн пюгпшб хкх мер 
				  IF(IIIFKL.EQ.0) THEN
			         ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			         ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			         ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			         IndexVilky=0 
                     INDEXDF=0 
					 Exf=ExfM(5) 
				  ENDIF

             ENDIF
			
		ENDIF 
        
        ! опнбепъел сярюмнбкем тюйр мюкхвхъ йнпмъ хкх мер
		IF(IKLPERZnak.EQ.1) THEN
           ! опнбепъел бшонкмхкняэ кх сякнбхе
           ! мнлеп мюидемнцн йнпмъ он онпъдйс (хмрепбюкю б йнрнпнл мюундхряъ йнпемэ) яннрберярбсер цкюбмнлс йбюмрнбнлс вхякс
           ISolution=0
           IF(IndexCorrect.EQ.NNEL) THEN
              ISolution=1
		   ENDIF
          		   
		   
		   IF(ISolution.EQ.1) THEN
              ! тхйяхпсел пефхл онхяйю йнпмъ б хмрепбюке
		      IndexVilky=2 
		      !WRITE(*,*) 'ZERO',IndexVilky
		      !WRITE(*,*) DetM(5),DetM(4),DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(5),ExfM(4),ExfM(2),ExfM(3)
		      !READ(*,*) 
		      ! гюохяшбюел цпюмхжш хмрепбюкю 
              ExfMA(2)=ExfM(4)
		      DetMA(2)=DetM(4)
		      ExfMA(3)=ExfM(3)
		      DetMA(3)=DetM(3)
		      !WRITE(*,*) DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(2),ExfM(3)
		      !READ(*,*)   
		      ! гюмскъел оепед онхяйнл йнпмъ
		      ExfMA(1)=0.D0
		      DetMA(1)=0.D0
		      ExfMA(4)=0.D0
              DetMA(4)=0.D0
		      ExfMA(5)=0.D0
		      DetMA(5)=0.D0
		      INDEXDF=0 
           ENDIF
        ENDIF 
	
		  		    
   ENDIF

   ! пефхл онхяйю йнпмъ б хмрепбюке
   IF(IndexVilky.EQ.2) THEN 
      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
      INDEXDF=INDEXDF+1
      ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	  IF(INDEXDF.EQ.1) THEN
           Exf=(DetMA(3)*ExfMA(2)-DetMA(2)*ExfMA(3))/(DetMA(3)-DetMA(2))
		   ExfMA(1)=ExfMA(2)
           ExfMA(2)=ExfMA(3)
		   ExfMA(3)=Exf
		 ELSE
           ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		   ! оепбши хмрепбюк
		   IF((DetMA(1)*DetMA(3)).LT.0.D0) THEN
               IIKK=1
		       Exf=(DetMA(3)*ExfMA(1)-DetMA(1)*ExfMA(3))/(DetMA(3)-DetMA(1))
			   ExfMA(2)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
		   ! брнпни хмрепбюк
           IF((DetMA(3)*DetMA(2)).LT.0.D0) THEN
		       IIKK=2
               Exf=(DetMA(2)*ExfMA(3)-DetMA(3)*ExfMA(2))/(DetMA(2)-DetMA(3))
		       ExfMA(1)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
      ENDIF		  		   
   

      ! опнбепъел мюидем йнпемэ хкх мер
	  ! опнбепъел бшонкмхкняэ кх сякнбхе
      ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
      ! мювхмюел опнбепърэ сякнбхе б яксвюе намюпсфемхъ хмрепбюкю цде мюундхряъ йнпемэ
      IF(INDEXDF.GT.1) THEN
         IF(DABS(DetMA(3)).LT.EPS) THEN
            ! йнпемэ мюидем 
			! тхйяхпсел мюидеммши йнпемэ
			IndexSolutions=IndexSolutions+1
			IndexSolutionsSET=IndexSolutions
			! гюохяшбюел мюидеммши йнпемэ
			ExfSolutions(IndexSolutions)=Exf
			! гюохяшбюел йнпемэ мюидеммши
			EnergyZeroZ(NMO,IndexCorrect,2)=ExfMA(3)  
			! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			IndexVilky=0 
            INDEXDF=0 
			Exf=ExfTemp 
			ExfSETT=ExfTemp 
			!write(6,*) 'kor',NNEL,IndexSolutions,ExfSolutions(IndexSolutions)
         ENDIF
      ENDIF 
	  
	  ! бшунд хг жхйкю он онхяйс йнпмеи опнхгнидер опх сякнбхх
	  IF(IndexSolutions.EQ.1) THEN
	     ! сякнбхе бшундю бшонкмемн
   	     IparametrDeterEnergy=2 
	  ENDIF 
	   
   ENDIF
    

  
   !WRITE(*,*) 'ITERASTION ZERO'
   !write(*,10685)  NEpsIter,IndexVilky,IndexCorrect,NNEL
   !write(*,19685)  DetM(1),DetM(3),DetM(2)
   !write(*,19685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,19685)  Exf
   !read(*,*)
   return
 end subroutine SearchSolutionAccelerated2



 ! ондопнцпюллю нясыеярбкъчыюъ онхяй хмрепбюкю мсфмнцн йнпмъ ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх дкъ яксвюъ онхяйю йнпмеи якебю
 ! NNEL- цкюбмне йбюмрнбне вхякн
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи дереплхмюмрю
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! DetMA-люяяхб гмювемхи дереплхмюмрю дкъ онхяйю йнпмъ
 ! ExfMA-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх дкъ онхяйю йнпмъ
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! AZ,BZ,E-бяонлнцюрекэмше люяяхбш
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ (оняке нйнмвюмхъ опнцпюллш дюер мскхбне опхакхфемхе дкъ йнпмъ)
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! IndexVilky=2-пефхл онхяйю йнпмъ 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя сйюгшбючыхи вхякн йнпмеи мюидеммшу мю бяел хмрепбюке  
 ! IndexSolutions-хмдейя вхякю пеьемхе, менаундхлн мюирх пеьемхъ дкъ NNEL-1, NNEL, NNEL+1, рн еярэ дкъ бшундю 
 ! менаундхлн, врнаш IndexSolutions=3 
 ! ExfSolutions(IndexSolutions)-щмепцхх мюидеммшу йнпмеи  NNEL-1, NNEL, NNEL+1
 ! ExfTemp-бяонлнцюрекэмюъ оепелеммюъ
 ! IKLKORX,IOPZ,IKLRXX-хмдейяш йнппейрхпнбйх дкъ сярюмюбкемхъ йнпемэ хкх пюгпшб
 ! INDEXAD-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу б яксвюе йнцдю хмрепбюк менопедекем 
 ! RcoffChengH-йнщттхжхемр хглемемхъ ьюцю
 ! ProDet(2)-люяяхб опнхгбндмшу
 ! INDEXRShag-пефхл хглемемхъ ьюцю
 ! IIHD-хмдейя слемэьемхъ ьюцю
 ! IIOXDPZ-йкчв сйюгшбючыхи врн лш оноюкх б пефхл слемэьемхъ ьюцю х йнмрпнкхпсел бнглнфмнярэ оепеяйнйю
 ! INDEXPeresc-хмдейя сйюгшбючыхи вюярх хрепюжхи б пефхле бнглнфмнцн оепеяйнйю
 ! DeltaEnerPers-хмрепбюк щмепцхи б накюярх оепеяйнйю
 ! DeltaFunPers-хмрепбюк тсмйжхи б накюярх оепеяйнйю 
 subroutine SearchSolutionAcceleratedLeft(NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::NNEL,NEpsIter,Kpoint,IndexM,IparametrDeterEnergy,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Exf,Henergy,EPS,ExfTemp,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::DetM,ExfM,DetMA,ExfMA,ExfSolutions,ProDet
   real(8),dimension(:,:)::AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,IIII,klkl,IKLPERZnak,ISolution,IIIFKL,IKLCOMBEC
   real(8)::DeltaEnergy,DDD,Kcoff,RcoffGran,DMtemp,Rdet
   !real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,ExfTemp,DMtemp,DDD
   ! нопедекъел дереплхмюмр яхярелш
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ) 
   10685 FORMAT(2X,100(I4,1X))
   19685 FORMAT(2X,100(F20.12,1X))

   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   !WRITE(100,*) 'AZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (AZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO 
   !WRITE(100,*) 'BZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (BZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   E=E-MATMUL(AZ,BZ)
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! гюохяшбюел пегскэрюр пюяверю
   ! нопедекъел дереплхмюмр яхярелш
   Rdet=Determenant(E)

   WRITE(17,*) Exf,Rdet,RcoffChengH
   
  
  
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(4)=Rdet
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(5)=Rdet
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(1)=DetM(2)
			DetM(2)=DetM(3)
            DetM(3)=DetM(4)
            DetM(4)=DetM(5)
            DetM(5)=Rdet
        ENDIF

		! хмдейя тхйяхпсчыхи вхякн хрепюжхи б яксвюе йнцдю хмрепбюк менопедекем
		INDEXAD=INDEXAD+1
		! б мювюке бнгбпюыюеляъ й опедшдсыелс ьюцс сахпюел йнщттхжхемр онксвеммши б опедедсыел хмрепбюке 
		IF(INDEXAD.EQ.4) THEN
		   Henergy=Henergy*RcoffChengH
        ENDIF 
		! мювхмюел йнппейрхпнбюрэ ьюц 
		IF(INDEXAD.EQ.5) THEN
		   ! бнгбпюыюел йнщттхжхемр б хяундмне янярнъмхе
		   RcoffChengH=1.D0
           ProDet(1)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
	    IF(INDEXAD.EQ.6) THEN
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
		IF(INDEXAD.GT.6) THEN
           ProDet(1)=ProDet(2)
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		   ! нясыеярбкъел пюявер йнщттхжхемрю х хглемемхе ьюцю
		   call ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)

		   ! опнбепъел оепеунд б пефхл слемэьемхъ ьюцю
		   IF(IKLCOMBEC.EQ.1) THEN
              ! бнгбпюыюеляъ мю ндмс рнвйс мюгюд
              Exf=ExfM(4)
			  ExfM(5)=ExfM(4)
              ExfM(4)=ExfM(3)
              ExfM(3)=ExfM(2)
              ExfM(2)=ExfM(1)
              DetM(5)=DetM(4)
              DetM(4)=DetM(3)
              DetM(3)=DetM(2)
              DetM(2)=DetM(1)   
		   ENDIF
		ENDIF
   ENDIF

   IF(IndexVilky.EQ.1) THEN	
      IF(INDEXDF.EQ.1) THEN 
         DetM(1)=DetM(2)
	     DetM(2)=DetM(3)
         DetM(3)=DetM(4)
         DetM(4)=DetM(5)
         DetM(5)=Rdet
      ENDIF  
	  IF(INDEXDF.EQ.2) THEN 
         ! опнбепъел йюйни хг пефхлнб бйкчвем
		 IF(IKLKORX.EQ.1) THEN
            ! йнппейрхпнбйю яопюбю
            ! опнбепъел врнаш гмювемхе яопюбю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(3).GT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(3).LT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
         IF(IKLKORX.EQ.2) THEN
            ! йнппейрхпнбйю якебю
            ! опнбепъел врнаш гмювемхе якебю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(4).GT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(4).LT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
		
      ENDIF    
   ENDIF 

   
   IF(IndexVilky.EQ.2) THEN
      IF(INDEXDF.EQ.1) THEN
           DetMA(4)=DetMA(1)
           DetMA(1)=DetMA(2)
           DetMA(2)=DetMA(3)
           DetMA(3)=Rdet
	     ELSE
	       ! йнпемэ оноюк б оепбши хрепбюк
           IF(IIKK.EQ.1) THEN
              DetMA(2)=DetMA(3)
              DetMA(3)=Rdet
	       ENDIF
           ! йнпемэ оноюк б брнпни хрепбюк
           IF(IIKK.EQ.2) THEN
              DetMA(1)=DetMA(3)
              DetMA(3)=Rdet
 	       ENDIF
      ENDIF
   ENDIF 
   		    
      
   ! опнбндхл онхяй йнпмеи мю бяел щмепцхрхвеяйнл хмрепбюке
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл оепеанп гмювемхи щмепцхх
      IF(NEpsIter.EQ.1) THEN
	        ExfM(4)=Exf
            Exf=Exf*(1.D0-Henergy)
            ExfM(5)=Exf  
          ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(4)*DetM(5)).LT.0.D0) THEN
	             IndexVilky=1 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи бме хмрепбюкю дкъ мнбнцн яксвюъ бме хмрепбюкю  
				 INDEXAD=0 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи я слемэьеммшл ьюцнл
				 INDEXRShag=0
				 ! нрйкчвюел йнмрпнкэ яйювйю
				 IIOXDPZ=0                            
			   ELSE
			     Exf=Exf*(1.D0-Henergy)
                 ExfM(1)=ExfM(2)
                 ExfM(2)=ExfM(3)
				 ExfM(3)=ExfM(4)
                 ExfM(4)=ExfM(5)
				 ExfM(5)=Exf  
           ENDIF      
	  ENDIF
   ENDIF
   
   ! сярюмнбкемю накюярэ оепелемш гмюйю 
   IF(IndexVilky.EQ.1) THEN
        ! опнбепъел оепелемю гмюйю ябъгюмю я бнгмхймнбемел пюгпшбю б дюммни накюярх хкх я мюкхвхел йнпмъ
        ! хмдейя сйюгшбюер, врн б дюммнл пефхле  опнцпюллю пюанрюер дкъ дюммнцн хмрепбюкю б оепбши пюг
		INDEXDF=INDEXDF+1
		! йкчв сйюгшбюер мю оепелемс гмюйю
		IKLPERZnak=0
        IF(INDEXDF.EQ.1) THEN
             Exf=Exf*(1.D0-Henergy/FLOAT(IndexCorrect+2))
             ExfM(1)=ExfM(2)
             ExfM(2)=ExfM(3)
		     ExfM(3)=ExfM(4)
             ExfM(4)=ExfM(5)
		     ExfM(5)=Exf  
			 ExfTemp=Exf 
			 ! тхйяхпсел IOPZ дкъ йнппейрхпнбйх йнпмеи
			 IOPZ=4  
			 ! тхйяхпсел вхякн рнвей йнрнпше ашкх свремш опх ясфемхх накюярх 
			 IKLRXX=0 
           ELSE   
             IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(2)).GT.DABS(DetM(3))) THEN
                 ! тхйяхпсел рнр тюйр, врн намюпсфем йнпемэ
				 IKLPERZnak=1 
				 ! намюпсфем хмрепбюк йнпмъ 
                 IndexCorrect=IndexCorrect+1
                 ! хглемемхе ьюцю 
				 Henergy=Henergy/FLOAT(IndexCorrect+1)  
                 ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			     ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			     ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			     IndexVilky=0 
                 INDEXDF=0
				 Exf=ExfM(5) 
				ELSE
                  ! сярюмюбкхбюел опхвхмю б пюгпшбе хкх мер
				  ! бшъямъел ясыеярбсер йнпемэ хкх хлеер леярн пюгпшб
				  ! опх щрнл, сярюмюбкхбюел рюйне онбедемхе ябъгюмн я йнкеаюмхълх б накюярх оепелемш гмюйю
				  ! нясыеярбкъел ясфемхе накюярх оепелемш гмюйю дкъ бмеяемхъ ъямнярх б онбедемхх нопедекхрекъ б накюярх мскъ
                   
				  IIIFKL=0
				  ! сярюмюбкхбюел опхвхмс нрясрярбхъ йнпмъ
				  ! опнбепъел бшонкмемхе оепбнцн сякнбхъ
			      IF(DABS(DetM(5)).GT.DABS(DetM(4))) THEN
				     ! брнпне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (яопюбю)     
				     IKLKORX=1
                     Exf=ExfM(3)-DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
				  ! опнбепъел бшонкмемхе брнпнцн сякнбхъ
				  IF(DABS(DetM(2)).GT.DABS(DetM(3))) THEN
				     ! оепбне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (якебю)     
				     IKLKORX=2
                     Exf=ExfM(4)+DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
                  
				  ! опнбепъел яйнкэйн пюг онксвемю мнбюъ рнвйю рн еярэ яйнкэйн пюг лш слемэьюкх накюярэ оепелемш гмюйю  
                  IF(IKLRXX.EQ.3) THEN
				     ! вхякн рнвей пюбмн рпел, опх щрнл ме опнхгнькн хглемемхи б онбедемхх нопедекхрекъ б акхгх рнвйх оепелемш гмюйю 
					 ! рюйхл напюгнл хлеер леярн пюгпшб 
					 IIIFKL=0 
				  ENDIF 
			
				  ! опнбепъел хлеер леярн пюгпшб хкх мер 
				  IF(IIIFKL.EQ.0) THEN
			         ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			         ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			         ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			         IndexVilky=0 
                     INDEXDF=0 
					 Exf=ExfM(5) 
				  ENDIF
			 ENDIF
		
		ENDIF 
        
        ! опнбепъел сярюмнбкем тюйр мюкхвхъ йнпмъ хкх мер
		IF(IKLPERZnak.EQ.1) THEN
           ! опнбепъел бшонкмхкняэ кх сякнбхе
           ! мнлеп мюидемнцн йнпмъ он онпъдйс(хмрепбюкю б йнрнпнл мюундхряъ йнпемэ) яннрберярбсер цкюбмнлс йбюмрнбнлс вхякс
           ISolution=0
           IF(IndexCorrect.EQ.(NNEL+1)) THEN
              ISolution=1
		   ENDIF
		   
		   
		   IF(ISolution.EQ.1) THEN
              ! тхйяхпсел пефхл онхяйю йнпмъ б хмрепбюке
		      IndexVilky=2 
		      !WRITE(*,*) 'ZERO',IndexVilky
		      !WRITE(*,*) DetM(5),DetM(4),DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(5),ExfM(4),ExfM(2),ExfM(3)
		      !READ(*,*) 
		      ! гюохяшбюел цпюмхжш хмрепбюкю 
              ExfMA(2)=ExfM(4)
		      DetMA(2)=DetM(4)
		      ExfMA(3)=ExfM(3)
		      DetMA(3)=DetM(3)
		      !WRITE(*,*) DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(2),ExfM(3)
		      !READ(*,*)   
		      ! гюмскъел оепед онхяйнл йнпмъ
		      ExfMA(1)=0.D0
		      DetMA(1)=0.D0
		      ExfMA(4)=0.D0
              DetMA(4)=0.D0
		      ExfMA(5)=0.D0
		      DetMA(5)=0.D0
		      INDEXDF=0 
           ENDIF
        ENDIF 
	
		  		    
   ENDIF

   ! пефхл онхяйю йнпмъ б хмрепбюке
   IF(IndexVilky.EQ.2) THEN 
      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
      INDEXDF=INDEXDF+1
      ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	  IF(INDEXDF.EQ.1) THEN
           Exf=(DetMA(3)*ExfMA(2)-DetMA(2)*ExfMA(3))/(DetMA(3)-DetMA(2))
		   ExfMA(1)=ExfMA(2)
           ExfMA(2)=ExfMA(3)
		   ExfMA(3)=Exf
		 ELSE
           ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		   ! оепбши хмрепбюк
		   IF((DetMA(1)*DetMA(3)).LT.0.D0) THEN
               IIKK=1
		       Exf=(DetMA(3)*ExfMA(1)-DetMA(1)*ExfMA(3))/(DetMA(3)-DetMA(1))
			   ExfMA(2)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
		   ! брнпни хмрепбюк
           IF((DetMA(3)*DetMA(2)).LT.0.D0) THEN
		       IIKK=2
               Exf=(DetMA(2)*ExfMA(3)-DetMA(3)*ExfMA(2))/(DetMA(2)-DetMA(3))
		       ExfMA(1)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
      ENDIF		  		   
   

      ! опнбепъел мюидем йнпемэ хкх мер
	  ! опнбепъел бшонкмхкняэ кх сякнбхе
      ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
      ! мювхмюел опнбепърэ сякнбхе б яксвюе намюпсфемхъ хмрепбюкю цде мюундхряъ йнпемэ
      IF(INDEXDF.GT.1) THEN
         IF(DABS(DetMA(3)).LT.EPS) THEN
            ! йнпемэ мюидем 
			! тхйяхпсел мюидеммши йнпемэ
			IndexSolutions=IndexSolutions+1
			! гюохяшбюел мюидеммши йнпемэ
			ExfSolutions(IndexSolutions)=Exf
			! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			IndexVilky=0 
            INDEXDF=0 
			Exf=ExfTemp 
			!write(6,*) 'kor',NNEL,IndexSolutions,ExfSolutions(IndexSolutions)
         ENDIF
      ENDIF 
	  
	  ! бшунд хг жхйкю он онхяйс йнпмеи опнхгнидер опх сякнбхх
	  IF(IndexSolutions.EQ.3) THEN
	     ! сякнбхе бшундю бшонкмемн
   	     IparametrDeterEnergy=2 
	  ENDIF 
	   
   ENDIF
    

  
   !WRITE(*,*) 'ITERASTION ZERO'
   !write(*,10685)  NEpsIter,IndexVilky,IndexCorrect,NNEL
   !write(*,19685)  DetM(1),DetM(3),DetM(2)
   !write(*,19685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,19685)  Exf
   !read(*,*)
   return
 end subroutine SearchSolutionAcceleratedLeft

 ! ондопнцпюллю нясыеярбкъчыюъ онхяй хмрепбюкю мсфмнцн йнпмъ ндмнщкейрпнммни щмепцхх уюпрпх-тнйю лнкейскъпмни нпахрюкх дкъ яксвюъ онхяйю йнпмеи якебю
 ! NNEL- цкюбмне йбюмрнбне вхякн
 ! NEpsIter-мнлеп хрепюжхх дкъ онксвемхъ ндмнщкейрпнммни щмепцхх
 ! Kpoint-рнвйю яьхбйх
 ! IndexM-пюглеп люяяхбю
 ! Vout(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)- люяяхб оюпюлерпнб опнцнмйх "бмсрпэ"
 ! DetM-люяяхб гмювемхи дереплхмюмрю
 ! ExfM-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх
 ! DetMA-люяяхб гмювемхи дереплхмюмрю дкъ онхяйю йнпмъ
 ! ExfMA-люяяхб гмювемхи ндмнщкейрпнммни тсмйжхх дкъ онхяйю йнпмъ
 ! Exf-ндмнщкейрпнммюъ щмепцхъ
 ! DeltaEnergy-пюгмхжю лефдс ндмнщкейрпнммшлх тсмйжхълх
 ! AZ,BZ,E-бяонлнцюрекэмше люяяхбш
 ! Henergy-ЬЮЦ ОН НДМНЩКЕЙРПНММНИ ЩМЕПЦХХ (оняке нйнмвюмхъ опнцпюллш дюер мскхбне опхакхфемхе дкъ йнпмъ)
 ! IparametrDeterEnergy-оюпюлерп бшундю хг жхйкю он срнвмемхч ндмнщкейрпнммни щмепцхх
 ! EPS-рнвмнярэ
 ! IndexVilky-хмдейя сйюгшбючыхи мю рн врн нопедекем хмрепбюк б йнрнпнл мюундхряъ йнпемэ
 ! IndexVilky=0-хмрепбюк ме нопедекем
 ! IndexVilky=1-хмрепбюк нопедекем 
 ! IndexVilky=2-пефхл онхяйю йнпмъ 
 ! INDEXDF-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу оняке оноюдюмхъ б бхкйс
 ! IndexCorrect-хмдейя сйюгшбючыхи вхякн йнпмеи мюидеммшу мю бяел хмрепбюке  
 ! IndexSolutions-хмдейя вхякю пеьемхе, менаундхлн мюирх пеьемхъ дкъ NNEL-1, NNEL, NNEL+1, рн еярэ дкъ бшундю 
 ! менаундхлн, врнаш IndexSolutions=3 
 ! ExfSolutions(IndexSolutions)-щмепцхх мюидеммшу йнпмеи  NNEL-1, NNEL, NNEL+1
 ! ExfTemp-бяонлнцюрекэмюъ оепелеммюъ
 ! IKLKORX,IOPZ,IKLRXX-хмдейяш йнппейрхпнбйх дкъ сярюмюбкемхъ йнпемэ хкх пюгпшб
 ! INDEXAD-хмдейя сйюгшбючыхи вхякн хрепюжхи ядекюммшу б яксвюе йнцдю хмрепбюк менопедекем 
 ! RcoffChengH-йнщттхжхемр хглемемхъ ьюцю
 ! ProDet(2)-люяяхб опнхгбндмшу
 ! INDEXRShag-пефхл хглемемхъ ьюцю
 ! INDEXGRNYT-йкчв бшундю
 ! INDEXSHAA-сйюгшбюер люйяхлюкэмне вхякн ьюцнб 
 ! IIHD-хмдейя слемэьемхъ ьюцю
 ! IIOXDPZ-йкчв сйюгшбючыхи врн лш оноюкх б пефхл слемэьемхъ ьюцю х йнмрпнкхпсел бнглнфмнярэ оепеяйнйю
 ! INDEXPeresc-хмдейя сйюгшбючыхи вюярх хрепюжхи б пефхле бнглнфмнцн оепеяйнйю
 ! DeltaEnerPers-хмрепбюк щмепцхи б накюярх оепеяйнйю
 ! DeltaFunPers-хмрепбюк тсмйжхи б накюярх оепеяйнйю 
 subroutine SearchSolutionAcceleratedLeft2(NNEL,NEpsIter,Kpoint,IndexM,Vout,Vin,DetM,ExfM,DetMA,ExfMA,Exf,AZ,BZ,E,Henergy,IparametrDeterEnergy,EPS,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,ExfSolutions,ExfTemp,IKLKORX,IOPZ,IKLRXX,INDEXAD,RcoffChengH,ProDet,INDEXRShag,INDEXGRNYT,INDEXSHAA,IIHD,IIOXDPZ,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::NNEL,NEpsIter,Kpoint,IndexM,IparametrDeterEnergy,IndexVilky,INDEXDF,IIKK,IndexCorrect,IndexSolutions,IKLKORX,IOPZ,IKLRXX,INDEXAD,INDEXRShag,INDEXGRNYT,INDEXSHAA,IIHD,IIOXDPZ,INDEXPeresc
   real(8)::Exf,Henergy,EPS,ExfTemp,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::DetM,ExfM,DetMA,ExfMA,ExfSolutions,ProDet
   real(8),dimension(:,:)::AZ,BZ,E
   real(8),dimension(:,:,:)::Vout,Vin
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   integer::IIDD,IBBB,IIII,klkl,IKLPERZnak,ISolution,IIIFKL,IKLCOMBEC
   real(8)::DeltaEnergy,DDD,Kcoff,RcoffGran,DMtemp,Rdet
   !real(8)::DeltaEnergy,RnormConstas,Kcoff,RcoffGran,ExfTemp,DMtemp,DDD
   ! нопедекъел дереплхмюмр яхярелш
   call ReadWriteMassiv2(2,Kpoint,IndexM,Vout,AZ) 
   call ReadWriteMassiv2(2,Kpoint+1,IndexM,Vin,BZ) 
   10685 FORMAT(2X,100(I4,1X))
   19685 FORMAT(2X,100(F20.12,1X))

   ! тнплхпсел едхмхвмсч люрпхжс
   E=0.D0
   DO IIDD=1,IndexM
      E(IIDD,IIDD)=1.D0
   ENDDO
   
   ! тнплхпсел люрпхжс яхярелш
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   !WRITE(100,*) 'AZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (AZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO 
   !WRITE(100,*) 'BZ'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (BZ(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   E=E-MATMUL(AZ,BZ)
   !WRITE(100,*) 'E'
   !DO IIDD=1,IndexM
   !   WRITE(100,19345 ) (E(IIDD,IBBB),IBBB=1,IndexM)
   !ENDDO
   ! гюохяшбюел пегскэрюр пюяверю
   ! нопедекъел дереплхмюмр яхярелш
   Rdet=Determenant(E)

   WRITE(17,*) Exf,Rdet,RcoffChengH
   
  
  
   
   ! гюохяшбюел пегскэрюр пюяверю
   IF(IndexVilky.EQ.0) THEN
        IF(NEpsIter.EQ.1) THEN
           DetM(4)=Rdet
		ENDIF
        IF(NEpsIter.EQ.2) THEN
           DetM(5)=Rdet
		ENDIF
        IF(NEpsIter.NE.1.AND.NEpsIter.NE.2) THEN
            DetM(1)=DetM(2)
			DetM(2)=DetM(3)
            DetM(3)=DetM(4)
            DetM(4)=DetM(5)
            DetM(5)=Rdet
        ENDIF

		! хмдейя тхйяхпсчыхи вхякн хрепюжхи б яксвюе йнцдю хмрепбюк менопедекем
		INDEXAD=INDEXAD+1
		! б мювюке бнгбпюыюеляъ й опедшдсыелс ьюцс сахпюел йнщттхжхемр онксвеммши б опедедсыел хмрепбюке 
		IF(INDEXAD.EQ.4) THEN
		   Henergy=Henergy*RcoffChengH
        ENDIF 
		! мювхмюел йнппейрхпнбюрэ ьюц 
		IF(INDEXAD.EQ.5) THEN
		   ! бнгбпюыюел йнщттхжхемр б хяундмне янярнъмхе
		   RcoffChengH=1.D0
           ProDet(1)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
	    IF(INDEXAD.EQ.6) THEN
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		ENDIF
		IF(INDEXAD.GT.6) THEN
           ProDet(1)=ProDet(2)
           ProDet(2)=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5))/(ExfM(5)*(1.D0-ExfM(4)/ExfM(5))))
		   ! нясыеярбкъел пюявер йнщттхжхемрю х хглемемхе ьюцю
		   call ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
		   ! опнбепъел оепеунд б пефхл слемэьемхъ ьюцю
		   IF(IKLCOMBEC.EQ.1) THEN
              ! бнгбпюыюеляъ мю ндмс рнвйс мюгюд
              Exf=ExfM(4)
			  ExfM(5)=ExfM(4)
              ExfM(4)=ExfM(3)
              ExfM(3)=ExfM(2)
              ExfM(2)=ExfM(1)
              DetM(5)=DetM(4)
              DetM(4)=DetM(3)
              DetM(3)=DetM(2)
              DetM(2)=DetM(1)   
		   ENDIF
		ENDIF

		! бшунд асдер оняке 150 ьюцнб б дюммнл пефхле
		INDEXSHAA=INDEXSHAA+1
		IF(INDEXSHAA.GT.150) THEN
		   ! нясыеярбхрэ бшунд
		   INDEXGRNYT=1
		   ! тхйяхпсел мюидеммши йнпемэ
		   IndexSolutions=IndexSolutions+1
		   ! гюохяшбюел мюидеммши йнпемэ
		   ExfSolutions(IndexSolutions)=Exf  
		   ! сякнбхе бшундю бшонкмемн
   	       IparametrDeterEnergy=2 
		ENDIF 
   ENDIF

   IF(IndexVilky.EQ.1) THEN	
      IF(INDEXDF.EQ.1) THEN 
         DetM(1)=DetM(2)
	     DetM(2)=DetM(3)
         DetM(3)=DetM(4)
         DetM(4)=DetM(5)
         DetM(5)=Rdet
      ENDIF  
	  IF(INDEXDF.EQ.2) THEN 
         ! опнбепъел йюйни хг пефхлнб бйкчвем
		 IF(IKLKORX.EQ.1) THEN
            ! йнппейрхпнбйю яопюбю
            ! опнбепъел врнаш гмювемхе яопюбю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(3).GT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(3).LT.0.D0) THEN
               DetM(1)=DetM(2)
	           DetM(2)=DetM(3)
			   DetM(3)=Rdet
               ExfM(1)=ExfM(2)
               ExfM(2)=ExfM(3)
			   ExfM(3)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
         IF(IKLKORX.EQ.2) THEN
            ! йнппейрхпнбйю якебю
            ! опнбепъел врнаш гмювемхе якебю ашкн рнцн фе гмюйю
			IIIFKL=0
			IF(Rdet.GT.0.D0.AND.DetM(4).GT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
            IF(Rdet.LT.0.D0.AND.DetM(4).LT.0.D0) THEN
               DetM(5)=DetM(4)
	           DetM(4)=Rdet
			   ExfM(5)=ExfM(4)
               ExfM(4)=Exf
			   IIIFKL=1
			   ! сбекхвхбюел вхякн рнвей онксвеммшу б акхгх рнвйх хглемемхъ гмюйю
			   IKLRXX=IKLRXX+1
			   ! рнвйю нопедекемю бнгбпюыюел хяундмне гмювемхе йнппейрпхпсчыецн лмнфхрекъ
			   IOPZ=4 
			ENDIF
			IF(IIIFKL.EQ.0) THEN
               IOPZ=IOPZ*2
			ENDIF
			! мнлеп хрепюжхх бнгбпюыюел мю 1 мюгюд
			INDEXDF=INDEXDF-1
		 ENDIF
		
      ENDIF    
   ENDIF 

   
   IF(IndexVilky.EQ.2) THEN
      IF(INDEXDF.EQ.1) THEN
           DetMA(4)=DetMA(1)
           DetMA(1)=DetMA(2)
           DetMA(2)=DetMA(3)
           DetMA(3)=Rdet
	     ELSE
	       ! йнпемэ оноюк б оепбши хрепбюк
           IF(IIKK.EQ.1) THEN
              DetMA(2)=DetMA(3)
              DetMA(3)=Rdet
	       ENDIF
           ! йнпемэ оноюк б брнпни хрепбюк
           IF(IIKK.EQ.2) THEN
              DetMA(1)=DetMA(3)
              DetMA(3)=Rdet
 	       ENDIF
      ENDIF
   ENDIF 
   		    
      
   ! опнбндхл онхяй йнпмеи мю бяел щмепцхрхвеяйнл хмрепбюке
   IF(IndexVilky.EQ.0) THEN
      ! опнбндхл оепеанп гмювемхи щмепцхх
      IF(NEpsIter.EQ.1) THEN
	        ExfM(4)=Exf
            Exf=Exf*(1.D0-Henergy)
            ExfM(5)=Exf  
          ELSE
	       ! опнбепъел оноюдюер ндмнщкейрпнммюъ щмепцхъ б бхкйс хкх мер (нопедекем кх хмрепбюк б йнрнпнл мюундхряъ йнпемэ)
	       IF((DetM(4)*DetM(5)).LT.0.D0) THEN
	             IndexVilky=1 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи бме хмрепбюкю дкъ мнбнцн яксвюъ бме хмрепбюкю  
				 INDEXAD=0 
				 ! оняйнкэйс хмрепбюк нопедекем гюмскъел вхякн хрепюжхи я слемэьеммшл ьюцнл
				 INDEXRShag=0  
				 ! нрйкчвюел йнмрпнкэ яйювйю
				 IIOXDPZ=0                          
			   ELSE
			     Exf=Exf*(1.D0-Henergy)
                 ExfM(1)=ExfM(2)
                 ExfM(2)=ExfM(3)
				 ExfM(3)=ExfM(4)
                 ExfM(4)=ExfM(5)
				 ExfM(5)=Exf  
           ENDIF      
	  ENDIF
   ENDIF
   
   ! сярюмнбкемю накюярэ оепелемш гмюйю 
   IF(IndexVilky.EQ.1) THEN
        ! опнбепъел оепелемю гмюйю ябъгюмю я бнгмхймнбемел пюгпшбю б дюммни накюярх хкх я мюкхвхел йнпмъ
        ! хмдейя сйюгшбюер, врн б дюммнл пефхле  опнцпюллю пюанрюер дкъ дюммнцн хмрепбюкю б оепбши пюг
		INDEXDF=INDEXDF+1
		! йкчв сйюгшбюер мю оепелемс гмюйю
		IKLPERZnak=0
        IF(INDEXDF.EQ.1) THEN
             Exf=Exf*(1.D0-Henergy/FLOAT(IndexCorrect+2))
             ExfM(1)=ExfM(2)
             ExfM(2)=ExfM(3)
		     ExfM(3)=ExfM(4)
             ExfM(4)=ExfM(5)
		     ExfM(5)=Exf  
			 ExfTemp=Exf 
			 ! тхйяхпсел IOPZ дкъ йнппейрхпнбйх йнпмеи
			 IOPZ=4  
			 ! тхйяхпсел вхякн рнвей йнрнпше ашкх свремш опх ясфемхх накюярх 
			 IKLRXX=0 
           ELSE   
             IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(2)).GT.DABS(DetM(3))) THEN
                 ! тхйяхпсел рнр тюйр, врн намюпсфем йнпемэ
				 IKLPERZnak=1 
				 ! намюпсфем хмрепбюк йнпмъ 
                 IndexCorrect=IndexCorrect+1
                 ! хглемемхе ьюцю 
				 Henergy=Henergy/FLOAT(IndexCorrect+1)  
                 ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			     ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			     ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			     IndexVilky=0 
                 INDEXDF=0
				 Exf=ExfM(5) 
				ELSE
                  ! сярюмюбкхбюел опхвхмю б пюгпшбе хкх мер
				  ! бшъямъел ясыеярбсер йнпемэ хкх хлеер леярн пюгпшб
				  ! опх щрнл, сярюмюбкхбюел рюйне онбедемхе ябъгюмн я йнкеаюмхълх б накюярх оепелемш гмюйю
				  ! нясыеярбкъел ясфемхе накюярх оепелемш гмюйю дкъ бмеяемхъ ъямнярх б онбедемхх нопедекхрекъ б накюярх мскъ
                   
				  IIIFKL=0
				  ! сярюмюбкхбюел опхвхмс нрясрярбхъ йнпмъ
				  ! опнбепъел бшонкмемхе оепбнцн сякнбхъ
			      IF(DABS(DetM(5)).GT.DABS(DetM(4))) THEN
				     ! брнпне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (яопюбю)     
				     IKLKORX=1
                     Exf=ExfM(3)-DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
				  ! опнбепъел бшонкмемхе брнпнцн сякнбхъ
				  IF(DABS(DetM(2)).GT.DABS(DetM(3))) THEN
				     ! оепбне сякнбхе ме бшонкмъеряъ 
					 ! опхакхфюер гмювемхе й йнпмъч (якебю)     
				     IKLKORX=2
                     Exf=ExfM(4)+DABS(ExfM(3)-ExfM(4))/DBLE(IOPZ)
				     IIIFKL=1
				  ENDIF 
                  
				  ! опнбепъел яйнкэйн пюг онксвемю мнбюъ рнвйю рн еярэ яйнкэйн пюг лш слемэьюкх накюярэ оепелемш гмюйю  
                  IF(IKLRXX.EQ.3) THEN
				     ! вхякн рнвей пюбмн рпел, опх щрнл ме опнхгнькн хглемемхи б онбедемхх нопедекхрекъ б акхгх рнвйх оепелемш гмюйю 
					 ! рюйхл напюгнл хлеер леярн пюгпшб 
					 IIIFKL=0 
				  ENDIF 
			
				  ! опнбепъел хлеер леярн пюгпшб хкх мер 
				  IF(IIIFKL.EQ.0) THEN
			         ! дюммюъ оепелемю гмюйнб ябъгюмю я мюкхвхел пюгпшбю
			         ! кхан йнпемэ мюидем опх щрнл бяе пюбмн мсфмн оепеирх б пефхл дюкэмеиьецн онхяйю
			         ! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			         IndexVilky=0 
                     INDEXDF=0 
					 Exf=ExfM(5) 
				  ENDIF
			 ENDIF
		
		ENDIF 
        
        ! опнбепъел сярюмнбкем тюйр мюкхвхъ йнпмъ хкх мер
		IF(IKLPERZnak.EQ.1) THEN
           ! опнбепъел бшонкмхкняэ кх сякнбхе
           ! мнлеп мюидемнцн йнпмъ он онпъдйс(хмрепбюкю б йнрнпнл мюундхряъ йнпемэ) яннрберярбсер цкюбмнлс йбюмрнбнлс вхякс
           ISolution=0
           IF(IndexCorrect.EQ.(NNEL+1)) THEN
              ISolution=1
		   ENDIF
		   
		   
		   IF(ISolution.EQ.1) THEN
              ! тхйяхпсел пефхл онхяйю йнпмъ б хмрепбюке
		      IndexVilky=2 
		      !WRITE(*,*) 'ZERO',IndexVilky
		      !WRITE(*,*) DetM(5),DetM(4),DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(5),ExfM(4),ExfM(2),ExfM(3)
		      !READ(*,*) 
		      ! гюохяшбюел цпюмхжш хмрепбюкю 
              ExfMA(2)=ExfM(4)
		      DetMA(2)=DetM(4)
		      ExfMA(3)=ExfM(3)
		      DetMA(3)=DetM(3)
		      !WRITE(*,*) DetM(2),DetM(3)
		      !WRITE(*,*) ExfM(2),ExfM(3)
		      !READ(*,*)   
		      ! гюмскъел оепед онхяйнл йнпмъ
		      ExfMA(1)=0.D0
		      DetMA(1)=0.D0
		      ExfMA(4)=0.D0
              DetMA(4)=0.D0
		      ExfMA(5)=0.D0
		      DetMA(5)=0.D0
		      INDEXDF=0 
           ENDIF
        ENDIF 
	
		  		    
   ENDIF

   ! пефхл онхяйю йнпмъ б хмрепбюке
   IF(IndexVilky.EQ.2) THEN 
      ! нясыеярбкъел онхяй йнпмъ б мюидемнл хмрепбюке
      INDEXDF=INDEXDF+1
      ! дкъ оепбнцн пюгю сйюгшбюел опхакхфеммши йнпемэ
	  IF(INDEXDF.EQ.1) THEN
           Exf=(DetMA(3)*ExfMA(2)-DetMA(2)*ExfMA(3))/(DetMA(3)-DetMA(2))
		   ExfMA(1)=ExfMA(2)
           ExfMA(2)=ExfMA(3)
		   ExfMA(3)=Exf
		 ELSE
           ! опнбепъел б йюйнл хг дбсу хмрепбюкюу йнпемэ
		   ! оепбши хмрепбюк
		   IF((DetMA(1)*DetMA(3)).LT.0.D0) THEN
               IIKK=1
		       Exf=(DetMA(3)*ExfMA(1)-DetMA(1)*ExfMA(3))/(DetMA(3)-DetMA(1))
			   ExfMA(2)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
		   ! брнпни хмрепбюк
           IF((DetMA(3)*DetMA(2)).LT.0.D0) THEN
		       IIKK=2
               Exf=(DetMA(2)*ExfMA(3)-DetMA(3)*ExfMA(2))/(DetMA(2)-DetMA(3))
		       ExfMA(1)=ExfMA(3)
               ExfMA(3)=Exf
		   ENDIF
      ENDIF		  		   
   

      ! опнбепъел мюидем йнпемэ хкх мер
	  ! опнбепъел бшонкмхкняэ кх сякнбхе
      ! опнбепъел гмювемхе нопедекхрекъ нмн днкфмн ашрэ лемэье рнвмнярх
      ! мювхмюел опнбепърэ сякнбхе б яксвюе намюпсфемхъ хмрепбюкю цде мюундхряъ йнпемэ
      IF(INDEXDF.GT.1) THEN
         IF(DABS(DetMA(3)).LT.EPS) THEN
            ! йнпемэ мюидем 
			! тхйяхпсел мюидеммши йнпемэ
			IndexSolutions=IndexSolutions+1
			! гюохяшбюел мюидеммши йнпемэ
			ExfSolutions(IndexSolutions)=Exf
			! оепеундхл б пефхл дюкэмеиьецн онхяйю йнпмеи
			IndexVilky=0 
            INDEXDF=0 
			Exf=ExfTemp 
			!write(6,*) 'kor',NNEL,IndexSolutions,ExfSolutions(IndexSolutions)
         ENDIF
      ENDIF 
	  
	  ! бшунд хг жхйкю он онхяйс йнпмеи опнхгнидер опх сякнбхх
	  IF(IndexSolutions.EQ.3) THEN
	     ! сякнбхе бшундю бшонкмемн
   	     IparametrDeterEnergy=2
		 ! ярюмдюпрмши бшунд 
		 INDEXGRNYT=0 
	  ENDIF 
	   
   ENDIF
    

  
   !WRITE(*,*) 'ITERASTION ZERO'
   !write(*,10685)  NEpsIter,IndexVilky,IndexCorrect,NNEL
   !write(*,19685)  DetM(1),DetM(3),DetM(2)
   !write(*,19685)  ExfM(1),ExfM(3),ExfM(2)
   !write(*,19685)  Exf
   !read(*,*)
   return
 end subroutine SearchSolutionAcceleratedLeft2




 ! ондопнцпюллю хглемемхъ ьюцю
 subroutine ChangeHR(INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,Henergy,RcoffChengH,ExfM,ProDet,DetM,INDEXPeresc,DeltaEnerPers,DeltaFunPers)
   implicit none
   integer::INDEXRShag,IKLCOMBEC,IIHD,IIOXDPZ,IndexCorrect,INDEXPeresc
   real(8)::Henergy,RcoffChengH,DeltaEnerPers,DeltaFunPers
   real(8),dimension(:)::ExfM,ProDet,DetM
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IndexZZ1,IndexZZ2
   real(8)::RcoffChengHOld,RDS,DeltaEnerPersOld,DeltaFunPersOld

   ! йкчв сйюгшбючыхи менаундхлнярэ бнгбпюыемхъ мю ндмс рнвйс мюгюд
   IKLCOMBEC=0 
   
   ! бнгбпюыюел б хяундмне янярнъмхе ьюц
   Henergy=Henergy*RcoffChengH 
   ! гюохяшбюел ярюпне гмювемхе 
   RcoffChengHOld=RcoffChengH
   
   ! пюявхршбюел мнбши йнщттхжхемр
   RcoffChengH=ProDet(2)/ProDet(1)
   ! опнбепъел йпхрепхи
   IF(RcoffChengH.GT.2.1D0) THEN
      RcoffChengH=RcoffChengH*2.D0
   ENDIF 
   ! нцпюмхвемхе он бекхвхме слемэьемхъ ьюцю
   ! мер пефхлю слемэьемхъ ьюцю б 10 пюг
   IF(RcoffChengH.GT.12.D0.AND.INDEXRShag.EQ.0) THEN
      RcoffChengH=12.D0
   ENDIF 
   ! еярэ пефхл слемэьемхъ ьюцю б 10 пюг 
   IF(RcoffChengH.GT.20.D0.AND.INDEXRShag.NE.0) THEN
      RcoffChengH=20.D0
   ENDIF 

   IF(RcoffChengH.LT.1.D0) THEN
      RcoffChengH=1.D0
   ENDIF

   ! бнглнфмше йнппейрхпнбйх б пефхле слемэьемхъ ьюцю б 10*IIHD пюг
   
   IF(INDEXRShag.NE.0) THEN
      INDEXRShag=INDEXRShag+1
   ENDIF 

   RDS=RcoffChengH/RcoffChengHOld
   IF(RDS.GT.2.D0.AND.INDEXRShag.EQ.0) THEN
      ! оепеундхл й пефхлс слемэьемхъ ьюцю
	  INDEXRShag=1 
	  ! нясыеярбкъел ялеыемхе мю ндмс рнвйс мюгюд
	  IKLCOMBEC=1
	  ! сйюгшбюел йнщттхжхемр слемэьемхъ ьюцю
	  IIHD=1
   ENDIF

   ! опнбепъел опнхяундхр опхакхфемхе й мскч хкх мер
   IF(INDEXRShag.NE.0) THEN
      ! опнбепъел мнлеп йнпмъ 
	  ! яйювей бнгмхйюер дкъ брнпнцн йнпмъ х бшье   
      IF(IndexCorrect.NE.0) THEN
         ! опнбепъел опнькх б пефхл слемэьемъ ьюцю 
	     IF(IIOXDPZ.EQ.1) THEN
	        ! опнбепъел врнаш тсмйжхх ашкх б ндмни накюярх
	        IndexZZ1=0 
	        IF(DetM(5).GT.0.D0) THEN
               IndexZZ1=IndexZZ1+1 
	        ENDIF
            IF(DetM(4).GT.0.D0) THEN
               IndexZZ1=IndexZZ1+1 
	        ENDIF
	        IF(DetM(3).GT.0.D0) THEN
               IndexZZ1=IndexZZ1+1 
	        ENDIF
		    IF(DetM(2).GT.0.D0) THEN
               IndexZZ1=IndexZZ1+1 
	        ENDIF
            IndexZZ2=0 
            IF(DetM(5).LT.0.D0) THEN
               IndexZZ2=IndexZZ2+1 
	        ENDIF
            IF(DetM(4).LT.0.D0) THEN
               IndexZZ2=IndexZZ2+1 
	        ENDIF
	        IF(DetM(3).LT.0.D0) THEN
               IndexZZ2=IndexZZ2+1 
	        ENDIF
            IF(DetM(2).LT.0.D0) THEN
               IndexZZ2=IndexZZ2+1 
	        ENDIF
	        ! опнбепъел тсмйжхъ мюундхряъ б ндмни накюярх
            IF(IndexZZ1.EQ.4.OR.IndexZZ2.EQ.4) THEN
               ! опнбепъел ме бнгмхйкю яхрсюжхъ оепеяйнйю нр опюбни вюярх пюгпшбю б кебсч вюярэ 
	           IF(DABS(DetM(3)).LT.DABS(DetM(2))) THEN
			      IF(DABS(DetM(5)).GT.DABS(DetM(4)).AND.DABS(DetM(4)).LT.DABS(DetM(3))) THEN
				     ! оноюкх б пефхл оепеяйнйю
					 INDEXPeresc=INDEXPeresc+1
                     ! пюгмхжю щмепцхи опедшдсыее гмювемхе
                     DeltaEnerPersOld=DeltaEnerPers
					 ! пюгмхжю тсмйжхи опедшдсыее гмювемхе 
					 DeltaFunPersOld=DeltaFunPers
					 ! пюгмхжю щмепцхи дюммши пюявер
                     DeltaEnerPers=DABS(ExfM(5)*(1.D0-ExfM(4)/ExfM(5)))
					 ! пюгмхжю тсмйжхи дюммши пюявер
                     DeltaFunPers=DABS(DetM(5)*(1.D0-DetM(4)/DetM(5)))
					 IF(INDEXPeresc.EQ.1) THEN
	                    ! нясыеярбкъел ялеыемхе мю ндмс рнвйс мюгюд
	                    IKLCOMBEC=1
                        ! сйюгшбюел йнщттхжхемр слемэьемхъ ьюцю
	                    IIHD=IIHD+2
		                ! опнбепъел бшунд хг пефхлю сфе лнфер ашрэ нясыеярбкем хкх мер
                        IF(INDEXRShag.GT.25) THEN
                           ! слемэьюел вхякн хрепюжхи б пефхле слемэьеммнцн ьюцю, врнаш
		                   ! нм ме гюйнмвхкяъ пюмэье рнцн йюй нопедекхряъ йнпемэ
			               INDEXRShag=INDEXRShag-1 
		                ENDIF
					 ENDIF
					 IF(INDEXPeresc.GT.1) THEN
                        ! пюяялюрпхбюел бшонкмемхъ йпхрепхъ 
						! хмрепбюк лефдс рнвюлх яйювйю слемэьюеряъ х сбекхвхбюеряъ пюгмхжю лефдс гмювемхълх тсмйжхх яйювйю
                        IF(DeltaEnerPers.LT.DeltaEnerPersOld.AND.(DeltaFunPers/DeltaEnerPersOld).GT.1.5D0) THEN
                             ! пюгмхжю лефдс гмювемхълх тсмйжхи я слемэьемхел хмрепбюкю бнгпюярюер якеднбюрекэмн щрн пюгпшб
						     ! еякх щрнцн ме опнхяундхр щрн яйювей х йнпмъ б акхгх яйювйю х бнгмхймер
						     ! нясыеярбкъел ялеыемхе мю ндмс рнвйс мюгюд
	                         IKLCOMBEC=1
                             ! сйюгшбюел йнщттхжхемр слемэьемхъ ьюцю
	                         IIHD=IIHD+2
		                     ! опнбепъел бшунд хг пефхлю сфе лнфер ашрэ нясыеярбкем хкх мер
                             IF(INDEXRShag.GT.25) THEN
                                ! слемэьюел вхякн хрепюжхи б пефхле слемэьеммнцн ьюцю, врнаш
		                        ! нм ме гюйнмвхкяъ пюмэье рнцн йюй нопедекхряъ йнпемэ
			                    INDEXRShag=INDEXRShag-1 
		                     ENDIF
						   ELSE
						     ! пюгпшбю мер бнгбпюыюеляъ б пефхл нрясрярбхъ оепеяйнйю 
							 INDEXPeresc=0
						ENDIF
					 ENDIF
	              ENDIF
			   ENDIF    
            ENDIF
         ENDIF
      ENDIF
	    
	  ! пюяявхршбюел йнщттхжхемр
      RcoffChengH=RcoffChengH*10.D0*DBLE(IIHD)
	  
	  ! опнбепъел опнхяундхр опхакхфемхе й мскч хкх мер
	  RDS=RcoffChengH/RcoffChengHOld
	  IF(RDS.GT.1.D0.AND.DABS(DetM(5)).LT.DABS(DetM(4))) THEN
	     ! тхйяхпсел, врн лш оноюкх б пефхл слемэьемхъ ьюцю
		 IIOXDPZ=1
	     ! слемэьюел вхякн хрепюжхи б пефхле слемэьеммнцн ьюцю, врнаш
		 ! нм ме гюйнмвхкяъ пюмэье рнцн йюй нопедекхряъ йнпемэ 
		 IF(INDEXRShag.NE.1) THEN
            INDEXRShag=INDEXRShag-1
		 ENDIF
	  ENDIF

   ENDIF
 
   
   IF(INDEXRShag.GT.25) THEN
      ! оняке 25-рх хрепюжхи
	  ! бнгбпюыюеляъ напюрмн
	  RcoffChengH=RcoffChengH/(10.D0*DBLE(IIHD)) 
   ENDIF
   
   Henergy=Henergy/RcoffChengH
   
   return
 end subroutine ChangeHR


 ! ондопнцпюллю нясыеярбкъер пюявер оюпюлерпнб опнцнмйх
 ! NKLZ-йкчв сйюгбючыхи рхо пюяверю
 ! NKLZ=1-пюявер оюпюлерпнб дкъ яксвюъ ндмнпндмни яхярелш спюбмемхи
 ! NKLZ=2-пюявер оюпюлерпнб дкъ яксвюъ мендмнпндмни яхярелш спюбмемхи
 ! IndexM-пюглеп люяяхбю
 ! Kpoint-мнлеп рнвйх яьхбйх
 ! Npoint-вхякн рнвей
 ! AA(IndexM,IndexM,Npoint)-люяяхб A 
 ! BB(IndexM,IndexM,Npoint)-люяяхб B
 ! FX(IndexM,Npoint)-люяяхб F
 ! Vout(IndexM,IndexM,Npoint)-люяяхб опнцнмйх он мюопюбкемхч "мюпсфс"
 ! Uout(IndexM,Npoint)-люяяхб опнцнмйх он мюопюбкемхч "мюпсфс"
 ! Vin(IndexM,IndexM,Npoint)-люяяхб опнцнмйх он мюопюбкемхч "бмсрпэ"
 ! Uin(IndexM,Npoint)-люяяхб опнцнмйх он мюопюбкемхч "бмсрпэ"
 ! VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX-бяонлнцюрекэмше люяяхбш дкъ пюяверю
 subroutine CalculationParameterProgonky(NKLZ,IndexM,Kpoint,Npoint,AA,BB,FX,Vout,Uout,Vin,Uin,VNout,UNout,VNin,UNin,AZ,BZ,FZ,DZ,E,AX,BX)
   implicit none
   integer::NKLZ,IndexM,NumeroMin,NumeroMax,Kpoint,Npoint
   real(8),dimension(:)::UNout,UNin,AX,BX
   real(8),dimension(:,:)::FX,Uout,Uin,VNout,VNin,AZ,BZ,FZ,DZ,E
   real(8),dimension(:,:,:)::AA,BB,Vout,Vin
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDF,IIXX,IIYY
   45352 FORMAT(2X,100(F15.10,1X))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!нясыеярбкъел опнцнмйс мюпсфс "OUT"!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call ReadWriteMassiv2(2,1,IndexM,AA,AZ)
   call ReadWriteMassiv2(2,3,IndexM,AA,FZ)
   call ReadWriteMassiv2(2,2,IndexM,BB,BZ) 
   ! пеьюел спюбмемхъ хрепюжхнммшл лернднл онксвюел оюпюлерпш опнцнмйх б оепбни рнвйх
   call IterationSolutionEquations(AZ,FZ,BZ,VNout) 
   ! ТНПЛХПСЕЛ ЛЮЯЯХБШ Б МЮВЮКЭМНИ РНВЙХ 1,2
   call ReadWriteMassiv2(1,1,IndexM,Vout,VNout)
   call ReadWriteMassiv2(1,2,IndexM,Vout,VNout)
   !WRITE(6,*) 'POINT',1
   !DO IIXX=1,IndexM
   !   WRITE(6,45352) (VNout(IIXX,IIYY),IIYY=1,IndexM)
   !ENDDO

   IF(NKLZ.EQ.2) THEN 
      AX=0.D0
      call ReadWriteMassiv1(1,1,IndexM,Uout,AX) 
      call ReadWriteMassiv1(1,2,IndexM,Uout,AX) 
   ENDIF
    
   DO IIDF=3,Kpoint+1
      ! явхршбюел люяяхбш
	  call ReadWriteMassiv2(2,IIDF-1,IndexM,AA,AZ)
      call ReadWriteMassiv2(2,IIDF+1,IndexM,AA,FZ)
      call ReadWriteMassiv2(2,IIDF,IndexM,BB,BZ) 
      call ReadWriteMassiv2(2,IIDF-1,IndexM,Vout,E)
      DZ=BZ-MATMUL(AZ,E)
      ! мюундхл напюрмсч люрпхжс
	  call InverseMatrix(DZ,DZ)
      VNout=MATMUL(DZ,FZ)
      ! гюохяшбюел пегскэрюр пюяверю дкъ дюммни рнвйх б люяяхб Vout
	  call ReadWriteMassiv2(1,IIDF,IndexM,Vout,VNout)
      
      !WRITE(6,*) 'POINT',IIDF
      !DO IIXX=1,IndexM
      !   WRITE(6,45352) (VNout(IIXX,IIYY),IIYY=1,IndexM)
      !ENDDO


	  ! сярюмюбкхбюел рхо пюяверю 
	  IF(NKLZ.EQ.2) THEN 
         call ReadWriteMassiv1(2,IIDF,IndexM,FX,AX)
		 call ReadWriteMassiv1(2,IIDF-1,IndexM,Uout,BX)
		 UNout=MATMUL(DZ,MATMUL(AZ,BX)-AX)  
		 ! гюохяшбюел пегскэрюр пюяверю
         call ReadWriteMassiv1(1,IIDF,IndexM,Uout,UNout) 
	  ENDIF 
   ENDDO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! нясыярбкъел опнцнмйс бмсрпэ "IN"!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call ReadWriteMassiv2(2,Npoint,IndexM,AA,AZ)
   call ReadWriteMassiv2(2,Npoint-2,IndexM,AA,FZ)
   call ReadWriteMassiv2(2,Npoint-1,IndexM,BB,BZ) 
   ! пеьюел спюбмемхъ хрепюжхнммшл лернднл онксвюел оюпюлерпш опнцнмйх б онякедмеи рнвйх
   call IterationSolutionEquations(AZ,FZ,BZ,VNin) 
   
   ! ТНПЛХПСЕЛ ЛЮЯЯХБШ Б ЙНМЕВМНИ РНВЙХ Npoint,Npoint-1
   call ReadWriteMassiv2(1,Npoint,IndexM,Vin,VNin)
   call ReadWriteMassiv2(1,Npoint-1,IndexM,Vin,VNin)
   
   IF(NKLZ.EQ.2) THEN 
      AX=0.D0
      call ReadWriteMassiv1(1,Npoint,IndexM,Uin,AX) 
      call ReadWriteMassiv1(1,Npoint-1,IndexM,Uin,AX) 
   ENDIF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!нясыеярбкъел опнцнмйс бмсрпэ "IN"!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO IIDF=Npoint-2,Kpoint+1,-1
      ! явхршбюел люяяхбш
	  call ReadWriteMassiv2(2,IIDF-1,IndexM,AA,AZ)
      call ReadWriteMassiv2(2,IIDF+1,IndexM,AA,FZ)
      call ReadWriteMassiv2(2,IIDF,IndexM,BB,BZ) 
      call ReadWriteMassiv2(2,IIDF+1,IndexM,Vin,E)
      DZ=BZ-MATMUL(FZ,E)
      ! мюундхл напюрмсч люрпхжс
      call InverseMatrix(DZ,DZ)
      VNin=MATMUL(DZ,AZ)
      ! гюохяшбюел пегскэрюр пюяверю дкъ дюммни рнвйх б люяяхб Vin
	  call ReadWriteMassiv2(1,IIDF,IndexM,Vin,VNin)
      
	  ! сярюмюбкхбюел рхо пюяверю 
	  IF(NKLZ.EQ.2) THEN 
         call ReadWriteMassiv1(2,IIDF,IndexM,FX,AX)
		 call ReadWriteMassiv1(2,IIDF+1,IndexM,Uin,BX)
		 UNin=MATMUL(DZ,MATMUL(FZ,BX)-AX)  
		 ! гюохяшбюел пегскэрюр пюяверю
         call ReadWriteMassiv1(1,IIDF,IndexM,Uin,UNin) 
	  ENDIF 
   ENDDO


   return
 end subroutine CalculationParameterProgonky


 ! ондопнцпюллю пеьемхъ люрпхвмнцн спюбмемхъ хрепюжхнммшл лернднл 
 subroutine IterationSolutionEquations(A1,A3,B2,V)
   implicit none
   real(8),dimension(:,:)::A1,A3,B2,V
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::n,ierr,IIKL,IOUP,IASS,IKLP
   real(8)::delta
   real(8),parameter::epsd=1.d-40
   real(8),allocatable,dimension(:,:)::V1,V2,DeltaV,C
   
   n=size(V,dim=1)

   ! бшдекъел оълърэ онд люяяхбш
   allocate(C(n,n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'MEMORY ON THE FILE "C" IS NOT SELECTED'
  	  read(*,*)
	  stop 
   endif  
   allocate(V1(n,n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'MEMORY ON THE FILE "V1" IS NOT SELECTED'
  	  read(*,*)
	  stop 
   endif   
   allocate(V2(n,n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'MEMORY ON THE FILE "V2" IS NOT SELECTED'
	  read(*,*)
	  stop 
   endif   
   allocate(DeltaV(n,n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'MEMORY ON THE FILE "DeltaV" IS NOT SELECTED'
	  read(*,*)
	  stop 
   endif

   59412 FORMAT(2X,100(F15.10,1X))
   !write(6,*) 'B2'
   !DO IASS=1,n
   !   WRITE(6,59412) (B2(IASS,IKLP),IKLP=1,n)
   !ENDDO
 
   !write(6,*) 'A1'
   !DO IASS=1,n
   !   WRITE(6,59412) (A1(IASS,IKLP),IKLP=1,n)
   !ENDDO
  
   !write(6,*) 'A3'
   !DO IASS=1,n
   !   WRITE(6,59412) (A3(IASS,IKLP),IKLP=1,n)
   !ENDDO
   !STOP
   ! мювюкэмне опхакхфемхе
   C=B2
   ! мюундхл нпаюрмсч люрпхжс
   call InverseMatrix(C,C)
   V1=MATMUL(C,A3)   
   V2=V1
   delta=1.D0
   DO WHILE(delta.GT.epsd) 
      C=B2-MATMUL(A1,V1)
	  ! мюундхл нпаюрмсч люрпхжс
      call InverseMatrix(C,C)
	  V1=MATMUL(C,A3)
      DeltaV=V1-V2
      ! гюохяшбюел пегскэрюр пюяверю дкъ мнбни хрепюжхх
	  V2=V1
      ! нопедекъел дереплхмюмр люрпхжш
	  delta=DABS(Determenant(DeltaV))
   ENDDO
   
   ! гюохяшбюел пегскэрюр пюяверю онксвеммнцн хрепюжхнммшл лернднл 
   V=V2 
 

   ! сдюкемхе люяяхбнб хг оълърх
   deallocate(C,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'THE FILE "C" IS NOT REMOVED FROM MEMORY'
	  read(*,*) 
	  stop 
   endif
   deallocate(V1,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'THE FILE "V1" IS NOT REMOVED FROM MEMORY'
	  read(*,*) 
	  stop 
   endif
   deallocate(V2,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'THE FILE "V2" IS NOT REMOVED FROM MEMORY'
	  read(*,*)
	  stop 
   endif
   deallocate(DeltaV,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'IterationSolutionEquations'
      write(*,*) 'THE FILE "DeltaV" IS NOT REMOVED FROM MEMORY'
	  read(*,*)
	  stop 
   endif

   return
 end subroutine IterationSolutionEquations


 ! ондопнцпюллю пюяверю нопедекхрекъ люрпхжш 
 ! n-ПЮГЛЕП ЛЮРПХЖШ
 ! a(n,n)-ЛЮРПХЖЮ	
 real(8) function Determenant(a)
   use msimsl,nouse=>fac
   implicit none
   real(8),dimension(:,:)::a
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::n,ierr
   real(8)::det1,det2
   integer,allocatable,dimension(:)::ipvt   
   real(8),allocatable,dimension(:,:)::fac
 
   n=size(a,dim=1)	
   ! бшдекъел оълърэ онд люяяхбш
   allocate(fac(n,n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Determenant'
      write(*,*) 'MEMORY ON THE FILE "fac" IS NOT SELECTED'
  	  read(*,*)
	  stop 
   endif   
   allocate(ipvt(n),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Determenant'
      write(*,*) 'MEMORY ON THE FILE "ipvt" IS NOT SELECTED'
	  read(*,*)
	  stop 
   endif   



   !LU-пюгкнфемхе
   call DLftrg(n,a,n,fac,n,ipvt)
   call DLfdrg(n,fac,n,ipvt,det1,det2)
   
   Determenant=det1*(10.D0)**det2 

   ! сдюкемхе люяяхбнб хг оълърх 
   deallocate(fac,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'Determenant'
      write(*,*) 'THE FILE "fac" IS NOT REMOVED FROM MEMORY'
	  read(*,*) 
	  stop 
	endif
      deallocate(ipvt,stat=ierr)
	if(ierr/=0) then
	  write(*,*) 'Determenant'
      write(*,*) 'THE FILE "ipvt" IS NOT REMOVED FROM MEMORY'
	  read(*,*)
	  stop 
	endif

	return    
  end function  Determenant

 ! ондопнцпюллю нопедекъчыюъ рнвйс яьхбйх пеьемхи OUT AND IN
 ! Npoint-вхякн рнвей
 ! NumeroMin-лхмхлюкэмши мнлеп цюплнмхйх
 ! Exf-ндмнщкейрпнммюъ щмепцхъ уюпрпх-тнйю (лнкейскъпмни нпахрюкх)
 ! Pz(NumbreGarmMO(NMO),NumbreGarmMO(NMO),Npoint) - люяяхб гмювемхи P(ro)-ВЮЯРХ КНЦЮКЭМНИ ТСМЙЖХХ P БУНДЪЫЕИ Б СПЮБМЕМХЕ Y''=P*Y+Q-кнйюкэмюъ вюярэ онремжхюкю
 ! RO1X(Npoint)- люяяхб гмювемхи оепбни опнхгбндмни он мнбни оепелеммни
 ! Kpoint-мнлеп рнвйх яьхбйх пеьемхи
 subroutine PointIntersection(Npoint,NumeroMin,Exf,Pz,RO1X,Kpoint)
   implicit none
   integer::Npoint,NumeroMin,Kpoint
   real(8)::Exf
   real(8),dimension(:)::RO1X
   real(8),dimension(:,:,:)::Pz
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDXK,IKLW,NumeroKoren,Isumk
   real(8)::F1,F2
   
   ! мнлеп йнпмъ опхмхлюецн гю рнвйс яьхбйх
   NumeroKoren=1
    
   IKLW=1
   Isumk=0
   DO WHILE(IKLW.EQ.1) 
	  F1=Pz(NumeroMin,NumeroMin,1)+Exf/RO1X(1)**2 
      DO IIDXK=2,Npoint
         F2=Pz(NumeroMin,NumeroMin,IIDXK)+Exf/RO1X(IIDXK)**2
		 IF(F1*F2.LT.0.D0) THEN
		    Isumk=Isumk+1
			! мюидем йнпемэ
            IF(Isumk.EQ.NumeroKoren) THEN 
	           Kpoint=IIDXK-1
	 	       IKLW=2
	           EXIT
		    ENDIF
		 ENDIF 
	     F1=F2
      ENDDO
      ! йнппейрхпсел гмювемхе щмепцхх
	  IF(IKLW.EQ.1) THEN
         Exf=Exf*0.8D0 
      ENDIF
   ENDDO 
   return
 end subroutine PointIntersection


 ! ондопнцпюллю времхъ хг люяяхбю б люяяхб
 ! NN-мнлеп рнвйх йнрнпни яннрберярбсер люрпхжш
 ! NumeroMin,NumeroMax-лхмхлюкэмне х люйяхлюкэмне гмювемхе хмдейяю
 ! D(N1,N2,N3)-рпеулепмши люяяхб
 ! A(N1,N2)-дбсулепмши люяяхб 
 subroutine ReadWriteMassiv2X(NN,NumeroMin,NumeroMax,D,A) 
   implicit none
   integer::KL,NN,NumeroMin,NumeroMax
   real(8),dimension(:,:)::A
   real(8),dimension(:,:,:)::D
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDX,IIDY,IXXK,YXXK

   ! времхе хг люяяхбю
   !IXXK=0
   !DO IIDX=NumeroMin,NumeroMax
   !   IXXK=IXXK+1
   !   YXXK=0
   !   DO IIDY=NumeroMin,NumeroMax
   !      YXXK=YXXK+1
   !	  A(IXXK,YXXK)=D(IIDY,IIDX,NN)
   !   ENDDO
   ! ENDDO


   FORALL (IIDX=NumeroMin:NumeroMax)
     FORALL (IIDY=NumeroMin:NumeroMax)
       A(IIDX-NumeroMin+1,IIDY-NumeroMin+1)=D(IIDY,IIDX,NN) 
	 END FORALL
   END FORALL

   return
 end subroutine ReadWriteMassiv2X

 ! ондопнцпюллю времхъ х гюохях  люяяхбю б люяяхб
 ! KL=1-йкчв времхъ
 ! KL=2-йкчв гюохях
 ! NN-мнлеп рнвйх йнрнпни яннрберярбсер люрпхжш
 ! NumeroMin,NumeroMax-лхмхлюкэмне х люйяхлюкэмне гмювемхе хмдейяю
 ! D(N1,N2)-дбсулепмши люяяхб
 ! A(N1)-ндмнлепмши люяяхб 
 subroutine ReadWriteMassiv1X(KL,NN,NumeroMin,NumeroMax,D,A) 
   implicit none
   integer::KL,NN,NumeroMin,NumeroMax
   real(8),dimension(:)::A
   real(8),dimension(:,:)::D
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDX,IXXD

   ! гюохях б люяяхб
   IF(KL.EQ.1) THEN  
      !IXXD=0
      !DO IIDX=NumeroMin,NumeroMax
      !   IXXD=IXXD+1
	  !   D(IIDX,NN)=A(IXXD)
      !ENDDO
	  FORALL (IIDX=NumeroMin:NumeroMax)
        D(IIDX,NN)=A(IIDX-NumeroMin+1)
	  END FORALL
   ENDIF
   
   ! времхе хг люяяхбю
   IF(KL.EQ.2) THEN  
	  !IXXD=0
      !DO IIDX=NumeroMin,NumeroMax
      !   IXXD=IXXD+1
	  !   A(IXXD)=D(IIDX,NN)
      !ENDDO
	  FORALL (IIDX=NumeroMin:NumeroMax) 
        A(IIDX-NumeroMin+1)=D(IIDX,NN) 
	  END FORALL
   ENDIF
   
   return
 end subroutine ReadWriteMassiv1X








 ! ондопнцпюллю времхъ х гюохях люяяхбю б люяяхб
 ! KL=1-йкчв времхъ
 ! KL=2-йкчв гюохях
 ! NN-мнлеп нвйх йнрнпни яннрберярбсер люрпхжш
 ! IndexM-люйяхлюкэмне гмювемхе хмдейяю
 ! D(N1,N2,N3)-рпеулепмши люяяхб
 ! A(N1,N2)-дбсулепмши люяяхб 
 subroutine ReadWriteMassiv2(KL,NN,IndexM,D,A) 
   implicit none
   integer::KL,NN,IndexM
   real(8),dimension(:,:)::A
   real(8),dimension(:,:,:)::D
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDX,IIDY

   ! гюохяэ б люяяхб
   IF(KL.EQ.1) THEN
      !DO IIDX=1,IndexM
      !   DO IIDY=1,IndexM
      !      D(IIDY,IIDX,NN)=A(IIDY,IIDX)
	  !	 ENDDO
	  !ENDDO
      FORALL (IIDX=1:IndexM)
        FORALL (IIDY=1:IndexM)
          D(IIDY,IIDX,NN)=A(IIDY,IIDX)
		END FORALL
	  END FORALL
   ENDIF 
   
   ! времхе хг люяяхбю
   IF(KL.EQ.2) THEN
      !DO IIDX=1,IndexM
      !   DO IIDY=1,IndexM
      !      A(IIDY,IIDX)=D(IIDY,IIDX,NN)
	  !	 ENDDO
	  !ENDDO
      FORALL (IIDX=1:IndexM)
	    FORALL (IIDY=1:IndexM)
          A(IIDY,IIDX)=D(IIDY,IIDX,NN)
		END FORALL
	  END FORALL 
   ENDIF 


   return
 end subroutine ReadWriteMassiv2

 ! ондопнцпюллю времхъ х гюохях люяяхбю б люяяхб
 ! KL=1-йкчв времхъ
 ! KL=2-йкчв гюохях
 ! NN-мнлеп нвйх йнрнпни яннрберярбсер люрпхжш
 ! IndexM-люйяхлюкэмне гмювемхе хмдейяю
 ! D(N1,N2)-дбсулепмши люяяхб
 ! A(N1)-ндмнлепмши люяяхб 
 subroutine ReadWriteMassiv1(KL,NN,IndexM,D,A) 
   implicit none
   integer::KL,NN,IndexM
   real(8),dimension(:)::A
   real(8),dimension(:,:)::D
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IIDX

   ! гюохяэ б люяяхб
   IF(KL.EQ.1) THEN
      !DO IIDX=1,IndexM
      !   D(IIDX,NN)=A(IIDX)
	  !ENDDO
	  FORALL (IIDX=1:IndexM) D(IIDX,NN)=A(IIDX)
   ENDIF 
   
   ! времхе хг люяяхбю
   IF(KL.EQ.2) THEN
      !DO IIDX=1,IndexM
      !   A(IIDX)=D(IIDX,NN)
	  !ENDDO
	  FORALL (IIDX=1:IndexM) A(IIDX)=D(IIDX,NN)
   ENDIF 


   return
 end subroutine ReadWriteMassiv1

 ! ондопнцпюллю пюяверю напюрмни люрпхжш
 ! N-пюглеп йбюдпюрмни люрпхжш
 ! A(N,N)-йбюдпюрмюъ люрпхжю
 ! B(N,N)-напюрмюъ йбюдпюрмюъ люрпхжю
 subroutine InverseMatrix(A,B)
   use msimsl
   implicit none
   real(8),dimension(:,:)::A,B
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::N
   
   N=size(A,dim=1)

   call DLINRG(N,A,N,B,N)

   return
 end subroutine InverseMatrix

end module msde
