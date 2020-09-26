       
module mmolecule
  implicit none 


  contains


!! SUBPROGRAMME OF THE SOLUTION OF THE HARTRY-FOCK EQUATION SYSTEM FOR A MOLECULAR WITH A SINGLE-CENTER METHOD
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! IndexHartreeFock-PARAMETER INDICATING WHAT EQUATION IS DECIDED
══!! IndexHartreeFock = 0-HARTRI EQUATION
══!! IndexHartreeFock = 1-HARTRI-FOCA EQUATION
══!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
══!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
══!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
══!! Nnuclei-NUMBER OF NUCLEI (DO NOT CONSIDER NUCLEAR AT THE BEGINNING OF THE COORDINATE)
══!! Z0-CHARGE OF THE NUCLEI AT THE BEGINNING OF THE COORDINATE
══!! Z (Nnuclei) -MASSIVE CHARGES OF NUCLEI
══!! COORRR (N) -MASSIVE RADIAL COORDINATES
══!! COORRR (N) -RADIAL COORDINATE OF THE NUCLEI (R-RADIUS)
══!! COOAAA (N, 2) -MASSIVE NORDER COORDINATE
══!! COOAAA (N, 1) -WIRTH COORDINATE -RTeta-ANGLE (0 = <RTeta <= PI)
══!! COOAAA (N, 2) - COORDINATE CO-ORDINATE -RFu-ANGLE (0 = <RFu <= 2 * PI)
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoffNucleusElectron (NumbreMO) -MASSIVE NUMBER OF TYPES OF ELECTRON INTERACTION COEFFICIENTS WITH NUCLEI OF MOLECULES
══!! NtypIndexCoffF (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF COEFFICIENTS OF DIRECT CULON INTERACTION
══!! NtypIndexCoffG (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF COEFFICIENTS OF EXCHANGE COULIN INTERACTION
══!! RcoffNucleusElectron (NtypIndexCoffNucleusElectron, NumbreMO) -MASIVE OF ELECTRON INTERACTION FACTORS WITH MOLECULES NUCLEI
══!! RcoffA (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION FACTORS
══!! RcoffB (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! MLNucleusElectron (2, NtypIndexCoffNucleusElectron, NumbreMO) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING INTERACTION OF ELECTRON WITH NUCLEI OF MOLECULES
══!! MassivMLA (NtypIndexCoffF, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE DIRECT PART OF COULON INTERACTION
══!! MassivMLB (NtypIndexCoffG, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! NN (NumbreMO) -MASSIVE OF MAIN QUANTUM NUMBERS
══!! ML (NumbreMO) -MASSIVE OF MOLECULAR ORBITAL PROJECTIONS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL
══!! NumbreGarmLMO (NumbreMO) - MASSIVE OF HARMONIC NUMBER IN MOLECULAR ORBITAL TO LgarmonicMO (3, NumbreMO) INCLUDING IT
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! LgarmonicMO (3, NumbreMO) - ORBITAL MOMENT IMPORTANCE AFTER WHICH CALCULATION IS DONE ON LEGAND FUNCTIONS (functions do not change when solving)
══!! EPS-accuracy of calculation (DEVIATION BETWEEN TWO SOLUTIONS)
══!! EPSXR-accuracy of calculation (DETERMINATION OF SINGLE-ELECTRON ENERGIES IN DECISION OF EQUATIONS)
  !! H-step
══!! N is the number of points
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RO2X (Npoint) - MASSIVE VALUES OF THE SECOND DERIVATIVE ON THE NEW VARIABLE
══!! RO3X (Npoint) - MASSIVE VALUES OF THIRD DERIVATIVE ON THE NEW VARIABLE
══!! RFun (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! NumbreGarmMOZero (NumbreMO) - HARMONIC NUMBER OF MASSES IN THE MOLECULAR ORBITAL OF THE NULL APPROXIMATION
══!! LgarmonicMOZero (3, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS OF THE ZERO APPROXIMATE
══!! RfunZero (NumbreMO, NumbreGarmMOZero (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the zeroth approximation configuration
══!! ENN0-LOWER ENERGY BORDER
══!! NumreISzam (NumbreMO) -MASSIVE NUMBER OF SHELLS OF EQUIVALENT DATA
══!! ISCalculZZ (NumbreMO, NumbreMO) -MASIVE OF EQUIVALENT DATA SHEET NUMBERS
══!! FUNBuffer-BUFFER FILE NAME (FOR THE INTERMEDIATE STORING OF RADIAL PARTS OF WAVE FUNCTIONS)
══!! IKLZERO-KEY APPROACHING THE TYPE OF THE NULL APPROACH
══!! IKLZERO = 0- CALCULATION WITHOUT EXCHANGE
══!! IKLZERO = 1- CALCULATION OF EXCHANGE
══!! NFunctionLigand (Nnuclei) -MASSIVE NUMBER OF FUNCTIONS FOR EACH LIGANDA
══!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
══!! NumbreGarmLigand (Nnuclei, NumbreFunctionLig) -MASSIVE NUMBER OF HARMONICS IN THE FUNCTION OF LIGANDA
══!! LgarmonicLigand (2, Nnuclei, NumbreFunctionLig) -MASSIVE OF THE MAXIMUM AND MINIMUM VALUE OF THE ORBITAL MOMENT OF THE LIGAN FUNCTION
══!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
══!! NomeroFPOT (NumbreMO) -NUMERA OF FILES IN WHICH POTENTIALS OF DIRECT INTERACTION WRITE
══!! IreshimSO-PARAMETER INDICATING THE TYPE OF HARMONIZATION OF MOLECULAR ORBITALS
══!! IreshimSO = 0-STANDARD MODE OF CODE
══!! IreshimSO = 1- "STRENGTHENING" HARMONIZATION MODE THE BINDING COEFFICIENTS ARE EXPRESSED
══!! FileParameterCorrection-FILE NAME IN WHICH PARAMETERS CORRECT THE WORK OF THE PROGRAM
  subroutine DHFMO(IndexHartreeFock,IndexExit,Nnuclei,Z0,Z,COORRR,COOAAA,NumbreMO,NtypIndexCoffNucleusElectron,NtypIndexCoffF,NtypIndexCoffG,IreshimSO,RcoffNucleusElectron,RcoffA,RcoffB,MLNucleusElectron,MassivMLA,MassivMLB,NN,ML,IQ,NumbreGarmMO,NumbreGarmLMO,LgarmonicMO,EPS,EPSXR,H,N,R,RO1X,RO2X,RO3X,RFun,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,ENN0,NumreISzam,ISCalculZZ,FUNBuffer,FUNBuffer1,FUNBuffer2,FUNBuffer3,FUNBuffer4,IKLZERO,NFunctionLigand,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,NumbreGarmLigand,LgarmonicLigand,RfunLigands,NomeroFPOT,FileParameterCorrection)  
	use mc3js,only:C3JS_VAR
	use msde,only:CalculationNormalizationConstaResult
	implicit none
    integer::IndexHartreeFock,IndexExit,Nnuclei,NumbreMO,N,IKLZERO,INDEXZ,IreshimSO                 
    real(8)::Z0,EPS,EPSXR,H,ENN0
    integer,dimension(:)::NN,ML,IQ,NtypIndexCoffNucleusElectron,NumbreGarmMO,NumbreGarmLMO,NumbreGarmMOZero,NumreISzam,NumbreLigand,NFunctionLigand,NomeroFPOT
    integer,dimension(:,:)::NtypIndexCoffF,NtypIndexCoffG,LgarmonicMO,LgarmonicMOZero,ISCalculZZ,NLigands,NFunLigands,NumbreGarmLigand
	integer,dimension(:,:,:)::MLNucleusElectron,NumbreFunctionLig,LgarmonicLigand
    integer,dimension(:,:,:,:)::MassivMLA,MassivMLB
    real(8),dimension(:)::Z,COORRR,R,RO1X,RO2X,RO3X
    real(8),dimension(:,:)::COOAAA,RcoffNucleusElectron
    real(8),dimension(:,:,:)::RcoffA,RcoffB,RFun,RfunZero
	real(8),dimension(:,:,:,:)::RfunLigands
	character(90)::FUNBuffer,FUNBuffer1,FUNBuffer2,FUNBuffer3,FUNBuffer4,FileParameterCorrection   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ЛЮЙЯХЛЮКЭМНЕ ВХЯКН ЩКЕЛЕМРНБ ЛЮЯЯХБЮ ДКЪ ПЮЯВЕРЮ ЯТЕПХВЕЯЙХУ ЦЮПЛНМХЙ  
    integer,parameter::NmaxIndexMassiv=400                     
    integer::MMAX,MBEG,NCLS,ierr,IXXX,IYYY,IZZZ,ILX1,ILX2,NtypIndexCoff,NmaxGarmonik,NmaxGarmonikL,Iterations,J,NMO,NumeroIter,NNNKL,IKLZExch,IRETURN,ICalFP,ISZETZK,ISUMZETZIK,IValueStronMixing
    integer::NtipMOZ,IKLXD,LmaxGarmonik
	real(4)::Time_start,Time_finish,TimeCalcul,Time_startX,Time_finishX,TimeCalculX,Time_startXX,Time_finishXX,TimeCalculXX1,TimeCalculXX2
	real(8)::VAL,FF,TETMAX,Rad,RcoffMixing,RRnormcoffZX,TBX,TEX,TBXX,TEXX
    integer,allocatable,dimension(:)::NumeroMinGarmon,MTipMOML
    integer,allocatable,dimension(:,:,:,:)::IndexMLA,IndexMLB
    real(8),allocatable,dimension(:)::F,TETA,E3,E4,E5,RnormMOXF,RparamterCalcul,RTimeCalcul
	real(8),allocatable,dimension(:,:)::DeltaEnegry
	real(8),allocatable,dimension(:,:,:)::NcrystallXX,FunMOOld,FunMONew,FunMOOldZZ,FunMONewZZ,FunMONewZZD,FunMONewZZDD,FunMOOldZZD,FunMOOldZZDD  
    real(8),allocatable,dimension(:,:,:,:)::Ncrystall
    real(8),allocatable,dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::NFunLigandsMax,NumbreLigandMax,ISUMINDEX,IXXXV,IXXXY,IXZZ1,IXZY1,IXZZ2,IXZY2,IXZZ3,IXZY3,IIL1,IIL2,IMO1,IMO2,IOX1,IOX2
	real(8)::RT1,RT2,RT3,RT4,EnergyXF,RnormMOXFII
	real(8),allocatable,dimension(:)::RPOTX
	real(8),allocatable,dimension(:,:,:)::AlfaCoffMO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::LmaxKGarm,LmaxGarm,IzonFmax,IzonGmax,IsumXX,IsumXXYY,NumbreSxemaCalculation
	integer,allocatable,dimension(:)::IZONAMASSIV
    integer,allocatable,dimension(:,:)::IZONAMASSIVG,NumbreGPOTMO,NumbreGPOTMORPOT
	integer(1),allocatable,dimension(:,:,:,:)::IFPOTIndex,IGPOTIndex
	integer,allocatable,dimension(:,:,:,:)::IFPOTZona,IGPOTZona
    real(8),allocatable,dimension(:,:)::RFILEPoltens
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::NNmax,NtypIndexCoffFMax,NtypIndexCoffGMax,IYPX,IXZX
    integer,allocatable,dimension(:)::NumbreRegion
    integer,allocatable,dimension(:,:)::IndexRegionHMIN,IndexRegionHMIND
	real(8),allocatable,dimension(:,:,:)::EnergyZeroZ,EnergyRegion,EnergyZeroZD,EnergyRegionD,EnergyRegionHMIN,EnergyRegionHMIND
    
	
	character(5),dimension(2)::LBZ  
    data LBZ/'sigma','pi'/

	3000  FORMAT(/'    n        Energy      Radius     Teta       Norm ')
	3100  FORMAT(I3,A5,2F12.6,2E11.3)
    3200  FORMAT(2X,'Niteration',1X,'NMO',5X,'Exf',6X,'Exf left',3X,'Exf zero',3X,'Exf right',1X,'Coefficient mixing',3X,'Time of iteration')            
    3300  FORMAT(5X,I3,5X,I2,1X,F10.5,1X,F10.5,1X,1X,F10.5,1X,F10.5,2X,F7.5,8X,5(F11.3,' min.',1X))
	3400  FORMAT(77X,5(F11.3,' min.',1X))
	3500  FORMAT(83X,'----------      -----------      -----------      -----------      -----------')
	3600  FORMAT(2X,'Time of calculation derect coulomb interaction full',F11.3,' min.')
	3700  FORMAT(2X,'Time of calculation exchang coulomb interaction full',F11.3,' min.') 
	3800  FORMAT(2X,'Niteration',1X,'NMO',2X,'Numbre Ligand',2X,'Numbre function',2X,'Alfa function Ligands')       
	3900  FORMAT(5X,I3,5X,I2,7X,I2,14X,I3,12X,F10.5)
	4000  FORMAT(83X,'  Total           Pot Full          Pot F	        Pot G           Equation')
	4100  FORMAT(2X,'Niteration',1X,'NMO',5X,'L<=Lk',10X,'L>Lk')
    4200  FORMAT(2X,'Niteration',1X,'NMO',5X,'Norma MO')
	4300  FORMAT(5X,I3,5X,I2,5X,F17.9)
	
	! оюпюлерпш янцкюянбюмхъ  
    ! лхмхлюкэмне х люйяхлюкэмне гмювемхе вхякю жхйкнб
	MBEG=1
	! ашк нясыеярбкем юбюпхимши яани хкх мер
	IF(IndexExit.EQ.1) THEN
       MBEG=2
	ENDIF
    MMAX=2500

	! йкчв сйюгшбючыхи мю опнжеяя ялеьхбюмхъ пеьемхъ онксвеммнцн мю опедшдсыел щрюое ян якедсчыел 
	NCLS=1 ! ялеьхбюмхе нясыеярбкъеряъ
    
	! НОПЕДЕКЪЕЛ ЛЮЙЯХЛЮКЭМНЕ ГМЮВЕМХЕ ЦЮПЛНМХЙ
	NmaxGarmonik=NumbreGarmMO(1)
	! НОНПЕДЕКЪЕЛ ЛЮЙЯХЛЮКЭМНЕ ГМЮВЕМХЕ ЦКЮБМНЦН ЙБЮМРНБНЦН ВХЯКЮ
	NNmax=NN(1)
    DO IXXX=1,NumbreMO
	   IF(NumbreGarmMO(IXXX).GT.NmaxGarmonik) THEN
          NmaxGarmonik=NumbreGarmMO(IXXX)
	   ENDIF
       IF(NN(IXXX).GT.NNmax) THEN
          NNmax=NN(IXXX)
	   ENDIF
    ENDDO
	! нопедекъел люйяхлюкэмне вхякн цюплнмхй дкъ люрпхжш
    NmaxGarmonikL=NumbreGarmLMO(1)
    DO IXXX=1,NumbreMO
	   IF(NumbreGarmLMO(IXXX).GT.NmaxGarmonikL) THEN
          NmaxGarmonikL=NumbreGarmLMO(IXXX)
	   ENDIF
    ENDDO
    
	! нопедекъел люйяхлюкэмне вхякн рхонб опълнцн х налеммнцн йскнмнбяйнцн бгюхлндеиярбхъ
    NtypIndexCoffFMax=0
	NtypIndexCoffGMax=0
    DO IXXX=1,NumbreMO
	   DO IYYY=1,NumbreMO
	      IF(NtypIndexCoffF(IXXX,IYYY).GT.NtypIndexCoffFMax) THEN
             NtypIndexCoffFMax=NtypIndexCoffF(IXXX,IYYY)
		  ENDIF
		  IF(NtypIndexCoffG(IXXX,IYYY).GT.NtypIndexCoffGMax) THEN
             NtypIndexCoffGMax=NtypIndexCoffG(IXXX,IYYY)
		  ENDIF
	   ENDDO
    ENDDO
    
	IF(NtypIndexCoffFMax.GT.NtypIndexCoffGMax) THEN
	    NtypIndexCoff=NtypIndexCoffFMax
      ELSE
        NtypIndexCoff=NtypIndexCoffGMax
    ENDIF
	
	! б яксвюе нрясрярбхъ опълшу хкх налеммшу йскнмнбяйху хмрецпюкнб 
	! сярюмюбкхбюел лхмхлюкэмне гмювемхе дкъ нрясрярбхъ ньханй я бшдекемхел оюлърх онд люяяхбш
	IF(NtypIndexCoffFMax.EQ.0) THEN
	   NtypIndexCoffFMax=1
	ENDIF 
    IF(NtypIndexCoffGMax.EQ.0) THEN
	   NtypIndexCoffGMax=1
	ENDIF
	IF(NtypIndexCoff.EQ.0) THEN
       NtypIndexCoff=1  
	ENDIF
	 

    !бшдекъел оюлърэ дкъ люяяхбнб 
	allocate(NumbreRegion(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "NumbreRegion" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(EnergyZeroZ(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyZeroZ" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(IndexRegionHMIN(NumbreMO,NNmax+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IndexRegionHMIN" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(EnergyRegionHMIN(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyRegionHMIN" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(IndexRegionHMIND(NumbreMO,NNmax+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IndexRegionHMIND" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(EnergyRegionHMIND(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyRegionHMIND" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(EnergyRegion(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyRegion" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(EnergyZeroZD(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyZeroZD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(EnergyRegionD(NumbreMO,NNmax+1,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "EnergyRegionD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(RTimeCalcul(5),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "RTimeCalcul" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(RparamterCalcul(5),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "RparamterCalcul" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(Ncrystall(NumbreMO,NmaxGarmonik,NmaxGarmonik,N),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "Ncrystall" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(NcrystallXX(NmaxGarmonik,NmaxGarmonik,N),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "NcrystallXX" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(F(NmaxIndexMassiv),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "F" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
    allocate(TETA(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "TETA" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif   
	allocate(E3(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "E3" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif   
    allocate(E4(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "E4" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(E5(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "E5" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(DeltaEnegry(NumbreMO,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "DeltaEnegry" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif    
    allocate(NumeroMinGarmon(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "NumeroMinGarmon" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	allocate(RnormMOXF(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "RnormMOXF" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	allocate(FunMOOld(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMOOld" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
    allocate(FunMONew(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMONew" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	allocate(FunMOOldZZ(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMOOldZZ" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	allocate(FunMOOldZZD(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMOOldZZD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
	allocate(FunMOOldZZDD(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMOOldZZDD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif   
    allocate(FunMONewZZ(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMONewZZ" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
    allocate(FunMONewZZD(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMONewZZD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
    allocate(FunMONewZZDD(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "FunMONewZZDD" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
	allocate(MTipMOML(2*NumbreMO*NtypIndexCoff),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "MTipMOML" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	allocate(IndexMLA(NtypIndexCoffFMax,NumbreMO,NumbreMO,4),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IndexMLA" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
	allocate(IndexMLB(NtypIndexCoffGMax,NumbreMO,NumbreMO,4),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IndexMLB" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif    
    allocate(RPOTX(N),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "RPOTX" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif    

    ! гюмскъел оепед пюявернл
    TETA=0.D0
    E3=0.D0
    E4=0.D0
	NumbreRegion=0
	EnergyZeroZ=0.D0
    EnergyRegion=0.D0
    EnergyZeroZD=0.D0
    EnergyRegionD=0.D0	
	IndexRegionHMIN=0
	EnergyRegionHMIN=0.D0
	IndexRegionHMIND=0
	EnergyRegionHMIND=0.D0
    

	! сЯРЮМЮБКХБЮЕЛ МНЛЕП ОЕПБНИ ЦЮПЛНМХЙХ НРКХВМНИ НР МСКЪ ДКЪ ЙЮФДНИ ЛНКЕЙСКЪПМНИ НПАХРЮКХ
	NumeroMinGarmon=1
	DO IXXX=1,NumbreMO
	   NNNKL=1
	   IYYY=0
	   DO ILX1=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	      IYYY=IYYY+1
		  IF(IABS(ML(IXXX)).GT.ILX1) THEN
             NNNKL=IYYY+1
		  ENDIF
	   ENDDO  
	   NumeroMinGarmon(IXXX)=NNNKL
    ENDDO

    ! щрюо 0. тнплхпсел бяонлнцюрекэмше люяяхбш х нясыеярбкъел пюявер ятепхвеяйху цюплнмхй
	call C3JS_VAR(NmaxIndexMassiv,FF,F)  
    
	! опнбндхл юмюкхг вхякю пюгкхвмшу рхонб опнейжхи лнкейскъпмшу нпахрюкеи
    
	NtipMOZ=0
	! опнейжхх опълни вюярх
    DO IXXX=1,NumbreMO
	   DO IYYY=1,NumbreMO
          ! опнбепъел нркхвхе вхякю йскнмнбяйху хмрецпюкнб нр мскъ
		  IF(NtypIndexCoffF(IXXX,IYYY).NE.0) THEN
             DO IXXXY=1,NtypIndexCoffF(IXXX,IYYY)
                DO IYPX=1,4
			       ! опнбепъел нрнапюмю дюммюъ опнейжхъ хкх мер
		           IKLXD=0
		           DO IZZZ=1,NtipMOZ
		              IF(MTipMOML(IZZZ).EQ.MassivMLA(IXXXY,IXXX,IYYY,IYPX)) THEN
                         ! дюммши рхо опнейжхх ясыеярбсер
				         IKLXD=1
			             EXIT
			          ENDIF
		           ENDDO
             
			       ! опнбепъел намюпсфемю дюммюъ опнейжхъ хкх мер
		           IF(IKLXD.EQ.0) THEN
		              ! дюммюъ опнейжхъ ме намюпсфемю, гюохяшбюел мнбши рхо опнейжхх 
                      NtipMOZ=NtipMOZ+1
                      MTipMOML(NtipMOZ)=MassivMLA(IXXXY,IXXX,IYYY,IYPX)
		           ENDIF
			    ENDDO  
		     ENDDO
          ENDIF 
	   ENDDO
    ENDDO


    ! опнейжхх налеммни вюярх
    DO IXXX=1,NumbreMO
	   DO IYYY=1,NumbreMO
	      ! опнбепъел нркхвхе вхякю йскнмнбяйху хмрецпюкнб нр мскъ
		  IF(NtypIndexCoffG(IXXX,IYYY).NE.0) THEN
             DO IXXXY=1,NtypIndexCoffG(IXXX,IYYY)
                DO IYPX=1,4
			       ! опнбепъел нрнапюмю дюммюъ опнейжхъ хкх мер
		           IKLXD=0
		           DO IZZZ=1,NtipMOZ
		              IF(MTipMOML(IZZZ).EQ.MassivMLB(IXXXY,IXXX,IYYY,IYPX)) THEN
                         ! дюммши рхо опнейжхх ясыеярбсер
				         IKLXD=1
			             EXIT
			          ENDIF
		           ENDDO
             
			       ! опнбепъел намюпсфемю дюммюъ опнейжхъ хкх мер
		           IF(IKLXD.EQ.0) THEN
		              ! дюммюъ опнейжхъ ме намюпсфемю, гюохяшбюел мнбши рхо опнейжхх 
                      NtipMOZ=NtipMOZ+1
                      MTipMOML(NtipMOZ)=MassivMLB(IXXXY,IXXX,IYYY,IYPX)
		           ENDIF
			    ENDDO  
    	     ENDDO
          ENDIF
	   ENDDO
    ENDDO


    ! тхйяхпсел рхо опнейжхх йюфдни лнкейскъпмни нпахрюкх
	DO IXXX=1,NumbreMO
	   DO IYYY=1,NumbreMO
	      ! опнбепъел нркхвхе вхякю йскнмнбяйху хмрецпюкнб нр мскъ
		  IF(NtypIndexCoffF(IXXX,IYYY).NE.0) THEN
	         DO IXXXY=1,NtypIndexCoffF(IXXX,IYYY)
		        DO IYPX=1,4
	               ! мюундхл дюммши рхо опнейжхх б наыел люяяхбе
		           IKLXD=0
				   DO IZZZ=1,NtipMOZ
		              IF(MTipMOML(IZZZ).EQ.MassivMLA(IXXXY,IXXX,IYYY,IYPX)) THEN
                         ! дюммши рхо опнейжхх ясыеярбсер
				         IKLXD=IZZZ
			             EXIT
			          ENDIF
		           ENDDO
				   ! опнбепъел намюпсфемю дюммюъ опнейжхъ хкх мер
		           IF(IKLXD.NE.0) THEN
		                ! дюммюъ опнейжхъ намюпсфемю, гюохяшбюел рхо опнейжхх 
                        IndexMLA(IXXXY,IXXX,IYYY,IYPX)=IKLXD
		              ELSE
			            WRITE(6,*) 'ERROR. search. not ML TYP',MassivMLA(IXXXY,IXXX,IYYY,IYPX)
			            STOP 
		           ENDIF  		 
                ENDDO
		     ENDDO
	      ENDIF
	   ENDDO
    ENDDO
	 
    
    ! тхйяхпсел рхо опнейжхх йюфдни лнкейскъпмни нпахрюкх
    DO IXXX=1,NumbreMO
	   DO IYYY=1,NumbreMO
	      ! опнбепъел нркхвхе вхякю йскнмнбяйху хмрецпюкнб нр мскъ
		  IF(NtypIndexCoffG(IXXX,IYYY).NE.0) THEN
	         DO IXXXY=1,NtypIndexCoffG(IXXX,IYYY)
		        DO IYPX=1,4
	               ! мюундхл дюммши рхо опнейжхх б наыел люяяхбе
		           IKLXD=0
				   DO IZZZ=1,NtipMOZ
		              IF(MTipMOML(IZZZ).EQ.MassivMLB(IXXXY,IXXX,IYYY,IYPX)) THEN
                         ! дюммши рхо опнейжхх ясыеярбсер
				         IKLXD=IZZZ
			             EXIT
			          ENDIF
		           ENDDO
				   ! опнбепъел намюпсфемю дюммюъ опнейжхъ хкх мер
		           IF(IKLXD.NE.0) THEN
		                ! дюммюъ опнейжхъ намюпсфемю, гюохяшбюел рхо опнейжхх 
                        IndexMLB(IXXXY,IXXX,IYYY,IYPX)=IKLXD
		              ELSE
			            WRITE(6,*) 'ERROR. search. not ML TYP',MassivMLB(IXXXY,IXXX,IYYY,IYPX)
			            STOP 
		           ENDIF  		 
                ENDDO
		     ENDDO
		  ENDIF
	   ENDDO
    ENDDO
	 

  
	! бшдекъел оюлърэ онд люяяхб гмювемхи ятепхвеяйху цюплнмхй
    allocate(CkkGarmonik(2*NmaxGarmonik,NmaxGarmonik,NmaxGarmonik,NtipMOZ,NtipMOZ),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "CkkGarmonik" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 

	
	
	! нопедекъел люйяхлюкэмне гмювемхе ятепхвеяйни цюплнмхйх
    LmaxGarmonik=LgarmonicMO(2,1)
	DO IXXX=2,NumbreMO
       IF(LgarmonicMO(2,IXXX).GT.LmaxGarmonik) THEN
          LmaxGarmonik=LgarmonicMO(2,IXXX)
	   ENDIF
	ENDDO
    
	  
  
	! нясыеярбкъел гюонкмемхе люяяхбю ятепхвеяйху цюплнмхй
    call CalculSphericalGarmonik(LmaxGarmonik,NtipMOZ,MTipMOML,CkkGarmonik,FF,F)
  
    	
    
	! щрюо 1. пюявер йпхярюккхвеяйнцн онкъ  
    ! нясыеярбкъел гюмскемхе оепед пюявернл
	Ncrystall=0.D0

	DO IXXX=1,NumbreMO
	   ! жхйк он йнщттхжхемрюл бгюхлндеиярбхъ щкейрпнмю я ъдпюлх лнкейскюлх
	   DO IXXXV=1,NtypIndexCoffNucleusElectron(IXXX)
	      NcrystallXX=0.D0
		  ! нясыеярбкъел пюявер дюммнцн якюцюелнцн
		  IYYY=0
	      DO ILX1=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	         IYYY=IYYY+1
             IZZZ=0
	         DO ILX2=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	            IZZZ=IZZZ+1
			    call MATRIX_CRYSTALLINE_FIELD_MM(IYYY,ILX1,MLNucleusElectron(1,IXXXV,IXXX),IZZZ,ILX2,MLNucleusElectron(2,IXXXV,IXXX),Nnuclei,Z,COORRR,COOAAA,N,R,NcrystallXX,FF,F)
			 ENDDO
	      ENDDO
          
		  ! ясллхпсел пегскэрюрш пюяверю
          IYYY=0
	      DO ILX1=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	         IYYY=IYYY+1
             IZZZ=0
	         DO ILX2=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	            IZZZ=IZZZ+1
			    Ncrystall(IXXX,IYYY,IZZZ,1:N)=Ncrystall(IXXX,IYYY,IZZZ,1:N)+RcoffNucleusElectron(IXXXV,IXXX)*NcrystallXX(IYYY,IZZZ,1:N)
			 ENDDO
	      ENDDO
       ENDDO
	ENDDO

    ! онксвюел люрпхжс ноепюрнпю йпхярюккхвеяйнцн онкъ 
    DO IXXX=1,NumbreMO
	   ! ясллхпсел пегскэрюрш пюяверю
       IYYY=0
	   DO ILX1=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	      IYYY=IYYY+1
          IZZZ=0
	      DO ILX2=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	         IZZZ=IZZZ+1
	      	 ! жхйк он рнвйюл
	         DO IXZX=1,N  
	            Ncrystall(IXXX,IYYY,IZZZ,IXZX)=Ncrystall(IXXX,IYYY,IZZZ,IXZX)/FLOAT(IQ(IXXX))
             ENDDO
		  ENDDO
	   ENDDO
	ENDDO

    ! щрюо 2. днаюбкъел онке жемрпюкэмнцн ъдпю
    ! опнбепъел мюкхвхе ъдпю б мювюке йннпдхмюр
	IF(Z0.NE.0.D0) THEN 
	   DO IXXX=1,NumbreMO
	      IYYY=0
	      DO ILX1=LgarmonicMO(1,IXXX),LgarmonicMO(2,IXXX)
	         IYYY=IYYY+1
             ! жхйк он рнвйюл
	         DO IZZZ=1,N 
                Ncrystall(IXXX,IYYY,IYYY,IZZZ)=Ncrystall(IXXX,IYYY,IYYY,IZZZ)-Z0/R(IZZZ)  
		     ENDDO
		  ENDDO
       ENDDO
    ENDIF



    

    	
    !55569 FORMAT(2X,100(F15.10,1X)) 	

    !ILX1=3
    ! жхйк он рнвйюл
    !DO IZZZ=1,N
    !   WRITE(6,*) 'POINTYY',IZZZ 
    !   DO ILX1=1,NumbreGarmMO(1) 
	!      WRITE(6,55569) (-Ncrystall(1,ILX1,IYYY,IZZZ)/RO1X(IZZZ)**2,IYYY=1,NumbreGarmMO(1))
    !   ENDDO 
    !ENDDO
    
	
   !STOP  




    
    ! щрюо 3. нясыеярбкъел пюявер ярпсйрспш бшяьху цюплнмхй 
	! пюявер нясыеярбкъеряъ мю тсмйжхъу кхцюмднб йнрнпше ундър б ярпсйрспс лнкейскъпмни нпахрюкх 
	! НОПЕДЕКЪЕЛ ЛЮЙЯХЛЮКЭМНЕ ВХЯКН ТСМЙЖХИ Х КХЦЮМДНБ
    NFunLigandsMax=1
    NumbreLigandMax=Nnuclei
	DO IXXX=1,Nnuclei
	   IF(NFunctionLigand(IXXX).GT.NFunLigandsMax) THEN
	      NFunLigandsMax=NFunctionLigand(IXXX)
	   ENDIF
	ENDDO
    
	
    ! БШДЕКЪЕЛ ОЮЛЪРЭ ДКЪ ЛЮЯЯХБНБ 
	allocate(AlfaCoffMO(NumbreMO,NumbreLigandMax,NFunLigandsMax),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "AlfaCoffMO" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif      


    ! опнбепъел опнцпюллю нрпюанрюкю аег яанеб хкх мер
    IF(IndexExit.EQ.0) THEN
       ! гюмскъел оепед пюявернл
	   AlfaCoffMO=0.D0
	   ! тхйяхпсел йнщттхжхемрш тсмйжхи кхцюмдю дкъ оепбни хрепюжхх нмх пюбмш 
       DO IXXX=1,NumbreMO
	      ! жхйк он кхцюмдюл
	      DO IXZZ1=1,NumbreLigand(IXXX)
             ! жхйк он тсмйжхъл кхцюмдю
	         DO IXZY1=1,NFunLigands(IXXX,IXZZ1)
	            AlfaCoffMO(IXXX,NLigands(IXXX,IXZZ1),NumbreFunctionLig(IXXX,IXZZ1,IXZY1))=1.D0
	         ENDDO
	      ENDDO
	      ! нясыеярбкъел оепеявер йнщттхжхемрнб юкэтю
	      !call CalculationCoffAlfaMOFull(IXXX,N,H,NumbreGarmLMO(IXXX),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunZero,RfunLigands,AlfaCoffMO,RO1X)             
       ENDDO
	 
       ! нясыеярбкъел онярпнемхе бшяьху цюплнмхй лнкейскъпмшу нпахрюкеи мскхбнцн опхакхфемхъ
       DO IXXX=1,NumbreMO
	      call FormStructureMOFull(IXXX,N,H,NumbreGarmLMO(IXXX),NumbreGarmMO(IXXX),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunZero,RfunLigands,AlfaCoffMO,RO1X) 
       ENDDO
    ENDIF 
       
   
   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! щрюо 3.5
    ! нясыеярбкъел юмюкхг дюммшу дкъ тнплхпнбюмхъ люяяхбнб гюохях онремжхюкнб!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! бшдекъел оюлърэ онд люяяхб 
	allocate(NumbreGPOTMO(NumbreMO,NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "NumbreGPOTMO" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif
	allocate(NumbreGPOTMORPOT(NumbreMO,NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "NumbreGPOTMORPOT" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif  
	
	! нопедекъел люйяхлюкэмше цюплнмхйх
	LmaxKGarm=0
	LmaxGarm=0
	DO IYYY=1,NumbreMO  
	   IF(LmaxKGarm.LT.LgarmonicMO(3,IYYY)) THEN
          LmaxKGarm=LgarmonicMO(3,IYYY) 
	   ENDIF
	   IF(LmaxGarm.LT.LgarmonicMO(2,IYYY)) THEN
          LmaxGarm=LgarmonicMO(2,IYYY)  
	   ENDIF
	ENDDO
	! тнплхпсел люяяхб хмдейя рхош бгюхлндеиярбхъ 
	NumbreGPOTMO=0
	IsumXX=0 
	DO IXXX=1,NumbreMO 
	    DO IYYY=1,NumbreMO 
	       IF(NtypIndexCoffG(IXXX,IYYY).NE.0) THEN
		      IsumXX=IsumXX+1
	       	  NumbreGPOTMO(IXXX,IYYY)=IsumXX 
	       ENDIF
	   ENDDO
	ENDDO
	! тнплхпсел люяяхб рхонб бгюхлндеиярбхи дкъ люяяхбю пюдхюкэмшу онремжхюкнб
	NumbreGPOTMORPOT=0
    DO IXXX=1,NumbreMO 
	   IsumXXYY=0 
	   DO IYYY=1,NumbreMO 
	      IF(NtypIndexCoffG(IXXX,IYYY).NE.0) THEN
		     IsumXXYY=IsumXXYY+1
	       	 NumbreGPOTMORPOT(IXXX,IYYY)=IsumXXYY
	      ENDIF
	   ENDDO
	ENDDO  	

	! бшдекъел люяяхбш гмювемхи дкъ тнплхпнбюмхъ тюикнб онремжхюкнб опълнцн бгюхлндеиярбхъ 
    allocate(IFPOTIndex(NumbreMO,LmaxGarm+1,LmaxGarm+1,2*LmaxGarm+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IFPOTIndex" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
    allocate(IFPOTZona(NumbreMO,LmaxGarm+1,LmaxGarm+1,LmaxGarm+LmaxGarm+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IFPOTZona" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
    allocate(IZONAMASSIV(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'MEMORY ON THE FILE "IZONAMASSIV" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
    
	! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
	  ! бшдекъел люяяхбш гмювемхи дкъ тнплхпнбюмхъ тюикнб онремжхюкнб налеммнцн бгюхлндеиярбхъ 
      allocate(IGPOTIndex(IsumXX,LmaxGarm+1,LmaxKGarm+1,LmaxGarm+LmaxKGarm+1),stat=ierr)
      if(ierr/=0) then
         write(*,*) 'DHFMO'
         write(*,*) 'MEMORY ON THE FILE "IGPOTIndex" IS NOT SELECTED'
         read(*,*)
	     stop 
      endif 
	  allocate(IGPOTZona(IsumXX,LmaxGarm+1,LmaxKGarm+1,LmaxGarm+LmaxKGarm+1),stat=ierr)
      if(ierr/=0) then
         write(*,*) 'DHFMO'
         write(*,*) 'MEMORY ON THE FILE "IGPOTZona" IS NOT SELECTED'
         read(*,*)
	     stop 
      endif 
 	  allocate(IZONAMASSIVG(NumbreMO,NumbreMO),stat=ierr)
      if(ierr/=0) then
         write(*,*) 'DHFMO'
         write(*,*) 'MEMORY ON THE FILE "IZONAMASSIVG" IS NOT SELECTED'
         read(*,*)
	     stop 
      endif 
    ENDIF

	! гюмскъел оепед пюявернл 
	IZONAMASSIV=0
    IFPOTIndex=0
    ! тнплхпсел люяяхб сйюгшбючыхи гнмш гюохях х цюплнмхйх опълнцн х налеммнцн бгюхлндеиярбхъ  
    DO IXXX=1,NumbreMO   
	   DO IIL1=LgarmonicMO(1,IXXX),LgarmonicMO(3,IXXX) 
          DO IIL2=IIL1,LgarmonicMO(3,IXXX) 
	         ! пюявер опълни вюярх
	         call FPOT__Direct_Coulomb_Interaction(IXXX,IIL1,IIL2,NumbreMO,LmaxGarm,NtypIndexCoffF,RcoffA,IZONAMASSIV,IFPOTIndex,IFPOTZona,LgarmonicMO,IndexMLA,CkkGarmonik)
		  ENDDO
       ENDDO
    ENDDO

    ! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
      IZONAMASSIVG=0
      IGPOTIndex=0
      DO IXXX=1,NumbreMO   
         DO IIL1=LgarmonicMO(1,IXXX),LgarmonicMO(3,IXXX) 
            DO IIL2=LgarmonicMO(1,IXXX),LgarmonicMO(3,IXXX) 
	           ! пюявер налеммни вюярх
	  	       call GPOT__Exchange_Coulomb_Interaction(IXXX,IIL1,IIL2,NumbreMO,LmaxGarm,LmaxKGarm,NtypIndexCoffG,RcoffB,IZONAMASSIVG,NumbreGPOTMO,IGPOTIndex,IGPOTZona,LgarmonicMO,IndexMLB,CkkGarmonik)
		    ENDDO
         ENDDO
	  ENDDO
    ENDIF

       
	! яуелю пюяверю опх йнрнпни хяонкэгсеряъ люяяхб
   	NumbreSxemaCalculation=1
	
    ! нопедекъел люйяхлюкэмне вхякн гюохяеи опълни х налеммни вюярх
	
	! опълюъ вюярэ 
	IzonFmax=1
	DO IXXX=1,NumbreMO 
       IF(IzonFmax.LT.IZONAMASSIV(IXXX)) THEN
          IzonFmax=IZONAMASSIV(IXXX) 
	   ENDIF
    ENDDO
    
	! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
       ! налеммюъ вюярэ
	   IzonGmax=1
	   DO IXXX=1,NumbreMO  
          DO IYYY=1,NumbreMO 
		     IF(IXXX.NE.IYYY) THEN
	            IF(IzonGmax.LT.IZONAMASSIVG(IXXX,IYYY)) THEN
			       IzonGmax=IZONAMASSIVG(IXXX,IYYY)
			    ENDIF
	         ENDIF
	      ENDDO
	   ENDDO
	ENDIF  	
	
    !STOP

    
	! опнбепъел йюйне спюбмемхе пеьюеряъ
	IF(IndexHartreeFock.EQ.0) THEN
	   ! бшдекъел оюлърэ онд люяяхб онремжхюкнб 
       allocate(RFILEPoltens(N,IzonFmax),stat=ierr) 
       if(ierr/=0) then
          write(6,*) 'DHFMO'
          write(6,*) 'MEMORY ON THE FILE "RFILEPoltens" IS NOT SELECTED'
          write(6,*) 'For allocation of the array',FLOAT(IzonFmax)*FLOAT(N)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	      ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б тюикш
	      NumbreSxemaCalculation=2 
       endif   
    ENDIF

    IF(IndexHartreeFock.EQ.1) THEN
	   IF(IzonGmax.GT.IzonFmax) THEN
	        ! бшдекъел оюлърэ онд люяяхб онремжхюкнб 
		    allocate(RFILEPoltens(N,IzonGmax),stat=ierr) 
            if(ierr/=0) then
               write(6,*) 'DHFMO'
               write(6,*) 'MEMORY ON THE FILE "RFILEPoltens" IS NOT SELECTED'
               write(6,*) 'For allocation of the array',FLOAT(IzonGmax)*FLOAT(N)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	           ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б тюикш
	           NumbreSxemaCalculation=2 
            endif   
          ELSE
            ! бшдекъел оюлърэ онд люяяхб онремжхюкнб 
		    allocate(RFILEPoltens(N,IzonFmax),stat=ierr) 
            if(ierr/=0) then
               write(6,*) 'DHFMO'
               write(6,*) 'MEMORY ON THE FILE "RFILEPoltens" IS NOT SELECTED'
               write(6,*) 'For allocation of the array',FLOAT(IzonFmax)*FLOAT(N)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	           ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б тюикш
	           NumbreSxemaCalculation=2 
            endif   
	   ENDIF
    ENDIF
    
    
	! опнбепъел опнхгньек юбюпхимеи бшунд хкх мер
	IF(IndexExit.EQ.1) THEN 
       call ReadBufferFunction(NumbreMO,N,NumbreGarmMO,RFun,FUNBuffer)
       call ReadBufferFunction(NumbreMO,N,NumbreGarmMO,FunMONewZZ,FUNBuffer1)
	   call ReadBufferFunction(NumbreMO,N,NumbreGarmMO,FunMOOldZZ,FUNBuffer2)
       call ReadBufferFunction(NumbreMO,N,NumbreGarmMO,FunMONewZZD,FUNBuffer3)
	   call ReadBufferFunction(NumbreMO,N,NumbreGarmMO,FunMOOldZZD,FUNBuffer4)
	   ! опх пюанре опнцпюллш опнхгньек яани х опнцпюлллю опейпюрхкю пюанрс 
       ! гюохяшбюел йпхрепхх нранпю ндмнщкейрпнммни щмепцхх уюпрпх-тнйю
       DO NMO=1,NumbreMO 
	      RnormMOXF(NMO)=RFun(NMO,NumeroMinGarmon(NMO),2)
	   ENDDO	   
    ENDIF


    

	WRITE(7,*)
    WRITE(7,3200)
	WRITE(7,4000)
	WRITE(25,3800)
    WRITE(27,4100)
	WRITE(29,4200)
	
	! йкчв сйюгшбючыхи мю пюявер б оепбни хрепюжхх аег налемю
	IKLZExch=0
	! опнбепъел опнхгньек юбюпхимеи бшунд хкх мер
	IF(IndexExit.EQ.1) THEN 
	   ! опнбепъел йюйне спюбмемхе пеьюеряъ
       IF(IndexHartreeFock.EQ.1) THEN
	      ! пюявер я налемнл 
		  IKLZExch=1
	   ENDIF
       IF(IndexHartreeFock.EQ.0) THEN
	      ! пюявер аег налемнл 
		  IKLZExch=0
	   ENDIF
	ENDIF
	! гюмскъел оепед пюявернл
	ICalFP=0
	ISZETZK=0
	IValueStronMixing=0
	INDEXZ=0
	DO NumeroIter=MBEG,MMAX
	   
	   ! нясыеярбкъел пюявер бпелемх
	   call CPU_TIME(Time_startX) 
	   
	   WRITE(7,*)
	   WRITE(25,*)
	   WRITE(27,*)
       WRITE(29,*)
	   ! люйяхлюкэмне нрйкнмемхе лефдс онякеднбюрекэмшлх хрепюжхълх опх пюявер бяеу лнкейскъпмшу нпахрюкеи 
	   TETMAX=0.D0
       ! бпелъ
       RT1=0.D0
	   RT2=0.D0
	   RT3=0.D0
	   RT4=0.D0
       IF(NumeroIter.EQ.3.AND.IndexExit.EQ.0) THEN
	      ! пюявер спюбмемхтъ уюпрпх опнхгбндхряъ еякх ме ашкн юбюпхимнцн бшундю 
          ! мю рперэеи хрепюжхх пеьюел спюбмемхе аег налемю (дкъ онксвемхъ мнбнцн мскхбнцн опхакхфемхъ) 
		  ! пюявер аег налемнл 
		  IKLZExch=0 
		  ! дюммне пеьемхе ме днкфмн ашрэ ябъгюмн я опедшдсыхлх
		  NCLS=0 
	   ENDIF
	   ! гюохяшбюел мнбне опхакхфемхе дкъ йнпмеи х хмрепбюкнб 
	   EnergyZeroZD=EnergyZeroZ
	   EnergyRegionD=EnergyRegion
	   ! гюохяшбюел мнбне опхакхфемхе дкъ ьюцю дкъ нопедекемхъ цпюмхж лхмхлслнб
	   IndexRegionHMIND=IndexRegionHMIN
	   EnergyRegionHMIND=EnergyRegionHMIN

	   ! напюыюеляъ й тюикс йнппейрхпсчыецн оюпюлерпш опнцпюллш
	   call ReadCorrectionParameters(NumeroIter,FileParameterCorrection,IreshimSO) 
	  
       ! жхйк он нанкнвйюл йнмтхцспюжхх
	   DO NMO=1,NumbreMO 
	      ! опнбепъел менаундхлн нясыеярбкъер пюявер дкъ дюммни нанкнвйх
	      IF(ISCalculZZ(NMO,1).EQ.0) THEN  
		     ! нясыеярбкъел пюявер бпелемх
			 call CPU_TIME(Time_start) 
			
			 call MAINMO(IndexHartreeFock,IndexExit,NN(NMO),NCLS,IreshimSO,NumeroIter,NMO,NumbreMO,N,NtypIndexCoffF,NtypIndexCoffG,RcoffA,RcoffB,IQ,NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon,LgarmonicMO,EPSXR,H,RFun,FunMOOld,FunMONew,R,RO1X,RO2X,RO3X,Ncrystall,TETMAX,TETA,E3,E4,E5,DeltaEnegry,RnormMOXFII,RnormMOXF,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,ENN0,NumreISzam,ISCalculZZ,IKLZExch,IRETURN,RparamterCalcul,RTimeCalcul,IndexMLA,IndexMLB,CkkGarmonik,AlfaCoffMO,RfunLigands,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,LmaxKGarm,LmaxGarm,NumbreSxemaCalculation,RFILEPoltens,NumbreRegion,EnergyZeroZD,EnergyRegionD,IndexRegionHMIND,EnergyRegionHMIND)
	         
			 ! нясыеярбкъел пюявер бпелемх
			 call CPU_TIME(Time_finish)
             TimeCalcul=(Time_finish-Time_start)/60.
			 
			 ! опнлефсрнвмше пегскэрюрш пюяверю
			 IF(IRETURN.NE.1) THEN
			    WRITE(7,3300) NumeroIter,NMO,RparamterCalcul(1),RparamterCalcul(2),RparamterCalcul(3),RparamterCalcul(4),RparamterCalcul(5),TimeCalcul+RTimeCalcul(5),RTimeCalcul(1),RTimeCalcul(3),RTimeCalcul(4),RTimeCalcul(2)
		        RT1=RT1+RTimeCalcul(1)
			    RT2=RT2+RTimeCalcul(3)
			    RT3=RT3+RTimeCalcul(4)
			    RT4=RT4+RTimeCalcul(2)
                ! жхйк он кхцюмдюл
	            DO IXZZ1=1,NumbreLigand(NMO)
                   ! жхйк он тсмйжхъл кхцюмдю
		           DO IXZY1=1,NFunLigands(NMO,IXZZ1)
	                  WRITE(25,3900) NumeroIter,NMO,NLigands(NMO,IXZZ1),NumbreFunctionLig(NMO,IXZZ1,IXZY1),AlfaCoffMO(NMO,NLigands(NMO,IXZZ1),NumbreFunctionLig(NMO,IXZZ1,IXZY1))
		           ENDDO
	            ENDDO
				! нясыеярбкъел бшдювс хмтнплюжхх на щкейрпнммни окнрмнярх
		        call AnalysMOConfig(NumeroIter,NMO,N,H,NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon,RO1X,RFun) 
				! бшдюел мнплс онксвеммнцн пеьемхъ
				WRITE(29,4300)  NumeroIter,NMO,RnormMOXFII
	         ENDIF
		  ENDIF
		  ! опнбепъел рхо бшундю
	      IF(IRETURN.EQ.1) THEN
             ! хмдей сйюгшбючыхи вхякн онбрнпемхи
             INDEXZ=INDEXZ+1
             
			 ! йкчв сйюгшбючыхи йнщттхжхемр опхлеьхбюмхъ
			 IValueStronMixing=INDEXZ
						
			 ! йнпемэ ме мюидем оепеундхл й пефхлс янцкюянбюмхъ 
			 ! я оепеопюбкеммшл онремжхюкнл 
			 ICalFP=1 
			 ! тхйяхпсел мювюкн нрверю дн нрйкчвемхъ пефхлю
			 ISZETZK=1
			 ISUMZETZIK=0 
			 EXIT 
		  ENDIF
		  ! нясыеяркъел пюявер  
          IF(ISZETZK.EQ.1.AND.IRETURN.EQ.0) THEN
		     ISUMZETZIK=ISUMZETZIK+1
		  ENDIF
	   ENDDO
       
	   ! опнбепъл вхякн онксвеммшу б дюммнл пефхле йнпмеи
	   IF(ISUMZETZIK.EQ.NumbreMO) THEN
	       ! бнгбпюыюеляъ й хяундмнлс пефхлс
		   ! пефхл янцкюянбюмхъ я лемъчыхляъ онремжхюкнл 
           ICalFP=0
		   INDEXZ=0
	   ENDIF
       
	   ! опнбепъел онксвемш бяе пеьемхъ 
	   IF(ICalFP.EQ.0) THEN
          ! гюохяшбюел мнбше пегскэрюрш он хмрепбюкюл
          EnergyZeroZ=EnergyZeroZD
		  EnergyRegion=EnergyRegionD
		  ! гюохяшбюел ! гюохяшбюел мнбне опхакхфемхе дкъ ьюцю дкъ нопедекемхъ цпюмхж лхмхлслнб 
	      IndexRegionHMIN=IndexRegionHMIND
	      EnergyRegionHMIN=EnergyRegionHMIND
       ENDIF   
	   ! опнбепъел хглемхкяъ пефхл хкх мер
       IF(ICalFP.EQ.1) THEN
	      ! йнщттхжхемр ялеьхбюмхъ
	      RcoffMixing=1.D0-(0.1D0)**(IValueStronMixing+1)
	      WRITE(7,*) ' Coff Mixing'
		  WRITE(7,*) IValueStronMixing,RcoffMixing
		  ! опюбхл онремжхюк
		  ! гюохяшбюел онксвеммсч лнкейскъпмсч нпахрюкэ
	      ! опнбепъел мнлеп йнщттхжхемрю опхлеьхбюмхъ
		  ! еякх анкэье IValueStronMixing>3 хяонкэгсел опедшдсысч (опедшдсыеи) хрепюжхх 
		  IF(IValueStronMixing.LT.3) THEN
		       DO NMO=1,NumbreMO 
		          DO IXXX=1,N+2  
	                 DO IYYY=NumeroMinGarmon(NMO),NumbreGarmMO(NMO)
	                    RFun(NMO,IYYY,IXXX)=RcoffMixing*FunMOOldZZ(NMO,IYYY,IXXX)+(1.D0-RcoffMixing)*FunMONewZZ(NMO,IYYY,IXXX)  
		             ENDDO
                  ENDDO
                  ! нясыеярбкъел оепеявер йнщттхжхемрнб юкэтю
				  call CalculationCoffAlfaMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
                  
				  ! тнплхпсел бшяьхе цюплнмхйх лнкейскъпмни нпахрюкх                   
                  call FormStructureMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreGarmMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
               ENDDO
	      ENDIF

		  IF(IValueStronMixing.GT.3.AND.IValueStronMixing.LT.6) THEN
		     DO NMO=1,NumbreMO 
		        DO IXXX=1,N+2  
	               DO IYYY=NumeroMinGarmon(NMO),NumbreGarmMO(NMO)
	                  RFun(NMO,IYYY,IXXX)=RcoffMixing*FunMOOldZZD(NMO,IYYY,IXXX)+(1.D0-RcoffMixing)*FunMONewZZD(NMO,IYYY,IXXX)  
		           ENDDO
                ENDDO
                
				! нясыеярбкъел оепеявер йнщттхжхемрнб юкэтю
				call CalculationCoffAlfaMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
                  
				! тнплхпсел бшяьхе цюплнмхйх лнкейскъпмни нпахрюкх                   
                call FormStructureMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreGarmMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
		     ENDDO

		  ENDIF
          IF(IValueStronMixing.GT.6) THEN
             DO NMO=1,NumbreMO 
		        DO IXXX=1,N+2  
	               DO IYYY=NumeroMinGarmon(NMO),NumbreGarmMO(NMO)
	                  RFun(NMO,IYYY,IXXX)=RcoffMixing*FunMOOldZZDD(NMO,IYYY,IXXX)+(1.D0-RcoffMixing)*FunMONewZZDD(NMO,IYYY,IXXX)  
		           ENDDO
                ENDDO
                
				! нясыеярбкъел оепеявер йнщттхжхемрнб юкэтю
				call CalculationCoffAlfaMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
                  
				! тнплхпсел бшяьхе цюплнмхйх лнкейскъпмни нпахрюкх                   
                call FormStructureMOFull(NMO,N,H,NumbreGarmLMO(NMO),NumbreGarmMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFun,RfunLigands,AlfaCoffMO,RO1X)  
		     ENDDO
		  ENDIF
	   ENDIF

	   IF(ICalFP.EQ.0) THEN
	      ! нясыеярбкъел гюохяэ б асттеп бнкмнбшу тсмйжхи онксвеммшу мю опнлефсрнвмнл щрюое
	      call WriteBufferFunction(NumbreMO,N,NumbreGarmMO,RFun,FUNBuffer)
		  
		  ! нясыеярбкъел пюявер бпелемх
		  call CPU_TIME(Time_finishX)
          TimeCalculX=(Time_finishX-Time_startX)/60.
          WRITE(7,3500)
          WRITE(7,3400) TimeCalculX,RT1,RT2,RT3,RT4
		  
		  Iterations=NumeroIter
          ! сякнбхе бшундю хг жхйкю  
	      IF(TETMAX/FLOAT(NumbreMO).LE.EPS) THEN 
	         EXIT
  	      ENDIF
          ! гюохяшбюел опшдедсысч хрепюжхч N-2
		  FunMOOldZZDD=FunMOOldZZD
		  FunMONewZZDD=FunMONewZZD
		  ! гюохяшбюел опшдедсысч хрепюжхч N-1
		  FunMOOldZZD=FunMOOldZZ
          FunMONewZZD=FunMONewZZ
		  ! гюохяшбюел онксвеммше тсмйжхх б дсакхпсчыхе люяяхбш
          FunMOOldZZ=FunMOOld
		  FunMONewZZ=FunMONew
          IF(NumeroIter.EQ.3.AND.IndexExit.EQ.0) THEN
             ! нясыеярбкъел оепемнплхпнбйс тсмйжхи
			 DO NMO=1,NumbreMO
			    DO IYYY=NumeroMinGarmon(NMO),NumbreGarmMO(NMO) 
		           DO IXXX=3,N+2  
	                  FunMONewZZ(NMO,IYYY,IXXX)=FunMONewZZ(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMONewZZ(NMO,IYYY,2))
					  FunMOOldZZ(NMO,IYYY,IXXX)=FunMOOldZZ(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMOOldZZ(NMO,IYYY,2))
                      FunMONewZZD(NMO,IYYY,IXXX)=FunMONewZZD(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMONewZZD(NMO,IYYY,2))
					  FunMOOldZZD(NMO,IYYY,IXXX)=FunMOOldZZD(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMOOldZZD(NMO,IYYY,2))
                      FunMONewZZDD(NMO,IYYY,IXXX)=FunMONewZZDD(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMONewZZDD(NMO,IYYY,2))
                      FunMOOldZZDD(NMO,IYYY,IXXX)=FunMOOldZZDD(NMO,IYYY,IXXX)*DSQRT(RFun(NMO,IYYY,2)/FunMOOldZZDD(NMO,IYYY,2))
		           ENDDO
                   FunMONewZZ(NMO,IYYY,2)=RFun(NMO,IYYY,2)
				   FunMOOldZZ(NMO,IYYY,2)=RFun(NMO,IYYY,2)
                   FunMONewZZD(NMO,IYYY,2)=RFun(NMO,IYYY,2)
				   FunMOOldZZD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
                   FunMONewZZDD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
				   FunMOOldZZDD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
                ENDDO
                DO IYYY=1,NumeroMinGarmon(NMO) 
			       FunMONewZZ(NMO,IYYY,2)=RFun(NMO,IYYY,2)
				   FunMOOldZZ(NMO,IYYY,2)=RFun(NMO,IYYY,2)
                   FunMONewZZD(NMO,IYYY,2)=RFun(NMO,IYYY,2)
				   FunMOOldZZD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
                   FunMONewZZDD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
				   FunMOOldZZDD(NMO,IYYY,2)=RFun(NMO,IYYY,2) 
				ENDDO
			 ENDDO   
		  ENDIF

          ! гюохяшбюел опнлефсрнвмше дюммше
		  ! ме гюохяшбюел мю оепбни хрепюжхх 
		  IF(NumeroIter.NE.1) THEN 
             call WriteBufferFunction(NumbreMO,N,NumbreGarmMO,FunMONewZZ,FUNBuffer1)
		     call WriteBufferFunction(NumbreMO,N,NumbreGarmMO,FunMOOldZZ,FUNBuffer2)
             call WriteBufferFunction(NumbreMO,N,NumbreGarmMO,FunMONewZZD,FUNBuffer3)
		     call WriteBufferFunction(NumbreMO,N,NumbreGarmMO,FunMOOldZZD,FUNBuffer4)
	      ENDIF
	   ENDIF
	   
	   ! опнбепъел йюйне спюбмемхе пеьюеряъ
       IF(IndexHartreeFock.EQ.1) THEN
	      ! пюявер я налемнл 
		  IKLZExch=1
	   ENDIF
       IF(IndexHartreeFock.EQ.0) THEN
	      ! пюявер аег налемнл 
		  IKLZExch=0
	   ENDIF
	   ! тхйяхпсел рнр тюйр, врн пеьемхъ днкфмш ашрэ ябъгюммшлх
	   NCLS=1 
	ENDDO
     
    
    
    WRITE(7,*)
	WRITE(25,*)
    ! ТХМЮКЭМШИ ПЮЯВЕР 
	! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
	   ! пюявер я налемнл 
	   IKLZExch=1
	ENDIF
    IF(IndexHartreeFock.EQ.0) THEN
	   ! пюявер аег налемнл 
	   IKLZExch=0
	ENDIF
	! пюявер я лемъчыхляъ онремжхюкюл
    NCLS=0 ! ялеьхбюмхе нрясрярбсер
    ! бпелъ
    RT1=0.D0
	RT2=0.D0
	RT3=0.D0
	RT4=0.D0
	! нясыеярбкъел пюявер бпелемх
	call CPU_TIME(Time_startX) 
    
	DO NMO=1,NumbreMO
	   ! опнбепъел менаундхлн нясыеярбкъер пюявер дкъ дюммни нанкнвйх
	   IF(ISCalculZZ(NMO,1).EQ.0) THEN  
	       ! нясыеярбкъел пюявер бпелемх
		   call CPU_TIME(Time_start) 
	       
		   call MAINMO(IndexHartreeFock,IndexExit,NN(NMO),NCLS,IreshimSO,Iterations+1,NMO,NumbreMO,N,NtypIndexCoffF,NtypIndexCoffG,RcoffA,RcoffB,IQ,NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon,LgarmonicMO,EPSXR,H,RFun,FunMOOld,FunMONew,R,RO1X,RO2X,RO3X,Ncrystall,TETMAX,TETA,E3,E4,E5,DeltaEnegry,RnormMOXFII,RnormMOXF,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,ENN0,NumreISzam,ISCalculZZ,IKLZExch,IRETURN,RparamterCalcul,RTimeCalcul,IndexMLA,IndexMLB,CkkGarmonik,AlfaCoffMO,RfunLigands,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,LmaxKGarm,LmaxGarm,NumbreSxemaCalculation,RFILEPoltens,NumbreRegion,EnergyZeroZ,EnergyRegion,IndexRegionHMIN,EnergyRegionHMIN)
	       
		   ! нясыеярбкъел пюявер бпелемх
		   call CPU_TIME(Time_finish)
           TimeCalcul=(Time_finish-Time_start)/60.
		   ! опнлефсрнвмше пегскэрюрш пюяверю
	   	   WRITE(7,3300) Iterations+1,NMO,RparamterCalcul(1),RparamterCalcul(2),RparamterCalcul(3),RparamterCalcul(4),RparamterCalcul(5),TimeCalcul+RTimeCalcul(5),RTimeCalcul(1),RTimeCalcul(3),RTimeCalcul(4),RTimeCalcul(2) 
	       RT1=RT1+RTimeCalcul(1)
		   RT2=RT2+RTimeCalcul(3)
		   RT3=RT3+RTimeCalcul(4)
		   RT4=RT4+RTimeCalcul(2)
           ! жхйк он кхцюмдюл
	       DO IXZZ1=1,NumbreLigand(NMO)
              ! жхйк он тсмйжхъл кхцюмдю
		      DO IXZY1=1,NFunLigands(NMO,IXZZ1)
	             WRITE(25,3900) Iterations+1,NMO,NLigands(NMO,IXZZ1),NumbreFunctionLig(NMO,IXZZ1,IXZY1),AlfaCoffMO(NMO,NLigands(NMO,IXZZ1),NumbreFunctionLig(NMO,IXZZ1,IXZY1))
		      ENDDO
	       ENDDO
	   ENDIF
	ENDDO
        
    ! нясыеярбкъел пюявер бпелемх
	call CPU_TIME(Time_finishX)
    TimeCalculX=(Time_finishX-Time_startX)/60.
	WRITE(7,3500)
    WRITE(7,3400) TimeCalculX,RT1,RT2,RT3,RT4

  

    ! сдюкъел мемсфмше люяяхбш 
    IF(NumbreSxemaCalculation.EQ.1) THEN
       deallocate(RFILEPoltens,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'DHFMO'
          write(*,*) 'THE FILE "RFILEPoltens" IS NOT REMOVED FROM MEMORY'
          read(*,*)
	      stop 
       endif 
	ENDIF

	deallocate(NumbreGPOTMO,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "NumbreGPOTMO" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(NumbreGPOTMORPOT,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "NumbreGPOTMORPOT" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(IFPOTIndex,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IFPOTIndex" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(IFPOTZona,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IFPOTZona" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
	   deallocate(IGPOTIndex,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'DHFMO'
          write(*,*) 'THE FILE "IGPOTIndex" IS NOT REMOVED FROM MEMORY'
          read(*,*)
	      stop 
       endif
	   deallocate(IGPOTZona,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'DHFMO'
          write(*,*) 'THE FILE "IGPOTZona" IS NOT REMOVED FROM MEMORY'
          read(*,*)
	      stop 
       endif  
	ENDIF
    deallocate(FunMOOld,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMOOld" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMONew,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMONew" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMONewZZD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMONewZZD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMONewZZDD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMONewZZDD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMOOldZZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMOOldZZ" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMOOldZZD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMOOldZZD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMOOldZZDD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMOOldZZDD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(FunMONewZZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "FunMONewZZ" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	! опнбепъел йюйне спюбмемхе пеьюеряъ
    IF(IndexHartreeFock.EQ.1) THEN
       deallocate(IZONAMASSIVG,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'DHFMO'
          write(*,*) 'THE FILE "IZONAMASSIVG" IS NOT REMOVED FROM MEMORY'
          read(*,*)
	      stop 
       endif 
	ENDIF
	deallocate(IZONAMASSIV,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IZONAMASSIV" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 



      
	! МНПЛХПНБМЙЮ ОНКСВЕММШУ БНКМНБШУ ТСМЙЖХИ Х ГЮОХЯЭ НДМНЩКЕЙРПНММШУ ЩМЕПЦХИ
    DO IXXX=1,NumbreMO 
       ! мнплхпнбнвмши йнщттхжхемр
	   VAL=1.D0/DSQRT(RFun(IXXX,NumeroMinGarmon(IXXX),2)) 
       EnergyXF=RFun(IXXX,NumeroMinGarmon(IXXX),1)
	   DO IYYY=1,NumbreGarmMO(IXXX)
          ! мнплхпнбйю тсмйжхх
		  DO J=1,N
             RFun(IXXX,IYYY,J+2)=RFun(IXXX,IYYY,J+2)*VAL
          ENDDO
		  ! гюохяэ ндмнщкейрпнммни щмепцхх
	      RFun(IXXX,IYYY,1)=-EnergyXF
       ENDDO
    ENDDO
	  
	
	! бшбнд опнлефсрнвмшу пегскэрюрнб
	! ндмнщкейрпнммшу щмепцхи х япедмху пюдхсянб лнкейскъпмшу нпахрюкеи  
    ! нясыеярбкъел пюявер япедмху пюдхсянб
    WRITE(6,3000)
    DO IXXX=1,NumbreMO
       Rad=MediumRadius(IXXX,NumbreGarmMO,H,N,RFun,R,RO1X)
	   WRITE(6,3100)  NN(IXXX),LBZ(IABS(ML(IXXX))+1),RFun(IXXX,NumeroMinGarmon(IXXX),1),Rad,TETA(IXXX),RFun(IXXX,NumeroMinGarmon(IXXX),2) 
	ENDDO 

   

	
	 
	
    ! нясыеярбкъел пюявер онкмни щмепцхх йнмтхцспюжхх
	call TOTALXF(IndexHartreeFock,Nnuclei,Z0,Z,COORRR,COOAAA,NumbreMO,NtypIndexCoffF,NtypIndexCoffG,NN,ML,RcoffA,RcoffB,IQ,NumbreGarmMO,LgarmonicMO,N,H,Ncrystall,RFun,R,RO1X,IndexMLA,IndexMLB,CkkGarmonik,NomeroFPOT)	 
   
      
    


    ! сдюкемхе люяяхбнб хг оюлърх 
	deallocate(NcrystallXX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "NcrystallXX" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(NumbreRegion,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "NumbreRegion" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
    deallocate(EnergyZeroZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyZeroZ" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(EnergyZeroZD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyZeroZD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(IndexRegionHMIN,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IndexRegionHMIN" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(EnergyRegionHMIN,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyRegionHMIN" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(IndexRegionHMIND,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IndexRegionHMIND" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(EnergyRegionHMIND,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyRegionHMIND" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif   
    deallocate(EnergyRegionD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyRegionD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
    deallocate(EnergyRegion,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "EnergyRegion" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(RPOTX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "RPOTX" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(AlfaCoffMO,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "AlfaCoffMO" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(IndexMLA,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IndexMLA" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
    deallocate(IndexMLB,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "IndexMLB" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(CkkGarmonik,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "CkkGarmonik" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(MTipMOML,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "MTipMOML" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(RTimeCalcul,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "RTimeCalcul" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(RparamterCalcul,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "RparamterCalcul" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	
	deallocate(Ncrystall,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "Ncrystall" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(F,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "F" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(TETA,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "TETA" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(E3,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "E3" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 
	deallocate(E4,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "E4" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(E5,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "E5" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(DeltaEnegry,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "DeltaEnegry" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
	deallocate(NumeroMinGarmon,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'DHFMO'
       write(*,*) 'THE FILE "NumeroMinGarmon" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif
    return
  end  subroutine DHFMO

  
!! PROGRAM OF SOLVING THE HATRY-FOCK EQUATION FOR A MOLECULAR SHELL FROM CONFIGURATION
══!! DESCRIPTION OF INTRODUCTORY PARAMETERS
══!! IndexHartreeFock-PARAMETER INDICATING WHAT EQUATION IS DECIDED
══!! IndexHartreeFock = 0-HARTRI EQUATION
══!! IndexHartreeFock = 1-HARTRI-FOCA EQUATION
══!! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
══!! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
══!! IndexExit = 1-HAS FAILED THROUGH READING INFORMATION FROM EMERGENCY FILES
══!! NNEL-MAIN QUANTUM NUMBER OF MOLECULAR ORBITALS FOR WHICH CALCULATION OF PARTIAL HARMONICS IS IMPLEMENTED
══!! NCLS-KEY INDICATORING THE PROCESS OF MIXING THE DECISION RECEIVED ON THE PREVIOUS STAGE WITH THE FOLLOWING
══!! IreshimSO-PARAMETER INDICATING THE TYPE OF HARMONIZATION OF MOLECULAR ORBITALS
══!! IreshimSO = 0-STANDARD MODE OF CODE
══!! IreshimSO = 1- "STRENGTHENING" HARMONIZATION MODE THE BINDING COEFFICIENTS ARE EXPRESSED
══!! NumeroIter-ITERRATION NUMBER
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH WE MAKE CALCULATION
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoffF (NumbreMO, NumbreMO) - NUMBER OF TYPES OF COEFFICIENTS OF DIRECT CULON INTERACTION
══!! NtypIndexCoffG (NumbreMO, NumbreMO) - NUMBER OF TYPES OF COEFFICIENTS OF EXCHANGE COULON INTERACTION
══!! RcoffA (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION FACTORS
══!! RcoffB (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL
══!! NumbreGarmLMO (NumbreMO) - MASSIVE OF HARMONIC NUMBER IN MOLECULAR ORBITAL TO LgarmonicMO (3, NumbreMO) INCLUDING IT
══!! NumeroMinGarmon (NumbreMO) - MASSIVE MINIMUM NUMBER OF HARMONICS EXCELLENT FROM ZERO IN MOLECULAR ORBITALS
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! EPS-accuracy of calculation (DETERMINATION OF SINGLE-ELECTRON ENERGIES IN SOLUTION OF EQUATIONS)
══!! H-step
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! FunMOOld (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) -mass of values ??of radial parts of harmonics of molecular orbitals configuration previous iteration
══!! FunMONew (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of harmonics of the molecular orbitals of the configuration not mixed with the solution of the previous iteration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RO2X (Npoint) - MASSIVE VALUES OF THE SECOND DERIVATIVE ON THE NEW VARIABLE
══!! RO3X (Npoint) - MASSIVE VALUES OF THIRD DERIVATIVE ON THE NEW VARIABLE
══!! FF, F-AUXILIARY PARAMETERS FOR CALCULATING ANGULAR HARMONICS
══!! Ncrystall (NumbreMO, NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -MASIVE of the CRYSTAL FIELD OPERATOR (FIELD CREATED BY NUCLEAR IN THE BEGINNING OF THE COORDINATE)
══!! TETMAX-MAXIMUM DEVIATION (TOTAL) OF TWO SEQUENTIAL ITERATIONS BETWEEN YOUR SELF (RECEIVED FUNCTIONS IN THESE ITERATIONS)
══!! TETA (NumbreMO) -MASSIVE OF DEVICES FOR THE FUNCTIONS OF EACH MOLECULAR SHELL
══!! E3 (NumbreMO) -MASSIVE OF SINGLE-ELECTRONIC ENERGIES RECEIVED ON D-1 STEP
══!! E4 (NumbreMO) -MASSIVE OF SINGLE-ELECTRONIC ENERGIES OBTAINED ON D-2 STEP
══!! E5 (NumbreMO) -MASSIVE OF SINGLE-ELECTRONIC ENERGIES RECEIVED AT D-3 STEP
══!! DeltaEnegry (NumbreMO, 2) -MASSIVE DISTORTION OF SINGLE-ELECTRONIC ENERGY UNDER THE SOLUTION OF THE HOMOGENEOUS AND INHOMOGENEOUS SYSTEM OF EQUATIONS
══!! RnormMOXFII-NORMAL FACTOR FOR THIS MOLECULAR ORBITAL FOR THIS ITERATION
══!! RnormMOXF (NumbreMO) -MASSIVE OF THE NORMATIVE COEFFICIENTS OF MOLECULAR ORBITALS (CRITERION OF ORDERING CORRESPONDENCE OF OBTAINED FUNCTIONS)
══!! NumbreGarmMOZero (NumbreMO) - HARMONIC NUMBER OF MASSES IN THE MOLECULAR ORBITAL OF THE NULL APPROXIMATION
══!! LgarmonicMOZero (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS OF THE NULE APPROXIMATE
══!! RfunZero (NumbreMO, NumbreGarmMOZero (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the zeroth approximation configuration
══!! ENN0 - LOWER ENERGY BORDER
══!! NumreISzam (NumbreMO) -MASSIVE NUMBER OF SHELLS OF EQUIVALENT DATA
══!! ISCalculZZ (NumbreMO, NumbreMO) -MASIVE OF EQUIVALENT DATA SHEET NUMBERS
══!! IKLZERO-KEY APPROACHING THE TYPE OF THE NULL APPROACH
══!! IKLZERO = 0- CALCULATION WITHOUT EXCHANGE
══!! IKLZERO = 1- CALCULATION OF EXCHANGE
  !! IRETURN-KEY APPROPRIATE TO THE RESULT OF THE WORK OF THE PROGRAM
══!! IRETURN = 0-REMOVAL SOLVED AND ROOF FOUND
══!! IRETURN = 1-EQUATION NOT SOLVED ROOTS NOT AVAILABLE
══!! RparamterCalcul () - OUTPUT CALCULATION PARAMETERS SINGLE-ELECTRONIC ENERGIES AND MIXING RATIO
══!! RTimeCalcul (5) -MASSIVE IN WHICH TIME OF CALCULATING THE MATRIX OF CULON INTERACTION AND SYSTEM SOLUTION
══!! IndexMLA (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE DIRECT PART OF CULON INTERACTION
══!! IndexMLB (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO)) - SPHERICAL HARMONIC array
══!! NumbreIndexLigandDirect (NumbreMO) -MASSIVE NUMBER OF INDICES (NUMBER OF COMPONENTS IN THE SUM) FOR DIRECT INTERACTION
══!! IndexNFunLigandDirect (NumbreMO, NumbreIndexLigandDirect (NumbreMO), 4) -MASSIVE INDICATOR WHICH SUGGESTED (COEFFICIENTS IN THE AMOUNT) SUMMIT FOR THIS MOLECULAR ORBITAL
══!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
══!! DCoulombFull (Npoint, IndexNFunLigandDirect, NtypIndexCoff, NumbreMO, Ngarmon1, Ngarmon2, NumbreMO) is the mass of the direct part of the Coulomb interaction of molecular orbitals (a part calculated on the functions of the ligand)
══!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
══!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
══!! NomeroFPOT (NumbreMO) -MASSIVE NUMBER OF FILES
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! NomeroGPOT (NumbreMO, NumbreMO) -MASSIVE IN WHICH THE FILE NUMBERS ARE FOR RECORDING EXCHANGE POTENTIALS
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! NumbreGPOTMORPOT (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION FOR THE MASSIVE OF RADIAL PARTS OF POTENTIALS
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! LmaxKGarm-MAXIMUM VALUE LgarmonicMO (3, NumbreMO)
══!! LmaxGarm-MAXIMUM VALUE LgarmonicMO (2, NumbreMO)
══!! NumbreSxemaCalculation-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculation = 1-CAPACITY RECORDING IS IN ARMIVE
══!! NumbreSxemaCalculation = 2-CAPACITY RECORDING IS IN FILE
══!! RFILEPoltens (Npoint, Izonmax) -FILE FOR CAPACITY BUILDING
══!! NumbreRegion (NumbreMO) -MASSIVE NUMBER OF REGIONS IN WHICH A ROOT REDUCTION IS CARRIED OUT
══!! EnergyZeroZ (NumbreMO, NNEL + 1,2) -MASSIVE OF THE ROOTS OF A HOMOGENEOUS EQUATION
══!! EnergyRegion (NumbreMO, NNEL + 1,2) -MASSIVE OF AREAS IN WHICH STEP DECREASES
══!! IndexRegionHMIN (NumbreMO, NNEL + 1) -MASSIVE INDEX OF INDICATORS INDEXING STEP FOR DETERMINING BORDERS FOR THIS MINIMUM
══!! EnergyRegionHMIN (NumbreMO, NNEL + 1,2) -MASSIVE OF INITIAL STEPS FOR DETERMINATION OF THE BOUNDARIES OF MINIMUMS 
  subroutine MAINMO(IndexHartreeFock,IndexExit,NNEL,NCLS,IreshimSO,NumeroIter,NMO,NumbreMO,Npoint,NtypIndexCoffF,NtypIndexCoffG,RcoffA,RcoffB,IQ,NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon,LgarmonicMO,EPS,H,RFunMO,FunMOOld,FunMONew,R,RO1X,RO2X,RO3X,Ncrystall,TETMAX,TETA,E3,E4,E5,DeltaEnegry,RnormMOXFII,RnormMOXF,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,ENN0,NumreISzam,ISCalculZZ,IKLZERO,IRETURN,RparamterCalcul,RTimeCalcul,IndexMLA,IndexMLB,CkkGarmonik,AlfaCoffMO,RfunLigands,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,LmaxKGarm,LmaxGarm,NumbreSxemaCalculation,RFILEPoltens,NumbreRegion,EnergyZeroZ,EnergyRegion,IndexRegionHMIN,EnergyRegionHMIN)
    use msde, only:SolutionSystemDifferentialEquations 
	implicit none
    integer::IndexHartreeFock,IndexExit,NNEL,NCLS,NumeroIter,NMO,NumbreMO,Npoint,IKLZERO,IRETURN,LmaxKGarm,LmaxGarm,NumbreSxemaCalculation,IreshimSO  
    real(8)::H,EPS,TETMAX,ENN0,RnormMOXFII
    integer,dimension(:)::IQ,NumeroMinGarmon,NumbreGarmMO,NumbreGarmLMO,NumbreGarmMOZero,NumreISzam,NumbreLigand,NomeroFPOT,NumbreRegion
    integer,dimension(:,:)::NtypIndexCoffF,NtypIndexCoffG,LgarmonicMO,LgarmonicMOZero,ISCalculZZ,NLigands,NFunLigands,NumbreGPOTMO,NumbreGPOTMORPOT,IndexRegionHMIN
    integer,dimension(:,:,:)::NumbreFunctionLig
	integer(1),dimension(:,:,:,:)::IFPOTIndex,IGPOTIndex
	integer,dimension(:,:,:,:)::IndexMLA,IndexMLB,IFPOTZona,IGPOTZona
	real(8),dimension(:)::R,RO1X,RO2X,RO3X,TETA,E3,E4,E5,RnormMOXF,RparamterCalcul,RTimeCalcul
	real(8),dimension(:,:)::DeltaEnegry,RFILEPoltens
    real(8),dimension(:,:,:)::RcoffA,RcoffB,RFunMO,FunMOOld,FunMONew,RfunZero,AlfaCoffMO,EnergyZeroZ,EnergyRegion,EnergyRegionHMIN
    real(8),dimension(:,:,:,:)::Ncrystall,RfunLigands
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::NI,LI,QI,IIKL,IIOP
    integer::ierr,IXXD,J,MAXCNV,IIXYT
	integer(2)::HZ,M,S,HS
	real(4)::TimeSXZ,TimeFXZ,TimeSXZX,TimeFXZX
    real(8)::EpsLagram,T,C,TXBX,TXEX,RR 
    real(8),allocatable,dimension(:)::RPOTXZ
	real(8),allocatable,dimension(:,:)::Qz
    real(8),allocatable,dimension(:,:,:)::Pz
	  
	556 FORMAT(' N=',I4,'    Shell=',I2,'    E=',F10.4,'    Teta=',D10.3) 
    
    ! WE ALLOW MEMORY MASSIVES
════ !! Pz (NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -LOCAL PART OF THE MATRIX DIFFERENTIAL EQUATION FOR THIS MOLECULAR SHELL
════ !! Qz (NumbreGarmMO (NumbreMO), Npoint + 2) -NONLOCAL PART OF THE MATRIX DIFFERENTIAL EQUATION FOR THIS MOLECULAR SHELL    allocate(Pz(NumbreGarmLMO(NMO),NumbreGarmLMO(NMO),Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'MEMORY ON THE FILE "Pz" IS NOT SELECTED'
       READ(*,*)
       stop 
    endif 
    allocate(Qz(NumbreGarmLMO(NMO),Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'MEMORY ON THE FILE "Qz" IS NOT SELECTED'
       READ(*,*)
       stop 
    endif 
	allocate(RPOTXZ(Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'MEMORY ON THE FILE "RPOTXZ" IS NOT SELECTED'
       READ(*,*)
       stop 
    endif 
    

  
    !! WE MAKE THE CALCULATION OF TIME
	call CPU_TIME(TimeSXZ) 
	      	  
    !! PROGRAM OF CALCULATION OF Pz AND Qz-FUNCTIONS NECESSARY FOR SOLVING DIFFERENTIAL EQUATION
    call Calcul_Pz_Qz_MO(IndexHartreeFock,NumeroIter,NMO,NumbreMO,Npoint,NtypIndexCoffF,NtypIndexCoffG,RcoffA,RcoffB,IQ,NumbreGarmMO,NumbreGarmLMO,LgarmonicMO,H,RFunMO,R,RO1X,RO2X,RO3X,Ncrystall,Pz,Qz,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,IKLZERO,IndexMLA,IndexMLB,CkkGarmonik,RTimeCalcul,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,NumbreSxemaCalculation,RFILEPoltens,LmaxKGarm,LmaxGarm) 
 	
    !! WE MAKE THE CALCULATION OF TIME
	call CPU_TIME(TimeFXZ) 
   
    RTimeCalcul(1)=DBLE((TimeFXZ-TimeSXZ)/60.)
	
    !! WE MAKE THE CALCULATION OF TIME
	call CPU_TIME(TimeSXZX) 
	  
    !! We Solve the Differential Equation by the Numeroff Vector Method
════!! Pz + EPS = R6-MASSIVE of VALUE FUNCTIONS P (local part) in the equation Y '' = P * Y + q
════!! Qz-MASSIVE of VALUES OF FUNCTIONS q (not local part) in the equation Y "= P * Y + q
    call SolutionSystemDifferentialEquations(IndexExit,NNEL,NumeroIter,NMO,NumeroMinGarmon(NMO),NumbreGarmLMO(NMO),NCLS,IreshimSO,Npoint,H,EPS,TETMAX,TETA,MAXCNV,E3,E4,E5,DeltaEnegry,RnormMOXFII,RnormMOXF,Pz,Qz,RFunMO,FunMOOld,FunMONew,RO1X,ENN0,NumreISzam,ISCalculZZ,IKLZERO,IRETURN,RparamterCalcul,NumbreGarmMO(NMO),NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RfunLigands,AlfaCoffMO,NumbreRegion,EnergyZeroZ,EnergyRegion,IndexRegionHMIN,EnergyRegionHMIN)
    
   
	!! WE MAKE THE CALCULATION OF TIME
	call CPU_TIME(TimeFXZX) 
    RTimeCalcul(2)=DBLE((TimeFXZX-TimeSXZX)/60.)

    !! CHECKING THE RESULT OF THE WORK
	IF(IRETURN.NE.1) THEN
       WRITE(*,556) NumeroIter,NMO,RFunMO(NMO,NumeroMinGarmon(NMO),1),TETA(NMO)
    ENDIF  
      


    ! сдюкемхе люяяхбнб хг оълърх 
	deallocate(RPOTXZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'THE FILE "RPOTXZ" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
       stop 
    endif
    deallocate(Pz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'THE FILE "Pz" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
       stop 
    endif
    deallocate(Qz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAINMO'
       write(*,*) 'THE FILE "Qz" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
       stop 
    endif
   
    return
  end subroutine MAINMO



!! The subgraph sets out the matrices Pz and Qz for solving the matrix differential equation for this shell of the molecule
══!! IndexHartreeFock-PARAMETER INDICATING WHAT EQUATION IS DECIDED
══!! IndexHartreeFock = 0-HARTRI EQUATION
══!! IndexHartreeFock = 1-HARTRI-FOCA EQUATION
══!! NumeroIter-ITERRATION NUMBER
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH WE MAKE CALCULATION
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoffF (NumbreMO, NumbreMO) - NUMBER OF TYPES OF COEFFICIENTS OF DIRECT CULON INTERACTION
══!! NtypIndexCoffG (NumbreMO, NumbreMO) - NUMBER OF TYPES OF COEFFICIENTS OF EXCHANGE COULON INTERACTION
══!! RcoffA (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION FACTORS
══!! RcoffB (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL
══!! NumbreGarmLMO (NumbreMO) - MASSIVE OF HARMONIC NUMBER IN MOLECULAR ORBITAL TO LgarmonicMO (3, NumbreMO) INCLUDING IT
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! H-step
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RO2X (Npoint) - MASSIVE VALUES OF THE SECOND DERIVATIVE ON THE NEW VARIABLE
══!! RO3X (Npoint) - MASSIVE VALUES OF THIRD DERIVATIVE ON THE NEW VARIABLE
══!! FF, F-AUXILIARY PARAMETERS FOR CALCULATING ANGULAR HARMONICS
══!! Ncrystall (NumbreMO, NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -MASIVE of the CRYSTAL FIELD OPERATOR (FIELD CREATED BY NUCLEAR IN THE BEGINNING OF THE COORDINATE)
══!! Pz (NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -LOCAL PART OF THE MATRIX DIFFERENTIAL EQUATION FOR THIS MOLECULAR SHELL
══!! Qz (NumbreGarmMO (NumbreMO), Npoint) -NONLOCAL PART OF THE MATRIX DIFFERENTIAL EQUATION FOR THIS MOLECULAR SHELL
══!! NumbreGarmMOZero (NumbreMO) - HARMONIC NUMBER OF MASSES IN THE MOLECULAR ORBITAL OF THE NULL APPROXIMATION
══!! LgarmonicMOZero (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS OF THE NULE APPROXIMATE
══!! RfunZero (NumbreMO, NumbreGarmMOZero (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the zeroth approximation configuration  ! IKLZERO-йкчв сйюгшбючыхи мю рхо мскхбнцн опхакхфемхъ 
  !! IKLZERO = 0- CALCULATION WITHOUT EXCHANGE
══!! IKLZERO = 1- CALCULATION OF EXCHANGE
══!! IndexMLA (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF TOMORROW PROJECTIONS (MOLECULAR ORBITALS) FOR CALCULATING THE DIRECT PART OF COULON INTERACTION
══!! IndexMLB (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO)) - SPHERICAL HARMONIC array
══!! NumbreIndexLigandDirect (NumbreMO) -MASSIVE NUMBER OF INDICES (NUMBER OF COMPONENTS IN THE SUM) FOR DIRECT INTERACTION
══!! IndexNFunLigandDirect (NumbreMO, NumbreIndexLigandDirect (NumbreMO), 4) -MASSIVE INDICATOR WHICH SUGGESTED (COEFFICIENTS IN THE AMOUNT) SUMMIT FOR THIS MOLECULAR ORBITAL
══!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
══!! DCoulombFull (Npoint, IndexNFunLigandDirect, NtypIndexCoff, NumbreMO, Ngarmon1, Ngarmon2, NumbreMO) is the mass of the direct part of the Coulomb interaction of molecular orbitals (a part calculated on the functions of the ligand)
══!! RTimeCalcul (4) -MASIVE TO WHICH TIME OF CALCULATING A DIRECT AND EXCHANGE PART OF CAPACITY IS RECORDED
══!! NomeroFPOT (NumbreMO) -MASSIVE NUMBER OF FILES
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! NomeroGPOT (NumbreMO, NumbreMO) -MASSIVE IN WHICH FILE NUMBERS ARE FOR RECORDING EXCHANGE POTENTIALS
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! NumbreGPOTMORPOT (NumbreMO, NumbreMO) - MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION FOR MASSIVE RADIAL PARTS OF POTENTIALS
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! NumbreSxemaCalculation-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculation = 1-CAPACITY RECORDING IS IN ARMIVE
══!! NumbreSxemaCalculation = 2-CAPACITY RECORDING IS IN FILE
══!! RFILEPoltens (Npoint, Izonmax) -MASSIVE FOR CAPACITY BUILDING
══!! LmaxKGarm-MAXIMUM VALUE LgarmonicMO (3, NumbreMO)
══!! LmaxGarm-MAXIMUM VALUE LgarmonicMO (2, NumbreMO)
  subroutine Calcul_Pz_Qz_MO(IndexHartreeFock,NumeroIter,NMO,NumbreMO,Npoint,NtypIndexCoffF,NtypIndexCoffG,RcoffA,RcoffB,IQ,NumbreGarmMO,NumbreGarmLMO,LgarmonicMO,H,RFunMO,R,RO1X,RO2X,RO3X,Ncrystall,Pz,Qz,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,IKLZERO,IndexMLA,IndexMLB,CkkGarmonik,RTimeCalcul,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,NumbreSxemaCalculation,RFILEPoltens,LmaxKGarm,LmaxGarm) 
    implicit none
    integer::IndexHartreeFock,NumeroIter,NMO,NumbreMO,Npoint,IKLZERO,NumbreSxemaCalculation,LmaxKGarm,LmaxGarm
    real(8)::H
    integer,dimension(:)::IQ,NumbreGarmMO,NumbreGarmLMO,NumbreGarmMOZero,NomeroFPOT
    integer,dimension(:,:)::NtypIndexCoffF,NtypIndexCoffG,LgarmonicMO,LgarmonicMOZero,NumbreGPOTMO,NumbreGPOTMORPOT
    integer(1),dimension(:,:,:,:)::IFPOTIndex,IGPOTIndex
	integer,dimension(:,:,:,:)::IndexMLA,IndexMLB,IFPOTZona,IGPOTZona
	real(8),dimension(:)::R,RO1X,RO2X,RO3X,RTimeCalcul
    real(8),dimension(:,:)::Qz,RFILEPoltens
    real(8),dimension(:,:,:)::RcoffA,RcoffB,RFunMO,Pz,RfunZero
    real(8),dimension(:,:,:,:)::Ncrystall
    real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::NI,LI,QI,IIL1,IIL2,IMO1,IMO2,INMOIDG,Indexcgkl
    integer::ierr,IISD,IXXD,J,MAXCNV,IOX
    real(4)::TimeSXZA,TimeFXZA
	real(8)::P,T,B,C,D,RO21,RO31,ROIF,Vmo,Vcfmo
    real(8),allocatable,dimension(:,:,:)::DCoulomb,XCoulomb
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8),allocatable,dimension(:)::RPOT,RPOTRez,RPOTDF
  
	 
    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(DCoulomb(NumbreGarmLMO(NMO),NumbreGarmLMO(NMO),Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'MEMORY ON THE FILE "DCoulomb" IS NOT SELECTED'
       READ(*,*)
	   stop 
    endif 
    allocate(XCoulomb(NumbreGarmLMO(NMO),NumbreGarmLMO(NMO),Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'MEMORY ON THE FILE "XCoulomb" IS NOT SELECTED'
       READ(*,*)
	   stop 
    endif 
    
	!бшдекъел оюлърэ дкъ люяяхбнб дкъ бяонлнцюрекэмшу люяяхбнб
    allocate(RPOT(Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'MEMORY ON THE FILE "RPOT" IS NOT SELECTED'
       read(*,*)
       stop 
    endif
	allocate(RPOTRez(Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'MEMORY ON THE FILE "RPOTRez" IS NOT SELECTED'
       read(*,*)
       stop 
    endif 
	allocate(RPOTDF(Npoint),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'MEMORY ON THE FILE "RPOTDF" IS NOT SELECTED'
       read(*,*)
       stop 
    endif   
      


   
   !! CALCULATION OF THE DIRECT AND EXCHANGE PART OF CAPACITY
    IF(NumeroIter.NE.1) THEN
        
		!! WE MAKE THE CALCULATION OF TIME
	    call CPU_TIME(TimeSXZA) 

		   !! CALCULATION OF THE DIRECT PART OF CAPACITY
           WRITE(*,*) 'DERECT Coulomb Interaction'
		   Indexcgkl=0
		   !! CYCLE ON MOLECULAR ORBITALS
		   DO INMOIDG=1,NumbreMO
			  !! CHECKING THERE IS INTERACTION BETWEEN THE DATA SHELLS (MOLECULAR ORBITALS)
              IF(NtypIndexCoffF(NMO,INMOIDG).NE.0) THEN
                 !! INDEX INDICATING ON THE NECESSITY OF ENERGY
                 Indexcgkl=Indexcgkl+1
			    !! WE MAKE CALCULATION OF POTENTIALS FOR THE DIRECT PART
			     WRITE(*,*) 'Matrix_Potential_Direct NMO=',INMOIDG 
	             call Matrix_Potential_Direct_Coulomb_Interaction(INMOIDG,Npoint,H,LmaxGarm,NomeroFPOT,LgarmonicMO,IFPOTIndex,IFPOTZona,RFunMO,R,RO1X,RPOTRez,NumbreSxemaCalculation,RFILEPoltens)
       
		         IMO1=0
                 DO IIL1=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
                    IMO1=IMO1+1
	                IMO2=IMO1-1 
	                DO IIL2=IIL1,LgarmonicMO(3,NMO) 
	                   IMO2=IMO2+1
			           WRITE(*,*) INMOIDG,IIL1,IIL2
                       ! опнбепъел менаундхлнярэ гюмскемхъ
					   IF(Indexcgkl.EQ.1) THEN
                          DCoulomb(IMO1,IMO2,1:Npoint)=0.D0
					   ENDIF
			           !! CALCULATION OF THE DIRECT PART
	                   call Matrix__Direct_Coulomb_Interaction_IJ(NMO,INMOIDG,IMO1,IIL1,IMO2,IIL2,NumbreMO,NtypIndexCoffF,H,RcoffA,IQ,NumbreGarmMO,LgarmonicMO,Npoint,RFunMO,DCoulomb,IndexMLA,CkkGarmonik,NomeroFPOT,IFPOTIndex,IFPOTZona,RPOT,RPOTRez,RPOTDF,NumbreSxemaCalculation,RFILEPoltens)
			        ENDDO
                 ENDDO
              ENDIF
		   ENDDO
 
           !! FILLING ELEMENTS BELOW DIAGONATED 
           IMO1=0
           DO IIL1=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
              IMO1=IMO1+1
	          IMO2=IMO1-1 
	          DO IIL2=IIL1,LgarmonicMO(3,NMO) 
	             IMO2=IMO2+1
		         WRITE(*,*) IIL1,IIL2
		   	     IF(IMO1.NE.IMO2) THEN
			        !! WE WRITE THE NON-IGONAL ELEMENT CORRESPONDING TO THIS DATA
				    DCoulomb(IMO2,IMO1,1:Npoint)=DCoulomb(IMO1,IMO2,1:Npoint)
			     ENDIF	
              ENDDO
           ENDDO  

		!! WE MAKE THE CALCULATION OF TIME
		call CPU_TIME(TimeFXZA) 
   
        RTimeCalcul(3)=DBLE((TimeFXZA-TimeSXZA)/60.)

        
		!! WE MAKE THE CALCULATION OF TIME
	    call CPU_TIME(TimeSXZA)
		 
	    !! CHECK THE NECESSITY OF THE CALCULATION OF THE EXCHANGE PART
		IF(IKLZERO.EQ.1) THEN 
           ! пюявер налеммни вюярх онремжхюкю
		   WRITE(*,*) 'EXCHANG Coulomb Interaction'
		   Indexcgkl=0
		   ! жхйк он лнкейскъпмшл нпахрюкъл 
		   DO INMOIDG=1,NumbreMO 
		      ! опнбепъел ясыеярбсер кх бгюхлндеиярбхе лефдс дюммшлх нанкнвйюлх (лнкейскъпмшлх нпахрюкълх)
              IF(NtypIndexCoffG(NMO,INMOIDG).NE.0) THEN 
		         ! хмдейя сйюгшбючыхи мю менаундхлнярэ гюмскемхъ
                 Indexcgkl=Indexcgkl+1
			     WRITE(*,*) 'Matrix_Potential_Exchange NMO=',NMO,INMOIDG
	             call Matrix_Potential_Exchange_Coulomb_Interaction(NMO,INMOIDG,Npoint,H,LmaxKGarm,LmaxGarm,NomeroFPOT,LgarmonicMO,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RFunMO,R,RO1X,RPOTRez,NumbreSxemaCalculation,RFILEPoltens)
	            
				 IMO1=0
                 DO IIL1=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
                    IMO1=IMO1+1
	                IMO2=0 
	                DO IIL2=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
	                   IMO2=IMO2+1
			           WRITE(*,*) INMOIDG,IIL1,IIL2
					   ! опнбепъел менаундхлнярэ гюмскемхъ
					   IF(Indexcgkl.EQ.1) THEN
                          XCoulomb(IMO1,IMO2,1:Npoint)=0.D0
					   ENDIF
                       ! пюявер налеммни вюярх
			           call Matrix__Exchange_Coulomb_Interaction_IJ(NMO,INMOIDG,IMO1,IIL1,IMO2,IIL2,NumbreMO,NtypIndexCoffG,H,RcoffB,IQ,NumbreGarmMO,LgarmonicMO,Npoint,RFunMO,XCoulomb,IndexMLB,CkkGarmonik,NomeroFPOT,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RPOT,RPOTRez,RPOTDF,NumbreSxemaCalculation,RFILEPoltens)
				    ENDDO
                 ENDDO
			  ENDIF
           ENDDO
		ENDIF

		! нясыеярбкъел пюявер бпелемх
		call CPU_TIME(TimeFXZA) 
   
        RTimeCalcul(4)=DBLE((TimeFXZA-TimeSXZA)/60.)
 	     
	   ELSE
	    ! пюявер йскнмнбяйнцн бгюхлндеиярбхъ дкъ оепбни хрепюжх    
        ! пюявер нясыеярбкъеряъ мю тсмйжхъу мскхбнцн опхакхфемхъ 
        ! нясыеярбкъел пюявер бпелемх
	    call CPU_TIME(TimeSXZA) 
		! опнбепъел рхо спюбмемхъ
		IF(IndexHartreeFock.EQ.1) THEN
		   ! пюявер опълни вюярх йскнмнбяйнцн онремжхюкю
		   WRITE(*,*) 'DERECT Coulomb Interaction zero approximation'
	       
		   Indexcgkl=0
		   ! жхйк он лнкейскъпмшл нпахрюкъл 
		   DO INMOIDG=1,NumbreMO 
		      ! опнбепъел ясыеярбсер кх бгюхлндеиярбхе лефдс дюммшлх нанкнвйюлх (лнкейскъпмшлх нпахрюкълх)
              IF(NtypIndexCoffF(NMO,INMOIDG).NE.0) THEN
                 ! хмдейя сйюгшбючыхи мю менаундхлнярэ гюмскемхъ
                 Indexcgkl=Indexcgkl+1
			     ! нясыеярбкъел пюявер онремжхюкнб дкъ опълни вюярх
			     WRITE(*,*) 'Matrix_Potential_Direct NMO=',INMOIDG 
	             call Matrix_Potential_Direct_Coulomb_Interaction(INMOIDG,Npoint,H,LmaxGarm,NomeroFPOT,LgarmonicMOZero,IFPOTIndex,IFPOTZona,RfunZero,R,RO1X,RPOTRez,NumbreSxemaCalculation,RFILEPoltens)
       		  
		         IMO1=0
                 DO IIL1=LgarmonicMOZero(1,NMO),LgarmonicMOZero(3,NMO) 
                    IMO1=IMO1+1
	                IMO2=IMO1-1
	                DO IIL2=IIL1,LgarmonicMOZero(3,NMO) 
	                   IMO2=IMO2+1
			           WRITE(*,*) INMOIDG,IIL1,IIL2
					   ! опнбепъел менаундхлнярэ гюмскемхъ
					   IF(Indexcgkl.EQ.1) THEN
                          DCoulomb(IMO1,IMO2,1:Npoint)=0.D0
					   ENDIF
			           ! пюявер опълни вюярх
	                   call Matrix__Direct_Coulomb_Interaction_IJ(NMO,INMOIDG,IMO1,IIL1,IMO2,IIL2,NumbreMO,NtypIndexCoffF,H,RcoffA,IQ,NumbreGarmMOZero,LgarmonicMOZero,Npoint,RfunZero,DCoulomb,IndexMLA,CkkGarmonik,NomeroFPOT,IFPOTIndex,IFPOTZona,RPOT,RPOTRez,RPOTDF,NumbreSxemaCalculation,RFILEPoltens)
			        ENDDO
                 ENDDO
              ENDIF
		   ENDDO
           
		   ! гюонкмъел щкелемрш мхфе дхюцнмюкх 
           IMO1=0
           DO IIL1=LgarmonicMOZero(1,NMO),LgarmonicMOZero(3,NMO) 
              IMO1=IMO1+1
	          IMO2=IMO1-1
	          DO IIL2=IIL1,LgarmonicMOZero(3,NMO) 
	             IMO2=IMO2+1
                 WRITE(*,*) IIL1,IIL2
                 IF(IMO1.NE.IMO2) THEN    
			        ! гюохяшбюел медхюцнмюкэмши щкелемр яннрберярбсчыхи дюммнлс
				    DCoulomb(IMO2,IMO1,1:Npoint)=DCoulomb(IMO1,IMO2,1:Npoint)
			     ENDIF
			  ENDDO
           ENDDO
		   
		   ! гюмскъел нярюкэмше щкелемрш люрпхжш опълнцн бгюхлндеиярбхъ
           ! гюмскъел оепбсч вюярэ люрпхжш
           DO IMO1=NumbreGarmMOZero(NMO)+1,NumbreGarmLMO(NMO)  
              DO IMO2=1,NumbreGarmLMO(NMO)  
	             WRITE(*,*) IMO1,IMO2
	             DO J=1,Npoint   
	                DCoulomb(IMO1,IMO2,J)=0.D0
	             ENDDO
		  	  ENDDO
           ENDDO
           DO IMO2=NumbreGarmMOZero(NMO)+1,NumbreGarmLMO(NMO)  
              DO IMO1=1,NumbreGarmMOZero(NMO)
	             WRITE(*,*) IMO1,IMO2
	             DO J=1,Npoint   
	                DCoulomb(IMO1,IMO2,J)=0.D0
	             ENDDO
		  	  ENDDO
           ENDDO
        ENDIF 
        ! б яксвюе спюбмемхъ уюрпх гюмскъел опълсч вюярэ йскнмнбяйнцн бгюхлндеирбхъ мю оепбни хрепюжхх
		IF(IndexHartreeFock.EQ.0) THEN 
           ! гюмскъел мю оепбни хрепюжхх
           DCoulomb=0.D0 
		ENDIF

        ! нясыеярбкъел пюявер бпелемх
		call CPU_TIME(TimeFXZA) 
   
        RTimeCalcul(3)=DBLE((TimeFXZA-TimeSXZA)/60.)

        
		! нясыеярбкъел пюявер бпелемх
	    call CPU_TIME(TimeSXZA) 
		! опнбепъел рхо спюбмемхъ
		IF(IndexHartreeFock.EQ.1) THEN
  		   ! пюявер налеммни вюярх йскнмнбяйнцн онремжхюкю 
           IF(IKLZERO.EQ.1) THEN         
              WRITE(*,*) 'EXCHANG Coulomb Interaction zero approximation'
			  Indexcgkl=0
			  DO INMOIDG=1,NumbreMO 
			     ! опнбепъел ясыеярбсер кх бгюхлндеиярбхе лефдс дюммшлх нанкнвйюлх (лнкейскъпмшлх нпахрюкълх)
                 IF(NtypIndexCoffG(NMO,INMOIDG).NE.0) THEN 
                    ! хмдейя сйюгшбючыхи мю менаундхлнярэ гюмскемхъ
                    Indexcgkl=Indexcgkl+1
				    WRITE(*,*) 'Matrix_Potential_Exchange NMO=',NMO,INMOIDG
	                call Matrix_Potential_Exchange_Coulomb_Interaction(NMO,INMOIDG,Npoint,H,LmaxKGarm,LmaxGarm,NomeroFPOT,LgarmonicMOZero,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RfunZero,R,RO1X,RPOTRez,NumbreSxemaCalculation,RFILEPoltens)
	                 

	                IMO1=0
                    DO IIL1=LgarmonicMOZero(1,NMO),LgarmonicMOZero(3,NMO) 
                       IMO1=IMO1+1
	                   IMO2=0
	                   DO IIL2=LgarmonicMOZero(1,NMO),LgarmonicMOZero(3,NMO) 
	                      IMO2=IMO2+1
					      ! опнбепъел менаундхлнярэ гюмскемхъ
					      IF(Indexcgkl.EQ.1) THEN
                             XCoulomb(IMO1,IMO2,1:Npoint)=0.D0
					      ENDIF
			              WRITE(*,*) INMOIDG,IIL1,IIL2
			              call Matrix__Exchange_Coulomb_Interaction_IJ(NMO,INMOIDG,IMO1,IIL1,IMO2,IIL2,NumbreMO,NtypIndexCoffG,H,RcoffB,IQ,NumbreGarmMOZero,LgarmonicMOZero,Npoint,RfunZero,XCoulomb,IndexMLB,CkkGarmonik,NomeroFPOT,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RPOT,RPOTRez,RPOTDF,NumbreSxemaCalculation,RFILEPoltens)
					   ENDDO
                    ENDDO
                 ENDIF
              ENDDO
		   ENDIF 
        ENDIF
        ! нясыеярбкъел пюявер бпелемх
		call CPU_TIME(TimeFXZA) 
   
        RTimeCalcul(4)=DBLE((TimeFXZA-TimeSXZA)/60.)
	ENDIF




     
    ! тнплхпсел кнйюкэмсч вюярэ люрпхжш ноепюрнпю тнйю  
    IMO1=0
    DO IIL1=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
       IMO1=IMO1+1
	   IMO2=0
	   DO IIL2=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
	      IMO2=IMO2+1
          IF(IIL1.EQ.IIL2) THEN
		       ! тнплхпсел дхюцнмюкэмше кнйюкэмше вюярх люрпхжш ноепюрнпю тнйю  
               D=FLOAT(IIL1*(IIL1+1)) 
               DO J=1,Npoint
                  RO21=(RO2X(J)/RO1X(J))**2
	              RO31=RO3X(J)/RO1X(J)
	              ROIF=0.5D0*(RO31-1.5D0*RO21)/RO1X(J)**2
			      ! опълюъ вюярэ онремжхюкю йскнмнбяйнцн бгюхлндеиярбхъ х йпхярюккхвеяйне онке
			      B=2.D0*DCoulomb(IMO1,IMO2,J)/(R(J)*RO1X(J)*RO1X(J))+2.D0*Ncrystall(NMO,IMO1,IMO2,J)/(RO1X(J)*RO1X(J))
	              C=D/(R(J)*R(J)*RO1X(J)*RO1X(J))
                  Pz(IMO1,IMO2,J)=B+C+ROIF
               ENDDO
			   !FORALL (J=1:Npoint)
			   !  Pz(IMO1,IMO2,J)=2.D0*DCoulomb(IMO1,IMO2,J)/(R(J)*RO1X(J)*RO1X(J))+2.D0*Ncrystall(NMO,IMO1,IMO2,J)/(RO1X(J)*RO1X(J))+D/(R(J)*R(J)*RO1X(J)*RO1X(J))+0.5D0*(RO3X(J)/RO1X(J)-1.5D0*(RO2X(J)/RO1X(J))**2)/RO1X(J)**2
			   !END FORALL  
             ELSE
			   ! тнплхпсел ме дхюцнмюкэмше кнйюкэмше вюярх люрпхжш ноепюрнпю тнйю  
               ! опълюъ вюярэ онремжхюкю йскнмнбяйнцн бгюхлндеиярбхъ х йпхярюккхвеяйне онке
			   DO J=1,Npoint
			      Pz(IMO1,IMO2,J)=2.D0*DCoulomb(IMO1,IMO2,J)/(R(J)*RO1X(J)*RO1X(J))+2.D0*Ncrystall(NMO,IMO1,IMO2,J)/(RO1X(J)*RO1X(J))
			   ENDDO 			 
		       !FORALL (J=1:Npoint)
		       !  Pz(IMO1,IMO2,J)=2.D0*DCoulomb(IMO1,IMO2,J)/(R(J)*RO1X(J)*RO1X(J))+2.D0*Ncrystall(NMO,IMO1,IMO2,J)/(RO1X(J)*RO1X(J)) 	
			   !END FORALL
		  ENDIF
       ENDDO
    ENDDO 
   
    
  
  
    ! тнплхпсел мекнйюкэмсч (налеммсч) вюярэ ноепюрнпю тнйю
    ! гюмскъел оепед пюявернл
    Qz=0.D0 
    ! опнбепъел рхо спюбмемхъ
	IF(IndexHartreeFock.EQ.1) THEN
	   ! пюявер нясыеярбкъеряъ б яксвюе еякх йкчв сйюгшбюер яннрберярбсчыхи рхо мскхбнцн опхакхфемхъ 
       IF(IKLZERO.EQ.1) THEN
	      IMO1=0
          DO IIL1=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
             IMO1=IMO1+1
	         IMO2=0
	         DO IIL2=LgarmonicMO(1,NMO),LgarmonicMO(3,NMO) 
	            IMO2=IMO2+1
	            DO J=1,Npoint
	  	           Qz(IMO1,J)=Qz(IMO1,J)+2.D0*XCoulomb(IMO1,IMO2,J)/(R(J)*RO1X(J)*RO1X(J))
			    ENDDO  
             ENDDO 
          ENDDO
       ENDIF 
	ENDIF
    54166 FORMAT(2X,100(F15.10,1X)) 	
  
   ! жхйк он рнвйюл
   !IF(NumeroIter.NE.1) THEN 
      !WRITE(6,*) 'QZ  QZQZQZQZQZQZQPQZQZQZQZ  '
      !DO J=1,Npoint
        ! WRITE(6,*) 'POINTXX',J
       !  IIL1=1 !DO IIL1=1,NumbreGarmMO(NMO)   
       !     WRITE(6,54166) R(J),(Pz(IIL1,IIL2,J),IIL2=1,NumbreGarmMO(NMO))
         !ENDDO
      !ENDDO
	  !DO J=1,Npoint
	 ! 	 WRITE(6,54166) R(J),(Qz(IIL2,J),IIL2=1,NumbreGarmMO(NMO))
	 ! ENDDO  
	 ! STOP
   !ENDIF 
   !STOP 
   
  !  IF(NumeroIter.NE.1) THEN 
      !WRITE(6,*) 'QZ  QZQZQZQZQZQZQPQZQZQZQZ  '
      !DO J=1,Npoint
        ! WRITE(6,*) 'POINTXX',J
       !  IIL1=1 !DO IIL1=1,NumbreGarmMO(NMO)   
       !     WRITE(6,54166) R(J),(Pz(IIL1,IIL2,J),IIL2=1,NumbreGarmMO(NMO))
         !ENDDO
      !ENDDO
!	  IF(NMO.EQ.3) THEN
!	     DO J=1,Npoint
!	   	    WRITE(40,54166) R(J),(Qz(IIL2,J),IIL2=1,NumbreGarmMO(NMO))
!	     ENDDO  
!	  ENDIF
!   ENDIF   


	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(DCoulomb,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'THE FILE "DCoulomb" IS NOT REMOVED FROM MEMORY'
       READ(*,*)
	   stop 
    endif
    deallocate(XCoulomb,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_MO'
       write(*,*) 'THE FILE "XCoulomb" IS NOT REMOVED FROM MEMORY'
       READ(*,*)
	   stop 
    endif
    
    return
  end subroutine Calcul_Pz_Qz_MO



  !! SUBPROGRAM OF OBTAINING THE DIRECT CULON INTERACTION MATRIX OF THE OPERATOR OPERATOR BETWEEN MOLECULAR SHELLS
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF DIRECT INTERACTION IS CONDUCTED
══!! NMOIIX-NUMBER OF MOLECULAR SHELL WITH WHICH THERE IS INTERACTION OF THIS NMO
══!! Ngarmon1-HARMONIC NUMBER OF MOLECULAR ORBITALS
══!! L1-ORBITAL MOMENT
══!! Ngarmon2-HARMONIC NUMBER OF MOLECULAR ORBITALS
══!! L2-ORBITAL MOMENT OF MOLECULAR ORBITAL
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoff (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF DIRECT INTERACTION COEFFICIENTS
══!! RcoffADI (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION RELATIONS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONIC MOLECULAR ORBITALS
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! DCoulomb (Ngarmon1, Ngarmon2, Npoint) is the mass of the operator of the direct part of the Coulomb interaction of molecular orbitals
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! IndexML (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE DIRECT PIECE OF COULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO, NumbreMO, 4), IndexMLAB (NtypIndexCoff, NumbreMO)) - MASIVE OF SPHERICAL HARMONICS
══!! NumbreIndexLigandDirect (NumbreMO) -MASSIVE NUMBER OF INDICES (NUMBER OF COMPONENTS IN THE SUM) FOR DIRECT INTERACTION
══!! IndexNFunLigandDirect (NumbreMO, NumbreIndexLigandDirect (NumbreMO), 4) -MASSIVE INDICATOR WHICH SUGGESTED (COEFFICIENTS IN THE AMOUNT) SUMMIT FOR THIS MOLECULAR ORBITAL
══!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
══!! NomeroFPOT (NumbreMO) -MASSIVE NUMBER OF FILES
══!! IndexFPOT (NumbreMO, LgarmonicMO (3, NumbreMO), LgarmonicMO (2, NumbreMO), KK, 2) -MASSIVE OF DIRECT INTERACTION POPULAR CAPACITY NUMBERS
══!! RPOT, RPOTRez, RPOTDF-AUXILIARY MASSIVE
══!! NumbreSxemaCalculationF-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationF = 1-DIRECT POTENTIALS RECORDED IN THE MASIVE
══!! NumbreSxemaCalculationF = 2-DIRECT POTENTIALS RECORDED IN FILE
══!! RFILEPoltensF (Npoint, IzonFmax) -MASSIVE FOR DIRECT POTENTIALS RECORDING
  subroutine Matrix__Direct_Coulomb_Interaction_IJ(NMO,NMOIIX,Ngarmon1,L1,Ngarmon2,L2,NumbreMO,NtypIndexCoff,H,RcoffADI,IQ,NumbreGarmMO,LgarmonicMO,Npoint,RFunMO,DCoulomb,IndexML,CkkGarmonik,NomeroFPOT,IFPOTIndex,IFPOTZona,RPOT,RPOTRez,RPOTDF,NumbreSxemaCalculationF,RFILEPoltensF)
    implicit none
    integer::NMO,NMOIIX,Ngarmon1,L1,Ngarmon2,L2,NumbreMO,Npoint,NumbreSxemaCalculationF
    real(8)::H
    integer,dimension(:)::NumbreGarmMO,IQ,NomeroFPOT
    integer,dimension(:,:)::LgarmonicMO,NtypIndexCoff
   	integer(1),dimension(:,:,:,:)::IFPOTIndex
	integer,dimension(:,:,:,:)::IndexML,IFPOTZona
	real(8),dimension(:,:)::RFILEPoltensF  
    real(8),dimension(:,:,:)::RcoffADI,RFunMO
    real(8),dimension(:,:,:)::DCoulomb
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    real(8),dimension(:)::RPOT,RPOTRez,RPOTDF
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IITyp,II1,II2,KK,MINN,MAXX,IIG1,IIG2,INCV,ierr,IndexML1,IndexML2,IndexML3,IndexML4,IICGF,IOIO,IDFSZ,IOSD
    real(8)::RcoffDrect,RGarm,RcoffMON,RcoffSimmetry,Rsum,RnormFun
    

    ! опнбепъел ясыеярбсер кх бгюхлндеиярбхе лефдс дюммшлх нанкнвйюлх (лнкейскъпмшлх нпахрюкълх)
    IF(NtypIndexCoff(NMO,NMOIIX).NE.0) THEN

	   ! пюявер опълни вюярх йскнмнбяйнцн бгюхлндеиярбхъ C NMOIIX лнкейскъпмни нанкнвйни
       ! мнлеп тюикю
	   IOSD=NomeroFPOT(1)
       IF(NMOIIX.EQ.NMO) THEN
           RcoffMON=2.D0
	      ELSE
           RcoffMON=1.D0
	   ENDIF
	  
	   ! гюмскъел оепед пюявернл
	   RPOTDF=0.D0
	
	   ! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх NMOIIX
       ! мнлеп цюплнмхйх
	   IIG1=0
	   DO II1=LgarmonicMO(1,NMOIIX),LgarmonicMO(2,NMOIIX)   
          IIG1=IIG1+1
          ! мнлеп цюплнмхйх
	      IIG2=IIG1-1
	      DO II2=II1,LgarmonicMO(2,NMOIIX)   
	         IIG2=IIG2+1
   	         ! свхршбюел яхллерпхч нрмняхрекэмн оепеярюмнбяйх хмдейянб 
             IF(IIG1.EQ.IIG2) THEN
                 RcoffSimmetry=1.D0
		        ELSE
                 RcoffSimmetry=2.D0 
             ENDIF 
		     MINN=MAX0(IABS(L1-L2),IABS(II1-II2))
             MAXX=MIN0(L1+L2,II1+II2)
		     DO KK=MINN,MAXX,2
		        ! жхйк он рхоюл йнщттхжхемрнб
		        RGarm=0.D0
	            DO IITyp=1,NtypIndexCoff(NMO,NMOIIX) 
                   ! хмдейяш опнейжхи  
		           IndexML1=IndexML(IITyp,NMO,NMOIIX,1)
                   IndexML2=IndexML(IITyp,NMO,NMOIIX,2)
			       IndexML3=IndexML(IITyp,NMO,NMOIIX,3)
			  	   IndexML4=IndexML(IITyp,NMO,NMOIIX,4)
                   ! йнщттхжхемр йскнмнбяйнцн бгюхлндеиярбхъ  
	               RcoffDrect=RcoffMON*RcoffADI(IITyp,NMO,NMOIIX)/FLOAT(IQ(NMO))
				   ! пюявер опнхгбедемхъ цюплнмхй 
				   RGarm=RGarm+RcoffDrect*CkkGarmonik(KK+1,II2+1,II1+1,IndexML4,IndexML2)*CkkGarmonik(KK+1,L1+1,L2+1,IndexML1,IndexML3) 
			    ENDDO
 
			      					 
			    IF(DABS(RGarm).GT.1.D-12) THEN
			       ! пюявер онремжхюкю
			       ! нясыеярбкъел явхршбюмхе онремжхюкю
			       IF(IFPOTIndex(NMOIIX,II1+1,II2+1,KK+1).EQ.0) THEN
                      ! дюммши онремжхюк ме гюохяюм
			          WRITE(6,*) 'FPOT NOT WRITING  NMO=',NMOIIX,' L1=',II1,' L2=',II2,' K=',KK
			  	      STOP
			       ENDIF
				   ! гнмю гюохях
                   IDFSZ=IFPOTZona(NMOIIX,II1+1,II2+1,KK+1)
				   ! явхршбюмхе хг люяяхб
				   IF(NumbreSxemaCalculationF.EQ.1) THEN
                      RPOT(1:Npoint)=RFILEPoltensF(1:Npoint,IDFSZ)
				   ENDIF
                   ! гюохяэ б тюик
				   IF(NumbreSxemaCalculationF.EQ.2) THEN
				      READ (IOSD,REC=IDFSZ) (RPOT(IOIO),IOIO=1,Npoint)  
			       ENDIF
						 
				   ! мюундхл ясллс (онремжхюк)
				   RPOTDF=RPOTDF+RcoffSimmetry*RGarm*RPOT
                ENDIF
		     ENDDO 
  	      ENDDO
       ENDDO
        
	   ! гюохяшбюел пегскэрюр пюяверю
	   DCoulomb(Ngarmon1,Ngarmon2,1:Npoint)=DCoulomb(Ngarmon1,Ngarmon2,1:Npoint)+RPOTDF(1:Npoint)
    ENDIF
	
    return
  end subroutine Matrix__Direct_Coulomb_Interaction_IJ



  !! SUBPROGRAMME OF OBTAINING THE MATRIX OF CAPACITY OF DIRECT CULONARY INTERACTION
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF DIRECT INTERACTION IS CONDUCTED
══!! LmaxKGarm-MAXIMUM VALUE OF HARMONIC LgarmonicMO (3, NumbreMO)
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RPOT-AUXILIARY MASSIVE
══!! NumbreSxemaCalculationF-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationF = 1-CAPACITY RECORDING IS MADE IN ARMY
══!! NumbreSxemaCalculationF = 2-CAPACITY RECORDING IS FILEED
══!! RFILEPoltensF (Npoint, IzonFmax) -MASSIVE FOR DIRECT POTENTIALS RECORDING
  subroutine Matrix_Potential_Direct_Coulomb_Interaction(NMO,Npoint,H,LmaxKGarm,NomeroFPOT,LgarmonicMO,IFPOTIndex,IFPOTZona,RFunMO,R,RO1X,RPOT,NumbreSxemaCalculationF,RFILEPoltensF)
    implicit none
    integer::NMO,Npoint,LmaxKGarm,NumbreSxemaCalculationF
    real(8)::H
    integer,dimension(:)::NomeroFPOT
    integer,dimension(:,:)::LgarmonicMO
    integer(1),dimension(:,:,:,:)::IFPOTIndex
	integer,dimension(:,:,:,:)::IFPOTZona
	real(8),dimension(:)::R,RO1X
    real(8),dimension(:,:)::RFILEPoltensF
    real(8),dimension(:,:,:)::RFunMO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    real(8),dimension(:)::RPOT
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,II1,II2,KK,IIG1,IIG2,IDXFZ,IFSD
   
    
    ! мнлеп тюикю гюохях
    IFSD=NomeroFPOT(1)

   	! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх NMO
    ! мнлеп цюплнмхйх
	IIG1=0
	DO II1=LgarmonicMO(1,NMO),LgarmonicMO(2,NMO) 
       IIG1=IIG1+1
       ! мнлеп цюплнмхйх
	   IIG2=IIG1-1
	   DO II2=II1,LgarmonicMO(2,NMO)   
	      IIG2=IIG2+1
   		  DO KK=0,2*LmaxKGarm
			 IF(IFPOTIndex(NMO,II1+1,II2+1,KK+1).EQ.1) THEN
			    ! пюявер онремжхюкю
			    call POTENS(DBLE(KK),NMO,IIG1,NMO,IIG2,RFunMO,Npoint,H,RPOT,R,RO1X)
			    ! нясыеярбкъел гюохяэ онксвеммнцн онремжхюкю
				! гнмю гюохях
				IDXFZ=IFPOTZona(NMO,II1+1,II2+1,KK+1) 
				! гюохяэ б люяяхб
				IF(NumbreSxemaCalculationF.EQ.1) THEN
                   RFILEPoltensF(1:Npoint,IDXFZ)=RPOT(1:Npoint)
				ENDIF
                ! гюохяэ б тюик
				IF(NumbreSxemaCalculationF.EQ.2) THEN
				   WRITE(IFSD,REC=IDXFZ) (RPOT(IIX),IIX=1,Npoint)  
			    ENDIF
			 ENDIF  
		  ENDDO 
	   ENDDO
	ENDDO

   

	
    return
  end subroutine Matrix_Potential_Direct_Coulomb_Interaction



  !! SUBPROGRAMME OF OBTAINING THE MATRIX OF CAPACITY OF DIRECT CULONARY INTERACTION
══!! FOR CALCULATION OF FULL ENERGY
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF DIRECT INTERACTION IS CONDUCTED
══!! LmaxGarm-MAXIMUM VALUE OF HARMONIC LgarmonicMO (3, NumbreMO)
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RPOT-AUXILIARY MASSIVE
══!! NumbreSxemaCalculationF-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationF = 1-POTENTIAL RECORDING IS INTO THE MASSIVE real (8)
══!! NumbreSxemaCalculationF = 2-CAPACITY BUILDING IS MADE IN THE MASSIVE real (4)
══!! NumbreSxemaCalculationF = 3-CAPACITY RECORDING IS FILEED
══!! RFILEPoltensF8 (Npoint, IzonFmax) -MASSIVE FOR DIRECT POTENTIALS RECORDING
══!! RFILEPoltensF4 (Npoint, IzonFmax) -MASSIVE FOR DIRECT POTENTIALS RECORDING
  subroutine Matrix_Energy_Potential_Direct_Coulomb_Interaction(NMO,Npoint,H,LmaxGarm,NomeroFPOT,LgarmonicMO,IFPOTIndex,IFPOTZona,RFunMO,R,RO1X,RPOT,NumbreSxemaCalculationF,RFILEPoltensF8,RFILEPoltensF4)
    implicit none
    integer::NMO,Npoint,LmaxGarm,NumbreSxemaCalculationF
    real(8)::H
    integer,dimension(:)::NomeroFPOT
    integer,dimension(:,:)::LgarmonicMO
    integer(1),dimension(:,:,:,:)::IFPOTIndex
	integer,dimension(:,:,:,:)::IFPOTZona
	real(8),dimension(:)::R,RO1X
    real(4),dimension(:,:)::RFILEPoltensF4
    real(8),dimension(:,:)::RFILEPoltensF8
    real(8),dimension(:,:,:)::RFunMO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    real(8),dimension(:)::RPOT
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,II1,II2,KK,IIG1,IIG2,IDXFZ,IFSD
   
    
    ! мнлеп тюикю гюохях
    IFSD=NomeroFPOT(1)

   	! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх NMO
    ! мнлеп цюплнмхйх
	IIG1=0
	DO II1=LgarmonicMO(1,NMO),LgarmonicMO(2,NMO)
       IIG1=IIG1+1
       ! мнлеп цюплнмхйх
	   IIG2=IIG1-1
	   DO II2=II1,LgarmonicMO(2,NMO)   
	      IIG2=IIG2+1
   		  DO KK=0,2*LmaxGarm
			 IF(IFPOTIndex(NMO,II1+1,II2+1,KK+1).EQ.1) THEN
			    ! пюявер онремжхюкю
			    call POTENSS(DBLE(KK),NMO,IIG1,NMO,IIG2,RFunMO,Npoint,H,RPOT,R,RO1X)
			    ! нясыеярбкъел гюохяэ онксвеммнцн онремжхюкю
				! гнмю гюохях
				IDXFZ=IFPOTZona(NMO,II1+1,II2+1,KK+1) 
				! гюохяэ б люяяхб
                IF(NumbreSxemaCalculationF.EQ.1) THEN
                   RFILEPoltensF8(1:Npoint,IDXFZ)=RPOT(1:Npoint)
				ENDIF
				IF(NumbreSxemaCalculationF.EQ.2) THEN
                   RFILEPoltensF4(1:Npoint,IDXFZ)=SNGL(RPOT(1:Npoint))
				ENDIF
                ! гюохяэ б тюик
				IF(NumbreSxemaCalculationF.EQ.3) THEN
				   WRITE(IFSD,REC=IDXFZ) (RPOT(IIX),IIX=1,Npoint)  
			    ENDIF
			 ENDIF  
		  ENDDO 
	   ENDDO
	ENDDO

   

	
    return
  end subroutine Matrix_Energy_Potential_Direct_Coulomb_Interaction






  !! SUBPROGRAMME OF RECEIVING THE MASSIVE OF THE INDICATING ZONE OF RECORDING AND MULTIPLETITY OF DIRECT INTERACTION POTENTIALS
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF DIRECT INTERACTION IS CONDUCTED
══!! L1-ORBITAL MOMENT
══!! L2-ORBITAL MOMENT OF MOLECULAR ORBITAL
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! LmaxKGarm-MAXIMUM VALUE LgarmonicMO (3, NumbreMO)
══!! NtypIndexCoff (NumbreMO, NumbreMO) -MASIVE OF NUMBERS OF TYPES OF DIRECT INTERACTION COEFFICIENTS
══!! RcoffADI (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION RELATIONS
══!! IZONAMASSIV (NumbreMO) - MASSIVE NUMBER OF CAPABILITIES OF POTENTIALS
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! IndexML (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE DIRECT PIECE OF COULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO)) - SPHERICAL HARMONIC array
  subroutine FPOT__Direct_Coulomb_Interaction(NMO,L1,L2,NumbreMO,LmaxKGarm,NtypIndexCoff,RcoffADI,IZONAMASSIV,IFPOTIndex,IFPOTZona,LgarmonicMO,IndexML,CkkGarmonik)
    implicit none
    integer::NMO,L1,L2,NumbreMO,LmaxKGarm
    integer,dimension(:)::IZONAMASSIV
	integer,dimension(:,:)::LgarmonicMO,NtypIndexCoff
	integer(1),dimension(:,:,:,:)::IFPOTIndex
    integer,dimension(:,:,:,:)::IFPOTZona,IndexML
    real(8),dimension(:,:,:)::RcoffADI
  	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IITyp,II1,II2,KK,MINN,MAXX,IIG1,IIG2,INCV,ierr,IndexML1,IndexML2,IndexML3,IndexML4,IICGF,IZONA
    real(8)::RGarm
    
    

    ! пюявер опълни вюярх йскнмнбяйнцн бгюхлндеиярбхъ
    DO IIX=1,NumbreMO 
	   ! опнбепъел ясыеярбсер бгюхлндеиярбхе лефдс дюммшлх лнкейскъпмшлх нпахрюкълх хкх мер
       IF(NtypIndexCoff(NMO,IIX).NE.0) THEN
          IZONA=IZONAMASSIV(IIX)
	      
		  ! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх IIX
          DO II1=LgarmonicMO(1,IIX),LgarmonicMO(2,IIX)
             DO II2=II1,LgarmonicMO(2,IIX)   
		            			       
		        MINN=MAX0(IABS(L1-L2),IABS(II1-II2))
                MAXX=MIN0(L1+L2,II1+II2)
			    DO KK=MINN,MAXX,2
			       ! жхйк он рхоюл йнщттхжхемрнб
			       RGarm=0.D0
	               DO IITyp=1,NtypIndexCoff(NMO,IIX) 
                      ! йнщттхжхемр йскнмнбяйнцн бгюхлндеиярбхъ  
		              IndexML1=IndexML(IITyp,NMO,IIX,1)
                      IndexML2=IndexML(IITyp,NMO,IIX,2)
			          IndexML3=IndexML(IITyp,NMO,IIX,3)
			  	      IndexML4=IndexML(IITyp,NMO,IIX,4)
					  ! пюявер опнхгбедемхъ цюплнмхй 
				      RGarm=RGarm+RcoffADI(IITyp,NMO,IIX)*CkkGarmonik(KK+1,II2+1,II1+1,IndexML4,IndexML2)*CkkGarmonik(KK+1,L1+1,L2+1,IndexML1,IndexML3) 
			       ENDDO

				   IF(DABS(RGarm).GT.1.D-12) THEN
                      ! опнбепъел бшундхр гю пюлйх хкх мер
				      IF(KK.GT.2*LmaxKGarm) THEN 
                         WRITE(6,*) 'FPOT__Direct_Coulomb_Interaction'
					     WRITE(6,*) 'Multy Garm KK>',2*LmaxKGarm,' KK= ',KK
					     STOP
				      ENDIF
				      ! опнбепъел тхйяхпнбюм дюммши онремжхюк хкх мер
				      IF(IFPOTIndex(IIX,II1+1,II2+1,KK+1).NE.1) THEN
                         ! тхйяхпсел рн, врн дюммши онремжхюк сфе свюбярбнбюк б пюявере
                         IFPOTIndex(IIX,II1+1,II2+1,KK+1)=1
				   	     ! гюохяшбюел мнлеп гнмш гюохях
					     IZONA=IZONA+1
					     IFPOTZona(IIX,II1+1,II2+1,KK+1)=IZONA 
				      ENDIF
                   ENDIF
		        ENDDO 
		     ENDDO
	      ENDDO       
	    
          IZONAMASSIV(IIX)=IZONA 
	   ENDIF
	ENDDO

    

	
    return
  end subroutine FPOT__Direct_Coulomb_Interaction


  !! SUBPROGRAMME OF RECEIVING THE MASSIVE OF THE INDICATING ZONE OF RECORDING AND MULTIPLETITY OF DIRECT INTERACTION POTENTIALS
══!! FOR CALCULATION OF CONFIGURATION ENERGY
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF DIRECT INTERACTION IS CONDUCTED
══!! L1-ORBITAL MOMENT
══!! L2-ORBITAL MOMENT OF MOLECULAR ORBITAL
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! LmaxGarm-MAXIMUM VALUE LgarmonicMO (2, NumbreMO)
══!! NtypIndexCoff (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF DIRECT INTERACTION COEFFICIENTS
══!! RcoffADI (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION RELATIONS
══!! IZONAMASSIV (NumbreMO) - MASSIVE NUMBER OF CAPABILITIES OF POTENTIALS
══!! IFPOTIndex (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE KEY OF CALCULATING DIRECT CAPACITY
══!! IFPOTZona (NumbreMO, LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK) -MASIVE IN WHICH THE DIRECT CAPACITY
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! IndexML (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF TOMORROW (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO, NumbreMO, 4), IndexMLAB (NtypIndexCoff, NumbreMO, NumbreMO, 4)) - MAP OF SPHERICAL HARMONIC
  subroutine FPOT_Energy_Direct_Coulomb_Interaction(NMO,L1,L2,NumbreMO,LmaxGarm,NtypIndexCoff,RcoffADI,IZONAMASSIV,IFPOTIndex,IFPOTZona,LgarmonicMO,IndexML,CkkGarmonik)
    implicit none
    integer::NMO,L1,L2,NumbreMO,LmaxGarm
    integer,dimension(:)::IZONAMASSIV
	integer,dimension(:,:)::NtypIndexCoff,LgarmonicMO
	integer(1),dimension(:,:,:,:)::IFPOTIndex
    integer,dimension(:,:,:,:)::IndexML,IFPOTZona
    real(8),dimension(:,:,:)::RcoffADI
  	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IITyp,II1,II2,KK,MINN,MAXX,IIG1,IIG2,INCV,ierr,IndexML1,IndexML2,IndexML3,IndexML4,IICGF,IZONA
    real(8)::RcoffDrect,RGarm
    
    

    ! пюявер опълни вюярх йскнмнбяйнцн бгюхлндеиярбхъ
    DO IIX=1,NumbreMO 
	   ! опнбепъел ясыеярбсер бгюхлндеиярбхе лефдс дюммшлх лнкейскъпмшлх нпахрюкълх хкх мер
       IF(NtypIndexCoff(NMO,IIX).NE.0) THEN
          IZONA=IZONAMASSIV(IIX)
	  
          ! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх IIX
          DO II1=LgarmonicMO(1,IIX),LgarmonicMO(2,IIX)
             DO II2=II1,LgarmonicMO(2,IIX)   
		          	       
     	        MINN=MAX0(IABS(L1-L2),IABS(II1-II2))
                MAXX=MIN0(L1+L2,II1+II2)
			    DO KK=MINN,MAXX,2
			       ! жхйк он рхоюл йнщттхжхемрнб
			       RGarm=0.D0
	               DO IITyp=1,NtypIndexCoff(NMO,IIX) 
                      ! йнщттхжхемр йскнмнбяйнцн бгюхлндеиярбхъ  
		              IndexML1=IndexML(IITyp,NMO,IIX,1)
                      IndexML2=IndexML(IITyp,NMO,IIX,2)
			          IndexML3=IndexML(IITyp,NMO,IIX,3)
			          IndexML4=IndexML(IITyp,NMO,IIX,4)
			          ! пюявер опнхгбедемхъ цюплнмхй 
			          RGarm=RGarm+RcoffADI(IITyp,NMO,IIX)*CkkGarmonik(KK+1,II2+1,II1+1,IndexML4,IndexML2)*CkkGarmonik(KK+1,L1+1,L2+1,IndexML1,IndexML3) 
			       ENDDO
				   
				   IF(DABS(RGarm).GT.1.D-12) THEN
                      ! опнбепъел бшундхр гю пюлйх хкх мер
				      IF(KK.GT.2*LmaxGarm) THEN 
                         WRITE(6,*) 'FPOT__Direct_Coulomb_Interaction'
				  	     WRITE(6,*) 'Multy Garm KK>',2*LmaxGarm,' KK= ',KK
					     STOP
				      ENDIF
				      ! опнбепъел тхйяхпнбюм дюммши онремжхюк хкх мер
				      IF(IFPOTIndex(IIX,II1+1,II2+1,KK+1).NE.1) THEN
                         ! тхйяхпсел рн, врн дюммши онремжхюк сфе свюбярбнбюк б пюявере
                         IFPOTIndex(IIX,II1+1,II2+1,KK+1)=1
				  	     ! гюохяшбюел мнлеп гнмш гюохях
					     IZONA=IZONA+1
					     IFPOTZona(IIX,II1+1,II2+1,KK+1)=IZONA 
				      ENDIF
                   ENDIF
		        ENDDO 
		     ENDDO
	      ENDDO
             
	      IZONAMASSIV(IIX)=IZONA 
	   ENDIF
	ENDDO

    

	
    return
  end subroutine FPOT_Energy_Direct_Coulomb_Interaction


  



  !! SUBPROGRAMME OF OBTAINING THE MATRIX OF THE OPERATOR OF THE EXCHANGE COULOMOUS INTERACTION FOR INTERACTION OF THIS SHELL WITH AN OTHER SHELL
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF EXCHANGE INTERACTION IS CONDUCTED
══!! NMOIIX-NUMBER OF THE MOLECULAR ORBITAL WITH WHICH EXCHANGE INTERACTION IS CARRIED OUT
══!! Ngarmon1-HARMONIC NUMBER OF MOLECULAR ORBITALS
══!! L1-ORBITAL MOMENT
══!! Ngarmon2-HARMONIC NUMBER OF MOLECULAR ORBITALS
══!! L2-ORBITAL MOMENT OF MOLECULAR ORBITAL
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoff (NumbreMO, NumbreMO) -MASIVE OF NUMBERS OF TYPES OF EXCHANGE INTERACTION COEFFICIENTS
══!! RcoffBEI (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONIC MOLECULAR ORBITALS
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! Npoint-number of points
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! XCoulomb (Ngarmon1, Ngarmon2, Npoint) is the mass of the operator of the direct part of the Coulomb interaction of molecular orbitals
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! IndexML (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF TOMORROW (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO)) - SPHERICAL HARMONIC array
══!! RPOT, RPOTRez, RPOTDG-SUPPORTING MASSIVE
══!! NomeroFPOT (NumbreMO) -MASSIVE IN WHICH FILE NUMBERS ARE FOR RECORDING EXCHANGE POTENTIALS
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! NumbreGPOTMORPOT (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION FOR THE MASSIVE OF RADIAL PARTS OF POTENTIALS
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! NumbreSxemaCalculationG-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationG = 1-EXCHANGE OF EXCHANGE POTENTIALS IS PROVIDED IN THE MASSIF
══!! NumbreSxemaCalculationG = 2-EXCHANGE OF EXCHANGE POTENTIALS FILES INTO FILE
══!! RFILEPoltensG (Npoint, IzonGmax) -FILE FOR RECORDING EXCHANGE POTENTIALS
  subroutine Matrix__Exchange_Coulomb_Interaction_IJ(NMO,NMOIIX,Ngarmon1,L1,Ngarmon2,L2,NumbreMO,NtypIndexCoff,H,RcoffBEI,IQ,NumbreGarmMO,LgarmonicMO,Npoint,RFunMO,XCoulomb,IndexML,CkkGarmonik,NomeroFPOT,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RPOT,RPOTRez,RPOTDG,NumbreSxemaCalculationG,RFILEPoltensG)
    implicit none
    integer::NMO,NMOIIX,Ngarmon1,L1,Ngarmon2,L2,NumbreMO,Npoint,NumbreSxemaCalculationG
    real(8)::H
    integer,dimension(:)::NumbreGarmMO,IQ,NomeroFPOT
    integer,dimension(:,:)::NtypIndexCoff,LgarmonicMO,NumbreGPOTMO,NumbreGPOTMORPOT
    integer(1),dimension(:,:,:,:)::IGPOTIndex
	integer,dimension(:,:,:,:)::IndexML,IGPOTZona
	real(8),dimension(:,:)::RFILEPoltensG
	real(8),dimension(:,:,:)::RcoffBEI,RFunMO
    real(8),dimension(:,:,:)::XCoulomb
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8),dimension(:)::RPOT,RPOTRez,RPOTDG
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IITyp,II1,II2,KK,MINN,MAXX,IIG1,IIG2,INCV,ierr,IndexML1,IndexML2,IndexML3,IndexML4,Lkk,IICGF,IOIU,ITYGR,IZOGR,IDSFZX,IDFHG
    real(8)::RcoffExchange,RGarm,RnormFun,RsumG
    
    
	! опнбепъел ясыеярбсер кх бгюхлндеиярбхе лефдс дюммшлх нанкнвйюлх (лнкейскъпмшлх нпахрюкълх)
    IF(NtypIndexCoff(NMO,NMOIIX).NE.0) THEN
    
       ! гюмскъел оепед пюявернл
 	   RPOTRez=0.D0
    
	   ! пюявер налеммни вюярх йскнмнбяйнцн бгюхлндеиярбхъ лефдс нанкнвйюлх NMO,NMOIIX
    
	   ! мнлеп рхою
       ITYGR=NumbreGPOTMO(NMO,NMOIIX)
       ! мнлеп тюикю
	   IDSFZX=NomeroFPOT(1) 
   
       ! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх NMOIIX
       ! мнлеп цюплнмхйх
       IIG1=0
       DO II1=LgarmonicMO(1,NMOIIX),LgarmonicMO(2,NMOIIX)
          IIG1=IIG1+1
          ! мнлеп цюплнмхйх
          IIG2=0
          DO II2=LgarmonicMO(1,NMOIIX),LgarmonicMO(2,NMOIIX)   
             IIG2=IIG2+1
             MINN=MAX0(IABS(L1-II2),IABS(II1-L2))
             MAXX=MIN0(L1+II2,II1+L2)
   	         DO KK=MINN,MAXX,2
		        ! пюявер опнхгбедемхъ цюплнмхй 
                RGarm=0.D0
		        ! жхйк он рхоюл йнщттхжхемрнб
	            DO IITyp=1,NtypIndexCoff(NMO,NMOIIX)
		           ! хмдейяш опнейжхи  
		           IndexML1=IndexML(IITyp,NMO,NMOIIX,1)
                   IndexML2=IndexML(IITyp,NMO,NMOIIX,2)
			       IndexML3=IndexML(IITyp,NMO,NMOIIX,3)
			       IndexML4=IndexML(IITyp,NMO,NMOIIX,4)
			       ! йнщттхжхемр йскнмнбяйнцн бгюхлндеиярбхъ  
	               RcoffExchange=RcoffBEI(IITyp,NMO,NMOIIX)/FLOAT(IQ(NMO))
			       ! пюявер опнхгбедемхъ цюплнмхй 
			       RGarm=RGarm+RcoffExchange*CkkGarmonik(KK+1,L2+1,II1+1,IndexML4,IndexML2)*CkkGarmonik(KK+1,L1+1,II2+1,IndexML1,IndexML3)
			    ENDDO

		        IF(DABS(RGarm).GT.1.D-12) THEN
		           ! пюявер онремжхюкю
		           IF(IGPOTIndex(ITYGR,II1+1,L2+1,KK+1).EQ.0) THEN
                      ! дюммши онремжхюк ме гюохяюм
		     	      WRITE(6,*) 'GPOT NOT WRITING  NMO=',NMO,NMOIIX,' L1=',L2,' L2=',II1,' K=',KK
				      STOP
			       ENDIF
              
			       ! гнмю гюохях
			       IZOGR=IGPOTZona(ITYGR,II1+1,L2+1,KK+1)
  			    
				   ! нясыеярбкъел явхршбюмхе налеммнцн онремжхюкю
				   ! явхршбюмхе хг люяяхбю 
				   IF(NumbreSxemaCalculationG.EQ.1) THEN
                      ! жхйк он рнвйюл
				      ! мнплхпнбйю тсмйжхх
                      RnormFun=1.D0/RFunMO(NMOIIX,IIG2,2)
				      RPOTRez(1:Npoint)=RPOTRez(1:Npoint)+RGarm*RnormFun*RFILEPoltensG(1:Npoint,IZOGR)*RFunMO(NMOIIX,IIG2,2+1:Npoint)
				   ENDIF
                
				   ! явхршбюмхе хг тюик
				   IF(NumbreSxemaCalculationG.EQ.2) THEN
				      READ(IDSFZX,REC=IZOGR) (RPOT(IOIU),IOIU=1,Npoint)
				      ! жхйк он рнвйюл
				      ! мнплхпнбйю тсмйжхх
                      RnormFun=1.D0/RFunMO(NMOIIX,IIG2,2)
				      RPOTRez(1:Npoint)=RPOTRez(1:Npoint)+RGarm*RnormFun*RPOT(1:Npoint)*RFunMO(NMOIIX,IIG2,2+1:Npoint)
				   ENDIF    
                ENDIF 
		     ENDDO 
	      ENDDO
	   ENDDO
       
       ! гюохяшбюел пегскэрюр пюяверю 
       XCoulomb(Ngarmon1,Ngarmon2,1:Npoint)=XCoulomb(Ngarmon1,Ngarmon2,1:Npoint)+RPOTRez(1:Npoint)
	ENDIF  
	
     	    
    return
  end subroutine Matrix__Exchange_Coulomb_Interaction_IJ



  

  !! SUBPROGRAMME OF OBTAINING THE MATRIX OF THE POTENTIAL OF EXCHANGE COULIN INTERACTION
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO1, NMO2-NUMBERS OF MOLECULAR ORBITALS FOR WHICH CALCULATION OF EXCHANGE INTERACTION IS CONDUCTED
══!! LmaxKGarm-MAXIMUM VALUE OF HARMONIC LgarmonicMO (3, NumbreMO)
══!! LmaxGarm-MAXIMUM VALUE OF HARMONIC LgarmonicMO (2, NumbreMO)
══!! NomeroFPOT (1) -MASSIVE NUMBERING FILE NUMBERS IN WHICH POTENTIALS ARE CARRIED OUT
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! NumbreGPOTMORPOT (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION FOR THE MASSIVE OF RADIAL POTENTIALS
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! Npoint-number of points
══!! H-STEP
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RPOT-AUXILIARY MASSIVE
══!! NumbreSxemaCalculationG-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationG = 1-CAPACITY BUILDING MADE IN ARMY
══!! NumbreSxemaCalculationG = 2-CAPACITY RECORDING IS IN THE FILE
══!! RFILEPoltensG (Npoint, IzonGmax) -FILE FOR RECORDING EXCHANGE POTENTIALS
  subroutine Matrix_Potential_Exchange_Coulomb_Interaction(NMO1,NMO2,Npoint,H,LmaxKGarm,LmaxGarm,NomeroFPOT,LgarmonicMO,NumbreGPOTMO,NumbreGPOTMORPOT,IGPOTIndex,IGPOTZona,RFunMO,R,RO1X,RPOT,NumbreSxemaCalculationG,RFILEPoltensG)
    implicit none
    integer::NMO1,NMO2,Npoint,LmaxKGarm,LmaxGarm,NumbreSxemaCalculationG
    real(8)::H
	integer,dimension(:)::NomeroFPOT 
    integer,dimension(:,:)::LgarmonicMO,NumbreGPOTMO,NumbreGPOTMORPOT
	integer(1),dimension(:,:,:,:)::IGPOTIndex
	integer,dimension(:,:,:,:)::IGPOTZona
    real(8),dimension(:)::R,RO1X
    real(8),dimension(:,:)::RFILEPoltensG
    real(8),dimension(:,:,:)::RFunMO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8),dimension(:)::RPOT
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,KK,IIG1,IIG2,II1,II2,IDGX,IXSG,ITYG,IIRTG
        
    
	! мнлеп тюикю
    IDGX=NomeroFPOT(1)
    ! мнлеп рхою бгюхлндеиярбхъ
    ITYG=NumbreGPOTMO(NMO1,NMO2)
	! мнлеп дкъ люяяхбю пюдхюкэмшу онремжхюкнб
    IIRTG=NumbreGPOTMORPOT(NMO1,NMO2)

    ! пюявер налеммни вюярх йскнмнбяйнцн бгюхлндеиярбхъ
    ! мнлеп цюплнмхйх
	IIG1=0
	DO II1=LgarmonicMO(1,NMO2),LgarmonicMO(2,NMO2)
       IIG1=IIG1+1
       ! мнлеп цюплнмхйх
	   IIG2=0
	   DO II2=LgarmonicMO(1,NMO1),LgarmonicMO(3,NMO1)   
	      IIG2=IIG2+1
	      DO KK=0,LmaxKGarm+LmaxGarm
	         IF(IGPOTIndex(ITYG,II1+1,II2+1,KK+1).EQ.1) THEN
			    ! пюявер онремжхюкю
			    call POTENG(DBLE(KK),NMO2,IIG1,NMO1,IIG2,RFunMO,Npoint,H,RPOT,R,RO1X)
				! нясыеярбкъел гюохяэ онксвеммнцн онремжхюкю
				! мнлеп гюохях
                IXSG=IGPOTZona(ITYG,II1+1,II2+1,KK+1)
                ! гюохяэ б люяяхб
				IF(NumbreSxemaCalculationG.EQ.1) THEN
                   RFILEPoltensG(1:Npoint,IXSG)=RPOT(1:Npoint)
				ENDIF
                ! гюохяэ б тюик
				IF(NumbreSxemaCalculationG.EQ.2) THEN
				   WRITE(IDGX,REC=IXSG) (RPOT(IIX),IIX=1,Npoint)  
			    ENDIF
			 ENDIF 
		  ENDDO 
	   ENDDO
	ENDDO
     
     	    
    return
  end subroutine Matrix_Potential_Exchange_Coulomb_Interaction


  !! SUBPROGRAMME OF OBTAINING THE MATRIX OF THE POTENTIAL OF EXCHANGE COULIN INTERACTION
══!! FOR CALCULATION OF FULL ENERGY CONFIGURATION
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO1, NMO2-NUMBERS OF MOLECULAR ORBITALS FOR WHICH CALCULATION OF EXCHANGE INTERACTION IS CONDUCTED
══!! LmaxGarm-MAXIMUM VALUE OF HARMONIC LgarmonicMO (2, NumbreMO)
══!! NomeroFPOT (1) -MASSIVE NUMBERING FILE NUMBERS IN WHICH POTENTIALS ARE CARRIED OUT
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! Npoint-number of points
══!! H-STEP
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! RPOT-AUXILIARY MASSIVE
══!! NumbreSxemaCalculationG-PARAMETER CALCULATION DIAGRAM
══!! NumbreSxemaCalculationG = 1-CAPACITY BUILDING IS REALIZED TO REAL (8)
══!! NumbreSxemaCalculationG = 2-CAPACITY BUILDING IS REALIZED TO REAL (4)
══!! NumbreSxemaCalculationG = 3-CAPACITY RECORDING IS IN FILE
══!! RFILEPoltensG8 (Npoint, IzonGmax) -FILE FOR RECORDING EXCHANGE POTENTIALS
══!! RFILEPoltensG4 (Npoint, IzonGmax) -FILE FOR RECORDING EXCHANGE POTENTIALS
  subroutine Matrix_Enegry_Potential_Exchange_Coulomb_Interaction(NMO1,NMO2,Npoint,H,LmaxGarm,NomeroFPOT,LgarmonicMO,NumbreGPOTMO,IGPOTIndex,IGPOTZona,RFunMO,R,RO1X,RPOT,NumbreSxemaCalculationG,RFILEPoltensG8,RFILEPoltensG4)
    implicit none
    integer::NMO1,NMO2,Npoint,LmaxGarm,NumbreSxemaCalculationG
    real(8)::H
	integer,dimension(:)::NomeroFPOT
    integer,dimension(:,:)::LgarmonicMO,NumbreGPOTMO
	integer(1),dimension(:,:,:,:)::IGPOTIndex
	integer,dimension(:,:,:,:)::IGPOTZona
    real(8),dimension(:)::R,RO1X
    real(4),dimension(:,:)::RFILEPoltensG4
    real(8),dimension(:,:)::RFILEPoltensG8
    real(8),dimension(:,:,:)::RFunMO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8),dimension(:)::RPOT
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,KK,IIG1,IIG2,II1,II2,IDGX,IXSG,ITYG
    	   
	! мнлеп тюикю
    IDGX=NomeroFPOT(1)
    ! мнлеп рхою бгюхлндеиярбхъ
    ITYG=NumbreGPOTMO(NMO1,NMO2)

    ! пюявер налеммни вюярх йскнмнбяйнцн бгюхлндеиярбхъ
    ! мнлеп цюплнмхйх
	IIG1=0
	DO II1=LgarmonicMO(1,NMO2),LgarmonicMO(2,NMO2)
       IIG1=IIG1+1
       ! мнлеп цюплнмхйх
	   IIG2=0
	   DO II2=LgarmonicMO(1,NMO1),LgarmonicMO(2,NMO1)   
	      IIG2=IIG2+1
	      DO KK=0,2*LmaxGarm
	         IF(IGPOTIndex(ITYG,II1+1,II2+1,KK+1).EQ.1) THEN
			    ! пюявер онремжхюкю
			   	call POTENSS(DBLE(KK),NMO2,IIG1,NMO1,IIG2,RFunMO,Npoint,H,RPOT,R,RO1X)
				! нясыеярбкъел гюохяэ онксвеммнцн онремжхюкю
				! мнлеп гюохях
                IXSG=IGPOTZona(ITYG,II1+1,II2+1,KK+1)
                ! гюохяэ б люяяхб
				IF(NumbreSxemaCalculationG.EQ.1) THEN
                   RFILEPoltensG8(1:Npoint,IXSG)=RPOT(1:Npoint)
				ENDIF
				IF(NumbreSxemaCalculationG.EQ.2) THEN
                   RFILEPoltensG4(1:Npoint,IXSG)=SNGL(RPOT(1:Npoint))
				ENDIF
                ! гюохяэ б тюик
				IF(NumbreSxemaCalculationG.EQ.3) THEN
				   WRITE(IDGX,REC=IXSG) (RPOT(IIX),IIX=1,Npoint)  
			    ENDIF
			 ENDIF 
		  ENDDO 
	   ENDDO
	ENDDO
     
     	    
    return
  end subroutine Matrix_Enegry_Potential_Exchange_Coulomb_Interaction





  
  !! SUBPROGRAMME FORMATES THE MASSIVE INDICATING RECORDING ZONE AND MULTI-FITNESS NUMBER OF EXCHANGE INTERACTION POTENTIALS
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! NMO-NUMBER OF MOLECULAR ORBITAL FOR WHICH CALCULATION OF EXCHANGE INTERACTION IS CONDUCTED
══!! L1-ORBITAL MOMENT
══!! L2-ORBITAL MOMENT OF MOLECULAR ORBITAL
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! LmaxGarm-MAXIMUM HARMONIC LgarmonicMO (3, NumbreMO)
══!! LmaxKGarm-MAXIMUM HARMONIC LgarmonicMO (2, NumbreMO)
══!! NtypIndexCoff (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF EXCHANGE INTERACTION COEFFICIENTS
══!! RcoffBEI (NtypIndexCoff, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! IZONAMASSIVG (NumbreMO, NumbreMO) -MASSIVE IN WHICH MAXIMUM NUMBER OF RECORDING ZONES IS RECORDED
══!! NumbreGPOTMO (NumbreMO, NumbreMO) -MASSIVE OF INDEXES INDICATED TYPE (NUMBER) OF INTERACTION
══!! IGPOTIndex (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE OF KEYS DESCRIBING THE NECESSITY OF EXCHANGE POTENTIAL CALCULATION
══!! IGPOTZona (NumbreGPOTMO (NumbreMO, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), KK + 1) -MASSIVE RECORDING ZONES
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! IndexML (NtypIndexCoff, NumbreMO, NumbreMO, 4) -MASSIVE OF TOMORROW (MOLECULAR ORBITAL) PROJECTIONS FOR CALCULATING THE EXCHANGE PART OF CULON INTERACTION
══!! CkkGarmonik (2 * LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), LgarmonicMO (2, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO), IndexMLAB (NtypIndexCoff, NumbreMO)) - SPHERICAL HARMONIC array
  subroutine GPOT__Exchange_Coulomb_Interaction(NMO,L1,L2,NumbreMO,LmaxGarm,LmaxKGarm,NtypIndexCoff,RcoffBEI,IZONAMASSIVG,NumbreGPOTMO,IGPOTIndex,IGPOTZona,LgarmonicMO,IndexML,CkkGarmonik)
    implicit none
    integer::NMO,L1,L2,NumbreMO,LmaxGarm,LmaxKGarm
    integer,dimension(:,:)::LgarmonicMO,IZONAMASSIVG,NumbreGPOTMO,NtypIndexCoff
    integer(1),dimension(:,:,:,:)::IGPOTIndex
	integer,dimension(:,:,:,:)::IGPOTZona,IndexML
	real(8),dimension(:,:,:)::RcoffBEI
   	real(8),dimension(:,:,:,:,:)::CkkGarmonik
   	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IITyp,II1,II2,KK,MINN,MAXX,IIG1,IIG2,INCV,ierr,IndexML1,IndexML2,IndexML3,IndexML4,Lkk,IICGF,IZONA,IXDG
    real(8)::RGarm
    
  
        
	  ! пюявер налеммни вюярх йскнмнбяйнцн бгюхлндеиярбхъ
      DO IIX=1,NumbreMO 
	     ! опнбепъел ясыеярбсчр дюммшу йскнмнбяйхи хмрецпюк хкх мер
         IF(NtypIndexCoff(NMO,IIX).NE.0) THEN
		    ! рхо бгюхлндеиярбхъ
		    IXDG=NumbreGPOTMO(NMO,IIX)
	        
		    IZONA=IZONAMASSIVG(NMO,IIX)
            ! жхйк он цюплнмхйюл лнкейскъпмни нпахрюкх IIX
            DO II1=LgarmonicMO(1,IIX),LgarmonicMO(2,IIX)
               DO II2=LgarmonicMO(1,IIX),LgarmonicMO(2,IIX)   
		          MINN=MAX0(IABS(L1-II2),IABS(II1-L2))
                  MAXX=MIN0(L1+II2,II1+L2)
		  	      DO KK=MINN,MAXX,2
				     ! пюявер опнхгбедемхъ цюплнмхй 
                     RGarm=0.D0
				 	 ! жхйк он рхоюл йнщттхжхемрнб
	                 DO IITyp=1,NtypIndexCoff(NMO,IIX)
					    ! йнщттхжхемр йскнмнбяйнцн бгюхлндеиярбхъ  
		                IndexML1=IndexML(IITyp,NMO,IIX,1)
                        IndexML2=IndexML(IITyp,NMO,IIX,2)
			            IndexML3=IndexML(IITyp,NMO,IIX,3)
			  	        IndexML4=IndexML(IITyp,NMO,IIX,4)
						
					    ! пюявер опнхгбедемхъ цюплнмхй 
				        RGarm=RGarm+RcoffBEI(IITyp,NMO,IIX)*CkkGarmonik(KK+1,L2+1,II1+1,IndexML4,IndexML2)*CkkGarmonik(KK+1,L1+1,II2+1,IndexML1,IndexML3)
					 ENDDO
					 
					 IF(DABS(RGarm).GT.1.D-12) THEN
                        ! опнбепъел бшундхр гю пюлйх хкх мер
					    IF(KK.GT.(LmaxGarm+LmaxKGarm)) THEN 
                           WRITE(6,*) 'GPOT__Exchange_Coulomb_Interaction'
					  	   WRITE(6,*) 'Multy Garm KK>',LmaxGarm+LmaxKGarm,' KK= ',KK
						   STOP
						ENDIF 

					    ! опнбепъел тхйяхпнбюм дюммши онремжхюк хкх мер
					    IF(IGPOTIndex(IXDG,II1+1,L2+1,KK+1).NE.1) THEN
                           ! тхйяхпсел рн, врн дюммши онремжхюк сфе свюбярбнбюк б пюявере
                           IGPOTIndex(IXDG,II1+1,L2+1,KK+1)=1
					       ! гюохяшбюел мнлеп гнмш гюохях
						   IZONA=IZONA+1
						   IGPOTZona(IXDG,II1+1,L2+1,KK+1)=IZONA 
						ENDIF	     			          
                     ENDIF 
				  ENDDO 
			   ENDDO
	   	    ENDDO
            IZONAMASSIVG(NMO,IIX)=IZONA
		 ENDIF
      ENDDO

    
	
     	    
    return
  end subroutine GPOT__Exchange_Coulomb_Interaction


  !! SUBPROGRAMME FOR OBTAINING THE MATRIX OF THE OPERATOR OF THE CRYSTALLINE FIELD
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! L1-ORBITAL MOMENT OF THE FIRST ELECTRON
══!! M1-ORBITAL MOMENT OF THE FIRST ELECTRON
══!! L2-ORBITAL MOMENT OF THE SECOND ELECTRON
══!! M2-ORBITAL MOMENT OF THE SECOND ELECTRON
══!! Nnuclei-NUMBER OF NUCLEI (DO NOT CONSIDER NUCLEAR AT THE BEGINNING OF THE COORDINATE)
══!! Z (Nnuclei) -MASSIVE CHARGES OF NUCLEI
══!! COORRR (N) -MASSIVE RADIAL COORDINATES
══!! COORRR (N) -RADIAL COORDINATE OF THE NUCLEI (R-RADIUS)
══!! COOAAA (N, 2) -MASSIVE NORDER COORDINATE
══!! COOAAA (N, 1) -WIRTH COORDINATE -RTeta-ANGLE (0 = <RTeta <= PI)
══!! COOAAA (N, 2) - COORDINATE CO-ORDINATE -RFu-ANGLE (0 = <RFu <= 2 * PI)
══!! Npoint-NUMBER OF POINTS
══!! NMO-NUMBER OF MOLECULAR ORBITALS
══!! Ngarmon1-HARMONIC FIRST NUMBER
══!! Ngarmon2-HARMONIC NUMBER SECOND
══!! R (Npoint) -MASSIVE ARGUMENT VALUES
══!! Ncrystall (NMO, Ngarmon1, Ngarmon2, Npoint) -MASIVE OF ELECTRON INTERACTION OPERATOR WITH NUCLEI MOLECULES
══!! FF, F-AUXILIARY MASSIVE FOR CORNER HARMONIC CALCULATION
  subroutine MATRIX_CRYSTALLINE_FIELD_M(NMO,Ngarmon1,L1,M1,Ngarmon2,L2,M2,Nnuclei,Z,COORRR,COOAAA,Npoint,R,Ncrystall,FF,F)
    use mcmav,only:CMAV_ANGULAR_STRUCTURE_CRYSTALLINE_FIELD
    implicit none
    integer::NMO,Ngarmon1,L1,M1,Ngarmon2,L2,M2,Nnuclei,Npoint
    real(8)::FF
    real(8),dimension(:)::Z,R,F,COORRR
    real(8),dimension(:,:)::COOAAA
    real(8),dimension(:,:,:,:)::Ncrystall
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::ierr,IKNumbre,ISCMME,IKK,ISA,I1
	real(8)::Ral,ZZ,RcoffAS
    real(8),allocatable,dimension(:,:,:)::RcoffACF


    ! нясыеярбкъел тнплхпнбюмхе ярпсйрспш ноепюрнпю бгюхлндеиярбхъ щкейрпнмю я ъдпюлх ме мюундъыхлхяъ б мювюке йннпдхмюр
    ! гюмскъел оепед пюявернл
	DO ISCMME=1,Npoint 
       Ncrystall(NMO,Ngarmon1,Ngarmon2,ISCMME)=0.D0
	ENDDO
	  
	! опнбепъел мюкхвхе ъдеп б лнкейске
	IF(Nnuclei.EQ.0) THEN
       return
	ENDIF
	  
	 
	! СЯРЮМЮБКХБЮЕЛ ПЮГЛЕП ЛЮЯЯХБЮ
    IKNumbre=0
    DO ISCMME=IABS(L1-L2),L1+L2,2
       IKNumbre=IKNumbre+1
	ENDDO
    
    !бшдекъел оюлърэ дкъ люяяхбнб 
    allocate(RcoffACF(Nnuclei,IKNumbre,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MATRIX_CRYSTALLINE_FIELD_M'
       write(*,*) 'MEMORY ON THE FILE "RcoffACF" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
     
    ! гюмскъел оепед хяонкэгнбюмхел
    RcoffACF=0.D0
    
	! пюявхршбюел йнщттхжхемрш йпхярюккхвеяйнцн онкъ  
    call CMAV_ANGULAR_STRUCTURE_CRYSTALLINE_FIELD(Nnuclei,L1,M1,L2,M2,COOAAA,RcoffACF,FF,F)


    
	! опнбепъел бнгмхйюер кх йнлокейярмнярэ йнщттхжхемрю йпхярюккхвеяйнцн онкъ
	ISA=0
	DO ISCMME=1,Nnuclei
       DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
	      IF(DABS(RcoffACF(ISCMME,ISA,2)).GT.1.D-8) THEN
             WRITE(*,*) 'ATTENTION COEFFICIENT OF THE FIELD COMPLEX'
             WRITE(*,*) RcoffACF(ISCMME,ISA,2)
	         READ(*,*)
	         STOP  
	      ENDIF
       ENDDO
	ENDDO
    
	! нясыеярбкъел пюявер ноепюрнпю йпхярюккхвеяйнцн онкъ 
    DO ISCMME=1,Nnuclei
       Ral=COORRR(ISCMME)
       ZZ=Z(ISCMME)
	   ISA=0
	   DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
	   	  RcoffAS=RcoffACF(ISCMME,ISA,1)
		  DO I1=1,Npoint
             IF(R(I1).LE.Ral) THEN
                 Ncrystall(NMO,Ngarmon1,Ngarmon2,I1)=Ncrystall(NMO,Ngarmon1,Ngarmon2,I1)-ZZ*RcoffAS*(R(I1)/Ral)**IKK/Ral
                ELSE 
                 Ncrystall(NMO,Ngarmon1,Ngarmon2,I1)=Ncrystall(NMO,Ngarmon1,Ngarmon2,I1)-ZZ*RcoffAS*(Ral/R(I1))**IKK/R(I1)
	         ENDIF
          ENDDO     
	   ENDDO
    ENDDO

   
	          
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(RcoffACF,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CMME_ME_CRYSTALLINE_FIELD_M'
       write(*,*) 'THE FILE "RcoffACF" IS NOT REMOVED FROM MEMORY'
	   read(*,*)
	   stop 
    endif
      
      
    return
  end subroutine MATRIX_CRYSTALLINE_FIELD_M







 
  !! SUBPROGRAMME FOR OBTAINING THE MATRIX OF THE OPERATOR OF THE CRYSTALLINE FIELD
══!! DESCRIPTION OF SUBPROGRAMME PARAMETERS
══!! L1-ORBITAL MOMENT OF THE FIRST ELECTRON
══!! M1-ORBITAL MOMENT OF THE FIRST ELECTRON
══!! L2-ORBITAL MOMENT OF THE SECOND ELECTRON
══!! M2-ORBITAL MOMENT OF THE SECOND ELECTRON
══!! Nnuclei-NUMBER OF NUCLEI (DO NOT CONSIDER NUCLEAR AT THE BEGINNING OF THE COORDINATE)
══!! Z (Nnuclei) -MASSIVE CHARGES OF NUCLEI
══!! COORRR (N) -MASSIVE RADIAL COORDINATES
══!! COORRR (N) -RADIAL COORDINATE OF THE NUCLEI (R-RADIUS)
══!! COOAAA (N, 2) -MASSIVE NORDER COORDINATE
══!! COOAAA (N, 1) -WIRTH COORDINATE -RTeta-ANGLE (0 = <RTeta <= PI)
══!! COOAAA (N, 2) - COORDINATE CO-ORDINATE -RFu-ANGLE (0 = <RFu <= 2 * PI)
══!! Npoint-NUMBER OF POINTS
══!! Ngarmon1-HARMONIC FIRST NUMBER
══!! Ngarmon2-HARMONIC NUMBER SECOND
══!! R (Npoint) -MASSIVE ARGUMENT VALUES
══!! Ncrystall (Ngarmon1, Ngarmon2, Npoint) -MASIVE OF ELECTRON INTERACTION OPERATOR WITH NUCLEI MOLECULES
══!! FF, F-AUXILIARY MASSIVE FOR CORNER HARMONIC CALCULATION
  subroutine MATRIX_CRYSTALLINE_FIELD_MM(Ngarmon1,L1,M1,Ngarmon2,L2,M2,Nnuclei,Z,COORRR,COOAAA,Npoint,R,Ncrystall,FF,F)
    use mcmav,only:CMAV_ANGULAR_STRUCTURE_CRYSTALLINE_FIELD
    implicit none
    integer::NMO,Ngarmon1,L1,M1,Ngarmon2,L2,M2,Nnuclei,Npoint
    real(8)::FF
    real(8),dimension(:)::Z,R,F,COORRR
    real(8),dimension(:,:)::COOAAA
    real(8),dimension(:,:,:)::Ncrystall
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::ierr,IKNumbre,ISCMME,IKK,ISA,I1
	real(8)::Ral,ZZ,RcoffAS
    real(8),allocatable,dimension(:,:,:)::RcoffACF


    ! нясыеярбкъел тнплхпнбюмхе ярпсйрспш ноепюрнпю бгюхлндеиярбхъ щкейрпнмю я ъдпюлх ме мюундъыхлхяъ б мювюке йннпдхмюр
    ! гюмскъел оепед пюявернл
	DO ISCMME=1,Npoint 
       Ncrystall(Ngarmon1,Ngarmon2,ISCMME)=0.D0
	ENDDO
	  
	! опнбепъел мюкхвхе ъдеп б лнкейске
	IF(Nnuclei.EQ.0) THEN
       return
	ENDIF
	  
	 
	! СЯРЮМЮБКХБЮЕЛ ПЮГЛЕП ЛЮЯЯХБЮ
    IKNumbre=0
    DO ISCMME=IABS(L1-L2),L1+L2,2
       IKNumbre=IKNumbre+1
	ENDDO
    
    !бшдекъел оюлърэ дкъ люяяхбнб 
    allocate(RcoffACF(Nnuclei,IKNumbre,2),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MATRIX_CRYSTALLINE_FIELD_M'
       write(*,*) 'MEMORY ON THE FILE "RcoffACF" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
     
    ! гюмскъел оепед хяонкэгнбюмхел
    RcoffACF=0.D0
    
	! пюявхршбюел йнщттхжхемрш йпхярюккхвеяйнцн онкъ  
    call CMAV_ANGULAR_STRUCTURE_CRYSTALLINE_FIELD(Nnuclei,L1,M1,L2,M2,COOAAA,RcoffACF,FF,F)


    
	! опнбепъел бнгмхйюер кх йнлокейярмнярэ йнщттхжхемрю йпхярюккхвеяйнцн онкъ
	ISA=0
	DO ISCMME=1,Nnuclei
       DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
	      IF(DABS(RcoffACF(ISCMME,ISA,2)).GT.1.D-8) THEN
             WRITE(*,*) 'ATTENTION COEFFICIENT OF THE FIELD COMPLEX'
             WRITE(*,*) RcoffACF(ISCMME,ISA,2)
	         READ(*,*)
	         STOP  
	      ENDIF
       ENDDO
	ENDDO
    
	! нясыеярбкъел пюявер ноепюрнпю йпхярюккхвеяйнцн онкъ 
    DO ISCMME=1,Nnuclei
       Ral=COORRR(ISCMME)
       ZZ=Z(ISCMME)
	   ISA=0
	   DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
	   	  RcoffAS=RcoffACF(ISCMME,ISA,1)
		  DO I1=1,Npoint
             IF(R(I1).LE.Ral) THEN
                 Ncrystall(Ngarmon1,Ngarmon2,I1)=Ncrystall(Ngarmon1,Ngarmon2,I1)-ZZ*RcoffAS*(R(I1)/Ral)**IKK/Ral
                ELSE 
                 Ncrystall(Ngarmon1,Ngarmon2,I1)=Ncrystall(Ngarmon1,Ngarmon2,I1)-ZZ*RcoffAS*(Ral/R(I1))**IKK/R(I1)
	         ENDIF
          ENDDO     
	   ENDDO
    ENDDO

   
	          
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(RcoffACF,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CMME_ME_CRYSTALLINE_FIELD_M'
       write(*,*) 'THE FILE "RcoffACF" IS NOT REMOVED FROM MEMORY'
	   read(*,*)
	   stop 
    endif
      
      
    return
  end subroutine MATRIX_CRYSTALLINE_FIELD_MM


   !! Subroutine for calculating spherical harmonics (the array is filled for subsequent calculation)
══ !! LmaxGarmonik-MAXIMUM VALUE OF HARMONIC
══ !! NtipMOZ-NUMBER OF PROJECTION TYPES
══ !! MTipMOML (NtipMOZ) -MASSIVE OF PROJECTIONS
══ !! CkkGarmonik (2 * LmaxGarmonik, LmaxGarmonik, LmaxGarmonik, NtipMOZ, NtipMOZ) -MASSIVE OF SPHERICAL HARMONICS
══ !! FF, F-AUXILIARY MASIVE
  subroutine CalculSphericalGarmonik(LmaxGarmonik,NtipMOZ,MTipMOML,CkkGarmonik,FF,F)
    use mc3js,only:C3JS_Ckq
	implicit none
    integer::LmaxGarmonik,NtipMOZ
    real(8)::FF
	integer,dimension(:)::MTipMOML
    real(8),dimension(:)::F
	real(8),dimension(:,:,:,:,:)::CkkGarmonik 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::IXKX,IYKY,IZKZL1,IZKZL2,IZKZKK

    ! жхйк он рхоюл опнейжхи
    DO IXKX=1,NtipMOZ
	   DO IYKY=1,NtipMOZ
	      ! жхйк он нпахрюкэмшл вхякюл
		  DO IZKZL1=0,LmaxGarmonik
		     DO IZKZL2=0,LmaxGarmonik
                ! жхйк он лскэрхокермнярх ятепхвеяйни цюплнмхйх
				DO IZKZKK=0,2*LmaxGarmonik !IZKZL1+IZKZL2
				   !WRITE(*,*) IZKZL1,IZKZL2,IZKZKK,CkkGarmonik(IZKZKK+1,IZKZL2+1,IZKZL1+1,IYKY,IXKX)
                   CkkGarmonik(IZKZKK+1,IZKZL2+1,IZKZL1+1,IYKY,IXKX)=C3JS_Ckq(IZKZL1,MTipMOML(IXKX),IZKZL2,MTipMOML(IYKY),IZKZKK,FF,F)  
				   !WRITE(*,*) IZKZL1,IZKZL2,IZKZKK,CkkGarmonik(IZKZKK+1,IZKZL2+1,IZKZL1+1,IYKY,IXKX)
				   !READ(*,*)
                   ! сярюмюбкхбюел нйнмвюмхе жхйкю
				   !IF(IZKZKK.GT.2) THEN
				   !  IF(DABS(CkkGarmonik(IZKZKK+1,IZKZL2+1,IZKZL1+1,IYKY,IXKX)).LT.10.D-10.AND.DABS(CkkGarmonik(IZKZKK,IZKZL2+1,IZKZL1+1,IYKY,IXKX)).LT.10.D-10) THEN
                   !     IF(DABS(CkkGarmonik(IZKZKK-1,IZKZL2+1,IZKZL1+1,IYKY,IXKX)).LT.10.D-10) THEN
                   !       EXIT
				   !	    ENDIF
				   !  ENDIF
                   !ENDIF 
				ENDDO 
			 ENDDO
		  ENDDO 
	   ENDDO
	ENDDO 

	return
  end subroutine CalculSphericalGarmonik



  !! SUBPROGRAMME IMPROVES THE FORMATION OF THE STRUCTURE OF THE HIGHER HARMONICS OF THE MOLECULAR ORBITAL
══!! NMO-NUMBER OF MOLECULAR ORBITALS
══!! Npoint-NUMBER OF POINTS
══!! H-STEP
══!! NumbreGarmMO-NUMBER HARMONIC OF MOLECULAR ORBITAL TO LK INCLUSIVE
══!! NumbreGarmMOLimid-NUMBER HARMONIC OF MOLECULAR ORBITAL AFTER LK TO MAXIMUM HARMONIC
══!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
══!! RFunMO (NMO, NumbreGarmMO, Npoint + 2) -mass of the values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
══!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
══!! RO1X (Npoint) - AUXILIARY MASSIVE FIRST DERIVATIVE NEW VARIABLE BY OLD
  subroutine FormStructureMOFull(NMO,Npoint,H,NumbreGarmMO,NumbreGarmMOLimid,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFunMO,RfunLigands,AlfaCoffMO,RO1X) 
    use msde,only:CalculationNormalizationConstaResult
	implicit none
    integer::NMO,Npoint,NumbreGarmMO,NumbreGarmMOLimid
	real(8)::H
    integer,dimension(:)::NumbreLigand
    integer,dimension(:,:)::NLigands,NFunLigands
    integer,dimension(:,:,:)::NumbreFunctionLig
	real(8),dimension(:)::RO1X
	real(8),dimension(:,:,:)::RFunMO,AlfaCoffMO
    real(8),dimension(:,:,:,:)::RfunLigands
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIXDGAR,IX1,IX2,IX3,NumbreLigandX,NumbreFunLigX,ierr
	real(8)::Rcoff,RRnormcoffZXZ
    real(8),allocatable,dimension(:)::RfunXD
    
	! бшдекъел оюлърэ онд люяяхб
	allocate(RfunXD(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*)'FormStructureMOZero'
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
       RFunMO(NMO,IIXDGAR,2+1:Npoint)=RfunXD(1:Npoint)
	ENDDO


    ! мнплхпсел онксвеммсч тсмйжхч
    RRnormcoffZXZ=CalculationNormalizationConstaResult(NMO,1,NumbreGarmMOLimid,Npoint,H,RO1X,RFunMO) 
    !WRITE(*,*) NMO,NumbreGarmMOLimid,RRnormcoffZXZ
	!READ(*,*) 
    ! гюохяшбюел  мнплхпнбнвмши йнщттхжхемр
    DO IX1=1,NumbreGarmMOLimid
       RFunMO(NMO,IX1,2)=RRnormcoffZXZ
    ENDDO 

	! сдюкъел люяяхбш хг оюлърх 
    deallocate(RfunXD,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'FormStructureMOZero'
       write(*,*) 'THE MASSIV "RfunXD" IS NOT REMOVED FROM MEMORY'
       read(*,*)
	   stop 
    endif 

    return
  end subroutine FormStructureMOFull


  !! SUBPROGRAMME DEFINING ALPHA COEFFICIENT INCLUDING MOLECULAR ORBITAL STRUCTURE
══!! NMO-NUMBER OF MOLECULAR ORBITALS
══!! Npoint-NUMBER OF POINTS
══!! H-STEP
══!! NumbreGarmMO-NUMBER HARMONIC OF MOLECULAR ORBITAL TO LK INCLUSIVE
══!! NumbreLigand (NumbreMO) -MASSIVE NUMBER OF LIGANDES WHICH FUNCTIONS DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NLigands (NumbreMO, NumbreLigand (NumbreMO)) - MASSIVE OF LIGAND NUMBERS FUNCTIONS WHICH DESCRIBE HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NFunLigands (NumbreMO, NumbreLigand (NumbreMO)) - A MASSIVE OF THE NUMBER OF FUNCTIONS OF EACH LIGAND DESCRIBING HIGHER HARMONICS OF THIS MOLECULAR ORBITAL
══!! NumbreFunctionLig (NumbreMO, NumbreLigand (NumbreMO), NFunLigands (NumbreMO, NumbreLigand (NumbreMO)))) - MASSIVE OF NUMBERS OF LIGAND FUNCTIONS THAT COMPLY WITH THIS MOLECULAR ORBITAL
══!! RFunMO (NMO, NumbreGarmMO, Npoint + 2) -mass of the values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! RfunLigands (Nnuclei, NumbreFunctionLig, NumbreGarmLigand, Npoint + 2) -MASSIVE OF VALUES OF LIGAND FUNCTIONS
══!! AlfaCoffMO (NumbreMO, NumbreLigand, NumbreFunLigand) -MASSIVE OF COEFFICIENTS WITH WHICH LIGANDA FUNCTIONS COME INTO THE STRUCTURE OF MOLECULAR ORBITAL
══!! RO1X (Npoint) - AUXILIARY MASSIVE FIRST DERIVATIVE NEW VARIABLE BY OLD
  subroutine CalculationCoffAlfaMOFull(NMO,Npoint,H,NumbreGarmMO,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,RFunMO,RfunLigands,AlfaCoffMO,RO1X) 
    implicit none
    integer::NMO,Npoint,NumbreGarmMO
	real(8)::H
    integer,dimension(:)::NumbreLigand
    integer,dimension(:,:)::NLigands,NFunLigands
    integer,dimension(:,:,:)::NumbreFunctionLig
	real(8),dimension(:)::RO1X 
	real(8),dimension(:,:,:)::RFunMO,AlfaCoffMO
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
	   RintMO=sum((RFunMO(NMO,NumbreGMOAS,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS,2) 
	   RintFL=sum(RFunMO(NMO,NumbreGMOAS,2+1:Npoint)*RfunLigands(NumbreLigandVV,NumbreFunLigZX,NumbreGMOAS,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS,2)*RfunLigands(NumbreLigandVV,NumbreFunLigZX,NumbreGMOAS,2)) 
	   ! пюявхршбюел йнщттхжхемр юкэтю
	   AlfaCoffMO(NMO,NumbreLigandVV,NumbreFunLigZX)=RintMO/RintFL
	   !WRITE(*,*) 'ALFA',NMO,NumbreLigandVV,NumbreFunLigZX,RintMO,RintFL,AlfaCoffMO(NMO,NumbreLigandVV,NumbreFunLigZX)
	   !READ(*,*)
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
       RintMO1=sum((RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS1,2)        
       RintMO2=sum((RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint))**2)*H !/RFunMO(NumbreGMOAS2,2)
       RintFL11=sum(RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS1,2)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2)) 
	   RintFL12=sum(RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS1,2)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2))
	   RintFL21=sum(RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS2,2)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2)) 
	   RintFL22=sum(RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H !/DSQRT(RFunMO(NumbreGMOAS2,2)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2))
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
	   RintMO1=sum((RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintMO2=sum((RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintMO3=sum((RFunMO(NMO,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint))**2)*H 
       RintFL11=sum(RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H  
	   RintFL12=sum(RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL13=sum(RFunMO(NMO,NumbreGMOAS1,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS1,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL21=sum(RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL22=sum(RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL23=sum(RFunMO(NMO,NumbreGMOAS2,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS2,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL31=sum(RFunMO(NMO,NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV1,NumbreFunLigZX1,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL32=sum(RFunMO(NMO,NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV2,NumbreFunLigZX2,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
	   RintFL33=sum(RFunMO(NMO,NumbreGMOAS3,2+1:Npoint)*RfunLigands(NumbreLigandVV3,NumbreFunLigZX3,NumbreGMOAS3,2+1:Npoint)/RO1X(1:Npoint)**2)*H 
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
  end subroutine CalculationCoffAlfaMOFull



  !! Subroutine for calculating the total energy of a molecule configuration
══!! IndexHartreeFock-PARAMETER INDICATING WHAT EQUATION IS DECIDED
══!! IndexHartreeFock = 0-HARTRI EQUATION
══!! IndexHartreeFock = 1-HARTRI-FOCA EQUATION
══!! Nnuclei-NUMBER OF NUCLEI (DO NOT CONSIDER NUCLEAR AT THE BEGINNING OF THE COORDINATE)
══!! Z0-CHARGE OF THE NUCLEI AT THE BEGINNING OF THE COORDINATE
══!! Z (Nnuclei) -MASSIVE CHARGES OF NUCLEI
══!! COORRR (N) -MASSIVE RADIAL COORDINATES
══!! COORRR (N) -RADIAL COORDINATE OF THE NUCLEI (R-RADIUS)
══!! COOAAA (N, 2) -MASSIVE NORDER COORDINATE
══!! COOAAA (N, 1) -WIRTH COORDINATE -RTeta-ANGLE (0 = <RTeta <= PI)
══!! COOAAA (N, 2) - COORDINATE CO-ORDINATE -RFu-ANGLE (0 = <RFu <= 2 * PI)
══!! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
══!! NtypIndexCoffF (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF DIRECT INTERACTION COEFFICIENTS
══!! NtypIndexCoffG (NumbreMO, NumbreMO) -MASSIVE NUMBER OF TYPES OF EXCHANGE INTERACTION COEFFICIENTS
══!! NN (NumbreMO) -MASSIVE OF MAIN QUANTUM NUMBERS
══!! ML (NumbreMO) -MASSIVE OF MOLECULAR ORBITAL PROJECTIONS
══!! RcoffADI (NtypIndexCoffF, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION FACTORS
══!! RcoffBEI (NtypIndexCoffG, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
══!! IQ (NumbreMO) -Chapter numbers (number of electrons on a molecular orbital)
══!! NumbreGarmMO (NumbreMO) -MASSIVE NUMBER OF HARMONIC MOLECULAR ORBITALS
══!! LgarmonicMO (2, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
══!! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
══!! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
══!! Npoint-NUMBER OF POINTS
══!! H-step
══!! Ncrystall (NumbreMO, NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -MASIVE OF ELECTRON INTERACTION OPERATOR WITH NUCLEI OF MOLECULES
══!! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values ??of the radial parts of the harmonics of the molecular orbitals of the configuration
══!! R (Npoint) -MASSIVE VALUES OF RADIUS
══!! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
══!! IndexMLA (NtypIndexCoffF, NumbreMO, NumbreMO, 4) -MASSIVE OF TYPES OF Moment (MOLECULAR ORBITAL) PROJECTIONS USED IN CALCULATING DIRECT INTERACTION
══!! IndexMLB (NtypIndexCoffG, NumbreMO, NumbreMO, 4) -MASSIVE OF TOMORROW TYPE (MOLECULAR ORBITAL) TYPES USED IN CALCULATING THE EXCHANGE INTERACTION
══!! CkkGarmonik (KK + 1, L2 + 1, L1 + 1, IML, IML1) -MASSIVE OF MATRIX ELEMENTS OF SPHERICAL HARMONIC
  subroutine TOTALXF(IndexHartreeFock,Nnuclei,Z0,Z,COORRR,COOAAA,NumbreMO,NtypIndexCoffF,NtypIndexCoffG,NN,ML,RcoffADI,RcoffBEI,IQ,NumbreGarmMO,LgarmonicMO,Npoint,H,Ncrystall,RFunMO,R,RO1X,IndexMLA,IndexMLB,CkkGarmonik,NomeroFPOT)
    implicit none
	integer::IndexHartreeFock,Nnuclei,NumbreMO,Npoint
	real(8)::Z0,H
    integer,dimension(:)::NN,ML,IQ,NumbreGarmMO,NomeroFPOT
    integer,dimension(:,:)::NtypIndexCoffF,NtypIndexCoffG,LgarmonicMO
    integer,dimension(:,:,:,:)::IndexMLA,IndexMLB
	real(8),dimension(:)::Z,COORRR,R,RO1X
	real(8),dimension(:,:)::COOAAA
	real(8),dimension(:,:,:)::RcoffADI,RcoffBEI,RFunMO
	real(8),dimension(:,:,:,:)::Ncrystall
	real(8),dimension(:,:,:,:,:)::CkkGarmonik 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIX,IIY,IIZ,ierr
	real(8)::EnergyNuclear,SumFCoulomb,SumGCoulomb,EPSFull,EnergyCrystallFieldX
    real(8),allocatable,dimension(:)::IntCrystallField
    real(8),allocatable,dimension(:,:,:)::FCoulombInt
    real(8),allocatable,dimension(:,:,:)::GCoulombInt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::LmaxGarmXX,IzonFmaxXX,IzonGmaxXX,IRsumXXX,NumbreSxemaCalculationFXX,NumbreSxemaCalculationGXX
	integer::IXXXI,IYYYI,IIL1I,IIL2I,IOX1X,IOX1I,IOX2I,NtypIndexCoffFmax,NtypIndexCoffGmax
	integer,allocatable,dimension(:)::IZONAMASSIVXX
    integer,allocatable,dimension(:,:)::IZONAMASSIVGXX,NumbreGPOTMOXX
	integer(1),allocatable,dimension(:,:,:,:)::IFPOTIndexXX,IGPOTIndexXX
	integer,allocatable,dimension(:,:,:,:)::IFPOTZonaXX,IGPOTZonaXX
	real(8),allocatable,dimension(:)::RPOTX
	real(4),allocatable,dimension(:,:)::RFILEPoltensGXX4 
    real(8),allocatable,dimension(:,:)::RFILEPoltensGXX8
   	real(4),allocatable,dimension(:,:)::RFILEPoltensFXX4
	real(8),allocatable,dimension(:,:)::RFILEPoltensFXX8
	character(5),dimension(2)::LB  
    data LB/'sigma','pi'/
	
    2000 FORMAT(2X,'  F',I2,'(',2(I2,A5),')=',F10.7)
    2100 FORMAT(2X,'  G',I2,'(',2(I2,A5),')=',F10.7)
	2200 FORMAT(2X,'Direct Coulomb interaction in a.u.') 
    2300 FORMAT(2X,'Exchange Coulomb interaction in a.u.')
	2400 FORMAT(2X,'Type of a Coulomb interaction N= ',I2) 
	2500 FORMAT(5X,'------------------------')
	2550 FORMAT(5X,'                         ------------------------')
	2555 FORMAT(5X,'                                   ------------------------')
    2600 FORMAT(5X,'ENERGY COMPONENTS IN A.U.')
    2700 FORMAT(5X,'ONE ELECTRON ENERGY      = ',F18.10)
    2800 FORMAT(5X,'TWO ELECTRON ENERGY      = ',F18.10) 
	2850 FORMAT(5X,'NUCLEAR REPULSION ENERGY = ',F18.10)
	2900 FORMAT(5X,'TOTAL ENERGY             =',F18.10)
	3000 FORMAT(5X,'ELECTRON-ELECTRON POTENTIAL ENERGY =',F18.10)
	3100 FORMAT(5X,'NUCLEUS-ELECTRON POTENTIAL ENERGY  =',F18.10)
    3200 FORMAT(5X,'NUCLEUS-NUCLEUS POTENTIAL ENERGY   =',F18.10)
    3300 FORMAT(5X,'TOTAL POTENTIAL ENERGY             =',F18.10)
    3400 FORMAT(5X,'TOTAL KINETIC ENERGY               =',F18.10)
    3500 FORMAT(5X,'VIRIAL RATIO (V/T)                 =',F18.10)

    
	
	
	
	
	! нопедекъел люйяхлюкэмне гмювемхе вхякю рхонб йскнмнбяйху хмрецпюкнб 
	NtypIndexCoffFmax=0
    NtypIndexCoffGmax=0 
	DO IXXXI=1,NumbreMO
	   DO IYYYI=1,NumbreMO
          IF(NtypIndexCoffF(IXXXI,IYYYI).GT.NtypIndexCoffFmax) THEN
             NtypIndexCoffFmax=NtypIndexCoffF(IXXXI,IYYYI)
		  ENDIF
          IF(NtypIndexCoffG(IXXXI,IYYYI).GT.NtypIndexCoffGmax) THEN
             NtypIndexCoffGmax=NtypIndexCoffG(IXXXI,IYYYI)
		  ENDIF
	   ENDDO
    ENDDO
    
	! опнбепъел мюкхвхе нркхвмнцн нр мскъ вхякю рхонб 
	! тхйяхпсел лхмхлюкэмне гмювемхе дкъ нрясрярбхъ ньханй опх бшдекемхх оюлърх дкъ люяяхбнб
	IF(NtypIndexCoffFmax.EQ.0) THEN
	   NtypIndexCoffFmax=1
	ENDIF 
    IF(NtypIndexCoffGmax.EQ.0) THEN
	   NtypIndexCoffGmax=1
	ENDIF 
	

	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RPOTX(Npoint),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'TOTALXF'
	   write(*,*) 'MEMORY ON THE FILE "RPOTX" IS NOT SELECTED'
	   read(*,*)
	   stop 
	endif  
	allocate(FCoulombInt(NtypIndexCoffFmax,NumbreMO,NumbreMO),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'TOTALXF'
	   write(*,*) 'MEMORY ON THE FILE "FCoulombInt" IS NOT SELECTED'
	   read(*,*)
	   stop 
	endif  
    allocate(GCoulombInt(NtypIndexCoffGmax,NumbreMO,NumbreMO),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'TOTALXF'
	   write(*,*) 'MEMORY ON THE FILE "GCoulombInt" IS NOT SELECTED'
	   read(*,*)
	   stop 
	endif 
	allocate(IntCrystallField(NumbreMO),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'TOTALXF'
	   write(*,*) 'MEMORY ON THE FILE "IntCrystallField" IS NOT SELECTED'
	   read(*,*)
	   stop 
	endif   


   
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! нясыеярбкъел юмюкхг дюммшу дкъ тнплхпнбюмхъ люяяхбнб гюохях онремжхюкнб! 
	! дкъ пюяверю онкмни щмепцхх йнмтхцспюжххх                                              !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!             нясыеярбкъел пюявер опълни вюярх йскнмнбяйнцн бгюхлндеиярбхъ	            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! нопедекъел люйяхлюкэмше цюплнмхйх
	LmaxGarmXX=0
	DO IYYYI=1,NumbreMO  
	   IF(LmaxGarmXX.LT.LgarmonicMO(2,IYYYI)) THEN
          LmaxGarmXX=LgarmonicMO(2,IYYYI)  
	   ENDIF
	ENDDO

	  	

	! бшдекъел люяяхбш гмювемхи дкъ тнплхпнбюмхъ тюикнб онремжхюкнб опълнцн бгюхлндеиярбхъ 
    allocate(IFPOTIndexXX(NumbreMO,LmaxGarmXX+1,LmaxGarmXX+1,2*LmaxGarmXX+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'MEMORY ON THE FILE "IFPOTIndexXX" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
	allocate(IFPOTZonaXX(NumbreMO,LmaxGarmXX+1,LmaxGarmXX+1,2*LmaxGarmXX+1),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'MEMORY ON THE FILE "IFPOTZonaXX" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
    allocate(IZONAMASSIVXX(NumbreMO),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'MEMORY ON THE FILE "IZONAMASSIVXX" IS NOT SELECTED'
       read(*,*)
	   stop 
    endif 
    

    
	

    ! гюмскъел оепед пюявернл
	IZONAMASSIVXX=0
    IFPOTIndexXX=0
    ! тнплхпсел люяяхб сйюгшбючыхи гнмш гюохях х цюплнмхйх опълнцн бгюхлндеиярбхъ х налеммнцн бгюхлндеиярбхъ  
    DO IXXXI=1,NumbreMO   
	   DO IIL1I=LgarmonicMO(1,IXXXI),LgarmonicMO(2,IXXXI) 
         DO IIL2I=IIL1I,LgarmonicMO(2,IXXXI) 
	        ! пюявер опълни вюярх
	         call FPOT_Energy_Direct_Coulomb_Interaction(IXXXI,IIL1I,IIL2I,NumbreMO,LmaxGarmXX,NtypIndexCoffF,RcoffADI,IZONAMASSIVXX,IFPOTIndexXX,IFPOTZonaXX,LgarmonicMO,IndexMLA,CkkGarmonik)
		  ENDDO
       ENDDO
    ENDDO

    ! яуелю пюяверю опх йнрнпни хяонкэгсеряъ люяяхб
   	NumbreSxemaCalculationFXX=1
    
    ! нопедекъел люйяхлюкэмне вхякн гюохяеи опълни х налеммни вюярх
	! опълюъ вюярэ 
	IzonFmaxXX=1
	DO IXXXI=1,NumbreMO  
       IF(IzonFmaxXX.LT.IZONAMASSIVXX(IXXXI)) THEN
          IzonFmaxXX=IZONAMASSIVXX(IXXXI) 
	   ENDIF
    ENDDO
    
	
    ! бшдекъел оюлърэ онд люяяхб опълшу онремжхюкнб 
    allocate(RFILEPoltensFXX8(Npoint,IzonFmaxXX),stat=ierr) 
    if(ierr/=0) then
       write(6,*) 'TOTALXF'
       write(6,*) 'MEMORY ON THE FILE "RFILEPoltensFXX8" IS NOT SELECTED'
       write(6,*) 'For allocation of the array',FLOAT(IzonFmaxXX)*FLOAT(Npoint)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	   ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю бшдекемхе оюлърх он люяяхб REAL(4)
	   NumbreSxemaCalculationFXX=2 
    endif
	IF(NumbreSxemaCalculationFXX.EQ.2) THEN   
       allocate(RFILEPoltensFXX4(Npoint,IzonFmaxXX),stat=ierr) 
       if(ierr/=0) then
          write(6,*) 'TOTALXF'
          write(6,*) 'MEMORY ON THE FILE "RFILEPoltensFXX4" IS NOT SELECTED'
          write(6,*) 'For allocation of the array',FLOAT(IzonFmaxXX)*FLOAT(Npoint)*4.D0/(1024.D0*1024.D0),'Mb are necessary'
	      ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б тюикш
	      NumbreSxemaCalculationFXX=3 
       endif
    ENDIF 
   


    
    ! гюмскъел оепед пюявернл
    FCoulombInt=0.D0
    GCoulombInt=0.D0
	IntCrystallField=0.D0

	! щрюо 1. пюявер йскнмнбяйху хмрецпюкнб
    ! пюявер опълшу хмрецпюкнб (б юрнлмшу едхмхжюу)
	DO IIY=1,NumbreMO 
       DO IIZ=IIY,NumbreMO
          DO IIX=1,NtypIndexCoffF(IIY,IIZ) 
		     IF(DABS(RcoffADI(IIX,IIY,IIZ)).GT.1.D-6) THEN
			    WRITE(*,*) 'Calculation FCoulombInt Ntype= ',IIX,' NMO1= ',IIY, 'NMO2= ',IIZ
				! нясыеярбкъел пюявер онремжхюкнб опълни вюярх 
                WRITE(*,*) 'Matrix_Potential_Direct FINAL NMO=',IIY
                call Matrix_Energy_Potential_Direct_Coulomb_Interaction(IIY,Npoint,H,LmaxGarmXX,NomeroFPOT,LgarmonicMO,IFPOTIndexXX,IFPOTZonaXX,RFunMO,R,RO1X,RPOTX,NumbreSxemaCalculationFXX,RFILEPoltensFXX8,RFILEPoltensFXX4)
                FCoulombInt(IIX,IIY,IIZ)=FCoulombMO(IIX,IIY,IIZ,LgarmonicMO,Npoint,H,RFunMO,R,RO1X,IndexMLA,CkkGarmonik,NomeroFPOT,IFPOTIndexXX,IFPOTZonaXX,NumbreSxemaCalculationFXX,RFILEPoltensFXX8,RFILEPoltensFXX4) 
				FCoulombInt(IIX,IIZ,IIY)=FCoulombInt(IIX,IIY,IIZ)
			 ENDIF
		  ENDDO  
	   ENDDO
	ENDDO
    
    ! сдюкъел мемсфмше люяяхбш  
	IF(NumbreSxemaCalculationFXX.EQ.1) THEN
	   deallocate(RFILEPoltensFXX8,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "RFILEPoltensFXX8" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
	ENDIF
	IF(NumbreSxemaCalculationFXX.EQ.2) THEN
	   deallocate(RFILEPoltensFXX4,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "RFILEPoltensFXX4" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
	ENDIF
    deallocate(IFPOTIndexXX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "IFPOTIndexXX" IS NOT REMOVED FROM MEMORY' 
       read(*,*)
   	   stop 
    endif
    deallocate(IFPOTZonaXX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "IFPOTZonaXX" IS NOT REMOVED FROM MEMORY' 
       read(*,*)
   	   stop 
    endif
    deallocate(IZONAMASSIVXX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "IZONAMASSIVXX" IS NOT REMOVED FROM MEMORY' 
       read(*,*)
   	   stop 
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            нясыеярбкъел пюявер налеммни вюярх йскнмнбяйнцн бгюхлндеиярбхъ	            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! опнбепъел рхо спюбмемхъ
	IF(IndexHartreeFock.EQ.1) THEN 
	   ! бшдекъел оюлърэ онд люяяхб
	   allocate(NumbreGPOTMOXX(NumbreMO,NumbreMO),stat=ierr)
       if(ierr/=0) then
          write(*,*) 'DHFMO'
          write(*,*) 'MEMORY ON THE FILE "NumbreGPOTMOXX" IS NOT SELECTED'
          read(*,*)
	      stop 
       endif 

       ! тнплхпсел люяяхб хмдейя рхош бгюхлндеиярбхъ 
	   NumbreGPOTMOXX=0
	   IRsumXXX=0
	   DO IXXXI=1,NumbreMO  
          DO IYYYI=1,NumbreMO 
	         IF(NtypIndexCoffG(IXXXI,IYYYI).NE.0) THEN
		        IRsumXXX=IRsumXXX+1
	       	    NumbreGPOTMOXX(IXXXI,IYYYI)=IRsumXXX
	         ENDIF
	      ENDDO
	   ENDDO

  
	   ! бшдекъел люяяхбш гмювемхи дкъ тнплхпнбюмхъ тюикнб онремжхюкнб налеммнцн бгюхлндеиярбхъ 
       allocate(IGPOTIndexXX(IRsumXXX,LmaxGarmXX+1,LmaxGarmXX+1,2*LmaxGarmXX+1),stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'MEMORY ON THE FILE "IGPOTIndexXX" IS NOT SELECTED'
          read(*,*)
	      stop 
       endif 
	   allocate(IGPOTZonaXX(IRsumXXX,LmaxGarmXX+1,LmaxGarmXX+1,2*LmaxGarmXX+1),stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'MEMORY ON THE FILE "IGPOTZonaXX" IS NOT SELECTED'
          read(*,*)
	      stop 
       endif 
	   allocate(IZONAMASSIVGXX(NumbreMO,NumbreMO),stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'MEMORY ON THE FILE "IZONAMASSIVGXX" IS NOT SELECTED'
          read(*,*)
	      stop 
       endif 
  

       IGPOTIndexXX=0
       IZONAMASSIVGXX=0

       DO IXXXI=1,NumbreMO   
          DO IIL1I=LgarmonicMO(1,IXXXI),LgarmonicMO(2,IXXXI) 
             DO IIL2I=LgarmonicMO(1,IXXXI),LgarmonicMO(2,IXXXI) 
	            ! пюявер налеммни вюярх
	  	        call GPOT__Exchange_Coulomb_Interaction(IXXXI,IIL1I,IIL2I,NumbreMO,LmaxGarmXX,LmaxGarmXX,NtypIndexCoffG,RcoffBEI,IZONAMASSIVGXX,NumbreGPOTMOXX,IGPOTIndexXX,IGPOTZonaXX,LgarmonicMO,IndexMLB,CkkGarmonik)
		     ENDDO
          ENDDO
 	   ENDDO


       
	   ! яуелю пюяверю опх йнрнпни хяонкэгсеряъ люяяхб
   	   NumbreSxemaCalculationGXX=1
   
       ! налеммюъ вюярэ
	   IzonGmaxXX=1
	   DO IXXXI=1,NumbreMO  
          DO IYYYI=1,NumbreMO 
	         IF(IXXXI.NE.IYYYI) THEN
	            IF(IzonGmaxXX.LT.IZONAMASSIVGXX(IXXXI,IYYYI)) THEN
			       IzonGmaxXX=IZONAMASSIVGXX(IXXXI,IYYYI)
			    ENDIF
	         ENDIF
	      ENDDO
	   ENDDO
	  	
    
	   ! бшдекъел оюлърэ онд люяяхб налеммшу онремжхюкнб 
       allocate(RFILEPoltensGXX8(Npoint,IzonGmaxXX),stat=ierr) 
       if(ierr/=0) then
          write(6,*) 'TOTALXF'
          write(6,*) 'MEMORY ON THE FILE "RFILEPoltensGXX8" IS NOT SELECTED'
          write(6,*) 'For allocation of the array',FLOAT(Npoint)*FLOAT(IzonGmaxXX)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	      ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б люяяхб REAL(4)
	      NumbreSxemaCalculationGXX=2 
       endif    
       IF(NumbreSxemaCalculationGXX.EQ.2) THEN
          allocate(RFILEPoltensGXX4(Npoint,IzonGmaxXX),stat=ierr) 
          if(ierr/=0) then
             write(6,*) 'TOTALXF'
             write(6,*) 'MEMORY ON THE FILE "RFILEPoltensGXX4" IS NOT SELECTED'
             write(6,*) 'For allocation of the array',FLOAT(Npoint)*FLOAT(IzonGmaxXX)*8.D0/(1024.D0*1024.D0),'Mb are necessary'
	         ! дкъ люяяхбю мебнглнфмн бшдекхрэ рюйне йнкхвеярбн оюлърх нясыеярбкъел оепеунд мю гюохяэ б тюикш
	         NumbreSxemaCalculationGXX=3 
          endif   
       ENDIF




	   ! пюявер налеммшу хмрецпюкнб (б юрнлмшу едхмхжюу)
       DO IIY=1,NumbreMO 
          DO IIZ=IIY,NumbreMO
	         DO IIX=1,NtypIndexCoffG(IIY,IIZ) 
		        IF(DABS(RcoffBEI(IIX,IIY,IIZ)).GT.1.D-6) THEN
			       WRITE(*,*) 'Matrix_Potential_Exchange FINAL NMO=',IIZ,IIY
    	           call Matrix_Enegry_Potential_Exchange_Coulomb_Interaction(IIZ,IIY,Npoint,H,LmaxGarmXX,NomeroFPOT,LgarmonicMO,NumbreGPOTMOXX,IGPOTIndexXX,IGPOTZonaXX,RFunMO,R,RO1X,RPOTX,NumbreSxemaCalculationGXX,RFILEPoltensGXX8,RFILEPoltensGXX4)
                   WRITE(*,*) 'Calculation GCoulombInt Ntype= ',IIX,' NMO1= ',IIY, 'NMO2= ',IIZ
				   GCoulombInt(IIX,IIY,IIZ)=GCoulombMO(IIX,IIY,IIZ,LgarmonicMO,Npoint,H,RFunMO,R,RO1X,IndexMLB,CkkGarmonik,NumbreGPOTMOXX,NomeroFPOT,IGPOTIndexXX,IGPOTZonaXX,NumbreSxemaCalculationGXX,RFILEPoltensGXX8,RFILEPoltensGXX4) 
				   GCoulombInt(IIX,IIZ,IIY)=GCoulombInt(IIX,IIY,IIZ)
			    ENDIF
		     ENDDO  
	      ENDDO
	   ENDDO

       ! сдюкъел мемсфмше люяяхбш 
       IF(NumbreSxemaCalculationGXX.EQ.1) THEN  
	      deallocate(RFILEPoltensGXX8,stat=ierr)
          if(ierr/=0) then
             write(*,*) 'TOTALXF'
             write(*,*) 'THE FILE "RFILEPoltensGXX8" IS NOT REMOVED FROM MEMORY' 
             read(*,*)
   	         stop 
          endif
       ENDIF
       IF(NumbreSxemaCalculationGXX.EQ.2) THEN  
	      deallocate(RFILEPoltensGXX4,stat=ierr)
          if(ierr/=0) then
             write(*,*) 'TOTALXF'
             write(*,*) 'THE FILE "RFILEPoltensGXX4" IS NOT REMOVED FROM MEMORY' 
             read(*,*)
   	         stop 
          endif
       ENDIF
	   deallocate(NumbreGPOTMOXX,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "NumbreGPOTMOXX" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
       deallocate(IGPOTIndexXX,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "IGPOTIndexXX" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
       deallocate(IGPOTZonaXX,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "IGPOTZonaXX" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
       deallocate(IZONAMASSIVGXX,stat=ierr)
       if(ierr/=0) then
          write(*,*) 'TOTALXF'
          write(*,*) 'THE FILE "IZONAMASSIVGXX" IS NOT REMOVED FROM MEMORY' 
          read(*,*)
   	      stop 
       endif
    ENDIF  

	! пюявер бгюхлндеиярбхъ щкейрпнмнб я йпхярюккхвеяйхл онкел 
    EnergyCrystallFieldX=0.D0
	DO IIX=1,NumbreMO
	   IntCrystallField(IIX)=EnergyCrystallField(IIX,Npoint,H,NumbreGarmMO,Ncrystall,RFunMO,RO1X)
       EnergyCrystallFieldX=EnergyCrystallFieldX+FLOAT(IQ(IIX))*IntCrystallField(IIX)
	ENDDO
	
	
	! пюявер щмепцхх лефъдепмнцн бгюхлндеиярбхъ
	EnergyNuclear=NuclearEnergy(Nnuclei,Z0,Z,COORRR,COOAAA)


    ! пюявер щмепцхх йскнмнбяйнцн бгюхлндеиярбхъ (нясыеярбкъеряъ б юрнлмшу едхмхжюу)
	! опълюъ вюярэ
	SumFCoulomb=0.D0
    DO IIY=1,NumbreMO 
       DO IIZ=1,NumbreMO
  	      DO IIX=1,NtypIndexCoffF(IIY,IIZ) 
		     IF(DABS(RcoffADI(IIX,IIY,IIZ)).GT.1.D-5) THEN
			    ! дкъ яксвюъ йскнмнбяйнцн бгюхлндеиярбхъ  бмсрпх нанкнвйх мер лмнфхрекъ 0.5 дкъ йнлохмяюжхх
				! слмнфюел мю дбю
				IF(IIY.EQ.IIZ) THEN
                    SumFCoulomb=SumFCoulomb+2.D0*RcoffADI(IIX,IIY,IIZ)*FCoulombInt(IIX,IIY,IIZ)
			       ELSE
                    SumFCoulomb=SumFCoulomb+RcoffADI(IIX,IIY,IIZ)*FCoulombInt(IIX,IIY,IIZ)
				ENDIF
			 ENDIF
		  ENDDO  
	   ENDDO
    ENDDO
    
	! пегскэрюр гюохяюм б юрнлмшу едемхжюу
    SumFCoulomb=SumFCoulomb*0.5D0

    
	! налеммюъ вюярэ (пюявер нясыеярбкъеряъ б юрнлмшу едхмхжюу)
    SumGCoulomb=0.D0
	! опнбепъел рхо спюбмемхъ
	IF(IndexHartreeFock.EQ.1) THEN 
	   DO IIY=1,NumbreMO 
          DO IIZ=1,NumbreMO
	         DO IIX=1,NtypIndexCoffG(IIY,IIZ) 
		        IF(DABS(RcoffBEI(IIX,IIY,IIZ)).GT.1.D-5) THEN
                   SumGCoulomb=SumGCoulomb+RcoffBEI(IIX,IIY,IIZ)*GCoulombInt(IIX,IIY,IIZ)
			    ENDIF
		     ENDDO  
	      ENDDO
       ENDDO
    
	   ! пегскэрюр гюохяюм б юрнлмшу едемхжюу
       SumGCoulomb=SumGCoulomb*0.5D0
    ENDIF 

	! пюявер щояхкнм
	EPSFull=0.D0
	DO IIY=1,NumbreMO 
       EPSFull=EPSFull+FLOAT(IQ(IIY))*RFunMO(IIY,1,1)
	ENDDO

	! пегскэрюр гюохяюм б юрнлмшу едемхжюу
    EPSFull=EPSFull*0.5D0

	! бшдювю пегскэрюрнб пюяверю
	! бшдювю йскнмнбяйху хмрецпюкнб
	! бяе бекхвхмш опхбндъряъ б юрнлмшу едхмхжюу 
	WRITE(6,*)
	WRITE(6,2200)
	WRITE(6,*)
     
    DO IIY=1,NumbreMO 
       DO IIZ=IIY,NumbreMO
          DO IIX=1,NtypIndexCoffF(IIY,IIZ) 
		     IF(DABS(RcoffADI(IIX,IIY,IIZ)).GT.1.D-5) THEN
                WRITE(6,2000) IIX,NN(IIY),LB(IABS(ML(IIY))+1),NN(IIZ),LB(IABS(ML(IIZ))+1),FCoulombInt(IIX,IIY,IIZ)
			 ENDIF
		  ENDDO  
	   ENDDO
    ENDDO  
	
	! опнбепъел рхо спюбмемхъ
	IF(IndexHartreeFock.EQ.1) THEN 
	   WRITE(6,*)
	   WRITE(6,2300)
	   DO IIY=1,NumbreMO 
          DO IIZ=1,NumbreMO
		     DO IIX=1,NtypIndexCoffG(IIY,IIZ) 
	            IF(DABS(RcoffBEI(IIX,IIY,IIZ)).GT.1.D-5) THEN
                   WRITE(6,2100) IIX,NN(IIY),LB(IABS(ML(IIY))+1),NN(IIZ),LB(IABS(ML(IIZ))+1),GCoulombInt(IIX,IIY,IIZ)
	  	        ENDIF
			 ENDDO
		  ENDDO  
	   ENDDO
	ENDIF
	

	WRITE(6,*)
	WRITE(6,2500) 
	WRITE(6,2600) 
    WRITE(6,2500) 
    WRITE(6,2700) EPSFull-2.D0*(SumFCoulomb+SumGCoulomb)
    WRITE(6,2800) SumFCoulomb+SumGCoulomb
	WRITE(6,2850) EnergyNuclear 
    WRITE(6,2550)
	WRITE(6,2900) EPSFull-(SumFCoulomb+SumGCoulomb)+EnergyNuclear
	WRITE(6,*)
    WRITE(6,3000) SumFCoulomb+SumGCoulomb
    WRITE(6,3100) EnergyCrystallFieldX
    WRITE(6,3200) EnergyNuclear 
    WRITE(6,2555) 
	WRITE(6,3300) SumFCoulomb+SumGCoulomb+EnergyCrystallFieldX+EnergyNuclear
    WRITE(6,3400) EPSFull-2.D0*(SumFCoulomb+SumGCoulomb)-EnergyCrystallFieldX
	WRITE(6,3500) (SumFCoulomb+SumGCoulomb+EnergyCrystallFieldX+EnergyNuclear)/(EPSFull-2.D0*(SumFCoulomb+SumGCoulomb)-EnergyCrystallFieldX)

	  
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(FCoulombInt,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "FCoulombInt" IS NOT REMOVED FROM MEMORY' 
       read(*,*)
   	   stop 
    endif
    deallocate(GCoulombInt,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "GCoulombInt" IS NOT REMOVED FROM MEMORY' 
       read(*,*)
       stop 
    endif
    deallocate(IntCrystallField,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'TOTALXF'
       write(*,*) 'THE FILE "IntCrystallField" IS NOT REMOVED FROM MEMORY' 
	   read(*,*)
	   stop 
    endif
    return
  end subroutine TOTALXF



    !! SUBPROGRAMME OF CALCULATION OF ENERGY INTERACTION OF ELECTRON WITH CRYSTALLINE FIELD
═══ !! NMO-NUMBER OF MOLECULAR ORBITALS
═══ !! Npoint-NUMBER OF POINTS
═══ !! H-step
═══ !! NumbreGarmMO (NumbreMO) -MASSIVE NUMBER OF HARMONIC MOLECULAR ORBITALS
═══ !! Ncrystall (NumbreMO, NumbreGarmMO (NumbreMO), NumbreGarmMO (NumbreMO), Npoint) -MASIVE OF ELECTRON INTERACTION OPERATOR WITH NUCLEI OF MOLECULES
═══ !! RFunMO (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values of the radial parts of the harmonics of the molecular orbitals of the configuration
═══ !! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
   real(8) function EnergyCrystallField(NMO,Npoint,H,NumbreGarmMO,Ncrystall,RFunMO,RO1X)
     implicit none 
     integer::NMO,Npoint
	 real(8)::H
	 integer,dimension(:)::NumbreGarmMO
	 real(8),dimension(:)::RO1X
     real(8),dimension(:,:,:)::RFunMO
     real(8),dimension(:,:,:,:)::Ncrystall
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IXX,IYY 
	 real(8)::SUMY

     EnergyCrystallField=0.D0
     
	 SUMY=0.D0

	 DO IXX=1,NumbreGarmMO(NMO)
        DO IYY=1,NumbreGarmMO(NMO)
           SUMY=SUMY+CrystallFieldGarmon(NMO,IXX,IYY,Npoint,H,Ncrystall,RFunMO,RO1X)
		ENDDO
     ENDDO

     EnergyCrystallField=SUMY

	 return
   end function EnergyCrystallField



   ! ондопнцпюллю пюяверю люрпхвмнцн щкелемрю йпхярюккхвеяйнцн онкъ
   ! NMO-мнлеп лнкейскъпмни нпахрюкх
   ! Npoint-вхякн рнвей
   ! H-ЬЮЦ
   ! Ncrystall(NumbreMO,NumbreGarmMO(NumbreMO),NumbreGarmMO(NumbreMO),Npoint)-люяяхб ноепюрнпю бгюхлндеиярбхъ щкейрпнмю я ъдпюлх лнкейскш
   ! RFunMO(NumbreMO,NumbreGarmMO(NumbreMO),Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
   ! RO1X(Npoint)- люяяхб гмювемхи оепбни опнхгбндмни он мнбни оепелеммни
   real(8) function CrystallFieldGarmon(NMO,IGAR1,IGAR2,Npoint,H,Ncrystall,RFunMO,RO1X)
     implicit none
	 integer::NMO,IGAR1,IGAR2,Npoint
	 real(8)::H
	 real(8),dimension(:)::RO1X
     real(8),dimension(:,:,:)::RFunMO
     real(8),dimension(:,:,:,:)::Ncrystall
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IXD
     real(8)::SUMM,SING2,SING1,R0,RN1,RN 

     !SUMM=0.D0
	 !DO IXD=1,Npoint
     !   SUMM=SUMM+RFunMO(NMO,IGAR1,IXD+2)*RFunMO(NMO,IGAR2,IXD+2)*Ncrystall(NMO,IGAR1,IGAR2,IXD)/RO1X(IXD)**2
     !ENDDO

     
	 SUMM=0.D0
	 ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
     IF(Npoint/2..EQ.Npoint/2) THEN
	     ! пюявер он тнплске яхлоянмю
         ! вермне вхякн рнвей
         SING2=0.D0
	     DO IXD=2,Npoint-2,2
            SING2=SING2+RFunMO(NMO,IGAR1,IXD+2)*RFunMO(NMO,IGAR2,IXD+2)*Ncrystall(NMO,IGAR1,IGAR2,IXD)/RO1X(IXD)**2
	     ENDDO
         SING1=0.D0
	     DO IXD=3,Npoint-3,2
            SING1=SING1+RFunMO(NMO,IGAR1,IXD+2)*RFunMO(NMO,IGAR2,IXD+2)*Ncrystall(NMO,IGAR1,IGAR2,IXD)/RO1X(IXD)**2
	     ENDDO
	     ! пегскэрхпсчыюъ ясллю 
	     R0=RFunMO(NMO,IGAR1,3)*RFunMO(NMO,IGAR2,3)*Ncrystall(NMO,IGAR1,IGAR2,1)/RO1X(1)**2
	     RN1=RFunMO(NMO,IGAR1,Npoint+1)*RFunMO(NMO,IGAR2,Npoint+1)*Ncrystall(NMO,IGAR1,IGAR2,Npoint-1)/RO1X(Npoint-1)**2   
	     RN=RFunMO(NMO,IGAR1,Npoint+2)*RFunMO(NMO,IGAR2,Npoint+2)*Ncrystall(NMO,IGAR1,IGAR2,Npoint)/RO1X(Npoint)**2   
         SUMM=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
       ELSE
         ! пювер он тнплске яхлоянмю
		 ! мевермне вхякн рнвей
         SING2=0.D0
		 DO IXD=2,Npoint-1,2
            SING2=SING2+RFunMO(NMO,IGAR1,IXD+2)*RFunMO(NMO,IGAR2,IXD+2)*Ncrystall(NMO,IGAR1,IGAR2,IXD)/RO1X(IXD)**2
	     ENDDO
         SING1=0.D0
		 DO IXD=3,Npoint-2,2
            SING1=SING1+RFunMO(NMO,IGAR1,IXD+2)*RFunMO(NMO,IGAR2,IXD+2)*Ncrystall(NMO,IGAR1,IGAR2,IXD)/RO1X(IXD)**2
	     ENDDO
         ! пегскэрхпсчыюъ ясллю 
		 R0=RFunMO(NMO,IGAR1,3)*RFunMO(NMO,IGAR2,3)*Ncrystall(NMO,IGAR1,IGAR2,1)/RO1X(1)**2
	     RN=RFunMO(NMO,IGAR1,Npoint+2)*RFunMO(NMO,IGAR2,Npoint+2)*Ncrystall(NMO,IGAR1,IGAR2,Npoint)/RO1X(Npoint)**2   
         SUMM=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
     ENDIF


     
	 CrystallFieldGarmon=SUMM*H

     return
   end function CrystallFieldGarmon

   ! ондопнцпюллю пюяверю щмепцхх лефъдепмнцн бгюхлндеиярбхъ 
   ! Nnuclei-вхякн ъдеп(ме свхрюбюел ъдпн б мювюке йннпдхмюр)
   ! Z0-гюпъд ъдпю б мювюке йннпдхмюр
   ! Z(Nnuclei)-люяяхб гюпъднб ъдеп
   ! COORRR(N)-люяяхб пюдхюкэмшу йнндхмюр 
   ! COORRR(N)-пюдхюкэмюъ йннпдхмюрю ъдпю (R-пюдхся)
   ! COOAAA(N,2)-люяяхб йннпдхмюр ъдеп
   ! COOAAA(N,1)-брнпюъ йннпдхмюрю -RTeta-сцнк (0=<RTeta<=PI)
   ! COOAAA(N,2)-рперэъ йннпдхмюрю -RFu-сцнк (0=<RFu<=2*PI) 
   real(8) function NuclearEnergy(Nnuclei,Z0,Z,COORRR,COOAAA)
     implicit none  
     integer::Nnuclei
	 real(8)::Z0
	 real(8),dimension(:)::Z,COORRR
     real(8),dimension(:,:)::COOAAA
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer::IIX,IIY
	 real(8)::SUME  
     
	 ! нясыеярбкъел пюявер щмепцхх бгюхлндеиярбхъ я ъдпнл б мювюке йннпдхмюр
     IF(Z0.NE.0.D0) THEN
	     SUME=0.D0
	     DO IIX=1,Nnuclei
            SUME=SUME+Z0*Z(IIX)/DistanceNuclear(0.D0,0.D0,0.D0,COORRR(IIX),COOAAA(IIX,1),COOAAA(IIX,2))
	     ENDDO
     ENDIF

     ! пюявер бгюхлндеиярбхъ ледкс ъдпюлх 
	 
	 DO IIX=1,Nnuclei
        DO IIY=1,Nnuclei
           IF(IIX.NE.IIY) THEN   
               SUME=SUME+0.5D0*Z(IIX)*Z(IIY)/DistanceNuclear(COORRR(IIX),COOAAA(IIX,1),COOAAA(IIX,2),COORRR(IIY),COOAAA(IIY,1),COOAAA(IIY,2))
		   ENDIF
		ENDDO
     ENDDO
     
     NuclearEnergy=SUME

	 return
   end function NuclearEnergy

   ! ондопнцпюллю пюяверю пюяярнъмхъ лефдс ъдпюлх
   ! R1,ATeta1,AFi1-яТЕПХВЕЯЙХЕ ЙННПДХМЮРШ ОЕПБНЦН ЪДПЮ
   ! R2,ATeta2,AFi2-яТЕПХВЕЯЙХЕ ЙННПДХМЮРШ БРНПНЦН ЪДПЮ
   real(8) function DistanceNuclear(R1,ATeta1,AFi1,R2,ATeta2,AFi2)
     implicit none  
     real(8)::R1,ATeta1,AFi1,R2,ATeta2,AFi2
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(8)::x1,y1,z1,x2,y2,z2  
     
	 ! ОЕПЕБНДХЛ ЯТЕПХВЕЯЙХЕ ЙННПДХМЮРШ Б ДЕЙЮПРНБШ
	 x1=R1*DSIN(ATeta1)*DCOS(AFi1)
     y1=R1*DSIN(ATeta1)*DSIN(AFi1)
	 z1=R1*DCOS(ATeta1)
     x2=R2*DSIN(ATeta2)*DCOS(AFi2)
     y2=R2*DSIN(ATeta2)*DSIN(AFi2)
	 z2=R2*DCOS(ATeta2)
     DistanceNuclear=DSQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	 return
   end function DistanceNuclear

   
  ! ондопнцпюллю пюяверю йскнмнбяйнцн хмрецпюкю лефдс лнкейскъпмшлх нпахрюкълх
  ! Itip-мнлеп рхою йскнмнбяйнцн бгюлндеиярбхъ
  ! IS1,IS2,IS3,IS4-мнлепю нанкнвей свюбярбсчыху б пюявере 
  ! MassivML(Itip,NumbreMO)-люяяхб опнейжхи лнлемрю (лнкейскъпмшу нпахрюкеи) хяонкэгселши опх пюявере йскнмняйнцн бгюхлндеиярбхъ
  ! LgarmonicMO(2,NumbreMO)-люяяхб гмювемхи нпахрюкэмшу лнлемрнб лнкейскъпмшу нпахрюкеи
  ! LgarmonicMO(1,NumbreMO)-лхмхлюкэмне гмювемхе
  ! LgarmonicMO(2,NumbreMO)-люйяхлюкэмне гмювемхе
  ! Npoint-ВХЯКН РНВЕЙ
  ! H-ЬЮЦ
  ! RFunMO(NumbreMO,NumbreGarmMO(NumbreMO),Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
  ! R(Npoint)-люяяхб гмювемхи пюдхсяю
  ! RO1X(Npoint)- люяяхб гмювемхи оепбни опнхгбндмни он мнбни оепелеммни
  ! IndexML(Itip,NumbreMO)-люяяхб рхонб опнейжхи лнлемрю (лнкейскъпмшу нпахрюкеи) хяонкэгселши опх пюявере йскнмняйнцн бгюхлндеиярбхъ
  ! CkkGarmonik(KK,L2,L1,ML2,ML1)-люяяяхб гмювемхи люрпхвмшу щкелемрнб ятепхвеяйху цюплнмхй
  real(8) function CoulombMO(Itip,IS1,IS2,IS3,IS4,MassivML,LgarmonicMO,Npoint,H,RFunMO,R,RO1X,IndexML,CkkGarmonik) 
    implicit none 
    integer::Itip,IS1,IS2,IS3,IS4,Npoint
	real(8)::H
	integer,dimension(:,:)::MassivML,IndexML,LgarmonicMO
    real(8),dimension(:)::R,RO1X
    real(8),dimension(:,:,:)::RFunMO
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::IL1,IL2,IL3,IL4,MINN,MAXX,KK
    integer::ML1,ML2,ML3,ML4,IML1,IML2,IML3,IML4,INX1,INX2,INX3,INX4,ierr
	real(8)::PID,SUMK
    real(8),allocatable,dimension(:)::R6

    
    
	CoulombMO=0.D0
	
	ML1=MassivML(Itip,IS1)
    ML2=MassivML(Itip,IS2)
    ML3=MassivML(Itip,IS3)
    ML4=MassivML(Itip,IS4)
    
	IML1=IndexML(Itip,IS1)
    IML2=IndexML(Itip,IS2)
    IML3=IndexML(Itip,IS3)
    IML4=IndexML(Itip,IS4)



	IF(ML1+ML2.NE.ML3+ML4) THEN 
	   return 
    ENDIF  

    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R6(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CoulombMO'
      write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	  read(*,*)
      stop 
    endif 



    
	SUMK=0.D0
	INX1=0
	DO IL1=LgarmonicMO(1,IS1),LgarmonicMO(2,IS1)
       INX1=INX1+1
	   INX2=0
	   DO IL2=LgarmonicMO(1,IS2),LgarmonicMO(2,IS2)    
          INX2=INX2+1
		  INX3=0
		  DO IL3=LgarmonicMO(1,IS3),LgarmonicMO(2,IS3)    
             INX3=INX3+1
			 INX4=0
			 DO IL4=LgarmonicMO(1,IS4),LgarmonicMO(2,IS4)  
		        INX4=INX4+1

				MINN=MAX0(IABS(IL1-IL3),IABS(IL2-IL4))
                MAXX=MIN0((IL1+IL3),(IL2+IL4)) 
                DO KK=MINN,MAXX,2 
				   PID=CkkGarmonik(KK+1,IL3+1,IL1+1,IML3,IML1)*CkkGarmonik(KK+1,IL2+1,IL4+1,IML2,IML4)
                   
                   IF(DABS(PID).LT.1.D-6) THEN 
				      CYCLE 
                   ENDIF
                   SUMK=SUMK+RRCoulomb(KK,IS1,INX1,IS2,INX2,IS3,INX3,IS4,INX4,R6,H,Npoint,RFunMO,R,RO1X)*PID
				ENDDO 
			 ENDDO 
		  ENDDO
       ENDDO  
    ENDDO

    CoulombMO=SUMK

    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R6,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CoulombMO'
       write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
    
	return   
  end function CoulombMO



 
  ! ондопнцпюллю пюяверю опълнцн йскнмнбяйнцн хмрецпюкю лефдс лнкейскъпмшлх нпахрюкълх
  ! Itip-мнлеп рхою йскнмнбяйнцн бгюлндеиярбхъ
  ! IS1,IS2,IS3,IS4-мнлепю нанкнвей свюбярбсчыху б пюявере 
  ! NomeroFPOT(NumbreMO)-люяяхб сйюгшбючыхи мнлепю тюикнб
  ! LgarmonicMO(2,NumbreMO)-люяяхб гмювемхи нпахрюкэмшу лнлемрнб лнкейскъпмшу нпахрюкеи
  ! LgarmonicMO(1,NumbreMO)-лхмхлюкэмне гмювемхе
  ! LgarmonicMO(2,NumbreMO)-люйяхлюкэмне гмювемхе
  ! Npoint-ВХЯКН РНВЕЙ
  ! H-ЬЮЦ
  ! RFunMO(NumbreMO,NumbreGarmMO(NumbreMO),Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
  ! R(Npoint)-люяяхб гмювемхи пюдхсяю
  ! RO1X(Npoint)- люяяхб гмювемхи оепбни опнхгбндмни он мнбни оепелеммни
  ! IndexML(Itip,NumbreMO)-люяяхб рхонб опнейжхи лнлемрю (лнкейскъпмшу нпахрюкеи) хяонкэгселши опх пюявере йскнмняйнцн бгюхлндеиярбхъ
  ! CkkGarmonik(KK,L2,L1,ML2,ML1)-люяяхб гмювемхи люрпхвмшу щкелемрнб ятепхвеяйху цюплнмхй
  real(8) function FCoulombMO(Itip,IS1,IS2,LgarmonicMO,Npoint,H,RFunMO,R,RO1X,IndexML,CkkGarmonik,NomeroFPOT,IFPOTIndex,IFPOTZona,NumbreSxemaCalculationF,RFILEPoltensF8,RFILEPoltensF4) 
    implicit none 
    integer::Itip,IS1,IS2,IS3,IS4,Npoint,NumbreSxemaCalculationF
	real(8)::H
	integer,dimension(:)::NomeroFPOT
	integer,dimension(:,:)::LgarmonicMO
	integer(1),dimension(:,:,:,:)::IFPOTIndex
	integer,dimension(:,:,:,:)::IndexML,IFPOTZona 
    real(8),dimension(:)::R,RO1X
    real(4),dimension(:,:)::RFILEPoltensF4
	real(8),dimension(:,:)::RFILEPoltensF8
    real(8),dimension(:,:,:)::RFunMO
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::IL1,IL2,IL3,IL4,MINN,MAXX,KK,IOIO,IDFSZ,IOSD,IDFZ
    integer::ML1,ML2,ML3,ML4,IML1,IML2,IML3,IML4,INX1,INX2,INX3,INX4,ierr
	integer::IndexML1,IndexML2,IndexML3,IndexML4
	real(8)::PID,SUMK,RcoffSimmetry,RcoffSimmetryC,SING,SING2,SING1,R0,RN1,RN
    real(8),allocatable,dimension(:)::RPOT,RPOTRez

    
    
	FCoulombMO=0.D0
	    
	
    ! хмдейяш опнейжхи   
    IndexML1=IndexML(Itip,IS1,IS2,1)
    IndexML2=IndexML(Itip,IS1,IS2,2)
	IndexML3=IndexML(Itip,IS1,IS2,3)
	IndexML4=IndexML(Itip,IS1,IS2,4)




    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(RPOT(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CoulombMO'
      write(*,*) 'MEMORY ON THE FILE "RPOT" IS NOT SELECTED'
	  read(*,*)
      stop 
    endif 
	allocate(RPOTRez(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CoulombMO'
      write(*,*) 'MEMORY ON THE FILE "RPOTRez" IS NOT SELECTED'
	  read(*,*)
      stop 
    endif

    ! мнлеп тюикю
	IOSD=NomeroFPOT(1)

    ! онремжхюк пюявхршбюеряъ мю IS1 х IS1 лнкейскъпмшу нпахрюкеи
	SUMK=0.D0
	INX2=0
	DO IL2=LgarmonicMO(1,IS2),LgarmonicMO(2,IS2)    
       INX2=INX2+1
	   INX4=INX2-1
	   DO IL4=IL2,LgarmonicMO(2,IS2)  
	      INX4=INX4+1
	      ! свхршбюел яхллерпхч нрмняхрекэмн оепеярюмнбяйх хмдейянб 
          IF(INX2.EQ.INX4) THEN
              RcoffSimmetryC=1.D0
		     ELSE
              RcoffSimmetryC=2.D0 
          ENDIF  

		  
		  ! пюявхршбюел онремжхюк 
	      RPOTRez=0.D0
		  INX1=0
		  DO IL1=LgarmonicMO(1,IS1),LgarmonicMO(2,IS1)
             INX1=INX1+1
	         INX3=INX1-1
		     DO IL3=IL1,LgarmonicMO(2,IS1)    
                INX3=INX3+1
                ! свхршбюел яхллерпхч нрмняхрекэмн оепеярюмнбяйх хмдейянб 
                IF(INX1.EQ.INX3) THEN
                     RcoffSimmetry=1.D0
		           ELSE
                     RcoffSimmetry=2.D0 
                ENDIF  
				MINN=MAX0(IABS(IL1-IL3),IABS(IL2-IL4))
                MAXX=MIN0((IL1+IL3),(IL2+IL4)) 
                DO KK=MINN,MAXX,2 
				   
				   PID=CkkGarmonik(KK+1,IL1+1,IL3+1,IndexML1,IndexML3)*CkkGarmonik(KK+1,IL4+1,IL2+1,IndexML4,IndexML2)
				   IF(DABS(PID).LT.1.D-6) THEN 
				      CYCLE 
                   ENDIF
				   ! пюявер онремжхюкю
				   ! нясыеярбкъел явхршбюмхе онремжхюкю
				   IF(IFPOTIndex(IS1,IL1+1,IL3+1,KK+1).EQ.0) THEN
                      ! дюммши онремжхюк ме гюохяюм
				      WRITE(6,*) 'FPOT NOT WRITING  NMO=',IS1,' L1=',IL1,' L2=',IL3,' K=',KK
					  WRITE(6,*) 'DATA',KK,IL1,IL3,IndexML1,IndexML3
                      WRITE(6,*) 'DATA',KK,IL4,IL2,IndexML4,IndexML2
					  WRITE(6,*) 'CKK',CkkGarmonik(KK+1,IL1+1,IL3+1,IndexML1,IndexML3),CkkGarmonik(KK+1,IL4+1,IL2+1,IndexML4,IndexML2)
					  STOP
				   ENDIF
				   ! гнмю гюохях
                   IDFSZ=IFPOTZona(IS1,IL1+1,IL3+1,KK+1)
				   ! явхршбюмхе хг люяяхб
				   IF(NumbreSxemaCalculationF.EQ.1) THEN
                      RPOT(1:Npoint)=RFILEPoltensF8(1:Npoint,IDFSZ)
				   ENDIF
				   ! явхршбюмхе хг люяяхб
				   IF(NumbreSxemaCalculationF.EQ.2) THEN
                      RPOT(1:Npoint)=DBLE(RFILEPoltensF4(1:Npoint,IDFSZ))
				   ENDIF
                   ! гюохяэ б тюик
				   IF(NumbreSxemaCalculationF.EQ.3) THEN
				      READ (IOSD,REC=IDFSZ) (RPOT(IOIO),IOIO=1,Npoint)  
			       ENDIF
						 
				   ! мюундхл ясллс
				   RPOTRez=RPOTRez+RcoffSimmetry*PID*RPOT
				ENDDO 
			 ENDDO 
		  ENDDO

          ! нясыеярбкъел пюявер йскнмнбяйнцн хмрецпюкю
		  SING=0.D0
		  ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
          IF(Npoint/2..EQ.Npoint/2) THEN
	          ! пюявер он тнплске яхлоянмю
              ! вермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,Npoint-2,2
                 SING2=SING2+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS2,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,Npoint-3,2
                 SING1=SING1+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS2,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
		      ! пегскэрхпсчыюъ ясллю 
		      R0=RPOTRez(1)*RFunMO(IS2,INX2,3)*RFunMO(IS2,INX4,3)/(R(1)*RO1X(1)**2)
		      RN1=RPOTRez(Npoint-1)*RFunMO(IS2,INX2,Npoint+1)*RFunMO(IS2,INX4,Npoint+1)/(R(Npoint-1)*RO1X(Npoint-1)**2)
		      RN=RPOTRez(Npoint)*RFunMO(IS2,INX2,Npoint+2)*RFunMO(IS2,INX4,Npoint+2)/(R(Npoint)*RO1X(Npoint)**2)
              SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
            ELSE
              ! пювер он тнплске яхлоянмю
			  ! мевермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,Npoint-1,2
                 SING2=SING2+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS2,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,Npoint-2,2
                 SING1=SING1+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS2,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              ! пегскэрхпсчыюъ ясллю 
		      R0=RPOTRez(1)*RFunMO(IS2,INX2,3)*RFunMO(IS2,INX4,3)/(R(1)*RO1X(1)**2)
		      RN=RPOTRez(Npoint)*RFunMO(IS2,INX2,Npoint+2)*RFunMO(IS2,INX4,Npoint+2)/(R(Npoint)*RO1X(Npoint)**2)
              SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
          ENDIF


          
		  ! пюяявхршбюел йскнмнбяйхи хмрецпюк я хяонкэгнбюмхел онксвеммнцн онремжхюкю
          SUMK=SUMK+RcoffSimmetryC*SING*H
       
	   ENDDO  
    ENDDO

    FCoulombMO=SUMK

    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(RPOT,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CoulombMO'
       write(*,*) 'THE FILE "RPOT" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
	deallocate(RPOTRez,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CoulombMO'
       write(*,*) 'THE FILE "RPOTRez" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
    
	return   
  end function FCoulombMO



  ! ондопнцпюллю пюяверю йскнмнбяйнцн хмрецпюкю лефдс лнкейскъпмшлх нпахрюкълх
  ! Itip-мнлеп рхою йскнмнбяйнцн бгюлндеиярбхъ
  ! IS1,IS2,IS3,IS4-мнлепю нанкнвей свюбярбсчыху б пюявере 
  ! MassivML(Itip,NumbreMO)-люяяхб опнейжхи лнлемрю (лнкейскъпмшу нпахрюкеи) хяонкэгселши опх пюявере йскнмняйнцн бгюхлндеиярбхъ
  ! LgarmonicMO(2,NumbreMO)-люяяхб гмювемхи нпахрюкэмшу лнлемрнб лнкейскъпмшу нпахрюкеи
  ! LgarmonicMO(1,NumbreMO)-лхмхлюкэмне гмювемхе
  ! LgarmonicMO(2,NumbreMO)-люйяхлюкэмне гмювемхе
  ! Npoint-ВХЯКН РНВЕЙ
  ! H-ЬЮЦ
  ! RFunMO(NumbreMO,NumbreGarmMO(NumbreMO),Npoint+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ  
  ! R(Npoint)-люяяхб гмювемхи пюдхсяю
  ! RO1X(Npoint)- люяяхб гмювемхи оепбни опнхгбндмни он мнбни оепелеммни
  ! IndexML(Itip,NumbreMO)-люяяхб рхонб опнейжхи лнлемрю (лнкейскъпмшу нпахрюкеи) хяонкэгселши опх пюявере йскнмняйнцн бгюхлндеиярбхъ
  ! CkkGarmonik(KK,L2,L1,ML2,ML1)-люяяяхб гмювемхи люрпхвмшу щкелемрнб ятепхвеяйху цюплнмхй
  real(8) function GCoulombMO(Itip,IS1,IS2,LgarmonicMO,Npoint,H,RFunMO,R,RO1X,IndexML,CkkGarmonik,NumbreGPOTMO,NomeroFPOT,IGPOTIndex,IGPOTZona,NumbreSxemaCalculationG,RFILEPoltensG8,RFILEPoltensG4) 
    implicit none 
    integer::Itip,IS1,IS2,Npoint,NumbreSxemaCalculationG
	real(8)::H
	integer,dimension(:)::NomeroFPOT 
	integer,dimension(:,:)::LgarmonicMO,NumbreGPOTMO
	integer(1),dimension(:,:,:,:)::IGPOTIndex
	integer,dimension(:,:,:,:)::IndexML,IGPOTZona
    real(8),dimension(:)::R,RO1X
    real(4),dimension(:,:)::RFILEPoltensG4
    real(8),dimension(:,:)::RFILEPoltensG8
    real(8),dimension(:,:,:)::RFunMO
	real(8),dimension(:,:,:,:,:)::CkkGarmonik
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::IL1,IL2,IL3,IL4,MINN,MAXX,KK,ITYGR,IDSFZX,IZOGR,IOIO,IDFZ
    integer::IndexML1,IndexML2,IndexML3,IndexML4,ML1,ML2,ML3,ML4,IML1,IML2,IML3,IML4,INX1,INX2,INX3,INX4,ierr
	real(8)::PID,SUMK,SING,SING2,SING1,R0,RN1,RN
    real(8),allocatable,dimension(:)::RPOT,RPOTRez

    
    
	GCoulombMO=0.D0

	! мнлеп рхою
    ITYGR=NumbreGPOTMO(IS2,IS1)
	! мнлеп тюикю
	IDSFZX=NomeroFPOT(1) 
    	
	! хмдейяш опнейжхи   
    IndexML1=IndexML(Itip,IS1,IS2,1)
    IndexML2=IndexML(Itip,IS1,IS2,2)
	IndexML3=IndexML(Itip,IS1,IS2,3)
	IndexML4=IndexML(Itip,IS1,IS2,4)

 

    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(RPOT(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'GCoulombMO'
      write(*,*) 'MEMORY ON THE FILE "RPOT" IS NOT SELECTED'
	  read(*,*)
      stop 
    endif 
	allocate(RPOTRez(Npoint),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'GCoulombMO'
      write(*,*) 'MEMORY ON THE FILE "RPOTRez" IS NOT SELECTED'
	  read(*,*)
      stop 
    endif



    
	SUMK=0.D0
    INX2=0
    DO IL2=LgarmonicMO(1,IS2),LgarmonicMO(2,IS2)    
       INX2=INX2+1
  	   INX4=0
	   DO IL4=LgarmonicMO(1,IS1),LgarmonicMO(2,IS1)  
	      INX4=INX4+1
          
		  ! пюявхршбюел онремжхюк 
	      RPOTRez=0.D0
		  INX1=0
	      DO IL1=LgarmonicMO(1,IS1),LgarmonicMO(2,IS1)
             INX1=INX1+1
             INX3=0
	         DO IL3=LgarmonicMO(1,IS2),LgarmonicMO(2,IS2)    
                INX3=INX3+1

				MINN=MAX0(IABS(IL1-IL3),IABS(IL2-IL4))
                MAXX=MIN0((IL1+IL3),(IL2+IL4)) 
                DO KK=MINN,MAXX,2 
				   PID=CkkGarmonik(KK+1,IL1+1,IL3+1,IndexML1,IndexML3)*CkkGarmonik(KK+1,IL4+1,IL2+1,IndexML4,IndexML2)
				       
                   IF(DABS(PID).LT.1.D-6) THEN 
				      CYCLE 
                   ENDIF
                   
				   ! пюявер онремжхюкю
				   IF(IGPOTIndex(ITYGR,IL1+1,IL3+1,KK+1).EQ.0) THEN
                      ! дюммши онремжхюк ме гюохяюм
					  WRITE(6,*) 'GPOT NOT WRITING  NMO=',IS2,IS1,' L1=',IL3,' L2=',IL1,' K=',KK
				      STOP
				   ENDIF
                   ! гнмю гюохях
				   IZOGR=IGPOTZona(ITYGR,IL1+1,IL3+1,KK+1)
                   
				   ! нясыеярбкъел явхршбюмхе налеммнцн онремжхюкю
				   ! явхршбюмхе хг люяяхбю 
				   IF(NumbreSxemaCalculationG.EQ.1) THEN
                      RPOT(1:Npoint)=RFILEPoltensG8(1:Npoint,IZOGR)
				   ENDIF 
                   IF(NumbreSxemaCalculationG.EQ.2) THEN
                      RPOT(1:Npoint)=DBLE(RFILEPoltensG4(1:Npoint,IZOGR))
				   ENDIF 
                   ! явхршбюмхе хг тюик
				   IF(NumbreSxemaCalculationG.EQ.3) THEN
				      READ (IDSFZX,REC=IZOGR) (RPOT(IOIO),IOIO=1,Npoint)  
			       ENDIF
				   
                    ! мюундхл ясллс
				   RPOTRez=RPOTRez+PID*RPOT
 
				ENDDO 
			 ENDDO 
		  ENDDO
          
		  ! нясыеярбкъел пюявер йскнмнбяйнцн хмрецпюкю
		  SING=0.D0
		  ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
          IF(Npoint/2..EQ.Npoint/2) THEN
	          ! пюявер он тнплске яхлоянмю
              ! вермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,Npoint-2,2
                 SING2=SING2+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS1,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,Npoint-3,2
                 SING1=SING1+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS1,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
		      ! пегскэрхпсчыюъ ясллю 
		      R0=RPOTRez(1)*RFunMO(IS2,INX2,3)*RFunMO(IS1,INX4,3)/(R(1)*RO1X(1)**2)
		      RN1=RPOTRez(Npoint-1)*RFunMO(IS2,INX2,Npoint+1)*RFunMO(IS1,INX4,Npoint+1)/(R(Npoint-1)*RO1X(Npoint-1)**2)
		      RN=RPOTRez(Npoint)*RFunMO(IS2,INX2,Npoint+2)*RFunMO(IS1,INX4,Npoint+2)/(R(Npoint)*RO1X(Npoint)**2)
              SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
            ELSE
              ! пювер он тнплске яхлоянмю
			  ! мевермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,Npoint-1,2
                 SING2=SING2+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS1,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,Npoint-2,2
                 SING1=SING1+RPOTRez(IDFZ)*RFunMO(IS2,INX2,2+IDFZ)*RFunMO(IS1,INX4,2+IDFZ)/(R(IDFZ)*RO1X(IDFZ)**2)
	          ENDDO
              ! пегскэрхпсчыюъ ясллю 
		      R0=RPOTRez(1)*RFunMO(IS2,INX2,3)*RFunMO(IS1,INX4,3)/(R(1)*RO1X(1)**2)
		      RN=RPOTRez(Npoint)*RFunMO(IS2,INX2,Npoint+2)*RFunMO(IS1,INX4,Npoint+2)/(R(Npoint)*RO1X(Npoint)**2)
              SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
          ENDIF

        
          ! пюяявхршбюел йскнмнбяйхи хмрецпюк я хяонкэгнбюмхел онксвеммнцн онремжхюкю
          SUMK=SUMK+SING*H
	   ENDDO  
    ENDDO

    GCoulombMO=SUMK
              
		  


    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(RPOT,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'GCoulombMO'
       write(*,*) 'THE FILE "RPOT" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
	deallocate(RPOTRez,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'GCoulombMO'
       write(*,*) 'THE FILE "RPOTRez" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
    
	return   
  end function GCoulombMO





  ! пюявер пюдхюкэмни вюярх йскнмняйнцн хмрецпюкю лефдс цюплнмхйюлх
  real(8) function RRCoulomb(K,NMO1,NGar1,NMO2,NGar2,NMO3,NGar3,NMO4,NGar4,R6,H,N,Rfun,R1,RO1X)
    implicit none 
	integer::K,NMO1,NGar1,NMO2,NGar2,NMO3,NGar3,NMO4,NGar4,N
	real(8)::H
	real(8),dimension(:)::R6,R1,RO1X
    real(8),dimension(:,:,:)::Rfun
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::I,IDFZ
	real(8)::SUMM,AK,SING,SING2,SING1,R0,RN1,RN
		 
    AK=DBLE(K)
    call POTENSS(AK,NMO1,NGar1,NMO3,NGar3,Rfun,N,H,R6,R1,RO1X)
    !SUMM=0.D0
    !DO I=1,N
    !   SUMM=SUMM+R6(I)*Rfun(NMO2,NGar2,I+2)*Rfun(NMO4,NGar4,I+2)/(R1(I)*RO1X(I)**2)   
    !ENDDO
	!RRCoulomb=sum(R6(1:N)*Rfun(NMO2,NGar2,2+1:N)*Rfun(NMO4,NGar4,2+1:N)/(R1(1:N)*RO1X(1:N)**2))*H
    
	SING=0.D0 
	! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
    IF(N/2..EQ.N/2) THEN
	     ! пюявер он тнплске яхлоянмю
         ! вермне вхякн рнвей
         SING2=0.D0
	     DO IDFZ=2,N-2,2
            SING2=SING2+R6(IDFZ)*Rfun(NMO2,NGar2,2+IDFZ)*Rfun(NMO4,NGar4,2+IDFZ)/(R1(IDFZ)*RO1X(IDFZ)**2)
	     ENDDO
         SING1=0.D0
	     DO IDFZ=3,N-3,2
            SING1=SING1+R6(IDFZ)*Rfun(NMO2,NGar2,2+IDFZ)*Rfun(NMO4,NGar4,2+IDFZ)/(R1(IDFZ)*RO1X(IDFZ)**2)
	     ENDDO
	     ! пегскэрхпсчыюъ ясллю 
	     R0=R6(1)*Rfun(NMO2,NGar2,3)*Rfun(NMO4,NGar4,3)/(R1(1)*RO1X(1)**2)
	     RN1=R6(N-1)*Rfun(NMO2,NGar2,N+1)*Rfun(NMO4,NGar4,N+1)/(R1(N-1)*RO1X(N-1)**2)
	     RN=R6(N)*Rfun(NMO2,NGar2,N+2)*Rfun(NMO4,NGar4,N+2)/(R1(N)*RO1X(N)**2)
         SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
       ELSE
         ! пювер он тнплске яхлоянмю
	     ! мевермне вхякн рнвей
         SING2=0.D0
	     DO IDFZ=2,N-1,2
            SING2=SING2+R6(IDFZ)*Rfun(NMO2,NGar2,2+IDFZ)*Rfun(NMO4,NGar4,2+IDFZ)/(R1(IDFZ)*RO1X(IDFZ)**2)
	     ENDDO
         SING1=0.D0
	     DO IDFZ=3,N-2,2
            SING1=SING1+R6(IDFZ)*Rfun(NMO2,NGar2,2+IDFZ)*Rfun(NMO4,NGar4,2+IDFZ)/(R1(IDFZ)*RO1X(IDFZ)**2)
	     ENDDO
         ! пегскэрхпсчыюъ ясллю 
	     R0=R6(1)*Rfun(NMO2,NGar2,3)*Rfun(NMO4,NGar4,3)/(R1(1)*RO1X(1)**2)
	     RN=R6(N)*Rfun(NMO2,NGar2,N+2)*Rfun(NMO4,NGar4,N+2)/(R1(N)*RO1X(N)**2)
         SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
    ENDIF
    
	RRCoulomb=SING*H
  
    return
  end function RRCoulomb



  ! ОНДОПНЦПЮЛЛЮ ЮМЮКХГЮ ЩКЕЙРПНММНИ ОКНРМНЯРХ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ
  ! NumeroIter-мнлеп хрепюжхх
  ! NumbreMO-вхякн лнкейскъпмшу нпахрюкеи
  ! N-вхякн рнвей
  ! H-ьюц
  ! NumbreGarmMO(NumbreMO)- люяяхб вхякю цюплнмхй б лнкейскъпмни нпахрюкх
  ! NumbreGarmLMO(NumbreMO)- люяяхб вхякю цюплнмхй б лнкейскъпмни нпахрюкх дн LgarmonicMO(3,NumbreMO) бйкчвюъ ее
  ! NumeroMinGarmon(NumbreMO)- люяяхб лхмхлюкэмнцн мнлепю цюплнмхйх нркхвмни нр мскъ б лнкейскъпмни нпахрюкх
  ! RO1X(N)-люяяхб оепбни опнхгбндмни мнбни оепелеммни он ярюпни
  ! RFun(NumbreMO,NumbreGarmMO(NumbreMO),N+2)-ЛЮЯЯХБ ГМЮВЕМХИ ПЮДХЮКЭМШУ ВЮЯРЕИ  ЦЮПЛНМХЙ ЛНКЕЙСКЪПМШУ НПАХРЮКЕИ ЙНМТХЦСПЮЖХХ 
  subroutine AnalysMOConfig(NumeroIter,NumbreMO,N,H,NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon,RO1X,RFun)  
    implicit none
    integer::NumeroIter,NumbreMO,N
    real(8)::H
    integer,dimension(:)::NumbreGarmMO,NumbreGarmLMO,NumeroMinGarmon
	real(8),dimension(:)::RO1X
	real(8),dimension(:,:,:)::RFun
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::INDEXS,INDXF,INDEKL,NpointMax,Lgarmonik,IDFZ
    real(8)::ROEDLK,ROEDL,SING,SING2,SING1,R0,RN1,RN

   
    5100 FORMAT(5X,I3,5X,I2,5X,F10.8,3X,D8.1)
    
	! мюундхл щкейрпнммсч окнрмнярэ дн Lk цюплнмхй     
	ROEDLK=0.D0
    DO INDXF=NumeroMinGarmon(NumbreMO),NumbreGarmLMO(NumbreMO) 
       !ROEDLK=ROEDLK+sum(RFun(NumbreMO,INDXF,2+1:N)**2/RO1X(1:N)**2)
       SING=0.D0 
	   ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
       IF(N/2..EQ.N/2) THEN
	        ! пюявер он тнплске яхлоянмю
            ! вермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-2,2
               SING2=SING2+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-3,2
               SING1=SING1+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
		    ! пегскэрхпсчыюъ ясллю 
		    R0=(RFun(NumbreMO,INDXF,3)/RO1X(1))**2
		    RN1=(RFun(NumbreMO,INDXF,N+1)/RO1X(N-1))**2 
		    RN=(RFun(NumbreMO,INDXF,N+2)/RO1X(N))**2 
            SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
          ELSE
            ! пювер он тнплске яхлоянмю
		    ! мевермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-1,2
               SING2=SING2+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-2,2
               SING1=SING1+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            ! пегскэрхпсчыюъ ясллю 
		    R0=(RFun(NumbreMO,INDXF,3)/RO1X(1))**2
		    RN=(RFun(NumbreMO,INDXF,N+2)/RO1X(N))**2 
            SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
       ENDIF
       ROEDLK=ROEDLK+SING
	ENDDO
    
	ROEDLK=ROEDLK*H/RFun(NumbreMO,1,2)
	
	! мюундхл щкейрпнммсч окнрмнярэ оняке Lk цюплнмхйх
	ROEDL=0.D0
    DO INDXF=NumbreGarmLMO(NumbreMO)+1,NumbreGarmMO(NumbreMO) 
       !ROEDL=ROEDL+sum(RFun(NumbreMO,INDXF,2+1:N)**2/RO1X(1:N)**2)
       SING=0.D0 
	   ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
       IF(N/2..EQ.N/2) THEN
	        ! пюявер он тнплске яхлоянмю
            ! вермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-2,2
               SING2=SING2+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-3,2
               SING1=SING1+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
		    ! пегскэрхпсчыюъ ясллю 
		    R0=(RFun(NumbreMO,INDXF,3)/RO1X(1))**2
		    RN1=(RFun(NumbreMO,INDXF,N+1)/RO1X(N-1))**2 
		    RN=(RFun(NumbreMO,INDXF,N+2)/RO1X(N))**2 
            SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
          ELSE
            ! пювер он тнплске яхлоянмю
		    ! мевермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-1,2
               SING2=SING2+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-2,2
               SING1=SING1+(RFun(NumbreMO,INDXF,IDFZ+2)/RO1X(IDFZ))**2
	        ENDDO
            ! пегскэрхпсчыюъ ясллю 
		    R0=(RFun(NumbreMO,INDXF,3)/RO1X(1))**2
		    RN=(RFun(NumbreMO,INDXF,N+2)/RO1X(N))**2 
            SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
       ENDIF
       ROEDL=ROEDL+SING
	ENDDO 
	
	ROEDL=ROEDL*H/RFun(NumbreMO,1,2)
	 

	WRITE(27,5100)  NumeroIter,NumbreMO,ROEDLK,ROEDL

   

    return
  end subroutine AnalysMOConfig

  ! ОНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ БЯОНЛЮЦЮРЕКЭМШУ ЛЮЯЯХБНБ ДКЪ ХМРЕЦППНБЮМХЪ
  subroutine VAR(N,H,ALFA,BET,GAMMA,Ral,RO1,RO2,R,RO1X,RO2X,RO3X)
    implicit none
    integer::I1,I,N
    real(8)::H,ALFA,BET,GAMMA,RO1,RO2,RT,ROT,A,Ral
    real(8),dimension(:)::R,RO1X,RO2X,RO3X
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8)::R2,R3,RR,RR2,RDELTA,RDELTA2,RDELTA3,XXDF
          
		       
	! лнкейскъпмюъ яерйю
    RO1=RO1 
    RT=(RO2-BET*DLOG(RO2)-DATAN2(RO2-Ral,GAMMA))/ALFA
    DO I1=1,N
       I=N-I1+1
       ROT=RO1+FLOAT(I)*H
    30 RR=RT-Ral
       RR2=RR*RR
       RDELTA=GAMMA**2+RR2
	   A=(ALFA*RT+BET*DLOG(RT)+DATAN2(RT-Ral,GAMMA)-ROT)*RDELTA/((ALFA*RT+BET)*RDELTA+GAMMA*RT)
       RT=RT*(1.D0-A)
       IF(DABS(A)-1.D-12) 20,20,30  
    20 R(I)=RT
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

        

    

    return 
  end subroutine VAR 

    ! ондопнцпюллю пюяверю япедмецн пюдхсяю лнкейскъпмни нпахрюкх
    ! NomerMO-мнлеп лнкейскъпмни нпахрюкх
    ! NumbreGarmMO(NomerMO)-вхякн цюплнмхй б лнкейскъпмни нпахрюкх
    ! H-ьюц хмрецпхпнбюмхъ
    ! N-вхякн рнвей
    ! RFun(NomerMO,NumbreGarmMO(NomerMO),N+2)-люяяхб гмювемхи пюдхюкэмшу вюяреи оюпжхюкэмшу цюплнмхй 
    ! R(N)-люяяхб гмювемхи пюдхсяю
    ! RO1X(N)-люяяхб оепбни опнхгбндмни мнбни оепелеммни он ярюпни 
    real(8) function MediumRadius(NomerMO,NumbreGarmMO,H,N,RFun,R,RO1X)
      implicit none
	  integer::NomerMO,N
	  integer,dimension(:)::NumbreGarmMO
	  real(8)::H
	  real(8),dimension(:)::R,RO1X
      real(8),dimension(:,:,:)::RFun
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::IX1,IX2,IDFZ
	  real(8)::SumM,SING,SING2,SING1,R0,RN1,RN

      SumM=0.D0
	  ! жхйк он оюпжхюкэмшл цюплнмхйюл
	  DO IX1=1,NumbreGarmMO(NomerMO)
	     !SumInt=0.D0
		 !DO IX2=1,N
         !   SumInt=SumInt+R(IX2)*RFun(NomerMO,IX1,IX2+2)**2/RO1X(IX2)**2
         !ENDDO      
         !SumM=SumM+sum(R(1:N)*RFun(NomerMO,IX1,2+1:N)**2/RO1X(1:N)**2) 
		 
		 SING=0.D0 
	     ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
         IF(N/2..EQ.N/2) THEN
	          ! пюявер он тнплске яхлоянмю
              ! вермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,N-2,2
                 SING2=SING2+R(IDFZ)*(RFun(NomerMO,IX1,IDFZ+2)/RO1X(IDFZ))**2
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,N-3,2
                 SING1=SING1+R(IDFZ)*(RFun(NomerMO,IX1,IDFZ+2)/RO1X(IDFZ))**2
	          ENDDO
		      ! пегскэрхпсчыюъ ясллю 
		      R0=R(1)*(RFun(NomerMO,IX1,3)/RO1X(1))**2
		      RN1=R(N-1)*(RFun(NomerMO,IX1,N+1)/RO1X(N-1))**2 
		      RN=R(N)*(RFun(NomerMO,IX1,N+2)/RO1X(N))**2 
              SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
            ELSE
              ! пювер он тнплске яхлоянмю
			  ! мевермне вхякн рнвей
              SING2=0.D0
		      DO IDFZ=2,N-1,2
                 SING2=SING2+R(IDFZ)*(RFun(NomerMO,IX1,IDFZ+2)/RO1X(IDFZ))**2
	          ENDDO
              SING1=0.D0
		      DO IDFZ=3,N-2,2
                 SING1=SING1+R(IDFZ)*(RFun(NomerMO,IX1,IDFZ+2)/RO1X(IDFZ))**2
	          ENDDO
              ! пегскэрхпсчыюъ ясллю 
		      R0=R(1)*(RFun(NomerMO,IX1,3)/RO1X(1))**2
		      RN=R(N)*(RFun(NomerMO,IX1,N+2)/RO1X(N))**2 
              SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
         ENDIF
         ! ясллхпсел хмрецпюкш
	     SumM=SumM+SING
      ENDDO
     
	  ! гмювемхе япедмецн пюдхсяю 
      MediumRadius=SumM*H
     
	  return
    end function MediumRadius


	! ондопнцпюллю хмрецпхпнбюмхъ
    real(8) function SIMPM(IN,N1,JN,N2,A,N,H,RO1X)
      implicit none
	  integer::IN,N1,JN,N2,N
	  real(8)::SUM,H
	  real(8),dimension(:)::RO1X
      real(8),dimension(:,:,:)::A
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::IDFZ
      real(8)::SING,SING2,SING1,R0,RN1,RN
	 
	  !SUM=0.D0
      !DO I=1,N
      !   SUM=SUM+A(IN,N1,I+2)*A(JN,N2,I+2)/RO1X(I)**2
      !ENDDO
       
      !SIMPM=sum(A(IN,N1,2+1:N)*A(JN,N2,2+1:N)/RO1X(1:N)**2)*H 
	  
	  SING=0.D0 
	  ! опнбепъел вермне вхякн рнвей хкх мевермне вхякн рнвей
      IF(N/2..EQ.N/2) THEN
	        ! пюявер он тнплске яхлоянмю
            ! вермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-2,2
               SING2=SING2+A(IN,N1,IDFZ+2)*A(JN,N2,IDFZ+2)/RO1X(IDFZ)**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-3,2
               SING1=SING1+A(IN,N1,IDFZ+2)*A(JN,N2,IDFZ+2)/RO1X(IDFZ)**2
	        ENDDO
		    ! пегскэрхпсчыюъ ясллю 
		    R0=A(IN,N1,3)*A(JN,N2,3)/RO1X(1)**2
		    RN1=A(IN,N1,N+1)*A(JN,N2,N+1)/RO1X(N-1)**2
		    RN=A(IN,N1,N+2)*A(JN,N2,N+2)/RO1X(N)**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN1)/3.D0+(RN1+RN)*0.5D0
         ELSE
            ! пювер он тнплске яхлоянмю
			! мевермне вхякн рнвей
            SING2=0.D0
		    DO IDFZ=2,N-1,2
               SING2=SING2+A(IN,N1,IDFZ+2)*A(JN,N2,IDFZ+2)/RO1X(IDFZ)**2
	        ENDDO
            SING1=0.D0
		    DO IDFZ=3,N-2,2
               SING1=SING1+A(IN,N1,IDFZ+2)*A(JN,N2,IDFZ+2)/RO1X(IDFZ)**2
	        ENDDO
            ! пегскэрхпсчыюъ ясллю 
		    R0=A(IN,N1,3)*A(JN,N2,3)/RO1X(1)**2
		    RN=A(IN,N1,N+2)*A(JN,N2,N+2)/RO1X(N)**2
            SING=(R0+4.D0*SING2+2.D0*SING1+RN)/3.D0
      ENDIF 

      SIMPM=SING*H
           
      return
    end function SIMPM



	! ондопнцпюллю явхршбюмхъ гмювемхи тсмйжхи мскхбнцн опхакхфемхъ 
    subroutine ReadDatIn(NumbreMO,Npoint,NumbreGarmMOZero,IZONGarmonZero,RfunZero)
      implicit none
      integer::NumbreMO,Npoint
      integer,dimension(:)::NumbreGarmMOZero
      integer,dimension(:,:)::IZONGarmonZero
      real(8),dimension(:,:,:)::RfunZero
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer::IIX,YYX,IZ,I2,IIZ,ierr
	  real(4),allocatable,dimension(:)::RfunIn
	  real(8),allocatable,dimension(:)::Rfun
      
	 
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  !allocate(RfunIn(Npoint+2),stat=ierr)
	  !if(ierr/=0) then
      !   write(*,*) 'ReadDatIn'
	  !   write(*,*) 'MEMORY ON THE FILE "RfunIn" IS NOT SELECTED'
	  !   stop 
	  !endif 
	  allocate(Rfun(Npoint+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'ReadDatIn'
	     write(*,*) 'MEMORY ON THE FILE "Rfun" IS NOT SELECTED'
	     stop 
	  endif 

      
	  DO IIX=1,NumbreMO
         DO YYX=1,NumbreGarmMOZero(IIX)
		    Rfun=0.D0
			! опх мске мнплю тхйяхпсел едхмхжс
			Rfun(2)=1.D0
            IF(IZONGarmonZero(IIX,YYX).NE.0) THEN
               REWIND 10
               IZ=IZONGarmonZero(IIX,YYX)-1
			   IF(IZ.GT.0) THEN 
                  DO I2=1,IZ
                     READ(10)
                  ENDDO
               ENDIF
			  
			   READ(10)(Rfun(I2),I2=1,Npoint+2)
               !DO I2=1,Npoint+2
               !   Rfun(I2)=DBLE(RfunIn(I2)) 
               !ENDDO  
            ENDIF 
            ! гюохяшбюел оюпжхюкэмсч цюплнмхйс 
			RfunZero(IIX,YYX,1)=Rfun(1)
            RfunZero(IIX,YYX,2)=Rfun(2)
            !WRITE(*,*) IIX,YYX,RfunZero(IIX,YYX,1)
			!WRITE(*,*) IIX,YYX,RfunZero(IIX,YYX,2)
			!READ(*,*)
			! оепебндхл тсмйжхч б мемнплхпнбюммне янярнъмхе
			DO IIZ=1,Npoint
			   RfunZero(IIX,YYX,IIZ+2)=Rfun(IIZ+2)*DSQRT(Rfun(2))
			   !WRITE(*,*) IIX,YYX,IIZ,RfunZero(IIX,YYX,IIZ+2),DABS(Rfun(2)),Rfun(2)
			   !READ(*,*) 
			   !WRITE(6,*) IIX,YYX,IIZ,RfunZero(IIX,YYX,IIZ+2)
			ENDDO 
		 ENDDO 
	  ENDDO
	  
	  ! сдюкемхе люяяхбнб хг оълърх 
      !deallocate(RfunIn,stat=ierr)
      !if(ierr/=0) then
      !   write(*,*) 'ReadDatIn'
      !   write(*,*) 'THE FILE "RfunIn" IS NOT REMOVED FROM MEMORY' 
 	  !   stop 
      !endif
      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(Rfun,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'ReadDatIn'
         write(*,*) 'THE FILE "Rfun" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif 

      return
    end subroutine ReadDatIn


    ! ондопнцпюллю явхршбюмхъ гмювемхи тсмйжхи кхцюмдю 
	! явхршбюел тсмйжхх кхцюмдю (тсмйжхх ъбкъчряъ мнплхпнбюммшлх мю едхмхжс)
    subroutine ReadDatInLigand(Nnuclei,Npoint,H,NFunctionLigand,NumbreGarmLigand,IZONGarmonLigand,RfunLigands,RO1X)
      implicit none
      integer::Nnuclei,Npoint
	  real(8)::H
      integer,dimension(:)::NFunctionLigand
      integer,dimension(:,:)::NumbreGarmLigand
	  integer,dimension(:,:,:)::IZONGarmonLigand
	  real(8),dimension(:)::RO1X
      real(8),dimension(:,:,:,:)::RfunLigands
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer::IIX,IIY,YYX,IZ,I2,IIZ,ierr,IIXS
	  real(8)::SUMA,SUMLI,CalculationAlfaNR
	  real(4),allocatable,dimension(:)::RfunIn
	  real(8),allocatable,dimension(:)::Rfun
     
	 
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  !allocate(RfunIn(Npoint+2),stat=ierr)
	  !if(ierr/=0) then
      !   write(*,*) 'ReadDatInLigand'
	  !   write(*,*) 'MEMORY ON THE FILE "RfunIn" IS NOT SELECTED'
	  !   stop 
	  !endif 
	  allocate(Rfun(Npoint+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'ReadDatInLigand'
	     write(*,*) 'MEMORY ON THE FILE "Rfun" IS NOT SELECTED'
	     stop 
	  endif 

      
	  DO IIX=1,Nnuclei
         DO IIY=1,NFunctionLigand(IIX) 
            DO YYX=1,NumbreGarmLigand(IIX,IIY)
		       Rfun=0.D0
               IF(IZONGarmonLigand(IIX,IIY,YYX).NE.0) THEN
                  REWIND 11
                  IZ=IZONGarmonLigand(IIX,IIY,YYX)-1
			      IF(IZ.GT.0) THEN 
                     DO I2=1,IZ
                        READ(11)
                     ENDDO
                  ENDIF 
			      READ(11) (Rfun(I2),I2=1,Npoint+2)   !(RfunIn(I2),I2=1,Npoint+2)
                  !DO I2=1,Npoint+2
                  !   Rfun(I2)=DBLE(RfunIn(I2)) 
                  !ENDDO  
               ENDIF 
               ! гюохяшбюел оюпжхюкэмсч цюплнмхйс 
		       RfunLigands(IIX,IIY,YYX,1)=Rfun(1)
               RfunLigands(IIX,IIY,YYX,2)=1.D0 
               DO IIZ=1,Npoint
			      RfunLigands(IIX,IIY,YYX,IIZ+2)=Rfun(IIZ+2)
			   ENDDO 
		    ENDDO 
            ! мнплхпсел онякедмсч цюплнмхйс дкъ сверю нярюбьеияъ щкейрпнммни окнрмнярх
            SUMA=0.D0
	        DO IIXS=1,NumbreGarmLigand(IIX,IIY)-1
               SUMA=SUMA+sum((RfunLigands(IIX,IIY,IIXS,2+1:Npoint)/RO1X(1:Npoint))**2)
	        ENDDO
            SUMLI=H*sum((RfunLigands(IIX,IIY,NumbreGarmLigand(IIX,IIY),2+1:Npoint)/RO1X(1:Npoint))**2)
	       
			CalculationAlfaNR=DSQRT(1.D0-SUMA*H)/DSQRT(SUMLI) 
			!сЛМЮФЮЕЛ ОНЯКЕДМСЧ ЦЮПЛНМХЙС МЮ ДЮММШИ ЙНЩТТХЖХЕМР
            RfunLigands(IIX,IIY,NumbreGarmLigand(IIX,IIY),2+1:Npoint)=CalculationAlfaNR*RfunLigands(IIX,IIY,NumbreGarmLigand(IIX,IIY),2+1:Npoint)
         ENDDO
	  ENDDO
	     
	  ! сдюкемхе люяяхбнб хг оълърх 
      !deallocate(RfunIn,stat=ierr)
      !if(ierr/=0) then
      !   write(*,*) 'ReadDatInLigand'
      !   write(*,*) 'THE FILE "RfunIn" IS NOT REMOVED FROM MEMORY' 
 	  !   stop 
      !endif
      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(Rfun,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'ReadDatInLigand'
         write(*,*) 'THE FILE "Rfun" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif 

      return
    end subroutine ReadDatInLigand

    
	! ондопнцпюллю нясыеярбкъер явхршбюмхе оюпюлерпнб йнппейрхпсчыху пюанрс опнцпюллш 
	subroutine ReadCorrectionParameters(NumbreIteration,FileParameterCorrection,IreshimSO) 
      implicit none
      integer::NumbreIteration,IreshimSO
      character(90)::FileParameterCorrection
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::NIreshimSO
      56547 FORMAT(2X,'***********************************************************************')
	  56548 FORMAT(2X,'NumbreIteraction=',I6,' The mode of linkage has changed IreshimSO= ',I3)
      
	  ! ябъгшбюел тюик я онрнйнл
	  ! нрйпшбюел тюик дкъ гюохях тсмйжхи
	  OPEN(4123,FILE=FileParameterCorrection)
      
      READ(4123,*)  NIreshimSO
	  
	  ! гюйпшбюел тюик
      CLOSE(4123)
	  
	  ! опнбепъел пефхл пюанрш рнр фе хкх мер
	  IF(IreshimSO.NE.NIreshimSO) THEN
	     IreshimSO=NIreshimSO
         WRITE(6,56547) 
		 WRITE(6,56548) NumbreIteration,IreshimSO
	     WRITE(6,56547)  
	  ENDIF 
	
      return
    end subroutine ReadCorrectionParameters
  

    
    ! ондопнцпюллю нясыеярбкъер гюохяэ опнлефсрнвмшу пегскэрюрнб б тюик 
	! NumbreMO-вхякн лнкейскъпмшу нпахрюкхи
	! N-вхякн рнвей
	! NumbreGarmMO(NumbreMO)-люяяхб вхякю цюплнмхй б лнкейскъпмни нпахрюкх
    ! RfuО(NumbreMO,NumbreGarmMO(NumbreMO),N+2)-люяяхб гмювемхи цюплнмхй лнкйскъпмни нпахрюкх
	! OUTFUNBuffer-хлъ асттепмнцн тюикю 
    subroutine WriteBufferFunction(NumbreMO,N,NumbreGarmMO,Rfun,FUNBuffer)
      implicit none
      integer::NumbreMO,N
	  integer,dimension(:)::NumbreGarmMO
      real(8),dimension(:,:,:)::Rfun
	  character(90)::FUNBuffer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::I,IW,IEND,IBGN,ierr,IZAP,INMO
	  real(4),allocatable,dimension(:)::ROUT
    
    
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(ROUT(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'WriteBufferFunction'
	     write(*,*) 'MEMORY ON THE FILE "ROUT" IS NOT SELECTED'
	     stop 
	  endif  
      
	  ! ябъгшбюел тюик я онрнйнл
	  ! нрйпшбюел тюик дкъ гюохях тсмйжхи
	  OPEN(UNIT=3123,FILE=FUNBuffer,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')  
	 
      DO INMO=1,NumbreMO
         DO I=1,NumbreGarmMO(INMO)
            ! гюохяшбюел ндмнщкейрпнммсч щмепцхч
			ROUT(1)=-SNGL(Rfun(INMO,I,1))
            ! гюохяшбюел мнплхпнбнвмсч йнмярюмрс 
			ROUT(2)=SNGL(Rfun(INMO,I,2))
	        DO IW=1,N
               ROUT(IW+2)=SNGL(Rfun(INMO,I,IW+2)/DSQRT(Rfun(INMO,I,2)))
            ENDDO
            WRITE(3123) (ROUT(IW),IW=1,N+2)  
	     ENDDO
	  ENDDO 

      ! гюйпшбюел тюик
      CLOSE(3123)

      
  
      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(ROUT,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'WriteBufferFunction'
         write(*,*) 'THE FILE "ROUT" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
	 
	 
	  return
    end subroutine WriteBufferFunction


	! ондопнцпюллю нясыеярбкъер	явхршбюмхе хг тюикнб опнлефсрнвмшу пегскэрюрнб 
	! NumbreMO-вхякн лнкейскъпмшу нпахрюкхи
	! N-вхякн рнвей
	! NumbreGarmMO(NumbreMO)-люяяхб вхякю цюплнмхй б лнкейскъпмни нпахрюкх
    ! RfuО(NumbreMO,NumbreGarmMO(NumbreMO),N+2)-люяяхб гмювемхи цюплнмхй лнкйскъпмни нпахрюкх
	! OUTFUNBuffer-хлъ асттепмнцн тюикю 
    subroutine ReadBufferFunction(NumbreMO,N,NumbreGarmMO,Rfun,FUNBuffer)
      implicit none
      integer::NumbreMO,N
	  integer,dimension(:)::NumbreGarmMO
      real(8),dimension(:,:,:)::Rfun
	  character(90)::FUNBuffer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::I,IW,IEND,IBGN,ierr,IZAP,INMO
	  real(4),allocatable,dimension(:)::ROUT
    
    
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(ROUT(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'ReadBufferFunction'
	     write(*,*) 'MEMORY ON THE FILE "ROUT" IS NOT SELECTED'
	     stop 
	  endif  
      
	  ! ябъгшбюел тюик я онрнйнл
	  ! нрйпшбюел тюик дкъ гюохях тсмйжхи
	  OPEN(UNIT=3123,FILE=FUNBuffer,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')  
	 
      DO INMO=1,NumbreMO
         DO I=1,NumbreGarmMO(INMO)
            READ(3123) (ROUT(IW),IW=1,N+2)  
	        ! гюохяшбюел ндмнщкейрпнммсч щмепцхч
			Rfun(INMO,I,1)=-DBLE(ROUT(1))
            ! гюохяшбюел мнплхпнбнвмсч йнмярюмрс 
			Rfun(INMO,I,2)=DBLE(ROUT(2))
	        ! оепебндхл тсмйжхч б мемнплхпнбюммне янярнъмхе
			DO IW=1,N
               Rfun(INMO,I,IW+2)=DBLE(ROUT(IW+2))*DSQRT(Rfun(INMO,I,2))
            ENDDO
		 ENDDO
	  ENDDO 

      ! гюйпшбюел тюик
      CLOSE(3123)

      
  
      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(ROUT,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'ReadBufferFunction'
         write(*,*) 'THE FILE "ROUT" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
	 
	 
	  return
    end subroutine ReadBufferFunction










    ! ондопнцпюллю нясыеярбкъер гюохяэ пегскэрюрнб б тюик 
    subroutine WriteDatOut(IZAP,NumbreMO,N,H,NN,ML,IQ,NumbreGarmMO,R,RO1X,Rfun,OUTFUN)
      implicit none
      integer::IZAP,NumbreMO,N
	  real(8)::H
	  integer,dimension(:)::NN,ML,IQ,NumbreGarmMO
	  real(8),dimension(:)::R,RO1X
      real(8),dimension(:,:,:)::Rfun
	  character(90),dimension(:)::OUTFUN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::I,J,IW,IEND,NU,IBGN,ND,K,QQ,IIIFD,ierr,INGmax,IXD,II,ICVBX
      real(8)::D,P,T,HH,B,C,SUMINT
	  character(5),dimension(2)::LBZZ  
	  real(4),allocatable,dimension(:)::R5OUT
      real(8),allocatable,dimension(:)::R5,R6
      real(8),allocatable,dimension(:,:)::RIntOrt
      
      data LBZZ/'sigma','pi'/	  
    
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(R5OUT(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'WriteDatOut'
	     write(*,*) 'MEMORY ON THE FILE "R5OUT" IS NOT SELECTED'
	     stop 
	  endif  
      allocate(RIntOrt(NumbreMO,NumbreMO),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'WriteDatOut'
	     write(*,*) 'MEMORY ON THE FILE "RIntOrt" IS NOT SELECTED'
	     stop 
	  endif  


       
      56  FORMAT(2X,A90)
      55  FORMAT(2X,I2,A5,' Written from record  ',I5,' Continue from record',I5)
      1   FORMAT(/' Delta RO = ',F7.5,'         F - Function  of P-type ')
      2   FORMAT(/'      R   ',4('  F(',I2,A1,')',1X),10X,4('  F(',I2,A1,')',1X)/)
      3   FORMAT(/' Total and partial charges in sphere of definite radius')
      4   FORMAT(/' Total Q  ',4('  Q(',I2,A1,')',1X),10X,4('  Q(',I2,A1,')',1X)/)
      5   FORMAT(5F9.3,10X,4F9.3)
      98  FORMAT(2X,F15.7,100(1X,D15.7))
      110 FORMAT(/12X,'r',14X,10(I2,A1,13X)/)
      120 FORMAT(/2X,'Radial part Int Ort',7X,10(I2,A5,9X)/)
      123 FORMAT(15X,I2,A5,10(F15.7,1X))
      
	  
	  ! пюявер хмрецпюкнб нпрнцнмюкэмнярх
	  DO I=1,NumbreMO
	     DO J=1,I
	        IF(NumbreGarmMO(I).GT.NumbreGarmMO(J)) THEN
                 INGmax=NumbreGarmMO(J)              
			   ELSE
                 INGmax=NumbreGarmMO(I)
			ENDIF
			IF(NumbreGarmMO(I).EQ.NumbreGarmMO(J)) THEN
               INGmax=NumbreGarmMO(I)
		    ENDIF
			SUMINT=0.D0
			DO IXD=1,INGmax
               SUMINT=SUMINT+SIMPM(I,IXD,J,IXD,Rfun,N,H,RO1X)       
            ENDDO 
            ! гюохяшбюел хмрецпюк нпрнцнмюкэмнярх
			RIntOrt(I,J)=SUMINT
		 ENDDO
	  ENDDO 
      
	  ! бшдювю хмрецпюкнб нпрнцнмюкэмнярх
	  WRITE(6,120)(NN(II),LBZZ(IABS(ML(II))+1),II=1,NumbreMO)
	  DO I=1,NumbreMO
	     WRITE(6,123) NN(I),LBZZ(IABS(ML(I))+1),(RIntOrt(I,J),J=1,I) 
	  ENDDO 
      
 
 
	 
      ! гЮОХЯЭ ТСМЙЖХИ Б ТЮИК
      IBGN=IZAP-1
      IF(IBGN) 50,51,52
  52  REWIND 1
  
      DO I=1,IBGN
         READ (1)
      ENDDO

  51  IEND=IZAP
          
	  WRITE(6,*)
	  DO I=1,NumbreMO
         ICVBX=IEND
	     DO J=1,NumbreGarmMO(I)
		    IEND=IEND+1  
            DO IW=1,N+2
               R5OUT(IW)=SNGL(Rfun(I,J,IW))
            ENDDO
            WRITE(1) (R5OUT(IW),IW=1,N+2)
         ENDDO
		 WRITE(6,56) OUTFUN(I)
         WRITE(6,55)  NN(I),LBZZ(IABS(ML(I))+1),ICVBX,IEND
	  ENDDO

      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R5OUT,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "R5OUT" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
      deallocate(RIntOrt,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "RIntOrt" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif

 
      
  50  CONTINUE

      ! бшдювю б тюик пюдхюкэмшу вюяреи бнкмнбшу тсмйжхи дкъ опнбепйх 
      DO II=1,NumbreMO
         DO J=1,N
            WRITE(70+II,98) R(J),(Rfun(II,I,J+2)/DSQRT(RO1X(J)),I=1,NumbreGarmMO(II))
         ENDDO
      ENDDO
      
    
  
      return
    end subroutine WriteDatOut

	
	! ондпнцпюллю пюяверю онремжхюкю 
	! пюявер нясыеярбкъеряъ я свернл мнплхпнбйх  б ъбмнл бхде
    subroutine POTENSS(AK,NMO1,NGar1,NMO2,NGar2,Rfun,N,H,R6,R1,RO1X)
      implicit none
	  integer::I,I1,N,NMO1,NGar1,NMO2,NGar2
	  real(8)::AK,D,P,T,B,C,H,RRT
	  real(8),dimension(:)::R6
      real(8),dimension(:)::R1,RO1X
      real(8),dimension(:,:,:)::Rfun
      real(8)::Rnorm
      D=0.D0
      P=0.D0
      T=0.D0
      B=H/12.D0
	       
             
	  do I=1,N
	     RRT=RO1X(I)*R1(I)
         C=Rfun(NMO1,NGar1,I+2)*Rfun(NMO2,NGar2,I+2)/(RO1X(I)*RO1X(I))
	     D=(8.D0*P-D+5.D0*C)*B
         T=T+D
         T=T/(1.D0+5.D0*AK*B/RRT)
         D=P
         P=C-T*AK/RRT
         R6(I)=T
      enddo


      D=0.D0
      P=0.D0

      do I1=1,N
         I=N-I1+1
	     RRT=RO1X(I)*R1(I)
         C=(2.D0*AK+1.D0)*R6(I)/RRT
         D=(8.D0*P-D-5.D0*C)*B
         T=T-D
         T=T/(1.D0+5.D0*(AK+1.D0)*B/RRT)
         D=P
         P=(AK+1.)*T/RRT-C
         R6(I)=T
      enddo

      return
    end subroutine POTENSS




    
    ! ондпнцпюллю пюяверю онремжхюкю 
	! пюявер нясыеярбкъеряъ я свернл мнплхпнбйх  б ъбмнл бхде
    subroutine POTENS(AK,NMO1,NGar1,NMO2,NGar2,Rfun,N,H,R6,R1,RO1X)
      implicit none
	  integer::I,I1,N,NMO1,NGar1,NMO2,NGar2
	  real(8)::AK,D,P,T,B,C,H,RRT
	  real(8),dimension(:)::R6
      real(8),dimension(:)::R1,RO1X
      real(8),dimension(:,:,:)::Rfun
      real(8)::Rnorm
      D=0.D0
      P=0.D0
      T=0.D0
      B=H/12.D0
	  ! МНПЛХПНБЙЮ ТСМЙЖХИ
      Rnorm=1.d0/DSQRT(Rfun(NMO1,NGar1,2)*Rfun(NMO2,NGar2,2))
             
	  do I=1,N
	     RRT=RO1X(I)*R1(I)
         C=Rfun(NMO1,NGar1,I+2)*Rfun(NMO2,NGar2,I+2)*Rnorm/(RO1X(I)*RO1X(I))
	     D=(8.D0*P-D+5.D0*C)*B
         T=T+D
         T=T/(1.D0+5.D0*AK*B/RRT)
         D=P
         P=C-T*AK/RRT
         R6(I)=T
      enddo


      D=0.D0
      P=0.D0

      do I1=1,N
         I=N-I1+1
	     RRT=RO1X(I)*R1(I)
         C=(2.D0*AK+1.D0)*R6(I)/RRT
         D=(8.D0*P-D-5.D0*C)*B
         T=T-D
         T=T/(1.D0+5.D0*(AK+1.D0)*B/RRT)
         D=P
         P=(AK+1.)*T/RRT-C
         R6(I)=T
      enddo

      return
    end subroutine POTENS  



    ! ондпнцпюллю пюяверю онремжхюкю мю тсмйжхъу кхцюмдю
	! пюявер нясыеярбкъеряъ я свернл мнплхпнбйх  б ъбмнл бхде
    subroutine POTENSLigand(AK,Nligand1,NFunLigand1,NGar1,Nligand2,NFunLigand2,NGar2,Rfun,N,H,R6,R1,RO1X)
      implicit none
	  integer::I,I1,N,Nligand1,NFunLigand1,NGar1,Nligand2,NFunLigand2,NGar2
	  real(8)::AK,D,P,T,B,C,H,RRT
	  real(8),dimension(:)::R6
      real(8),dimension(:)::R1,RO1X
      real(8),dimension(:,:,:,:)::Rfun
      real(8)::Rnorm
      D=0.D0
      P=0.D0
      T=0.D0
      B=H/12.D0
	  ! МНПЛХПНБЙЮ ТСМЙЖХИ
      Rnorm=1.d0/DSQRT(Rfun(Nligand1,NFunLigand1,NGar1,2)*Rfun(Nligand2,NFunLigand2,NGar2,2))
             
	  do I=1,N
	     RRT=RO1X(I)*R1(I)
         C=Rfun(Nligand1,NFunLigand1,NGar1,I+2)*Rfun(Nligand2,NFunLigand2,NGar2,I+2)*Rnorm/(RO1X(I)*RO1X(I))
	     D=(8.D0*P-D+5.D0*C)*B
         T=T+D
         T=T/(1.D0+5.D0*AK*B/RRT)
         D=P
         P=C-T*AK/RRT
         R6(I)=T
      enddo


      D=0.D0
      P=0.D0

      do I1=1,N
         I=N-I1+1
	     RRT=RO1X(I)*R1(I)
         C=(2.D0*AK+1.D0)*R6(I)/RRT
         D=(8.D0*P-D-5.D0*C)*B
         T=T-D
         T=T/(1.D0+5.D0*(AK+1.D0)*B/RRT)
         D=P
         P=(AK+1.)*T/RRT-C
         R6(I)=T
      enddo

      return
    end subroutine POTENSLigand 



	! ондпнцпюллю пюяверю онремжхюкю дкъ налеммни вюярх бгюхлндеиярбхъ мю тсмйжхъу кхцюмдю
	! ме мнплхпсеряъ лнкейскъпмюъ нпахрюкэ дкъ йнрнпни нясыеярбкъеряъ пеьемхе
	! пюявер нясыеярбкъеряъ я свернл мнплхпнбйх  б ъбмнл бхде
    subroutine POTENGLigand(AK,Nligand1,NFunLigand1,NGar1,Nligand2,NFunLigand2,NGar2,Rfun,N,H,R6,R1,RO1X)
      implicit none
	  integer::I,I1,N,Nligand1,NFunLigand1,NGar1,Nligand2,NFunLigand2,NGar2
	  real(8)::AK,D,P,T,B,C,H,RRT
	  real(8),dimension(:)::R6
      real(8),dimension(:)::R1,RO1X
      real(8),dimension(:,:,:,:)::Rfun
      real(8)::Rnorm
      D=0.D0
      P=0.D0
      T=0.D0
      B=H/12.D0
	  ! МНПЛХПНБЙЮ ТСМЙЖХИ
      Rnorm=1.d0/DSQRT(Rfun(Nligand1,NFunLigand1,NGar1,2))
             
	  do I=1,N
	     RRT=RO1X(I)*R1(I)
         C=Rfun(Nligand1,NFunLigand1,NGar1,I+2)*Rfun(Nligand2,NFunLigand2,NGar2,I+2)*Rnorm/(RO1X(I)*RO1X(I))
	     D=(8.D0*P-D+5.D0*C)*B
         T=T+D
         T=T/(1.D0+5.D0*AK*B/RRT)
         D=P
         P=C-T*AK/RRT
         R6(I)=T
      enddo


      D=0.D0
      P=0.D0

      do I1=1,N
         I=N-I1+1
	     RRT=RO1X(I)*R1(I)
         C=(2.D0*AK+1.D0)*R6(I)/RRT
         D=(8.D0*P-D-5.D0*C)*B
         T=T-D
         T=T/(1.D0+5.D0*(AK+1.D0)*B/RRT)
         D=P
         P=(AK+1.)*T/RRT-C
         R6(I)=T
      enddo

      return
    end subroutine POTENGLigand 




	
	! ондпнцпюллю пюяверю онремжхюкю дкъ налеммни вюярх бгюхлндеиярбхъ
	! ме мнплхпсеряъ лнкейскъпмюъ нпахрюкэ дкъ йнрнпни нясыеярбкъеряъ пеьемхе
	! пюявер нясыеярбкъеряъ я свернл мнплхпнбйх  б ъбмнл бхде
    subroutine POTENG(AK,NMO1,NGar1,NMO2,NGar2,Rfun,N,H,R6,R1,RO1X)
      implicit none
	  integer::I,I1,N,NMO1,NGar1,NMO2,NGar2
	  real(8)::AK,D,P,T,B,C,H,RRT
	  real(8),dimension(:)::R6
      real(8),dimension(:)::R1,RO1X
      real(8),dimension(:,:,:)::Rfun
      real(8)::Rnorm
      D=0.D0
      P=0.D0
      T=0.D0
      B=H/12.D0
	  ! МНПЛХПНБЙЮ ТСМЙЖХИ
      !Rnorm=1.d0/DSQRT(Rfun(NMO1,NGar1,2))
             
	  do I=1,N
	     RRT=RO1X(I)*R1(I)
         C=Rfun(NMO1,NGar1,I+2)*Rfun(NMO2,NGar2,I+2)/(RO1X(I)*RO1X(I))
	     D=(8.D0*P-D+5.D0*C)*B
         T=T+D
         T=T/(1.D0+5.D0*AK*B/RRT)
         D=P
         P=C-T*AK/RRT
         R6(I)=T
      enddo


      D=0.D0
      P=0.D0

      do I1=1,N
         I=N-I1+1
	     RRT=RO1X(I)*R1(I)
         C=(2.D0*AK+1.D0)*R6(I)/RRT
         D=(8.D0*P-D-5.D0*C)*B
         T=T-D
         T=T/(1.D0+5.D0*(AK+1.D0)*B/RRT)
         D=P
         P=(AK+1.)*T/RRT-C
         R6(I)=T
      enddo

      return
    end subroutine POTENG    
	 
	 
end module mmolecule
