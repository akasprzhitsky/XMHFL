!! PROGRAM MAKES SOLUTIONS OF THE SYSTEM OF INTEGRO-DIFFERENTIAL HARTRI-FOCAL EQUATIONS FOR MOLECULES
!! IN THE APPROXIMATION OF THE BORN-OPPENGEIMER. In the case of a single-centered representation of molecular orbits.
!! In the present version of the program, one of the variants of the calculation of the completeness of the decomposition of the molecular orbital
!! AT THIS, THE SOLUTION IS CREATED FOR THE FINAL NUMBER OF HARMONICS, THE HIGHER HARMONICS ARE FORMATED BY THE BAND WITH THE LIGANDA V FUNCTIONS
!! LgarmonicMO (3, I). LIGANDA FUNCTION - NUCLEAR FUNCTIONS DECLINED ABOUT THE COORDINATE OF THE MOLECULAR SYSTEM.

program XMoleculeHFL 
   real(4) Times,Timef
 
  !! WE MAKE THE CALCULATION OF TIME
   call CPU_TIME(Times) 
 
   !! main subroutine
   call Molecule
      
   !! WE MAKE THE CALCULATION OF TIME  
   call CPU_TIME(Timef) 
   
   WRITE(6,1) (Timef-Times)/60.
 
 1 FORMAT(' Commercial    time :',F10.2,' min.')
   stop
end program XMoleculeHFL 

 

!! PROGRAM MAKES WAVE FUNCTIONS RECEIVED BY SOLUTION
!! Systems of the Hartree-Fock equations for a molecule in a single-center representation of molecular orbits
subroutine Molecule
  use mmolecule, only:DHFMO,VAR,ReadDatIn,ReadDatInLigand,WriteDatOut
  implicit real(8)(A-H,O,P,R-Z)
  
      
!! MASSIVE FOR CALCULATION 
  integer,allocatable,dimension(:)::NN,ML,IQ,NumbreGarmMO,NumbreGarmLMO
  integer,allocatable,dimension(:,:)::LgarmonicMO
  integer,allocatable,dimension(:,:,:,:)::MassivMLaf,MassivMLbg
  real(8),allocatable,dimension(:)::Z,CoorR,R,RO1X,RO2X,RO3X
  real(8),allocatable,dimension(:,:)::CoorA,EEX
  real(8),allocatable,dimension(:,:,:)::RcoffAF,RcoffBG
  real(8),allocatable,dimension(:,:,:)::RFun
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,allocatable,dimension(:)::NumbreGarmMOZero,NumbreGarmLMOZero,NumreISzam
  integer,allocatable,dimension(:,:)::LgarmonicMOZero,IZONGarmonZero,ISCalculZZ
  real(8),allocatable,dimension(:,:,:)::RfunZero
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,allocatable,dimension(:)::NFunctionLigand,NumbreLigand
  integer,allocatable,dimension(:,:)::MLligands,NnLigands,NLigands,Lligands,NFunLigands,NumbreGarmLigand,NomeroGPOT
  integer,allocatable,dimension(:,:,:)::NumbreFunctionLig,LgarmonicLigand,IZONGarmonLigand 
  real(8),allocatable,dimension(:,:,:,:)::RfunLigands
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,allocatable,dimension(:)::NtypIndexCoffNucleusElectron
  integer,allocatable,dimension(:,:)::NtypIndexCoffF,NtypIndexCoffG
  integer,allocatable,dimension(:,:,:)::MLNucleusElectron 
  real(8),allocatable,dimension(:,:)::RcoffNucleusElectron
  integer,dimension(1)::NomeroFPOT
  character(5),dimension(2)::LB 
  character(1),dimension(20)::LB1   
  character(90),dimension(100)::OUTFUN
  character(90)::DATFIL,DAFILES,FILE1,FILE2,OUTFIL,TITLE,TITWRK,OUTFTOTALE,FILEZERO,OUTFILEOE,FUNBuffer,FUNBuffer1,FUNBuffer2,FUNBuffer3,FUNBuffer4,OUTFILSEARCHE,OUTFILSPECTRUME,FILELigands,OUTFILAlfaLigands,FilePOT 
  character(90)::FileElectronicDensity,FileMORnorm,FileParameterCorrection       
   
                
  data LB/'sigma','pi'/
  data LB1/'s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u','v','w','x','y'/   
  1500 FORMAT(5X,' The report of program XMHFL ver 1.14.1.5 March, 2020  Build 22.03.2020')
  1000 FORMAT(1X,'Configuration nucleus')
  1010 FORMAT(2X,'Z= ',F3.0,' Rmol= ',F10.5,' a.e. ','Angle Teta= ',F10.5,' radian ',' Angle Fi= ',F10.5,' radian ') 
  1020 FORMAT(1X,'IndexHartreeFock=',I2,1X,'NumbreMO=',I2,3X,'N=',I4,2X,'H=',F5.3,3X,'Rmax=',F10.3,'  Eps=',D7.1,2X,'Izap=',I3,' IreshimSO=',I2,' IndexExit=',I2)	
  1030 FORMAT(1X,'Parameters setky')
  1040 FORMAT(2X,'  Alfa=',F14.8,2X,' Beta=',F6.3,'  Gamma=',F14.8)  
  1050 FORMAT(/' Electron configuration :   ',20(I3,A5,I2)/)
  1060 FORMAT(1X,'Parameters one-centered approximation')
  1070 FORMAT(2X,I3,A5,I2,' Ngarmonik= ',I3,' Lmin= ',I3,' Lkgran= ',I3,' Lmax= ',I3)
  1080 FORMAT(/' Electroctatic coefficients ',I4,'(',I4,' Exchange )')
  1090 FORMAT(1X,'Coefficients of a direct Coulomb interaction') 
  1100 FORMAT(1X,'Coefficients of an exchange Coulomb interaction')
  1110 FORMAT(13X,20(I2,A5,I2,1X))
  1120 FORMAT(2X,I2,A5,I2,1X,20(F9.5,1X))
  1130 FORMAT(2X,'Type of a Coulomb interaction N= ',I2) 
  1140 FORMAT(2X,'       ML',4X,20(I3,8X)) 
  1150 FORMAT(2X,'Low bound for one-electron energy (Ry)')
  1160 FORMAT(2X,'EnergyXF= ',F10.5,' Ry')  
  1170 FORMAT(1X,'Parameters one-centered approximation for zero approximation')
  2000 FORMAT(1X,'Function Ligands') 
  2100 FORMAT(2X,'Ligand Numbre=',I2, ' Z=',F3.0,' Functions of Ligands: ',50(A5,'(',I2,A1,')',2X))
  2200 FORMAT(1X,'Structure Molecular Orbital') 	
  2300 FORMAT(/' Molecular Orbital :   ',20(I3,A5,I2)/)
  2400 FORMAT(1X,'Parameters functions of ligands')
  2500 FORMAT(2X,'Ligand Numbre=',I2, ' Z=',F3.0,' Function: ',A5,'(',I2,A1,')',' Ngarmonik= ',I3,' Lmin= ',I3,' Lmax= ',I3)
  2600 FORMAT(3X,' R(1)=',E8.3,'      R(N)=',F10.3)
  2700 FORMAT(2X,'  IS1     IS2 ','  ML1 ML2 ML3 ML4', '   Fcoff')
  2800 FORMAT(2X,'  IS1     IS2 ','  ML1 ML2 ML3 ML4', '   Gcoff')          
  2900 FORMAT(2X,I2,A5,1X,I2,A5,1X,I2,2X,I2,2X,I2,2X,I2,2X,F9.5)
  3000 FORMAT(1X,'Coefficients of interaction of an electron with nucleus of a molecule') 
  3100 FORMAT(2X,'  IS  ','  ML1 ML2 ', '   Rcoff')
  3200 FORMAT(2X,I2,A5,1X,I2,2X,I2,2X,F9.5)

  WRITE(*,'(A\)') ' Input TITLE file name or * (MTTITLE.DAT) or END--->'
  
  READ(*,'(A90)') TITLE
     
  IF(TITLE.EQ.'*') THEN
     TITLE='MTTITLE.DAT'
  ENDIF

  IF(TITLE.EQ.'END'.OR.TITLE.EQ.'end') THEN
     RETURN
  ENDIF
  
 
  OPEN (UNIT = 5, FILE = TITLE)
  !! DATA FILE NAME
  READ (5, '(A90)') DATFIL
  !! FILE FILE NAME IN WHICH CALCULATION RESULT WILL BE WRITTEN
  READ (5, '(A90)') FILE1
  !! FILE NAME OF CALCULATION RESULTS
  READ (5, '(A90)') OUTFIL
  !! FILE NAME OF INTERMEDIATE RESULTS OF HARMONIZATION (SINGLE-ELECTRONIC ENERGIES OF MOLECULAR ORBITALS)
  READ (5, '(A90)') OUTFILEOE
  !! FILE NAME IN WHICH COEFFICIENTS OF MIXING LIGANDA FUNCTIONS AT EVERY ITERATION WILL BE RECORDED
  READ (5, '(A90)') OUTFILAlfaLigands
  !! FILE NAME IN WHICH THE SEARCH OF THE SOLUTION OF THE SYSTEM OF EQUATIONS IS WRITTEN
  READ (5, '(A90)') OUTFILSEARCHE
  !! FILE NAME IN WHICH DEPENDENCE OF THE INTEGRAL OF ORTHOGONALITY FROM ENERGY IS RECORDED
  !! IN THE EVENT OF ROOT ABSENCE
  READ (5, '(A90)') OUTFILSPECTRUME
  !! FILE NAME FILE CONTAINING ZERO APPROXIMATION
  READ (5, '(A90)') FILEZERO
  !! FILE NAME FOR RECORDING INTERIM RESULTS
  READ (5, '(A90)') FUNBuffer
  READ (5, '(A90)') FUNBuffer1
  READ (5, '(A90)') FUNBuffer2
  READ (5, '(A90)') FUNBuffer3
  READ (5, '(A90)') FUNBuffer4
  !! FILE NAME IN WHICH FUNCTIONS OF LIGANDS RECORDED
  READ (5, '(A90)') FILELigands
  !! FILE NAME IN WHICH DATA ABOUT ELECTRONIC DENSITY IS RECORDED
  READ (5, '(A90)') FileElectronicDensity
  !! FILE NAME IN WHICH NORMAL PERMANENT RECORDED
  READ (5, '(A90)') FileMORnorm
  !! FILE NAME IN WHICH POTENTIALS WILL BE WRITTEN
  READ (5, '(A90)') FilePOT
  !! FILE NAME IN WHICH PARAMETERS CORRECT THE WORK OF THE PROGRAM
  READ (5, '(A90)') FileParameterCorrection
  
  write(*,*) DATFIL
  write(*,*) FILE1
  write(*,*) OUTFIL
  write(*,*) OUTFILEOE
  write(*,*) FILEZERO
  write(*,*) FILELigands 
  write(*,*) FilePOT
  write(*,*) FileParameterCorrection 
 
  !! BINDING WITH THE FLOW  	
  OPEN(3,FILE=DATFIL)
  OPEN(6,FILE=OUTFIL)
  OPEN(7,FILE=OUTFILEOE)
  OPEN(17,FILE=OUTFILSEARCHE)
  OPEN(18,FILE=OUTFILSPECTRUME)
  OPEN(25,FILE=OUTFILAlfaLigands)
  OPEN(27,FILE=FileElectronicDensity)
  OPEN(29,FILE=FileMORnorm)   
 
  
 
  REWIND 3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
   !!!!!!!!! BLOCK OF DATA INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !! INTRODUCING INFORMATION ABOUT CONFIGURATION OF MOLECULE NUCLEI
   !! Nnuclei-NUMBER OF NUCLEI (DO NOT CONSIDER NUCLEAR AT THE BEGINNING OF THE COORDINATE)
  
  77 READ(3,*) Nnuclei
 
  !! STOPPING PARAMETER OF THE PROGRAM
  IF(Nnuclei.EQ.0) THEN 
     RETURN
  ENDIF  
  
  !! We select the memory for the mass
  allocate(Z(Nnuclei),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "Z" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(CoorR(Nnuclei),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "CoorR" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(CoorA(Nnuclei,2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "CoorA" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NFunctionLigand(Nnuclei),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NFunctionLigand" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
    

   !! INTRODUCE THE CONFIGURATION OF MOLECULE NUCLEI
   !! Z0-CHARGE OF THE NUCLEI AT THE BEGINNING OF THE COORDINATE
   !! Z (Nnuclei) -MASSIVE CHARGES OF NUCLEI
   !! CoorR (Nnuclei) -MASSIVE RADIAL COORDINATES
   !! CoorR (Nnuclei) -RADIAL NUCLEAR COORDINATE (R-RADIUS)
   !! CoorA (Nnuclei, 2) -MASSIVE NORDIC COORDINATE
   !! CoorA (Nnuclei, 1) -THE COORDINATE -RTeta-ANGLE (0 = <RTeta <= PI)
   !! CoorA (Nnuclei, 2) - COORDINATE CO-ORDINATE -RFu-ANGLE (0 = <RFu <= 2 * PI)
   !! NFunctionLigand (Nnuclei) -Number of Ligand functions (functions used to describe the higher harmonics of the representation of a molecular orbital)  READ(3,*) Z0   
  DO IIZ=1,Nnuclei
     READ(3,*)  Z(IIZ),CoorR(IIZ),CoorA(IIZ,1),CoorA(IIZ,2),NFunctionLigand(IIZ)
  ENDDO

  !! FIND THE MAXIMUM NUMBER OF FUNCTIONS
  NFunctionLigandMAX=NFunctionLigand(1) 
  DO IIZ=2,Nnuclei
     IF(NFunctionLigand(IIZ).GT.NFunctionLigandMAX) THEN
        NFunctionLigandMAX=NFunctionLigand(IIZ)
	 ENDIF
  ENDDO
  
  !We select the memory for the Array  
  allocate(MLligands(Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "MLligands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NnLigands(Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NLigands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(Lligands(Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "Lligands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif

   !! CONSIDER THE TYPES OF LIGAND FUNCTIONS
   !! MLligands (IIZ, IXID) -MASSIVE OF THE PROJECTS OF LIGANDA FUNCTIONS
   !! NLigands (IIZ, IXID) -MASSIVE OF THE MAIN QUANTUM NUMERALS OF LIGAND FUNCTIONS
   !! Lligands (IIZ, IXID) -MASSIVE OF ORBITAL QUANTUM NUMBERS
  DO IIZ=1,Nnuclei
     READ(3,*) (MLligands(IIZ,IXID),NnLigands(IIZ,IXID),Lligands(IIZ,IXID),IXID=1,NFunctionLigand(IIZ))
  ENDDO
  
  WRITE(6,1500)
  WRITE(6,*) 
  ! ÂÛÂÎÄÈÌ ÍÀ ÏÅ×ÀÒÜ ÊÎÍÔÈÃÓÐÀÖÈÞ ßÄÅÐ
  WRITE(6,1000)
  
  IF(Z0.NE.0.D0) THEN
    WRITE(6,1010) Z0,0.D0,0.D0,0.D0
  ENDIF

  DO IIZ=1,Nnuclei
      WRITE(6,1010) Z(IIZ),CoorR(IIZ),CoorA(IIZ,1),CoorA(IIZ,2)
  ENDDO
  WRITE(6,*)
  WRITE(6,2000)
  DO IIZ=1,Nnuclei
     WRITE(6,2100)  IIZ,Z(IIZ),(LB(IABS(MLligands(IIZ,IXID))+1),NnLigands(IIZ,IXID),LB1(Lligands(IIZ,IXID)+1),IXID=1,NFunctionLigand(IIZ))
  ENDDO

      
  !! INTRODUCING CONFIGURATION INFORMATION
  !! IndexHartreeFock-PARAMETER INDICATING WHAT EQUATION IS DECIDED
  !! IndexHartreeFock = 0-HARTRI EQUATION
  !! IndexHartreeFock = 1-HARTRI-FOCA EQUATION
  !! NumbreMO-NUMBER OF MOLECULAR ORBITALS IN CONFIGURATION
  !! PARAMETER PARAMETERS
  !! ALFA
  !! BET
  !! GAMMA
  !! N is the number of points
  !! H-step
  !! Rmax-The maximum value of the integration radius
  !! EPS-accuracy of calculation (DEVIATION BETWEEN TWO SOLUTIONS)
  !! EPSXR-accuracy of calculation (DETERMINATION OF SINGLE-ELECTRON ENERGIES IN DECISION OF EQUATIONS)
  !! IZAP-record number for the first harmonic of the first molecular orbital
  !! IKLZERO-KEY APPROACHING THE TYPE OF THE NULL APPROACH
  !! AT IKLZERO = 0-NULL APPROXIMATION IS "BAD" CALCULATION FOR THE FIRST ITERATION IS IMPLEMENTED WITHOUT EXCHANGE
  !! AT IKLZERO = 1-NULL APPROACH IS A "GOOD" CALCULATION ON THE FIRST ITERATION IS IMPLEMENTED WITH THE EXCHANGE
  !! IreshimSO-PARAMETER INDICATING THE TYPE OF HARMONIZATION OF MOLECULAR ORBITALS
  !! IreshimSO = 0-STANDARD MODE OF CODE
  !! IreshimSO = 1- "STRENGTHENING" HARMONIZATION MODE THE BINDING COEFFICIENTS ARE EXPRESSED
  !! IndexExit-PARAMETER INDICATING EMERGENCY STOPPING OF THE PROGRAM
  !! IndexExit = 0-PROGRAM WORKS WITHOUT FAILURES
  !! IndexExit = 1-HAS FAILED TO REALIZE INFORMATION FROM EMERGENCY FILES
  READ(3,*) IndexHartreeFock,NumbreMO,BET,GAMMA,N,H,Rmax,EPS,EPSXR,IZAP,IKLZERO,IreshimSO,IndexExit

  WRITE(6,*) 
  WRITE(6,1020) IndexHartreeFock,NumbreMO,N,H,Rmax,EPS,IZAP,IreshimSO,IndexExit
  
  
  !! WE READ THE FILE NAMES WHICH THE RADIAL PARTS OF MOLECULAR ORBITALS WILL BE WRITTEN
  DO IIZ = 1, NumbreMO
      READ (5, '(A90)') OUTFUN (IIZ)
      !! BINDING WITH THE FLOW
      OPEN (70 + IIZ, FILE = OUTFUN (IIZ))
   ENDDO
  
   !! BINDING WITH THE FLOW
   OPEN (147, FILE = FilePOT, ACCESS = 'DIRECT', RECL = 8 * N)
   NomeroFPOT (1) = 147

   !! Calculation of the integration parameters of ALFA, BET
   !! INTERVAL FOR PO
   RO1 = -15.D0
   RO2 = RO1 + FLOAT (N) * H
  
  ! ÎÏÐÅÄÅËßÅÌ ALFA 
  ALFA=(RO2-BET*DLOG(Rmax))/Rmax
  RTEMP=ALFA*CoorR(1)+BET*DLOG(CoorR(1))-RO1
  ALFA=(RO1+H*DINT(RTEMP/H)-BET*DLOG(CoorR(1)))/CoorR(1)

 !! WE MAKE THE CALCULATION OF PARAMETERS OF THE MESH
  allocate(R(N),stat=ierr)
  if(ierr/=0) then
      write(*,*) 'Molecule MEMORY ON THE FILE "R" IS NOT SELECTED'
      read(*,*)
	  stop 
  endif 
  allocate(RO1X(N),stat=ierr)
  if(ierr/=0) then
      write(*,*) 'Molecule MEMORY ON THE FILE "RO1X" IS NOT SELECTED'
      read(*,*)
	  stop 
  endif 
  allocate(RO2X(N),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule MEMORY ON THE FILE "RO2X" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 
  allocate(RO3X(N),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule MEMORY ON THE FILE "RO3X" IS NOT SELECTED'
      read(*,*)
	 stop 
  endif 
 
  ! FORMING WEIGHTS FOR THIS NET
  call VAR(N,H,ALFA,BET,GAMMA,CoorR(1),RO1,RO2,R,RO1X,RO2X,RO3X)
  

    
  WRITE(6,1030) 
  WRITE(6,1040) ALFA,BET,GAMMA 
  WRITE(6,2600) R(1),R(N)
 

  ! We highlight the memory for the massifs
  allocate(NN(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NN" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(ML(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "ML" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(IQ(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "IQ" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NumbreGarmMO(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NumbreGarmMO" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NumbreGarmLMO(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NumbreGarmLMO" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(LgarmonicMO(3,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "LgarmonicMO" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif



   !! INTRODUCE THE ELECTRONIC CONFIGURATION OF MOLECULES
   !! NN (NumbreMO) -MASSIVE OF MAIN QUANTUM NUMBERS
   !! ML (NumbreMO) -MASSIVE OF MOLECULAR ORBITAL PROJECTIONS
   !! IQ (NumbreMO) -MASSIVE NUMBERS OF FILLING (NUMBER OF ELECTRONS ON MOLECULAR ORBITAL)
  READ(3,*)(NN(I),ML(I),IQ(I),I=1,NumbreMO)

  WRITE(6,1050) (NN(I),LB(IABS(ML(I))+1),IQ(I),I=1,NumbreMO)
  WRITE(6,*)

  !! We introduce the parameters of the single-center decomposition of molecular orbits
  !! NumbreGarmMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL
  !! NumbreGarmLMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL TO LgarmonicMO (3, NumbreMO) INCLUDING THIS HARMONIC
  !! LgarmonicMO (3, NumbreMO) -MASSIVE VALUES OF ORBITAL MOMENTS OF MOLECULAR ORBITALS
  !! LgarmonicMO (1, NumbreMO) - MINIMUM VALUE
  !! LgarmonicMO (2, NumbreMO) -MAXIMUM VALUE
  !! LgarmonicMO (3, NumbreMO) - ORBITAL MOMENT IMPORTANCE AFTER WHICH CALCULATION IS DONE ON LEGAND FUNCTIONS (functions do not change when solving)
  DO I=1,NumbreMO
     READ(3,*) NumbreGarmMO(I),LgarmonicMO(1,I),LgarmonicMO(2,I),LgarmonicMO(3,I)
  ENDDO 

 

  !! We select the memory for the array NFunctionLigandMAX
  allocate(NumbreLigand(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NumbreLigand" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NLigands(NumbreMO,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NLigands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NFunLigands(NumbreMO,Nnuclei),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NFunLigands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NumbreFunctionLig(NumbreMO,Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NFunLigands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(NumbreGarmLigand(Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NumbreGarmLigand" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(LgarmonicLigand(2,Nnuclei,NFunctionLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "LgarmonicLigand" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif


  
   !! INTRODUCE THE STRUCTURE OF A MOLECULAR ORBITAL CONSTRUCTED FROM THE FUNCTIONS OF LIGANDS
   !! NumbreLigand (I) - NUMBER OF LIGANDS OF FUNCTIONS USED FOR A MOLECULAR ORBITAL
   !! NLigands (I, II) -MASSIVE OF LIGAND NUMBERS USED FOR MOLECULAR ORBITALS
   !! NFunLigands (I, IIXZ) -MASSIVE NUMBER OF FUNCTIONS OF THIS LIGAND USED FOR THIS MOLECULAR ORBITAL
  DO I=1,NumbreMO 
     READ(3,*) NumbreLigand(I),(NLigands(I,II),II=1,NumbreLigand(I))
     !! WE WRITE THE NUMBERS OF FUNCTIONS 
	 DO IIXZ=1,NumbreLigand(I)
	    READ(3,*) NFunLigands(I,IIXZ),(NumbreFunctionLig(I,IIXZ,IIYZ),IIYZ=1,NFunLigands(I,IIXZ))
     ENDDO
  ENDDO
  
  NfunLigMax=0 
  !! DETERMINATE THE MAXIMAL NUMBER OF FUNCTIONS OF THE ORIGINAL ORBITAL
  DO I=1,NumbreMO 
     NfunLigI=0
	 DO IIXZ=1,NumbreLigand(I)
        NfunLigI=NfunLigI+NFunLigands(I,IIXZ) 
	 ENDDO
     IF(NfunLigI.GT.NfunLigMax) THEN
        NfunLigMax=NfunLigI 
	 ENDIF
  ENDDO

  !! DEFINITION Lk CORRECT FOR THE IMPLEMENTATION OF THE ÌÎ
  IF(NfunLigMax.NE.0) THEN
     DO I=1,NumbreMO
        LgarmonicMO(3,I)=LgarmonicMO(3,I)+NfunLigMax-1
	 ENDDO
  ENDIF

  !! WRITE THE NUMBER OF HARMONICS
  !! NumbreGarmLMO (NumbreMO) - MASSIVE NUMBER OF HARMONICS IN THE MOLECULAR ORBITAL TO LgarmonicMO (3, NumbreMO) INCLUDING THIS HARMONIC
   DO I=1,NumbreMO
      NumbreGarmLMO(I)=LgarmonicMO(3,I)-LgarmonicMO(1,I)+1
   ENDDO

!! INFORMATION ABOUT MOLECULAR ORBITALS
  write(6,1060)
  DO I=1,NumbreMO
     WRITE(6,1070) NN(I),LB(IABS(ML(I))+1),IQ(I),NumbreGarmMO(I),LgarmonicMO(1,I),LgarmonicMO(3,I),LgarmonicMO(2,I)
  ENDDO
  
  WRITE(6,*)
  WRITE(6,2200)
  DO I=1,NumbreMO  
     WRITE(6,2300) NN(I),LB(IABS(ML(I))+1),IQ(I)
	!! CYCLE ON LIGANDS
	 DO IJ=1,NumbreLigand(I)
        WRITE(6,2100) NLigands(I,IJ),Z(NLigands(I,IJ)),(LB(IABS(MLligands(NLigands(I,IJ),NumbreFunctionLig(I,IJ,IXID)))+1),NnLigands(NLigands(I,IJ),NumbreFunctionLig(I,IJ,IXID)),LB1(Lligands(NLigands(I,IJ),NumbreFunctionLig(I,IJ,IXID))+1),IXID=1,NFunLigands(I,IJ))
	 ENDDO
	 WRITE(6,*)
  ENDDO
  
   !! ENTER THE PARAMETERS OF THE LIGANDA FUNCTIONS
   !! NumbreGarmLigand (Nnuclei, NFunctionLigand (IIZ)) - HARMONIC NUMBER OF MASSES IN THE FUNCTION OF LIGANDA
   !! LgarmonicLigand (2, Nnuclei, NFunctionLigand (IIZ)) - MASSIVE OF ORBITAL MOMENTS VALUES OF LIGAND FUNCTION
   !! LgarmonicLigand (1, Nnuclei, NFunctionLigand (IIZ)) - MINIMUM VALUE
   !! LgarmonicLigand (2, Nnuclei, NFunctionLigand (IIZ)) - MAXIMUM VALUE
   !! CYCLE ON LIGANDS
  DO I=1,Nnuclei
     DO IDC=1,NFunctionLigand(I)  
        READ(3,*) NumbreGarmLigand(I,IDC),LgarmonicLigand(1,I,IDC),LgarmonicLigand(2,I,IDC)
     ENDDO
  ENDDO 

!! DETERMINATE THE MAXIMUM NUMBER OF HARMONICS
  NumbreGarmLigandMAX=NumbreGarmLigand(1,1)
  DO I=1,Nnuclei
     DO IDC=1,NFunctionLigand(I)  
        IF(NumbreGarmLigand(I,IDC).GT.NumbreGarmLigandMAX) THEN
           NumbreGarmLigandMAX=NumbreGarmLigand(I,IDC) 
		ENDIF
     ENDDO
  ENDDO 
  
  
  WRITE(6,2400)
  DO I=1,Nnuclei
     DO IXID=1,NFunctionLigand(I)   
       WRITE(6,2500) I,Z(I),LB(IABS(MLligands(I,IXID))+1),NnLigands(I,IXID),LB1(Lligands(I,IXID)+1),NumbreGarmLigand(I,IXID),LgarmonicLigand(1,I,IXID),LgarmonicLigand(2,I,IXID)
     ENDDO
  ENDDO
  


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! INTRODUCE THE NUMBER OF TYPES OF ELECTRON INTERACTION COEFFICIENTS WITH THE MOLECULES OF THE MOLECULES !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(NtypIndexCoffNucleusElectron(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NtypIndexCoffNucleusElectron" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 

  !! Let us introduce the number of types of electron interaction coefficients with nuclei of a molecule
  READ(3,*) (NtypIndexCoffNucleusElectron(I),I=1,NumbreMO)

 !! We determine the maximum number of types of electron interaction coefficients with nuclei of a molecule
  NtypIndexCoffNucleusElectronMax=0
  DO I=1,NumbreMO
     IF(NtypIndexCoffNucleusElectronMax.LT.NtypIndexCoffNucleusElectron(I)) THEN
        NtypIndexCoffNucleusElectronMax=NtypIndexCoffNucleusElectron(I)
	 ENDIF
  ENDDO
  
  allocate(MLNucleusElectron(2,NtypIndexCoffNucleusElectronMax,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "MLNucleusElectron" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 
  allocate(RcoffNucleusElectron(NtypIndexCoffNucleusElectronMax,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "RcoffNucleusElectron" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 

   !! INTRODUCE THE PROJECTIONS OF MOLECULAR ORBITALS AND COEFFICIENTS OF INTERACTION
   !! CYCLE ON MOLECULAR ORBITALS
  DO I=1,NumbreMO
     !! CYCLE FOR TYPES OF COEFFICIENTS
     DO J=1,NtypIndexCoffNucleusElectron(I)
        READ(3,*) MLNucleusElectron(1,J,I),MLNucleusElectron(2,J,I),RcoffNucleusElectron(J,I)
     ENDDO
  ENDDO

  WRITE(6,*)  
  WRITE(6,3000)
  WRITE(6,3100)
  DO I=1,NumbreMO
     !! CYCLE FOR TYPES OF COEFFICIENTS
     DO J=1,NtypIndexCoffNucleusElectron(I)
        WRITE(6,3200) NN(I),LB(IABS(ML(I))+1),MLNucleusElectron(1,J,I),MLNucleusElectron(2,J,I),RcoffNucleusElectron(J,I)
     ENDDO
  ENDDO 
  
  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!
   !! ENTER THE NUMBER OF TYPES OF COEFFICIENTS OF CULON INTERACTION!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!
  allocate(NtypIndexCoffF(NumbreMO,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NtypIndexCoffF" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 
  allocate(NtypIndexCoffG(NumbreMO,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "NtypIndexCoffG" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif

  !! ENTER THE NUMBER OF TYPES OF DIRECT CUULON INTEGRALS
  DO I=1,NumbreMO
     READ(3,*) (NtypIndexCoffF(I,J),J=1,NumbreMO)
  ENDDO 
  !! Let us introduce the number of types of exchangeable Coulomb integrals
  DO I=1,NumbreMO
     READ(3,*) (NtypIndexCoffG(I,J),J=1,NumbreMO)
  ENDDO 

  !! We determine the maximum number of types of direct Coulomb integrals
  NtypIndexCoffFmax=0
  NtypIndexCoffGmax=0
  DO I=1,NumbreMO
     DO J=1,NumbreMO
        IF(NtypIndexCoffF(I,J).GT.NtypIndexCoffFmax) THEN 
		   NtypIndexCoffFmax=NtypIndexCoffF(I,J) 
        ENDIF
		IF(NtypIndexCoffG(I,J).GT.NtypIndexCoffGmax) THEN 
		   NtypIndexCoffGmax=NtypIndexCoffG(I,J) 
        ENDIF
	 ENDDO
  ENDDO 
  
 
  ! INTRODUCE COEFFICIENTS OF CULON INTERACTION
  allocate(RcoffAF(NtypIndexCoffFmax,NumbreMO,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "RcoffAF" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif 
  allocate(RcoffBG(NtypIndexCoffGmax,NumbreMO,NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "RcoffBG" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(MassivMLaf(NtypIndexCoffFmax,NumbreMO,NumbreMO,4),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "MassivMLaf" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(MassivMLbg(NtypIndexCoffGmax,NumbreMO,NumbreMO,4),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "MassivMLbg" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif  
  
   !! INTRODUCTION OF DIRECT CULON INTERACTION COEFFICIENTS
   !! MassivMLaf (NtypIndexCoffF, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS WHEN CALCULATING DIRECT INTERACTION
   !! RcoffAF (NtypIndexCoffF, NumbreMO, NumbreMO) -MASIVE OF DIRECT INTERACTION FACTORS
  DO I=1,NumbreMO
     DO J=1,NumbreMO
	    IF(NtypIndexCoffF(I,J).NE.0) THEN
           DO II=1,NtypIndexCoffF(I,J)
              READ(3,*) MassivMLaf(II,I,J,1),MassivMLaf(II,I,J,2),MassivMLaf(II,I,J,3),MassivMLaf(II,I,J,4),RcoffAF(II,I,J) 
		   ENDDO
	    ENDIF
	 ENDDO
  ENDDO

 
   !! INTRODUCTION OF EXCHANGE INTERACTION COEFFICIENTS
   !! MassivMLbg (NtypIndexCoffG, NumbreMO, NumbreMO, 4) -MASSIVE OF Moment (MOLECULAR ORBITAL) PROJECTIONS WHEN CALCULATING EXCHANGE INTERACTION
   !! RcoffBG (NtypIndexCoffG, NumbreMO, NumbreMO) -MASIVE OF EXCHANGE INTERACTION COEFFICIENTS
  DO I=1,NumbreMO
     DO J=1,NumbreMO
	    IF(NtypIndexCoffG(I,J).NE.0) THEN
           DO II=1,NtypIndexCoffG(I,J)
              READ(3,*) MassivMLbg(II,I,J,1),MassivMLbg(II,I,J,2),MassivMLbg(II,I,J,3),MassivMLbg(II,I,J,4),RcoffBG(II,I,J) 
		   ENDDO
	    ENDIF
	 ENDDO
  ENDDO
  
  !! ISSUE OF COEFFICIENTS OF CULON INTERACTION! We set the number of direct and exchangeable coefficients
  NFcoff=0
  DO I=1,NumbreMO
     DO J=1,NumbreMO
        DO II=1,NtypIndexCoffF(I,J)
		    IF(RcoffAF(II,I,J).NE.0.D0) THEN
			   NFcoff=NFcoff+1
			ENDIF 
	    ENDDO
	 ENDDO 
  ENDDO
  NGcoff=0
  DO I=1,NumbreMO
     DO J=1,NumbreMO
        DO II=1,NtypIndexCoffG(I,J)
		   IF(RcoffBG(II,I,J).NE.0.D0) THEN
		      NGcoff=NGcoff+1
		   ENDIF 
	    ENDDO
	 ENDDO 
  ENDDO           
  

  WRITE(6,1080)  NFcoff+NGcoff,NGcoff
  WRITE(6,1090)
  
  WRITE(6,2700)
  DO I=1,NumbreMO
     DO J=1,NumbreMO
        DO II=1,NtypIndexCoffF(I,J)
           WRITE(6,2900) NN(I),LB(IABS(ML(I))+1),NN(J),LB(IABS(ML(J))+1),MassivMLaf(II,I,J,1),MassivMLaf(II,I,J,2),MassivMLaf(II,I,J,3),MassivMLaf(II,I,J,4),RcoffAF(II,I,J) 
		ENDDO
     ENDDO
  ENDDO
 
  WRITE(6,*)
  WRITE(6,1100)
  WRITE(6,2800)
  DO I=1,NumbreMO
     DO J=1,NumbreMO
        DO II=1,NtypIndexCoffG(I,J)
           WRITE(6,2900) NN(I),LB(IABS(ML(I))+1),NN(J),LB(IABS(ML(J))+1),MassivMLbg(II,I,J,1),MassivMLbg(II,I,J,2),MassivMLbg(II,I,J,3),MassivMLbg(II,I,J,4),RcoffBG(II,I,J) 
		ENDDO
	 ENDDO
  ENDDO
   
  

  !! WE WRITE THE MODULE OF THE LOWER ENERGY BORDER
  !! MODULE OF THE ENERGY OF THE POINT OF THE POTENTIAL PIT (CONDITIONAL POINT OF REPORT) IN RIDBERGS (WITHOUT MINUS)
  READ(3,*) EzeroXF 

!! CHECK OUT OF OR FAILURE OR NOT
  IF(IndexHartreeFock.EQ.1) THEN
     IF(IndexExit.EQ.0) THEN
        !! We highlight the memory for the massifs
        allocate(NumbreGarmMOZero(NumbreMO),stat=ierr)
        if(ierr/=0) then
           write(*,*)'Molecule MEMORY ON THE FILE "NumbreGarmMOZero" IS NOT SELECTED'
           read(*,*)
	       stop 
        endif 
        allocate(NumbreGarmLMOZero(NumbreMO),stat=ierr)
        if(ierr/=0) then
           write(*,*)'Molecule MEMORY ON THE FILE "NumbreGarmLMOZero" IS NOT SELECTED'
           read(*,*)
	       stop 
        endif 
        allocate(LgarmonicMOZero(3,NumbreMO),stat=ierr)
        if(ierr/=0) then
           write(*,*)'Molecule MEMORY ON THE FILE "LgarmonicMOZero" IS NOT SELECTED'
           read(*,*)
	       stop 
        endif


        !! We introduce the parameters of the single-center decomposition of the Molecular Orbital of a non-zero approximation
        !! NumbreGarmMOZero (NumbreMO) - HARMONIC NUMBER OF MASSES IN THE MOLECULAR ORBITAL OF THE NULL APPROXIMATION
        !! LgarmonicMOZero (3, NumbreMO) -MASSIVE OF THE ORBITAL MOMENT VALUES OF MOLECULAR ORBITALS OF THE ZERO APPROXIMATION
        !! LgarmonicMOZero (1, NumbreMO) - MINIMUM VALUE
        !! LgarmonicMOZero (2, NumbreMO) -MAXIMUM VALUE
        !! LgarmonicMOZero (2, NumbreMO) -MAXIMUM VALUE
        !! LgarmonicMOZero (3, NumbreMO) - ORBITAL MOMENT IMPORTANCE AFTER WHICH CALCULATION IS DONE ON LEGAND FUNCTIONS (functions do not change when solving)       
        DO I=1,NumbreMO
           READ(3,*) NumbreGarmMOZero(I),LgarmonicMOZero(1,I),LgarmonicMOZero(2,I),LgarmonicMOZero(3,I)
        ENDDO 

        !! DEFINITION Lk CORRECT FOR THE IMPLEMENTATION OF THE ÌÎ
        IF(NfunLigMax.NE.0) THEN
           DO I=1,NumbreMO
              LgarmonicMOZero(3,I)=LgarmonicMOZero(3,I)+NfunLigMax-1
	       ENDDO
        ENDIF
        !! DETERMINATE THE NUMBER OF HARMONICS TO LgarmonicMOZero (3, I) INCLUSIVE
        DO I=1,NumbreMO
           NumbreGarmLMOZero(I)=LgarmonicMOZero(3,I)-LgarmonicMOZero(1,I)+1
        ENDDO
   
        !! DETERMINATE THE MAXIMUM NUMBER OF HARMONICS
        NgarmZeroMax=NumbreGarmLMOZero(1)
        DO I=1,NumbreMO
           IF(NgarmZeroMax.LT.NumbreGarmLMOZero(I)) THEN
              NgarmZeroMax=NumbreGarmLMOZero(I)   
	       ENDIF
        ENDDO

        !! We select the memory for the array
        allocate(IZONGarmonZero(NumbreMO,NgarmZeroMax),stat=ierr)
        if(ierr/=0) then
           write(*,*)'Molecule MEMORY ON THE FILE "IZONGarmonZero" IS NOT SELECTED'
           read(*,*)
	       stop 
        endif


        !! INTRODUCTION ZONES OF RECORDING PARTIAL HARMONICS OF THE NULL APPROACH
        !! IZONGarmonZero (NumbreMO, NumbreGarmLMOZero (NumbreMO)) - MASSIVE of ZONES OF RECORDING PARTIAL HARMONICS OF THE ZERO ACCELERATION TO LgarmonicMOZero (3, NumbreMO) INCLUSIVE
        !! IZONGarmonZero (NumbreMO, NumbreGarmLMOZero (NumbreMO)) = 0-MEANS THAT THE PARTIAL HARMONIC IS EQUAL TO ZERO
        DO I=1,NumbreMO
           READ(3,*) (IZONGarmonZero(I,J),J=1,NumbreGarmLMOZero(I))
        ENDDO 
    ENDIF
  ENDIF 

  WRITE(6,*)
  WRITE(6,1150) 
  WRITE(6,1160) EzeroXF
  !! CHECK OUT OF OR FAILURE OR NOT
  IF(IndexHartreeFock.EQ.1) THEN
     IF(IndexExit.EQ.0) THEN 
        WRITE(6,*)
        write(6,1170)
        DO I=1,NumbreMO
            WRITE(6,1070) NN(I),LB(IABS(ML(I))+1),IQ(I),NumbreGarmMOZero(I),LgarmonicMOZero(1,I),LgarmonicMOZero(3,I),LgarmonicMOZero(2,I)
        ENDDO
     ENDIF
  ENDIF
   
  !! We select the memory for the array
  allocate(IZONGarmonLigand(Nnuclei,NFunctionLigandMAX,NumbreGarmLigandMAX),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "IZONGarmonLigand" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif

  !! ENTER THE ZONES OF RECORDING THE FUNCTIONS OF LIGANDS
  !! IZONGarmonLigand (Nnuclei, NFunctionLigandMAX, NumbreGarmLigandMAX) -MASSIVE ZONES OF RECORDING PARTIAL HARMONICS OF LIGANDE FUNCTIONS
  DO I=1,Nnuclei
     DO IIX=1,NFunctionLigand(I)
        READ(3,*) (IZONGarmonLigand(I,IIX,IIO),IIO=1,NumbreGarmLigand(I,IIX))
	 ENDDO
  ENDDO

  
  !!!!!!!!!!! END OF READING THE PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BINDING FILE FUNCTIONS WITH FLOW
  IF(IZAP.NE.0) THEN 
     OPEN(UNIT=1,FILE=FILE1,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  ENDIF
  !! BINDING THE FILE FUNCTION FUNCTIONS OF A NULL-FREQUENCY APPROXIMATION WITH THE FLOW
  OPEN(UNIT=10,FILE=FILEZERO,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  !! BINDING THE FILE FUNCTION OF LIGANDS WITH FLOW
  OPEN(UNIT=11,FILE=FILELigands,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! R (Npoint) -MASSIVE VALUES OF RADIUS
  !! RO1X (Npoint) - MASSIVE VALUES OF THE FIRST DERIVATIVE ON THE NEW VARIABLE
  !! RO2X (Npoint) - MASSIVE VALUES OF THE SECOND DERIVATIVE ON THE NEW VARIABLE
  !! RO3X (Npoint) - MASSIVE VALUES OF THIRD DERIVATIVE ON THE NEW VARIABLE
  !! RFun (NumbreMO, NumbreGarmMO (NumbreMO), Npoint + 2) is the array of values of the radial parts of the harmonics of the molecular orbitals of the configuration  
  ! îïðåäåëÿåì ìàêñèìàëüíîå çíà÷åíèå ãàðìîíèê
  NmaxGarmonik=NumbreGarmMO(1)
      
  DO IXXX=1,NumbreMO
     IF(NumbreGarmMO(IXXX).GT.NmaxGarmonik) THEN
        NmaxGarmonik=NumbreGarmMO(IXXX)
	 ENDIF
  ENDDO
    
   
  ! WE ALLOW MEMORY MASSIVES 
  allocate(Rfun(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "Rfun" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(RfunZero(NumbreMO,NmaxGarmonik,N+2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "RfunZero" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
  allocate(RfunLigands(Nnuclei,NFunctionLigandMAX,NumbreGarmLigandMAX,N+2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'Molecule MEMORY ON THE FILE "RfunLigands" IS NOT SELECTED'
     read(*,*)
	 stop 
  endif
 
  
  allocate(ISCalculZZ(NumbreMO,NumbreMO+1),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule MEMORY ON THE FILE "ISCalculZZ" IS NOT SELECTED'
      read(*,*)
	 stop 
  endif 
  allocate(NumreISzam(NumbreMO),stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule MEMORY ON THE FILE "NumreISzam" IS NOT SELECTED'
      read(*,*)
	 stop 
  endif 
  
  
  !! WRITE THE MOLECULAR ORBITALS OF THE ZERO APPROXIMATION
  !! CHECK OUT OF OR FAILURE OR NOT
  IF(IndexHartreeFock.EQ.1) THEN
     IF(IndexExit.EQ.0) THEN
        call ReadDatIn(NumbreMO,N,NumbreGarmLMOZero,IZONGarmonZero,RfunZero)
     ENDIF
  ENDIF

  
  !! WRITE THE FUNCTIONS OF LIGANDS
  call ReadDatInLigand(Nnuclei,N,H,NFunctionLigand,NumbreGarmLigand,IZONGarmonLigand,RfunLigands,RO1X)
  
  !! ZANULYAYEM BEFORE CALCULATION OF MASSIVE FUNCTIONS
  Rfun=0.D0 

  
  !! ZANULYAYEM BEFORE PAYMENT
  ISCalculZZ=0 
  NumreISzam=0
  !! FORM THE ARRAY IN WHICH WE INDICATE THE WAVE FUNCTIONS OF WHICH SHELLS WE RECEIVE THE SOLUTION OF THE SYSTEM OF DIFFERENTIAL EQUATIONS,
  !! AND WHAT FUNCTIONS RECORD FROM RECEIVED
  DO I=1,NumbreMO
     !! CHECK NOT EXCLUDED THIS SHELL FROM THE CONSIDERATION
	 IF(ISCalculZZ(I,1).NE.1) THEN 
       !! EQUIVOLENT SHELL INDEX 
       NumbreISEQ=0
       DO J=1,NumbreMO  
	      !! EXCLUDING THIS SHELL I
		  IF(I.NE.J) THEN
		    !! CHECK NOT EXCLUDED THIS SHELL FROM THE CONSIDERATION
		     IF(ISCalculZZ(J,1).NE.1) THEN
                !! DETECT THE SAME SHELLS
		        IF(NN(J).EQ.NN(I).AND.IABS(ML(J)).EQ.IABS(ML(I))) THEN
                   NumbreISEQ=NumbreISEQ+1
                   !! WRITE THE SHELL WHICH IS EQUIVOLENT AND WHICH DOES NOT NEED TO DECIDE THE EQUATION
			       ISCalculZZ(I,NumbreISEQ+1)=J
			       !! WE INDICATE THAT THIS SHELL IS NOT PARTICIPATED IN DIRECT CALCULATION
                   ISCalculZZ(J,1)=1 
		        ENDIF
		     ENDIF
		  ENDIF
	   ENDDO 
	   !! WE WRITE THE NUMBER OF EXCLUDED SHELLS
	   NumreISzam(I)=NumbreISEQ 
     ENDIF
   ENDDO

   
  
  !! Subprogram of the solution of the system of Hartree-Fock equations for this congruence
  call DHFMO(IndexHartreeFock,IndexExit,Nnuclei,Z0,Z,CoorR,CoorA,NumbreMO,NtypIndexCoffNucleusElectron,NtypIndexCoffF,NtypIndexCoffG,IreshimSO,RcoffNucleusElectron,RcoffAF,RcoffBG,MLNucleusElectron,MassivMLaf,MassivMLbg,NN,ML,IQ,NumbreGarmMO,NumbreGarmLMO,LgarmonicMO,EPS,EPSXR,H,N,R,RO1X,RO2X,RO3X,RFun,NumbreGarmMOZero,LgarmonicMOZero,RfunZero,EzeroXF,NumreISzam,ISCalculZZ,FUNBuffer,FUNBuffer1,FUNBuffer2,FUNBuffer3,FUNBuffer4,IKLZERO,NFunctionLigand,NumbreLigand,NLigands,NFunLigands,NumbreFunctionLig,NumbreGarmLigand,LgarmonicLigand,RfunLigands,NomeroFPOT,FileParameterCorrection) 
    
 
   
  !! SUBPROGRAMME OF RECORDING RESULTS
  call WriteDatOut(IZAP,NumbreMO,N,H,NN,ML,IQ,NumbreGarmMO,R,RO1X,Rfun,OUTFUN)
  !-----------------------------------------------------
 
 


  
  !! REMOVING FILM MASSIVES
  deallocate(MLNucleusElectron,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "MLNucleusElectron" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(RcoffNucleusElectron,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "RcoffNucleusElectron" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif   
  deallocate(NtypIndexCoffNucleusElectron,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NtypIndexCoffNucleusElectron" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(NtypIndexCoffF,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NtypIndexCoffF" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(NtypIndexCoffG,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NtypIndexCoffG" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif   
  deallocate(RfunLigands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "RfunLigands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(IZONGarmonLigand,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "IZONGarmonLigand" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif   
  deallocate(LgarmonicLigand,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "LgarmonicLigand" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif   
  deallocate(NumbreGarmLigand,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NumbreGarmLigand" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(NumbreFunctionLig,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NumbreFunctionLig" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(NFunLigands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NFunLigands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(NumbreLigand,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NumbreLigand" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(Lligands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "Lligands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(NnLigands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NnLigands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(NLigands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NLigands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(MLligands,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "MLligands" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(NFunctionLigand,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NFunctionLigand" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(NumreISzam,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "NumreISzam" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(ISCalculZZ,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "ISCalculZZ" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  ! ÏÐÎÂÅÐßÅÌ ÏÐÎÈÇÎØÅË ÑÁÎÉ ÈËÈ ÍÅÒ
  IF(IndexHartreeFock.EQ.1) THEN
     IF(IndexExit.EQ.0) THEN 
        deallocate(NumbreGarmMOZero,stat=ierr)
        if(ierr/=0) then
           write(*,*) 'Molecule'
           write(*,*) 'THE MASSIV "NumbreGarmMOZero" IS NOT REMOVED FROM MEMORY'
           read(*,*)
	       stop 
        endif 
        deallocate(NumbreGarmLMOZero,stat=ierr)
        if(ierr/=0) then
           write(*,*) 'Molecule'
           write(*,*) 'THE MASSIV "NumbreGarmLMOZero" IS NOT REMOVED FROM MEMORY'
           read(*,*)
	       stop 
        endif 
        deallocate(LgarmonicMOZero,stat=ierr)
        if(ierr/=0) then
           write(*,*) 'Molecule'
           write(*,*) 'THE FILE "LgarmonicMOZero" IS NOT REMOVED FROM MEMORY'
           read(*,*)
	       stop 
        endif
        deallocate(IZONGarmonZero,stat=ierr)
        if(ierr/=0) then
           write(*,*) 'Molecule'
           write(*,*) 'THE FILE "IZONGarmonZero" IS NOT REMOVED FROM MEMORY'
           read(*,*)
	       stop 
        endif
     ENDIF
  ENDIF
  deallocate(Rfun,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "Rfun" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(RfunZero,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE MASSIV "RfunZero" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif  
  deallocate(R,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "R" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(RO1X,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "RO1X" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(RO2X,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "RO2X" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(RO3X,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "RO3X" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(NN,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "NN" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(ML,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "ML" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(IQ,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "IQ" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(NumbreGarmMO,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "NumbreGarmMO" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(NumbreGarmLMO,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "NumbreGarmLMO" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(MassivMLaf,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "MassivMLaf" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(MassivMLbg,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "MassivMLbg" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(LgarmonicMO,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "LgarmonicMO" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(Z,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "Z" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(CoorR,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "CoorR" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  deallocate(CoorA,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "CoorA" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif

  deallocate(RcoffAF,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "RcoffAF" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif 
  deallocate(RcoffBG,stat=ierr)
  if(ierr/=0) then
     write(*,*) 'Molecule'
     write(*,*) 'THE FILE "RcoffBG" IS NOT REMOVED FROM MEMORY'
     read(*,*)
	 stop 
  endif
  
  GOTO 77

  return
end subroutine Molecule


 
