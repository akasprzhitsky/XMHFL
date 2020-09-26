      
	
  
program ExpansionFunSH
	  
 integer(2)  Hz1,Mz1,Sz1,HSz1,Hz2,Mz2,Sz2,HSz2
   call GETTIM(Hz1,Mz1,Sz1,HSz1)
   TB=60.*Hz1+Mz1+Sz1/60. 
   ! DATA ENTRY SUB-PROGRAM
   call readfile
   ! MATRIX ELEMENT CALCULATION SUBPROGRAM
   call CalculExpansionFunctionSphericalHarmonics
   ! RESULTS SUB-PROGRAM
   call writefile(6)
   call GETTIM(Hz2,Mz2,Sz2,HSz2)
   TE=60.*Hz2+Mz2+Sz2/60.
   T=TE-TB
   WRITE(6,*)
   WRITE(6,1) T
 1 FORMAT(' Commercial    time :',F10.2,' min.')

   STOP
end program ExpansionFunSH

       
! DATA ENTRY SUB-PROGRAM
subroutine readfile
  implicit real(8) (A-H,O,P,R-Z)
  
  common /C1/IS,M,M2,Mnew,Lmin,Lmax,NtipSetkyatom,NtipSetky,NKie,IZFUN
  common /C2/Z,H,BET,ALFA,RalFL,GAMMA,Ral,BETnew,ALFAnew,GAMMAnew,Hnew
  common /C3/Nn(10),L(10),LML(10)
  common /C4/IZON(10)
  common /CTX/OUTF(10)     
  character(90) TITLE,DATF,OUTFF,OUTF,FILE,OUTF0,OUTENERGY,OutEnergyTab,FILEfun

    
  WRITE(*,'(A\)') ' ‚‚…„ˆ’… TITLE- ”€‰‹ ˆ‹ˆ "*"-->'
  READ(*,'(A90)') TITLE
      
  IF(TITLE.EQ.'*') THEN
     TITLE='TITLE.DAT'
  ENDIF
	
  OPEN(4,FILE=TITLE)
      
  READ(4,'(A90)') DATF
  READ(4,'(A90)') OUTFF
  READ(4,'(A90)') FILE
  READ(4,'(A90)') OUTF0
  READ(4,'(A90)') OUTENERGY
  READ(4,'(A90)') OutEnergyTab
  READ(4,'(A90)') FILEfun
      
  OPEN(5,FILE=DATF)
  OPEN(6,FILE=OUTFF)
  OPEN(7,FILE=OUTF0)
  OPEN(8,FILE=OUTENERGY)
  OPEN(18,FILE=OutEnergyTab)

    
  WRITE(*,*) TITLE
  WRITE(*,*) DATF
  WRITE(*,*) OUTFF
  WRITE(*,*) FILE
  WRITE(*,*) OUTF0
  WRITE(*,*) OUTENERGY
  WRITE(*,*) OutEnergyTab 
  ! Numerical data input block
  ! Ral-distance from the atomic coordinate system to the decomposition center coordinate system
  ! Ral> 0 - displacement occurs along the OZ axis (IN POSITIVE DIRECTION)
  ! Ral <0 - offset occurs against the OZ axis (IN NEGATIVE DIRECTION)
  ! Lmax-orbital moment of the last term in which the expansion is performed
  ! NtipSetkyatom-mesh type for expandable function
  ! NtipSetky-type mesh
  ! NtipSetky = 1-ATOMIC GRID
  ! NtipSetky = 2-MOLECULAR GRID (GRID FOR THE CASE OF NON-COINCIDENCE OF THE BEGINNING OF THE SYSTEM OF COORDINATES WITH THE ATOM)
  ! NKie-CONDITION FOR RECORDING DECOMPOSITION COEFFICIENTS
  ! NKie = 0-DECOMPOSITION RATES ARE NOT WRITTEN TO FILE
  ! NKie = 1-DECOMPOSITION COEFFICIENTS WRITE TO FILE
  ! IZFUN-NUMBER OF THE FIRST RECORD INTO FILE DECOMPOSITION COEFFICIENT
  READ(5,*) Z,IS,NtipSetkyatom,NtipSetky,NKie

  IF(NtipSetkyatom.EQ.1) THEN
     READ(5,*) ALFA,BET,H,M
  ENDIF
  IF(NtipSetkyatom.EQ.2) THEN
     READ(5,*) RalFL,ALFA,BET,GAMMA,H,M
  ENDIF
	
  IF(NtipSetky.EQ.1) THEN
    READ(5,*) Ral,ALFAnew,BETnew,Hnew,Mnew,Lmin,Lmax
  ENDIF
  IF(NtipSetky.EQ.2) THEN
    READ(5,*) Ral,ALFAnew,BETnew,GAMMAnew,Hnew,Mnew,Lmin,Lmax
  ENDIF

  READ(5,*) (Nn(K),L(K),LML(K),K=1,IS)
  READ(5,*) (IZON(K),K=1,IS)
  IF(NKie.EQ.1) THEN
	READ(5,*) IZFUN 
  ENDIF

  M2=M+2
    
  
      
 ! WE READ THE NAMES OF THE FILES IN WHICH THE EXPANSIONS WILL BE WRITTEN
  DO K=1,IS 
     READ(4,'(A90)') OUTF(K)
     WRITE(*,*) OUTF(K)
  ENDDO
	
 ! LINKING FILES TO CHANNELS
  DO K=1,IS 
     OPEN(100+K,FILE=OUTF(K))
  ENDDO
  ! LINKING THE FILE TO THE CHANNEL (DECOMPOSABLE FUNCTIONS)
  OPEN(UNIT=1,FILE=FILE,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  
  IF(NKie.EQ.1) THEN 
    ! LINKING THE FILE TO THE CHANNEL (DECOMPOSITION RATES)
     OPEN(UNIT=2,FILE=FILEfun,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  ENDIF
             
  return
end subroutine readfile 
     

! SUBPROGRAM FOR CALCULATING THE EXPANSION COEFFICIENTS
subroutine  CalculExpansionFunctionSphericalHarmonics
  use mefsh,only:EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS,EFSH_RW
  use mcmri,only:CMRI_VAR,CMRI_ROFUN,CMRI_INT_ORT
  use mcmme,only:CMME_ONE_ELECTRONIC_INTEGRAL,CMME_TWO_ELECTRONIC_INTEGRAL
  implicit real(8)(A-H,O,P,R-Z)
  
  common /C1/IS,M,M2,Mnew,Lmin,Lmax,NtipSetkyatom,NtipSetky,NKie,IZFUN
  common /C2/Z,H,BET,ALFA,RalFL,GAMMA,Ral,BETnew,ALFAnew,GAMMAnew,Hnew
  common /C3/Nn(10),L(10),LML(10)
  common /C4/IZON(10)
      
  real(4),allocatable,dimension(:)::REW
  integer,allocatable,dimension(:)::NRabParametrs
  real(8),allocatable,dimension(:)::Rtemp,DM,R,Rnew
  real(8),allocatable,dimension(:)::RFUN,Acoff,RFUNAPROC,RO1,RO2,RO3
  real(8),allocatable,dimension(:)::RFUNN,RFUNN1,RFUNN2,RO1N,RO2N,RO3N
  real(8),allocatable,dimension(:)::RORTGarmon
  real(8),allocatable,dimension(:,:)::REnergyGarmon
  real(8),allocatable,dimension(:,:)::REnergyIntGarmon
  real(8),allocatable,dimension(:,:)::RKinEnergyGarmon  
  real(8),allocatable,dimension(:,:)::DMApro,RcoffSHE
  real(8),allocatable,dimension(:,:,:)::RCoffSH
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  dimension LB(10),ZZX(1),RCOR(1,2),ACOR(1,2)
  character(1)  LB
  data LB/'s','p','d','f','g','h','i','j','k','l'/
 

  ! WE ALWAYS MEMORY FOR ARRAYS
  allocate(NRabParametrs(2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "NRabParametrs" IS NOT SELECTED'
	 stop 
  endif
 
  allocate(Rtemp(M2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "Rtemp" IS NOT SELECTED'
     stop 
  endif
  allocate(DM(M2*IS),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "DM" IS NOT SELECTED'
  	 stop 
  endif
  allocate(DMApro(IS,M2),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "DMApro" IS NOT SELECTED'
	 stop 
  endif
  allocate(R(M),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "R" IS NOT SELECTED'
	 stop 
  endif
  allocate(Rnew(Mnew),stat=ierr)
  if(ierr/=0) then
     write(*,*)'MEMORY ON THE FILE "Mnew" IS NOT SELECTED'
	 stop 
  endif
  IF(Mnew.GT.M) THEN 
	 allocate(REW(Mnew+2),stat=ierr)
     if(ierr/=0) then
        write(*,*)'MEMORY ON THE FILE "REW" IS NOT SELECTED'
	    stop 
     endif
        allocate(RFUN(Mnew+2),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUN" IS NOT SELECTED'
          stop 
        endif
	    allocate(RFUNN(Mnew),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN" IS NOT SELECTED'
          stop 
        endif
        allocate(RFUNN1(Mnew),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN1" IS NOT SELECTED'
          stop 
        endif
		allocate(RFUNN2(Mnew),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN2" IS NOT SELECTED'
          stop 
        endif
        allocate(RFUNAPROC(Mnew),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNAPROC" IS NOT SELECTED'
          stop 
        endif
       ELSE 
	    allocate(REW(M2),stat=ierr)
        if(ierr/=0) then
            write(*,*)'MEMORY ON THE FILE "REW" IS NOT SELECTED'
	        stop 
        endif
        allocate(RFUN(M+2),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUN" IS NOT SELECTED'
          stop 
        endif
	    allocate(RFUNN(M),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN" IS NOT SELECTED'
          stop 
        endif
		allocate(RFUNN1(M),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN1" IS NOT SELECTED'
          stop 
        endif
		allocate(RFUNN2(M),stat=ierr)
        if(ierr/=0) then
          write(*,*)'MEMORY ON THE FILE "RFUNN2" IS NOT SELECTED'
          stop 
        endif
        allocate(RFUNAPROC(M),stat=ierr)
        if(ierr/=0) then
           write(*,*)'MEMORY ON THE FILE "RFUNAPROC" IS NOT SELECTED'
           stop 
        endif
      ENDIF 
      allocate(Acoff(Mnew),stat=ierr)
      if(ierr/=0) then
        write(*,*)'MEMORY ON THE FILE "Acoff" IS NOT SELECTED'
        stop 
      endif  
	  allocate(RCoffSHE(Lmax-Lmin+1,Mnew),stat=ierr)
      if(ierr/=0) then
         write(*,*)'MEMORY ON THE FILE "RCoffSHE" IS NOT SELECTED'
         stop 
      endif
      allocate(RCoffSH(IS,Lmax-Lmin+1,Mnew),stat=ierr)
      if(ierr/=0) then
         write(*,*)'MEMORY ON THE FILE "RCoffSH" IS NOT SELECTED'
         stop 
      endif
      allocate(RORTGarmon(Lmax-Lmin+1),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RORTGarmon" IS NOT SELECTED'
	stop 
	endif    
	allocate(REnergyGarmon(Lmax-Lmin+1,Lmax-Lmin+1),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "REnergyGarmon" IS NOT SELECTED'
	stop 
	endif
	allocate(REnergyIntGarmon(Lmax-Lmin+1,Lmax-Lmin+1),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "REnergyIntGarmon" IS NOT SELECTED'
	stop 
	endif
	allocate(RKinEnergyGarmon(Lmax-Lmin+1,Lmax-Lmin+1),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RKinEnergyGarmon" IS NOT SELECTED'
	stop 
	endif
    allocate(RO1(M),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO1" IS NOT SELECTED'
	stop 
	endif
	allocate(RO2(M),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO2" IS NOT SELECTED'
	stop 
	endif
	allocate(RO3(M),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO3" IS NOT SELECTED'
	stop 
	endif
      allocate(RO1N(Mnew),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO1N" IS NOT SELECTED'
	stop 
	endif
	allocate(RO2N(Mnew),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO2N" IS NOT SELECTED'
	stop 
	endif
	allocate(RO3N(Mnew),stat=ierr)
	if(ierr/=0) then
      write(*,*)'MEMORY ON THE FILE "RO3N" IS NOT SELECTED'
	stop 
	endif

  

      
	! ZERO BEFORE WORK
    REW=0.D0
	Rtemp=0.D0
	DM=0.D0
	R=0.D0
	Rpar=0.D0
	Rnew=0.D0
	Rparnew=0.D0
    RFUN=0.D0
	Acoff=0.D0
	RFUNAPROC=0.D0
    DMApro=0.D0
    RCoffSH=0.D0
	RcoffSHE=0.D0
    RO1=0.D0
	RO2=0.D0
	RO3=0.D0
	RO1N=0.D0
	RO2N=0.D0
	RO3N=0.D0
	RORTGarmon=0.D0
    REnergyGarmon=0.D0
	REnergyIntGarmon=0.D0
	RKinEnergyGarmon=0.D0
	!RCoulomb=0.D0  
    ! STEP 0 Reading data from file
      
	do  I=1,IS
	    ! ZERO FOR THIS FUNCTION
	    Rtemp=0.D0
          REWIND 1
          IZ=IZON(I)-1
         IF(IZ.LE.0) THEN 
		   READ (1) (REW(J1),J1=1,M2)
             do J1=1,M2
                Rtemp(J1)=DBLE(REW(J1))
			 enddo
		 ELSE
             do J=1,IZ
                READ(1)
             enddo
             READ (1) (REW(J1),J1=1,M2)
             do J1=1,M2
		      Rtemp(J1)=DBLE(REW(J1))
	          
		   enddo
	     ENDIF
	     call EFSH_RW(1,I,Rtemp,DM,M2)
      enddo
	
	
	CLOSE(1)



      ! STEP 1 DETERMINATION OF THE VALUES OF THE RADIUS AND THE CONVERSION FACTOR
      ! FOR ATOMIC FUNCTIONS AND DECOMPOSITION COEFFICIENTS
      call CMRI_VAR(NtipSetkyatom,M,H,ALFA,BET,GAMMA,DABS(RalFL),R,RO1,RO2,RO3)
      call CMRI_VAR(NtipSetky,Mnew,Hnew,ALFAnew,BETnew,GAMMAnew,DABS(Ral),Rnew,RO1N,RO2N,RO3N)
      ! DETERMINING THE LIGHT BOUNDARY OF THE INTERVAL BY THE TYPE OF THE NET
	  IF(NtipSetkyatom.EQ.1) THEN
          RO0=-14.D0
      ENDIF
      IF(NtipSetkyatom.EQ.2) THEN
          RO0=-15.D0
      ENDIF
	  IF(NtipSetky.EQ.1) THEN
          RO0new=-14.D0
      ENDIF
      IF(NtipSetky.EQ.2) THEN
          RO0new=-15.D0
      ENDIF
      
	  ! ZERO BEFORE CALCULATION
	  NRabParametrs=0

	
     
    ! STEP 2 CALCULATION OF THE EXPANSION COEFFICIENTS INTO SPHERICAL HARMONICS
    ! CONFIGURATION SHELL CYCLE
	DO IJ=1,IS
	   ! STEP 0. RECORDING THE RADIAL PART OF THE WAVE FUNCTION OF THE IJ SHELL
	   call EFSH_RW(2,IJ,Rtemp,DM,M2)
	   DO IZX=1,M
	      RFUN(IZX)=Rtemp(IZX+2)/DSQRT(RO1(IZX))
	   ENDDO
      ! STEP 1. SUBPROGRAM OF THE EXPANSION OF A FUNCTION INTO A SERIES IN SPHERICAL HARMONICS
	   call EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS(Nn(IJ),L(IJ),LML(IJ),M,R,RFUN,RFUNAPROC,Lmin,Lmax,Ral,Mnew,Rnew,RcoffSHE,NRabParametrs)   
       ! STEP 3. WRITING THE VALUES OF THE APPROXIMATED RADIAL PART OF THE WAVE FUNCTION
	   DO IZX=1,M
	      DMApro(IJ,IZX)=RFUNAPROC(IZX)
       ENDDO
       ! STEP 4. WRITE DOWN THE OBTAINED COEFFICIENT
	   IndexGarmons=0
	   DO ILS=Lmin,Lmax  
	      IndexGarmons=IndexGarmons+1
		  DO IZXC=1,Mnew
		     ! ZERO AT LARGE VALUES TO ELIMINATE THE CONSUMPTION
	         IF(Rnew(IZXC).LT.35.D0) THEN
		          RCoffSH(IJ,IndexGarmons,IZXC)=RCoffSHE(IndexGarmons,IZXC)
		        ELSE
                  RCoffSH(IJ,IndexGarmons,IZXC)=0.D0
			 ENDIF
		  ENDDO  
       ENDDO
     ENDDO


      ! REMOVING ARRAYS FROM MEMORY 
	 deallocate(NRabParametrs,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'CalculExpansionFunctionSphericalHarmonics'
      write(*,*) 'THE FILE "NRabParametrs" IS NOT REMOVED FROM MEMORY'
	stop 
	endif




   
  50  FORMAT(2X,'START',I4,1X,'Npoints= ',I5)
 100  FORMAT(4X,I3,A1,'(m=',I3,')') 
 200  FORMAT(6X,'Distance',9X,200('Y(l=',I4,')',7X))
 300  FORMAT(8X,'Distance',9X,10(I3,A1,12X))
 400  FORMAT(4X,E14.6,2X,10(E14.6,2X))
 500  FORMAT(4X,E14.6,2X,200(E14.6,2X))   
 3000 FORMAT(2X,'EnergyXF('I3,A1,'(m= ',I2,'))= ',F15.10,1X,F20.10,1X,F20.10 ,' Int_Ort= ',F15.5)   
 3100 FORMAT(2X,'  Energy(L=',I3 ',ML=',I2,')= ',F15.10,1X,F20.10,1X,F20.10,' Int_Ort= ',F15.13)
 3200 FORMAT(2X,I3,F20.10,2X,F15.13)
 3500 FORMAT(2X,200(E14.6,2X))         
 4000 FORMAT(2X,'RCulombXF('I3,A1,I3,A1,')= ',F20.10)
 4100 FORMAT(2X,'RCulomb(L=',I3 ',ML=',I2,' ',I3,A1,I3,A1')= ',F20.10)
 4200 FORMAT(2X,I3,2X,F20.10,2X,F20.10,2X,F20.10)  
      ! RECORDING THE INITIAL FUNCTION
      WRITE(7,300) (Nn(I),LB(L(I)+1),Nn(I),LB(L(I)+1),I=1,IS)
	DO IZX=1,M
      WRITE(7,400)  R(IZX),(DM((NF-1)*M2+IZX+2)/DSQRT(RO1(IZX)),DMApro(NF,IZX),NF=1,IS)   
	ENDDO 
	 
	 ! RECORDING THE DECOMPOSITION COEFFICIENTS
     	DO IC=1,IS 
         WRITE(100+IC,100) Nn(IC),LB(L(IC)+1),LML(IC)
	   WRITE(100+IC,200) (ISDA,ISDA=Lmin,Lmax)
     ! EXPANSION COEFFICIENTS CYCLE
	   DO IDFX=1,Mnew
      WRITE(100+IC,500) Rnew(IDFX),(RCoffSH(IC,ILS-Lmin+1,IDFX)*Rnew(IDFX),ILS=Lmin,Lmax)
	   ENDDO
      ENDDO

     
	
    
	

      ! CALCULATION OF SINGLE-ELECTRONIC ENERGY (FOR ESTABLISHING CONVERGENCE)
	DO IJ=1,IS
      ! ZERO BEFORE CALCULATION
       RORTGarmon=0.D0
       REnergyGarmon=0.D0
	   
	  ! RECORDING THE RADIAL PART OF THE WAVE FUNCTION OF THE IJ SHELL
	   call EFSH_RW(2,IJ,Rtemp,DM,M2)
	   DO IZX=1,M
          RFUN(IZX)=Rtemp(IZX+2)
	   ENDDO
      ! CALCULATING THE ONE-ELECTRONIC ENERGY OF THE DECOMPOSABLE FUNCTION
	   ZZX=0.D0
	   RCOR=0.D0
	   ACOR=0.D0
	   EnergyXFnew=CMME_ONE_ELECTRONIC_INTEGRAL(L(IJ),LML(IJ),L(IJ),LML(IJ),0,Z,ZZX,RCOR,ACOR,-14.D0,H,M,R,RFUN,RFUN,RO1,RO2,RO3,EkinXF,EpotXF) 
	   ! WE CALCULATE THE INTEGRAL OF ORTHOGONALITY
	   call EFSH_RW(2,IJ,Rtemp,DM,M2)
	   DO IZX=1,M
          RFUN(IZX)=Rtemp(IZX+2)
	   ENDDO
       RORTXFnew=CMRI_INT_ORT(M,RO0,H,RFUN,RFUN,RO1)
	   !WRITE(8,3000) Nn(IJ),LB(L(IJ)+1),LML(IJ),EnergyXF,RORTXF
       WRITE(8,*)
	   WRITE(8,3000) Nn(IJ),LB(L(IJ)+1),LML(IJ),EnergyXFnew,EkinXF,EpotXF,RORTXFnew
       WRITE(18,4200) 0,EnergyXFnew,0.D0,0.D0
	   WRITE(18,*)
       ! THE ORBITAL DECOMPOSITION CYCLE
       ! WRITE DATA ABOUT THE CONFIGURATION OF NUCLEARS
	   ZZX(1)=Z
	   ! RADIAL COORDINATES
	   RCOR(1,1)=DABS(Ral)
	   RCOR(1,2)=CMRI_ROFUN(ALFAnew,BETnew,GAMMAnew,DABS(Ral),DABS(Ral),0.D0)
       ! CORNER COORDINATES
       ACOR(1,2)=0.D0
	   IF(Ral.GT.0.D0) THEN
          ACOR(1,1)=0.D0
	     ELSE
          ACOR(1,2)=3.14159265358979D0
	   ENDIF
	   
	   Ngarmonik1=0
       DO ILS=Lmin,Lmax
	      ! harmonic index
	      Ngarmonik1=Ngarmonik1+1
	      ! WRITE FACTOR AS F-type
	      DO IZXC=1,Mnew
	         RRRFGH=Rnew(IZXC)*DSQRT(RO1N(IZXC))
	         RFUN(IZXC)=RCoffSH(IJ,Ngarmonik1,IZXC)*RRRFGH
    	  ENDDO
          Ngarmonik2=0
	      DO ILS1=Lmin,Lmax 
             Ngarmonik2=Ngarmonik2+1
	         ! WRITE FACTOR AS F-type
	         DO IZXC=1,Mnew
	            RRRFGH=Rnew(IZXC)*DSQRT(RO1N(IZXC))
	            RFUNN(IZXC)=RCoffSH(IJ,Ngarmonik2,IZXC)*RRRFGH
    	     ENDDO
	         ! calculate the two-electron harmonic ILS 
	         REnergyGarmon(Ngarmonik1,Ngarmonik2)=CMME_ONE_ELECTRONIC_INTEGRAL(ILS,LML(IJ),ILS1,LML(IJ),1,0.D0,ZZX,RCOR,ACOR,RO0new,Hnew,Mnew,Rnew,RFUN,RFUNN,RO1N,RO2N,RO3N,RENERGYKIN,RENERGYINTZ) 
	         RKinEnergyGarmon(Ngarmonik1,Ngarmonik2)=RENERGYKIN
	         REnergyIntGarmon(Ngarmonik1,Ngarmonik2)=RENERGYINTZ
          ENDDO 
	      ! calculate the harmonic orthogonality integral
          RORTGarmon(Ngarmonik1)=CMRI_INT_ORT(Mnew,RO0new,Hnew,RFUN,RFUN,RO1N)
    	  !WRITE(8,*) 'RORT',RIOSDF,RORTGarmon(Ngarmonik1)
	   ENDDO
         
	   ! WRITE(8,*) 'ENERGY'
	   !  DO IIIX=1,Ngarmonik1
       !     WRITE(8,3500) (REnergyGarmon(IIIX,JJJJX),JJJJX=1,Ngarmonik2)
	   !  ENDDO
       !WRITE(8,*) 'EPOT'
       !DO IIIX=1,Ngarmonik1
       !   WRITE(8,3500) (REnergyIntGarmon(IIIX,JJJJX),JJJJX=1,Ngarmonik2)
       !ENDDO
	   ! ÐÀÑ×ÅÒ ÝÍÅÐÃÈÈ
	   LIIO=Lmin-1  
	   DO ILS=1,Ngarmonik1  
          LIIO=LIIO+1
		  ! CALCULATING THE RATING MULTIPLIER 
		  SUM=0.D0
	      DO IIIX=1,ILS
             SUM=SUM+RORTGarmon(IIIX)
	      ENDDO
	      RALFA=1.D0/SUM
          ! WE CALCULATE SINGLE ELECTRONIC ENERGY
	      SUMENERGY=0.D0
          SUMENERGYKIN=0.D0
          SUMENERGYPOT=0.D0
          DO JJJX=1,ILS
	         DO JJJY=1,ILS
                SUMENERGY=SUMENERGY+REnergyGarmon(JJJX,JJJY)
                SUMENERGYKIN=SUMENERGYKIN+RKinEnergyGarmon(JJJX,JJJY)
                SUMENERGYPOT=SUMENERGYPOT+REnergyIntGarmon(JJJX,JJJY)
	         ENDDO
		  ENDDO
	     ! WE WRITE THE RESULT OF THE CALCULATION
	      ENERGY=RALFA*SUMENERGY
		  EKIN=RALFA*SUMENERGYKIN
		  EPOT=RALFA*SUMENERGYPOT 
		  WRITE(8,3100) LIIO,LML(IJ),ENERGY,EKIN,EPOT,SUM
	      WRITE(18,4200) LIIO,ENERGY,EkinXF-EKIN,EPOT-EpotXF 
		  !WRITE(8,3200) LIIO,ENERGY,SUM
	   ENDDO
      ENDDO

    ! WE WRITE DECOMPOSITION COEFFICIENTS INTO A FILE FROM A RECORD IZFUN NUMBER
	IF(NKie.EQ.1) THEN
	   NfunSD=0
	   DO IJ=1,IS
	      IndexGarmons=0
	      DO ILS=Lmin,Lmax
	         IndexGarmons=IndexGarmons+1
	   	     NfunSD=NfunSD+1
		     REWIND 2
             IZ=(IZFUN+NfunSD-1)-1
             WRITE(6,*) 'ZONA',IJ,IZ+1
			 IF(IZ.LE.0) THEN
		         REW=0.D0
	             REW(2)=SNGL(RORTGarmon(IndexGarmons))
			     DO IZXC=1,Mnew
	                RRRFGH=Rnew(IZXC)*DSQRT(RO1N(IZXC))
	                REW(IZXC+2)=SNGL(RCoffSH(IJ,IndexGarmons,IZXC)*RRRFGH)
			     ENDDO
			     WRITE(2) (REW(J1),J1=1,Mnew+2)
				
               ELSE
              
			     DO J=1,IZ
                    READ(2)
                 ENDDO

                 REW=0.D0
			     REW(2)=SNGL(RORTGarmon(IndexGarmons))
			     DO IZXC=1,Mnew
	                RRRFGH=Rnew(IZXC)*DSQRT(RO1N(IZXC))
	                REW(IZXC+2)=SNGL(RCoffSH(IJ,IndexGarmons,IZXC)*RRRFGH)
    	            ! WRITE(6,*) REW(IZXC+2)
			     ENDDO
                 WRITE(2) (REW(J1),J1=1,Mnew+2)
	         ENDIF
	      ENDDO
	   ENDDO
  	   CLOSE(2)
	ENDIF

     

   
	! REMOVING ARRAYS FROM MEMORY 
	deallocate(RORTGarmon,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RORTGarmon" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(REnergyGarmon,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "REnergyGarmon" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(REnergyIntGarmon,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "REnergyIntGarmon" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RKinEnergyGarmon,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RKinEnergyGarmon" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
    deallocate(RO1,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO1" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RO2,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO2" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
   	deallocate(RO3,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO3" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(RO1N,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO1N" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RO2N,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO2N" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
   	deallocate(RO3N,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RO3N" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
   	deallocate(REW,stat=ierr)
	if(ierr/=0) then
	  write(*,*) 'THE FILE "DM" IS NOT REMOVED FROM MEMORY'
	stop 
	endif  
	deallocate(DM,stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'THE FILE "DM" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(DMApro,stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'THE FILE "DMApro" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Rtemp,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "Rtemp" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(R,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "R" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Rnew,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "Rnew" IS NOT REMOVED FROM MEMORY'
	stop 
	endif   
	deallocate(RFUNN,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RFUNN" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RFUNN1,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RFUNN1" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RFUNN2,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RFUNN2" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
    deallocate(RFUN,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RFUN" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RFUNAPROC,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RFUNAPROC" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Acoff,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "Acoff" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(RCoffSH,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RCoffSH" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RCoffSHE,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'THE FILE "RCoffSHE" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      return
	end subroutine  CalculExpansionFunctionSphericalHarmonics



     ! results output subroutine (report file)
	subroutine writefile(n)
	
	implicit real(8)(A-H,O,P,R-Z)
      
    common /C1/IS,M,M2,Mnew,Lmin,Lmax,NtipSetkyatom,NtipSetky,NKie,IZFUN
    common /C2/Z,H,BET,ALFA,RalFL,GAMMA,Ral,BETnew,ALFAnew,GAMMAnew,Hnew
    common /C3/Nn(10),L(10),LML(10)
    common /C4/IZON(10)
    common /CTX/OUTF(10)   
      
	dimension LB(10)
	character (90) OUTF
    character (1)  LB
    data LB/'s','p','d','f','g','h','i','j','k','l'/
	


700   FORMAT(2X,'ALFA=',F8.6,2X,'BETA =',F8.6,2X,'M =',I6,2X,'H =',F6.4,2X)
750   FORMAT(2X,'ALFA=',F8.6,2X,'BETA =',F8.6,2X,'GAMMA =',F8.6,2X,'M =',I6,2X,'H =',F6.4,2X,'Ral= ',F6.3)
800   FORMAT(2X,'ALFAnew=',F8.6,2X,'BETAnew =',F8.6,2X,'Mnew =',I6,2X,'Hnew =',F6.4,2X,'Ral= ',F6.3,' a.e',2X,'Lmin= ',I2,2X,'Lmax= ',I2)  
850   FORMAT(2X,'ALFAnew=',F8.6,2X,'BETAnew =',F8.6,2X,'GAMMAnew =',F8.6,2X,'Mnew =',I6,2X,'Hnew =',F6.4,2X,'Ral= ',F6.3,' a.e',2X,'Lmin= ',I3,2X,'Lmax= ',I3)    
900   FORMAT(/,'  Configuration :   ',20(I3,A1,'m= 'I2)/)
1000  FORMAT(2X,I3,A1,'(m= ',I2,') ',A90)
1100  FORMAT(2X,'MOLECULAR GRID')
1200  FORMAT(2X,'NUCLEAR GRID')


	WRITE(n,*) 'The report of the program ExpansionFunSH ver 1.4.0 21.06.'
    WRITE(n,*)
	IF(NtipSetkyatom.EQ.1) THEN
	   WRITE(n,700) ALFA,BET,M,H
    ENDIF
	
	IF(NtipSetkyatom.EQ.2) THEN
	   WRITE(n,750) ALFA,BET,GAMMA,M,H,RalFL
    ENDIF
    

	IF(NtipSetky.EQ.1) THEN
	  WRITE(n,1200)
	  WRITE(n,800) ALFAnew,BETnew,Mnew,Hnew,Ral,Lmin,Lmax
      ENDIF
	
	IF(NtipSetky.EQ.2) THEN
	  WRITE(n,1100)
	  WRITE(n,850) ALFAnew,BETnew,GAMMAnew,Mnew,Hnew,Ral,Lmin,Lmax
      ENDIF
	
	WRITE(n,900) (Nn(I),LB(L(I)+1),LML(I),I=1,IS)    
      WRITE(n,*) 'Files writing expansion function' 	 
	DO IXC=1,IS
         WRITE(n,1000) Nn(IXC),LB(L(IXC)+1),LML(IXC),OUTF(IXC)   
	ENDDO




	return
	end subroutine writefile

   


    


    

    



	 
