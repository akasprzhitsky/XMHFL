
program GeneratorAF
   integer(2) H1,M1,S1,HS1,H2,M2,S2,HS2
 
   call GETTIM(H1,M1,S1,HS1)
   TB=60.*H1+M1+S1/60.
 
  ! main subroutine
   call AF
      
   call GETTIM(H2,M2,S2,HS2)
   TE=60.*H2+M2+S2/60.
   T=TE-TB
   WRITE(6,1) T
 1 FORMAT(' Commercial    time :',F10.2,' min.')
   stop
end program GeneratorAF 

 

! SUBPROGRAM IMPLEMENTS OBTAINING WAVE FUNCTIONS BY THE SOLUTION
! SYSTEMS OF HARTRE-FOCK EQUATIONS
subroutine AF
  use mAF
  implicit real(8)(A-H,O,P,R-Z)
  implicit integer(Q)
      
  common /SHORT/ IDAM,N,NDU,IB,M,NIT,NPS,MF,NZ,NWR,NLST,NN(20),L(20),Q(20),NCN,IW,I30,NST(20),IS,IZAP,IOUT,IFUNC,IPOT,NU,NG,NFK,NGK,IPKOEF,NONREL,J,H,Z,RAL,ACY,VAL,CFS2,E3(20),E4(20),TETA(20),A1,B1,A2,OMEGA,RO1,RO2,TETMAX,RAN,EPS,BET,ALFA
  common /FUNPOT/IS1(5000),IS2(5000),KK(5000),ISD(50),KD(50),R1(6000),RAB1(6000),RAB2(6000),P4(6000),GAM(5000),DEL(50)
  common /FUNC/R8(120400)
  
  character(90) DATFIL,DAFILES,FILE1,FILE2,OUTFIL,TITLE,TITWRK,OUTFTOTALE,OUTFUN
  
  ! Ã¿——»¬€ ƒÀﬂ «¿Ã≈Õ€  Œ›‘‘»÷»≈Õ“Œ¬  ”ÀŒÕŒ¬— Œ√Œ ¬«¿»ÃŒƒ≈…—“¬»ﬂ    
  integer,dimension(2,1000)::IX,IY,K
  real(8),dimension(2,1000)::DCOFF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Ã¿——»¬ ƒÀﬂ «¿œ»—» ¬Õ≈ÿÕ≈√Œ œŒÀﬂ
  real(8),allocatable,dimension(:,:)::RPOT,EepsLagran
  ! M¿——»¬ ƒÀﬂ «¿œ»—» –¿ƒ»¿À‹Õ€’ ◊¿—“≈… ¬ŒÀÕŒ¬€’ ‘”Õ ÷»»
  !real(8),allocatable,dimension(:)::R8
  real(8),allocatable,dimension(:,:)::RFun
  ! ¬—œŒÃŒ√¿“≈À‹Õ€≈ Ã¿——»¬€
  real(8),allocatable,dimension(:)::RO1X,RO2X,RO3X 
  character(1),dimension(10)::LB  
      
  data LB/'s','p','d','f','g','h','i','j','k','l'/
	
	
	
	
	
  WRITE(*,'(A\)') ' Input TITLE file name or * (ATTITLE.DAT) or END--->'
  READ(*,'(A90)') TITLE
     
  IF(TITLE.EQ.'*') TITLE='ATTITLE.DAT'
  IF(TITLE.EQ.'END'.OR.TITLE.EQ.'end') RETURN

  OPEN(UNIT=5,FILE=TITLE)
 
  READ(5,'(A90)') DATFIL
  READ(5,'(A90)') FILE1
  READ(5,'(A90)') FILE2
  READ(5,'(A90)') OUTFIL
  READ(5,'(A90)') OUTFTOTALE
  READ(5,'(A90)') DAFILES
  READ(5,'(A90)') OUTFUN



  write(*,*) DATFIL
  write(*,*) FILE1
  write(*,*) FILE2
  write(*,*) OUTFIL
  write(*,*) OUTFTOTALE
  write(*,*) DAFILES
  write(*,*) OUTFUN



    	
  OPEN(3,FILE=DATFIL)
  OPEN(6,FILE=OUTFIL)
  OPEN(10,FILE=OUTFTOTALE)
  OPEN(70,FILE=OUTFUN)   
 
  REWIND 3
! NcoffCalculLagran-parameter indicating the addition of off-diagonal factors
    ! NcoffCalculLagran = 0 - there are no off-diagonal factors
    ! NcoffCalculLagran = 1 - calculation taking into account off-diagonal factors
    ! NtipSetky-MESH TYPE
    ! NtipSetky = 1-ATOMIC GRID
    ! NtipSetky = 2-MOLECULAR GRID
77 READ(3,*) Z,IS,NtipSetky,BET,GAMMA,RAL,N,H,R,EPS,IZAP,IOUT,IFUNC,IPOT,NFK,NGK,NDP,NZ,NONREL,IPKOEF,NcoffCalculLagran
   
   
   ! œ¿–¿Ã≈“– Œ—“¿ÕŒ¬ » œ–Œ√–¿ÃÃ€
   IF(N.EQ.0) RETURN
   
  
    
   WRITE(6,10) Z,IS,RAL,N,H,R,EPS,IZAP,IOUT,IFUNC,IPOT

10 FORMAT(' Z=',F5.2,7X,'Is=',I2,3X,' Ral=',F6.3,3X,'N=',I4,2X,'H=',F5.3,3X,'R=',F10.3,'  Eps=',D7.1,2X,'Izap=',I3,1X,'Iout=',I2,4X,' Ifunc=',I1,2X,'Ipot=',I3)

   IF(IZAP.NE.0) OPEN(UNIT=1,FILE=FILE1,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
   IF(NDP.NE.0) OPEN(UNIT=2,FILE=FILE2,STATUS='OLD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      
   ! Calculation of integration parameters ALFA, BET
   ! INTERVAL FOR RO
   IF(NtipSetky.EQ.1) THEN
      RO1=-14.D0
	  RO2=RO1+FLOAT(N)*H
      ALFA=(RO2-BET*DLOG(R))/R
   ENDIF
   IF(NtipSetky.EQ.2) THEN
      RO1=-15.D0
      RO2=RO1+FLOAT(N)*H
      ALFA=(RO2-BET*DLOG(R))/R	  
   ENDIF
   
   IF(RAL.GT.1.D-3) THEN
      RTEMP=ALFA*RAL+BET*DLOG(RAL)-RO1
      ALFA=(RO1+H*DINT(RTEMP/H)-BET*DLOG(RAL))/RAL
   ENDIF
   
   
   IF(NtipSetky.EQ.1) THEN
      WRITE(6,5) NZ,NONREL,IPKOEF,ALFA,BET
   ENDIF
   IF(NtipSetky.EQ.2) THEN
      WRITE(6,35) NZ,NONREL,IPKOEF,ALFA,BET,GAMMA
   ENDIF
 5  FORMAT(1X,'Nz=',I3,2X,'Nonrel=',I1,2X,'Ipcoef=',I2,'  Alfa=',F14.8,2X,' Beta=',F6.3)
 35 FORMAT(1X,'Nz=',I3,2X,'Nonrel=',I1,2X,'Ipcoef=',I2,'  Alfa=',F14.8,2X,' Beta=',F6.3,'  Gamma=',F14.8)
   
   NWR=0
   IF(NZ.EQ.-1) NWR=1
   IF(NZ.EQ.-1) NZ=0
   NZ=NZ+1

   ! ¬¬Œƒ»Ã  ŒÕ‘»√”–¿÷»ﬁ 
   READ(3,*)(NN(I),L(I),Q(I),NST(I),I=1,IS)
   WRITE(6,50) Z,(NN(I),LB(L(I)+1),Q(I),I=1,IS)
50 FORMAT(/' Z=',F3.0,'  Configuration :   ',20(I3,A1,I2)/)

   !  Œ›‘‘»÷»≈Õ“€ ¬Õ≈ÿÕ≈√Œ œŒ“≈Õ÷»¿À¿
   NUF=1
   ACY=0
   VAL=0
   IF(RAL.GT.1.D-3) THEN
      READ(3,*) NUF,ACY,VAL
   	  READ(3,*)(ISD(I),KD(I),DEL(I),I=1,NUF)
      WRITE(6,9) RAL,ACY,(NN(ISD(I)),LB(L(ISD(I))+1),KD(I),DEL(I),I=1,NUF)
9     FORMAT(/' Molecule radius  ',F6.3,6X,' Cryst.field coeff.  Y00=',F12.4/5('       K   Delta    ')/20(5(I3,A1,I4,F8.4,4X)/)) 
   END IF
  
   ! «¿Õ”Àﬂ≈Ã œ≈–≈ƒ ¬¬ŒƒŒÃ
   IX=0
   IY=0
   K=0
   DCOFF=0.D0 

   ! ¬¬Œƒ  Œ›‘‘»÷»≈Õ“Œ¬  ”ÀŒÕŒ¬— Œ√Œ ¬«¿»ÃŒƒ≈…—“¬»ﬂ 
   IF(NGK.NE.0) THEN
      DO I=1,NGK
         READ(3,*) IX(1,I),IY(1,I),K(1,I),DCOFF(1,I)   
      ENDDO
   ENDIF

   IF(NFK.NE.0) THEN
      DO I=1,NFK
         READ(3,*) IX(2,I),IY(2,I),K(2,I),DCOFF(2,I) 
      ENDDO
   ENDIF



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ‚˚‰ÂÎÂÌËÂ Ô‡ÏˇÚË ÔÓ‰ Ï‡ÒÒË‚˚
   allocate(Rfun(IS,N+2),stat=ierr)
   if(ierr/=0) then
      write(*,*)'AF MEMORY ON THE FILE "Rfun" IS NOT SELECTED'
	  stop 
   endif
   !¬€ƒ≈Àﬂ≈Ã œ¿Ãﬂ“‹ ƒÀﬂ Ã¿——»¬Œ¬
   allocate(RO1X(N),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'MEMORY ON THE FILE "RO1X" IS NOT SELECTED'
      stop 
   endif 
   allocate(RO2X(N),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'MEMORY ON THE FILE "RO2X" IS NOT SELECTED'
      stop 
   endif 
   allocate(RO3X(N),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'MEMORY ON THE FILE "RO3X" IS NOT SELECTED'
      stop 
   endif 
   allocate(EepsLagran(IS,IS),stat=ierr)
    if(ierr/=0) then
       write(6,*) 'MEMORY ON THE FILE "EepsLagran" IS NOT SELECTED'
      stop 
    endif


   ! «¿Õ”Àﬂ≈Ã œ≈–≈ƒ –¿—◊≈“ŒÃ
   Rfun=0.D0
      


   call ENTR(NF,NG,NU,IS,L)


  ! find auxiliary arrays for integration
   call VAR(NtipSetky,N,H,ALFA,BET,GAMMA,RAL,RO1,RO2,R1,RAB1,RAB2,RO1X,RO2X,RO3X)

   ! WE ALWAYS MEMORY FOR THE EXTERNAL FIELD RECORDING ARRAY
   IF(RAL.GT.1.D-3) THEN
      allocate(RPOT(IS,N),stat=ierr)
	  if(ierr/=0) then
         write(*,*)'AF MEMORY ON THE FILE "RPOT" IS NOT SELECTED'
	     stop 
	  endif
	  RPOT=0.D0
   ENDIF


   ! SUBPROGRAM READS THE EXTERNAL POTENTIAL, WRITES FOR EACH SHELL
   call KOEFRAL(NUF,NDP,N,IS,RAL,ACY,ISD,KD,DEL,R1,RPOT)
   ! SUB-PROGRAM FINDS ANGULAR COEFFICIENTS FOR THE CULONIC INTERACTION
   ! IN THE "MIDDLE FIELD"
   call KOEF(NU,NG,IS,L,Q,NST,IS1,IS2,KK,GAM,NFK,NGK,IX,IY,K,DCOFF)

   ! ISSUANCE OF COULON INTERACTION COEFFICIENTS   
19 FORMAT(/' Electroctatic coefficients ',I4,'(',I4,' Exchange )')
20 FORMAT(2X,I2,A1,I2,A1,I3,F9.4,4X)
22 FORMAT('     G    K    Gamma ')
25 FORMAT('     F    K    Gamma ')
   IF(IPKOEF.NE.-1) THEN
      WRITE(6,19) NU,NG
      WRITE(6,22)
	  DO I=1,NG
         WRITE(6,20) NN(IS1(I)),LB(L(IS1(I))+1),NN(IS2(I)),LB(L(IS2(I))+1),KK(I),GAM(I)
	  ENDDO
	  WRITE(6,*)
	  WRITE(6,25)
      DO I=NG+1,NU
         WRITE(6,20) NN(IS1(I)),LB(L(IS1(I))+1),NN(IS2(I)),LB(L(IS2(I))+1),KK(I),GAM(I)
	  ENDDO
	  WRITE(6,*)
   ENDIF 
   
   
   
   
   
   
     

  
   ! SUBROGRAM GETS A ZERO APPROXIMATION
   call ZERO(Z,IS,NN,Q,L,NZ,NWR,N,RAL,ACY,VAL,R1,R8,P4,E3,Rfun)
    
   

   ! SUBPROGRAM FOR SOLVING THE SYSTEM OF HartreeñFock  EQUATIONS FOR THIS CONFIGURATION
   call DHF(N,NG,NU,NLST,IPOT,NZ,IS,NIT,NWR,NONREL,NcoffCalculLagran,Z,H,EPS,ALFA,BET,RO1,TETMAX,A1,B1,A2,OMEGA,GAMA,RAL,VAL,RAN,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X,RPOT,Rfun,EepsLagran) 
   
   ! WE CARRY OUT UNITARY CONVERSION
   !call UnitaryTransformation(IS,N,R8,EepsLagran) 

   ! SUB-PROGRAM FOR CALCULATING THE TOTAL ENERGY OF AN ATOM WITH RELATIVISTIC CORRECTIONS
   IF(NtipSetky.EQ.1) THEN
      IF(NONREL.NE.0) THEN 
         call TOTAL(Z,N,H,NU,NG,IS,NCN,RAN,VAL,TETA,Q,L,KK,IS1,IS2,GAM,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X)	      
      ENDIF
   ENDIF
   ! SUBPROGRAM FOR CALCULATING THE TOTAL ENERGY OF AN ATOM IN THE HartreeñFock  APPROXIMATION
   IF(NONREL.EQ.0) THEN 
      call TOTALXF(Z,N,H,NU,NG,IS,Q,L,KK,IS1,IS2,GAM,NN,NST,R1,R8,RO1X)
   ENDIF
   
   ! œŒƒœ–Œ√–¿ÃÃ¿ «¿œ»—» –≈«”À‹“¿“Œ¬
   call INTER(IZAP,IS,N,IOUT,IPOT,IFUNC,H,NN,L,Q,R1,RO1X,R8)
   !-----------------------------------------------------

15 FORMAT(/)	
   WRITE(6,15)

   ! REMOVING ARRAYS FROM MEMORY
   deallocate(Rfun,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'AF THE MASSIV "Rfun" IS NOT REMOVED FROM MEMORY'
 	  stop 
   endif 
   IF(RAL.GT.1.D-3) THEN
      deallocate(RPOT,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'AF THE MASSIV "RPOT" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
   ENDIF

   ! REMOVING ARRAYS FROM MEMORY
    deallocate(RO1X,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "RO1X" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    deallocate(RO2X,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "RO2X" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
	deallocate(RO3X,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "RO3X" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    deallocate(EepsLagran,stat=ierr)
    if(ierr/=0) then
       write(6,*) 'THE FILE "EepsLagran" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
  
   GOTO 77

  return
end subroutine AF



