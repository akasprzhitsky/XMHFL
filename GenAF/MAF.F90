       
module mAF
  implicit none 


  contains

  
  ! ондопнцпюллю пеьемхъ яхярелш спюбмемхи уюпрпх-тнйю
  ! нохяюмхе оюпюлерпнб ондопнцпюллш 


  subroutine DHF(N,NG,NU,NLST,IPOT,NZ,IS,NIT,NWR,NONREL,NcoffCalculLagran,Z,H,EPS,ALFA,BET,RO1,TETMAX,A1,B1,A2,OMEGA,GAMA,RAL,VAL,RAN,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X,RPOT,Rfun,EepsLagran) 
    implicit none
	integer::N,NG,NU,NLST,IPOT,NZ,IS,NIT,NWR,NONREL,NcoffCalculLagran
	real(8)::Z,H,EPS,ALFA,BET,RO1,TETMAX,A1,B1,A2,OMEGA,GAMA,RAL
	real(8)::VAL,RAN
	integer,dimension(:)::NN,L,Q,KK,IS1,IS2
    real(8),dimension(:)::E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X 
    real(8),dimension(:,:)::RPOT,Rfun,EepsLagran
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
	integer::MMAX,MBEG,KKK,M,IB,I,J,K,ierr,IXXX,IYYY,IZZZ,NonEdf,IDISXC,NfactorLagrang,NfindD,NnumeroM,Nclusd
	integer,allocatable,dimension(:)::Liqw,NumeroIS
    real(8),allocatable,dimension(:)::R2,Rfun12
	real(8),allocatable,dimension(:,:,:)::RintOrt,EespLagr                                                                                                                                                                                                                      
    character(1),dimension(10)::LB 
	     
	DATA LB/'s','p','d','f','g','h','i','j','k','l'/
	

     
             

    MMAX=2500
    MBEG=1
    KKK=(N-2)*IS
    NLST=0
   	IF(NZ.NE.1) MBEG=2
    IF(IPOT.GT.5) MMAX=IPOT

    !бшдекъел оюлърэ дкъ люяяхбнб 
    allocate(Liqw(IS),stat=ierr)
    if(ierr/=0) then
      write(6,*) 'DHF'
      write(6,*) 'MEMORY ON THE FILE "Liqw" IS NOT SELECTED'
      stop 
    endif
	allocate(NumeroIS(IS),stat=ierr)
    if(ierr/=0) then
      write(6,*) 'DHF'
      write(6,*) 'MEMORY ON THE FILE "NumeroIS" IS NOT SELECTED'
      stop 
    endif  
    allocate(RintOrt(3,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(6,*) 'DHF'
      write(6,*) 'MEMORY ON THE FILE "RintOrt" IS NOT SELECTED'
      stop 
    endif 
    allocate(EespLagr(3,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(6,*) 'DHF'
      write(6,*) 'MEMORY ON THE FILE "EespLagr" IS NOT SELECTED'
      stop 
    endif 
    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R2(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'MAING'
      write(*,*) 'MEMORY ON THE FILE "R2" IS NOT SELECTED'
      stop 
    endif 
	allocate(Rfun12(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'MAING'
      write(*,*) 'MEMORY ON THE FILE "Rfun2" IS NOT SELECTED'
      stop 
    endif 
	!бшдекъел оюлърэ дкъ люяяхбнб
  	RintOrt=0.D0
	EepsLagran=0.D0
    EespLagr=0.D0
    
	
	! нопедекъел вхякн медхюцнмюкэмшу лмнфхрекеи кюцпюмфю
    Liqw=0
	NumeroIS=0
	IDISXC=0
	DO IXXX=1,IS
	   ! опнбепъел еярэ кх б люяяхбе нанкнвйх я L(IXXX)
	   NfindD=0
	   DO IZZZ=1,IDISXC
          IF(L(IXXX).EQ.Liqw(IZZZ)) THEN
             NfindD=1
			 EXIT
		  ENDIF
	   ENDDO
       IF(Q(IXXX).NE.(4*L(IXXX)+2).AND.NfindD.EQ.0) THEN
          IDISXC=IDISXC+1
		  Liqw(IDISXC)=L(IXXX)
	   ENDIF
    ENDDO  
    
	! тнплхпсел люяяхб нанкнвей дкъ йнрнпшу менаундхлн б спюбмемхх хяонкэгнбюрэ медхюцнмюкэмше лмнфхрекх кюцпюмфю
    NnumeroM=0
	DO IXXX=1,IDISXC
       DO IYYY=1,IS
          IF(Liqw(IXXX).EQ.L(IYYY)) THEN
             NnumeroM=NnumeroM+1
			 ! тнплхпсел люяяхб нанкнвей йнмтхцспюжхх дкъ йнрнпшу менаундхлн хяонкэгнбюрэ медхюцнмюкэмше лмнфхрекх кюцпюмфю
			 NumeroIS(NnumeroM)=IYYY 
		  ENDIF
	   ENDDO
    ENDDO

    NfactorLagrang=0
	DO IYYY=1,IS
       ! опнбепъел еярэ кх дюммюъ нанкнвйю б яохяйе дкъ йнрнпшу менаундхлн хяонкэгнбюрэ лмнфхрекх кюцпюмфю
       Nclusd=0
	   DO IXXX=1,NnumeroM
	      IF(IYYY.EQ.NumeroIS(IXXX)) THEN
             Nclusd=1
			 EXIT
		  ENDIF 
	   ENDDO 
       ! пюявхршбюел вхякн медхюцнмюкэмшу лмнфхрекеи кюцпюмфю 
       IF(Nclusd.EQ.1) THEN
	      DO IXXX=1,IS 
	         IF(NN(IYYY).NE.NN(IXXX).AND.L(IYYY).EQ.L(IXXX)) THEN
	             NfactorLagrang=NfactorLagrang+1
			 ENDIF
		  ENDDO
	   ENDIF 
	ENDDO

    
    
    

	! жхйкш он янцкнянбюмхч бнкмнбшу тсмйжхи 
	DO M=MBEG,MMAX
       TETMAX=0.D0
	   NonEdf=0
	   ! пюявер медхюцнмюкэмшу лмнфхрекеи кюцпюмфю
	   IF(M.NE.1.AND.NcoffCalculLagran.EQ.1) THEN
	       IF(NfactorLagrang.NE.0) THEN
		      call CalculationEpsLagrangALFA(NfactorLagrang,NnumeroM,NumeroIS,M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran,EespLagr,RintOrt)
	          NonEdf=1 
		   ENDIF
         ! ELSE 
         !  EepsLagran=0.D0
		 !  WRITE(*,*) 'DDDD' 
	   ENDIF 
       ! жхйк он нанкнвйюл йнмтхцспюжхх
	   DO IB=NZ,IS
	    !  WRITE(6,*) 'ETAP 1',IB 
          call MAING(M,IB,IS,N,NG,NU,NLST,IPOT,NonEdf,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran)
	   ENDDO
       IF(IPOT.EQ.3.AND.NcoffCalculLagran.EQ.1) THEN 
	      write(6,98) M
		   DO IB=NZ,IS
		      write(6,198) (EepsLagran(IB,IXXX),IXXX=1,IS)     
           ENDDO
	   ENDIF
	  ! IF(M.NE.1) THEN
	  !    call UnitaryTransformation(IS,N,H,R8,RO1X,EepsLagran) 
	  ! ENDIF 
	   NIT=M
	   IF(TETMAX/(IS-NZ+1).LE.EPS) THEN 
	      EXIT
  	   ENDIF
	ENDDO



    IF(NZ.NE.1) THEN
	   NonEdf=0 
       NLST=1
	   ! пюявер медхюцнмюкэмшу лмнфхрекеи кюцпюмфю
	   IF(NcoffCalculLagran.EQ.1) THEN
	      IF(NfactorLagrang.NE.0) THEN
	         call CalculationEpsLagrangALFA(NfactorLagrang,NnumeroM,NumeroIS,M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran,EespLagr,RintOrt)
	         NonEdf=1
		  ENDIF	  
	   ENDIF 
	   ! жхйк он нанкнвйюл йнмтхцспюжхх 
       DO IB=1,IS
          call MAING(NIT,IB,IS,N,NG,NU,NLST,IPOT,NonEdf,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran)
       ENDDO
	   IF(IPOT.EQ.3.AND.NcoffCalculLagran.EQ.1) THEN 
	      write(6,98) M
		   DO IB=1,IS
		      write(6,198) (EepsLagran(IB,IXXX),IXXX=1,IS)     
           ENDDO
	   ENDIF 
    ENDIF
    
    
     




    IF(NWR.NE.0) THEN 
       call REWRIT(NWR,IS,N,R2,R8,Rfun)
    ENDIF
      
	! МНПЛХПНБМЙЮ ОНКСВЕММШУ БНКМНБШУ ТСМЙЖХИ Х ГЮОХЯЭ НДМНЩКЕЙРПНММШУ ЩМЕПЦХИ
    DO I=1,IS
       call MDRUM(I,2,N,R2,R8)
       VAL=1./DSQRT(R2(2))
       DO J=1,N
          R2(J+2)=R2(J+2)*VAL
       ENDDO
	   R2(1)=-R2(1)
       call MDRUM(I,7,N,R2,R8)
    ENDDO


    RAN=(2.2677D-5)*Z**0.3333D0
      
    
	! бшдювю пегскэрюрнб

8   FORMAT(' R(1)=',E8.3,'      R(N)=',F10.3,'      Number of iterations = ',I3,'    Nucl. R=',E15.7)
88  FORMAT(/'  nl     Energy      Radius      Teta     Norm ')
9   FORMAT(I3,A1,2F12.6,2E11.3)
98  FORMAT(' N=',I2,' One Electron Energy')
198 FORMAT(4X,20(1X,F14.7))   
     
	WRITE(6,8) R1(1),R1(N),NIT,RAN
    IF(NONREL.EQ.0) WRITE(6,88)
    

   

    DO I=1,IS
       K=(I-1)*(N+2)+2
       DO J=1,N
          R2(J)=R8(J+K)**2*R1(J)
       ENDDO
       VAL=SIMP(R2,N,0,H,RO1X)
       IF(NONREL.EQ.0) THEN
	      WRITE(6,9) NN(I),LB(L(I)+1),R8(K-1),VAL,TETA(I),R8(K)
       ENDIF
    ENDDO

      
	IF(NONREL.EQ.0) WRITE(6,9)

    ! сдюкемхе люяяхбнб хг оълърх 
	deallocate(EespLagr,stat=ierr)
    if(ierr/=0) then
       write(6,*) 'MAING'
       write(6,*) 'THE FILE "EespLagr" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
	deallocate(NumeroIS,stat=ierr)
    if(ierr/=0) then
       write(6,*) 'MAING'
       write(6,*) 'THE FILE "NumeroIS" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
	deallocate(Liqw,stat=ierr)
    if(ierr/=0) then
       write(6,*) 'MAING'
       write(6,*) 'THE FILE "Liqw" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    deallocate(RintOrt,stat=ierr)
    if(ierr/=0) then
       write(6,*) 'MAING'
       write(6,*) 'THE FILE "RintOrt" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    deallocate(R2,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "Rfun1" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    deallocate(Rfun12,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "Rfun2" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    return
  end  subroutine DHF


  ! ондопнцпюллю пюяверю медхюцнмюкэмшу лмнфхрекеи кюцпюмфю
  subroutine CalculationEpsLagrang(NfactorLagrang,NnumeroM,NumeroIS,M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran,Eesp,Rint)
    implicit none
    integer::NfactorLagrang,NnumeroM,M,IB,IS,N,NG,NU,NLST,IPOT 
   	real(8)::Z,H,EPS,ALFA,BET,TETMAX,RAL
    integer,dimension(:)::NN,L,Q,KK,IS1,IS2,NumeroIS
    real(8),dimension(:)::E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X 
    real(8),dimension(:,:)::RPOT,EepsLagran
    real(8),dimension(:,:,:)::Eesp,Rint
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IXKLS,IDPS,IDFGH,YDFGH,Nsum,ITIT,Nx,Lx,Qx,NcofLagER,ierr,IIYU,MAXCNVx,Nparam,INDEXEXIT
	! люйяхлюкэмне гмювемхъ вхякю жхйкнб опх онхяйх лмнфхрекеи кюцпюмфю
	integer,parameter::NMaxIteration=200
	real(8)::RGGJ,EcoffLag,Tx,Cx,TETMAXx
	! рнвмнярэ нпрнцнмюкхгхжхх тсмйжхи
	real(8),parameter::EpsInt=0.00000000001D0
	integer,allocatable,dimension(:,:)::INDEXKLUZ
    real(8),allocatable,dimension(:)::Pznl,Qznl,RfunNX,RfunNY,E3x,E4x,TETAx,R8x    
	real(8),allocatable,dimension(:,:)::Pz,Qz
    real(8),allocatable,dimension(:,:,:)::RIntort,EepsL
	real(8),dimension(6002)::R2x,R5x  
    

    !бшдекъел оюлърэ дкъ люяяхбнб 
   	allocate(RfunNX(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RfunNX" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(RfunNY(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RfunNY" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(INDEXKLUZ(IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "INDEXKLUZ" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(E3x(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "E3x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(E4x(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "E4x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif
	allocate(TETAx(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "TETAx" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif  
    allocate(Pz(IS,N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Pz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(Qz(IS,N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Qz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(Pznl(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Pznl" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(Qznl(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Qznl" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 

    allocate(RIntort(4,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RIntort" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(EepsL(4,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "EepsL" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(R8x((N+2)*IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "R8x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
     
    ! гюмскъел оепед пюявернл
    Pz=0.D0
	Qz=0.D0
	Pznl=0.D0
	Qznl=0.D0
    RIntort=0.D0
    EepsL=0.D0
    RfunNX=0.D0
	RfunNY=0.D0
	R8x=0.D0
    
	! нясыеярбкъел пюявер хмрецпюкнб нпрнцнмюкэмнярх х лмнфхрекеи кюцпюмфю дкъ мювюкэмнцн щрюою
    IF(M.EQ.2) THEN
	    DO IDFGH=1,IS
	       call MDRUM(IDFGH,2,N,RfunNX,R8)
	       DO YDFGH=1,IS
              IF(NN(IDFGH).NE.NN(YDFGH).AND.L(IDFGH).EQ.L(YDFGH)) THEN
                 call MDRUM(YDFGH,2,N,RfunNY,R8)
	             RGGJ=RIntOrtXX(RfunNX,RfunNY,N,H,RO1X) 
		         RIntort(1,IDFGH,YDFGH)=RGGJ
                 RIntort(2,IDFGH,YDFGH)=RGGJ
                 RIntort(3,IDFGH,YDFGH)=RGGJ
                 EepsL(1,IDFGH,YDFGH)=0.D0 
                 EepsL(2,IDFGH,YDFGH)=0.D0 
                 EepsL(3,IDFGH,YDFGH)=RGGJ
                 ! наыхи оюпюлерп
				 Rint(1,IDFGH,YDFGH)=RGGJ
                 Rint(2,IDFGH,YDFGH)=RGGJ
                 Rint(3,IDFGH,YDFGH)=RGGJ
	             Eesp(1,IDFGH,YDFGH)=0.D0
                 Eesp(2,IDFGH,YDFGH)=0.D0
				 Eesp(3,IDFGH,YDFGH)=RGGJ
			  ENDIF
	       ENDDO
	    ENDDO 
	   ELSE
        DO IDFGH=1,IS
	       call MDRUM(IDFGH,2,N,RfunNX,R8)
	       DO YDFGH=1,IS
              IF(NN(IDFGH).NE.NN(YDFGH).AND.L(IDFGH).EQ.L(YDFGH)) THEN
                 call MDRUM(YDFGH,2,N,RfunNY,R8)
	             RGGJ=RIntOrtXX(RfunNX,RfunNY,N,H,RO1X) 
		         RIntort(1,IDFGH,YDFGH)=Rint(1,IDFGH,YDFGH)
                 RIntort(2,IDFGH,YDFGH)=Rint(2,IDFGH,YDFGH)
                 RIntort(3,IDFGH,YDFGH)=RGGJ
                 EepsL(1,IDFGH,YDFGH)=Eesp(1,IDFGH,YDFGH)
                 EepsL(2,IDFGH,YDFGH)=Eesp(2,IDFGH,YDFGH)
                 EepsL(3,IDFGH,YDFGH)=Eesp(3,IDFGH,YDFGH)
                 WRITE(6,*) Eesp(1,IDFGH,YDFGH),Eesp(2,IDFGH,YDFGH),Eesp(3,IDFGH,YDFGH)
                 WRITE(6,*) Rint(1,IDFGH,YDFGH),Rint(2,IDFGH,YDFGH),RGGJ
                 WRITE(6,*)
 				 ! наыхи оюпюлерп
				 Eesp(1,IDFGH,YDFGH)=Eesp(2,IDFGH,YDFGH)
                 Eesp(2,IDFGH,YDFGH)=Eesp(3,IDFGH,YDFGH)
				 Rint(1,IDFGH,YDFGH)=Rint(2,IDFGH,YDFGH)
                 Rint(2,IDFGH,YDFGH)=RGGJ
                 
				

	          ENDIF
	       ENDDO
	    ENDDO 
	ENDIF  
	
	
     



	
	  
			
    ! тнплхпсел люяяхбш тсмйжхи Pz х Qz дкъ бяеу нанкнвей
	DO IXKLS=1,IS   
	   call Calcul_Pz_Qz_nl(M,IXKLS,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,Pznl,Qznl) 
       DO IDPS=1,N
          Pz(IXKLS,IDPS)=Pznl(IDPS)
	      Qz(IXKLS,IDPS+1)=Qznl(IDPS+1)
	   ENDDO
    ENDDO
    
    ! жхйк дкъ онксвемхе лмнфхрекеи кюцпюмфю 
    Nsum=0
	INDEXEXIT=0
	DO WHILE(Nsum.NE.NfactorLagrang) 
	   ! опнбепъел опхбшьем опедек вхякю хрепюжхи  
       IF(INDEXEXIT.EQ.NMaxIteration) THEN 
	      EXIT
	   ENDIF
	   !IF(M.EQ.25) THEN
	   !   WRITE(6,*) 'param',Nsum,NfactorLagrang
	   !ENDIF
	   DO IXKLS=1,IS 
          E3x(IXKLS)=E3(IXKLS)
		  E4x(IXKLS)=E4(IXKLS) 
		  TETAx(IXKLS)=TETA(IXKLS)    
       ENDDO    
       
       ! щрюо 1. пеьемхе яхярелш спюбмемхи уюпрпх-тнйю дкъ дюммшу гмювемхи медхюцнмюкэмшу лмнфхрекеи кюцпюмфю 
       DO IXKLS=1,IS 
          Nx=NN(IXKLS)
          Lx=L(IXKLS)
          Qx=Q(IXKLS)
          ! гюохяшбюел тсмйжхх Pz х Qz дкъ дюммни нанкнвйх 
          DO ITIT=1,N   
             Pznl(ITIT)=Pz(IXKLS,ITIT)
			 Qznl(ITIT+1)=Qz(IXKLS,ITIT+1)
          ENDDO
          ! бшъямъел менаундхлнярэ днаюбйх лмнфхрекеи кюцпюмфю
          NcofLagER=0
		  DO ITIT=1,NnumeroM
             IF(IXKLS.EQ.NumeroIS(ITIT)) THEN
                NcofLagER=1
				EXIT
			 ENDIF
          ENDDO
          ! опнбепъел менаундхлнярэ днаюбйх якюцюелшу я медхюцнмюкэмшлх лмнфхрекълх кюцпюмфю
          IF(NcofLagER.EQ.1) THEN
             DO ITIT=1,IS
                IF(NN(ITIT).NE.NN(IXKLS).AND.L(ITIT).EQ.L(IXKLS)) THEN
				   call MDRUM(ITIT,2,N,RfunNX,R8)
			       EcoffLag=EepsL(3,IXKLS,ITIT)/FLOAT(Q(IXKLS))
                   DO IIYU=1,N
                      Qznl(IIYU+1)=Qznl(IIYU+1)+EcoffLag*RfunNX(IIYU+2)/(DSQRT(RfunNX(2))*RO1X(IIYU)**2)
                   ENDDO 
		        ENDIF 
             ENDDO 
          ENDIF
          
		  ! яВХРШБЮЕЛ ТСМЙЖХЧ ОПЕДШДСЫЕИ ХРЕПЮЖХХ
          call MDRUM(IXKLS,2,N,R2x,R8)
		  call MDRUM(IXKLS,2,N,R5x,R8)
		  ! пеьюел спюбмемхе дкъ дюммни нанкнвйх 
          call PrimitiveDifferentialEquation(IXKLS,Nx,Lx,N,M,NLST,MAXCNVx,H,TETMAXx,EPS,E3x,E4x,TETAx,R2x,Qznl,R5x,Pznl,RO1X,Tx,Cx)     
		  ! ГЮОХЯШБЮЕЛ ОНКСВЕММСЧ ТСМЙЖХЧ
	      call MDRUM(IXKLS,7,N,R2x,R8x)               	    
       ENDDO
       !  щрюо 2. бшвхякъел хмрецпюкш нпрнцнмюкэмнярх мю онквеммшу тсмйжхъу 
       !гюмскъел оепед пюявернл
	   INDEXKLUZ=0
	   DO IDFGH=1,IS
	      call MDRUM(IDFGH,2,N,RfunNX,R8x)
	      DO YDFGH=1,IS
             call MDRUM(YDFGH,2,N,RfunNY,R8x)
	         RGGJ=RIntOrtXX(RfunNX,RfunNY,N,H,RO1X) 
			 IF(DABS(RGGJ).GT.EpsInt) THEN
		         ! ХМРЕЦПЮКШ НПРНЦНМЮКЭМНЯРХ
			     ! АСТЕП УПЮМЕМХЪ ЯЮЛНЦН "ЯРЮПНЦН" ГМЮВЕМХЪ
                 RIntort(4,IDFGH,YDFGH)=RIntort(1,IDFGH,YDFGH)
                 ! ЯЛЕЫЮЕЛ ГМЮВЕМХЪ
			     RIntort(1,IDFGH,YDFGH)=RIntort(2,IDFGH,YDFGH)
                 RIntort(2,IDFGH,YDFGH)=RGGJ
			     ! ЯЛЕЫЮЕЛ ОПХАКХФМХЪ ДКЪ ГЮОХЯХ ОНЯКЕДМЕЦН
                 ! АСТЕП УПЮМЕМХЪ ЯЮЛНЦН "ЯРЮПНЦН" ГМЮВЕМХЪ
			     EepsL(4,IDFGH,YDFGH)=EepsL(1,IDFGH,YDFGH)
                 ! ЯЛЕЫЮЕЛ ГМЮВЕМХЪ
                 EepsL(1,IDFGH,YDFGH)=EepsL(2,IDFGH,YDFGH)
                 EepsL(2,IDFGH,YDFGH)=EepsL(3,IDFGH,YDFGH)
	            ELSE
                 ! сйюгшбюел врн дюммши хмрецпюк нпрнцнмюкэмнярх яннрберярбсер рнвмнярх пюяверю 
				 ! х лмнфхрекэ кюцпюмфю дкъ мецн ме мсфмн онйю лемърэ
				 INDEXKLUZ(IDFGH,YDFGH)=1
			  ENDIF
		  ENDDO
	   ENDDO
	   ! щрюо 3. оНКСВЮЕЛ "МНБНЕ" ОПХАКХФЕМХЕ ДКЪ ЛМНФХРЕКЕИ КЮЦПЮМФЮ
	   Nsum=0
	   DO IXKLS=1,IS  
          ! бшъямъел мсфмн кх дкъ дюммни тсмйжхх нанкнвйх пюявхршбюрэ лнфхрекеи кюцпюмфю
          NcofLagER=0
		  DO ITIT=1,NnumeroM
             IF(IXKLS.EQ.NumeroIS(ITIT)) THEN
                NcofLagER=1
				EXIT
			 ENDIF
          ENDDO
          ! опнбепъел менаундхлнярэ пюяверю медхюцнмюкэмшлх лмнфхрекълх кюцпюмфю
          IF(NcofLagER.EQ.1) THEN
             DO ITIT=1,IS
                IF(NN(ITIT).NE.NN(IXKLS).AND.L(ITIT).EQ.L(IXKLS)) THEN
				   IF(INDEXKLUZ(IXKLS,ITIT).EQ.0) THEN
				       ! пЮЯВЕР МНБНЦН ОПХАКХФЕМХЪ ДКЪ МЕДХЮЦНМЮКЭМНЦН ЛМНФХРЕКЪ КЮЦПЮМФЮ
  	                   !write(6,*) 'IS',IXKLS,ITIT
				       call CalculNewAproxsimationEpsLagran(Nparam,EepsL(1,IXKLS,ITIT),RIntort(1,IXKLS,ITIT),EepsL(2,IXKLS,ITIT),RIntort(2,IXKLS,ITIT),EepsL(3,IXKLS,ITIT))
				       !IF(M.EQ.25) THEN
				       !   write(6,*) IXKLS,ITIT
				       !   write(6,*) EepsL(1,IXKLS,ITIT),RIntort(1,IXKLS,ITIT)
				       !   write(6,*) EepsL(2,IXKLS,ITIT),RIntort(2,IXKLS,ITIT)
                       !   write(6,*) Nparam,EepsL(3,IXKLS,ITIT)
				       !ENDIF
                      ELSE
                       ! дюммши хмрецпюк нпрнцнмюкэмнярх сдбнкербнпъер сякнбхч нпрнцнмюкэмнярх
					   Nparam=1 
				   ENDIF
				   ! опнбепъел бнгмхй йнпемэ хкх мер (мюидем кх лмнфхрекэ кюцпюмфю йнрнпши опхбндхр й гюмскемхч хмрецпюкю нпрнцнмюкэмнярх)
				   IF(Nparam.EQ.1) THEN
                      Nsum=Nsum+1
				   ENDIF
				ENDIF 
             ENDDO 
          ENDIF
       ENDDO
       
      !хмдей бшундю хг жхйкю  
       INDEXEXIT=INDEXEXIT+1 

    ENDDO
   !WRITE(6,*) 'EXIT',INDEXEXIT,Nsum
   IF(INDEXEXIT.NE.1) THEN
      ! онксвеммше лмнфхрекх кюцпюмфю
      DO IDFGH=1,IS
         DO YDFGH=1,IS 
            EepsLagran(IDFGH,YDFGH)=EepsL(3,IDFGH,YDFGH)
            Eesp(3,IDFGH,YDFGH)=EepsL(3,IDFGH,YDFGH)
         ENDDO
      ENDDO
   ENDIF  


  
    ! сдюкъел хг оюлърх люяяхбш	
	deallocate(INDEXKLUZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "INDEXKLUZ" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(R8x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "R8x" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(TETAx,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "TETAx" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(E3x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "E3x" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(E4x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "E4x" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
	deallocate(Pz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Pz" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(Qz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Qz" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
    deallocate(Pznl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Pznl" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(Qznl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Qznl" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
	deallocate(RIntort,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RIntort" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(EepsL,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "EepsL" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
    deallocate(RfunNX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RfunNX" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(RfunNY,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RfunNY" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif

    return
  end  subroutine CalculationEpsLagrang




  



! ондопнцпюллю пюяверю мнбнцн опхакхфемхъ дкъ лмнфхрекъ кюцпюмфю 
 subroutine CalculNewAproxsimationEpsLagran(NNZX,X1,Y1,X2,Y2,X3) 
    implicit none
    integer::NNZX
   	real(8)::X1,Y1,X2,Y2,X3

        
    NNZX=0
    ! ОПНБЕПЪЕЛ ЯКСВЮИ НДМНЦН ГМЮЙЮ ТСМЙЖХИ
    IF(ISIGNUM(Y1).EQ.ISIGNUM(Y2)) THEN 
       ! НОПЕДЕКЪЕЛ Б ЙЮЙНИ НАКЮЯРХ МЮУНДХРЯЪ РНВЙХ 
  	   IF(Y1.GE.0.D0) THEN
         ! нопедекъел сцнк мюйкнмю цпютхйю
	  	 IF((Y2-Y1)/(X2-X1).GT.0.D0) THEN
             IF(X2.GT.X1) THEN
	   	         X3=X1-0.01D0
	 		    ELSE
                 X3=X2-0.01D0
	  	     ENDIF    
		    ELSE
		     IF(X2.GT.X1) THEN
	   	         X3=X2+0.01D0
			    ELSE
                 X3=X1+0.01D0
	  	     ENDIF    
	 	ENDIF       
	   ELSE
        ! нопедекъел сцнк мюйкнмю цпютхйю
		IF((Y2-Y1)/(X2-X1).GT.0.D0) THEN
            IF(X2.GT.X1) THEN
	 	         X3=X2+0.01D0
		       ELSE
                 X3=X1+0.01D0
	  	    ENDIF    
	 	   ELSE
	        IF(X2.GT.X1) THEN
		        X3=X1-0.01D0
	 	       ELSE
                X3=X2-0.01D0
	  	    ENDIF    
		ENDIF       
     ENDIF
  ENDIF 
  
  ! яксвюи пюгмшу гмюйнб (мюидем йнпемэ б дюммнл хмрепбюке)
  IF(ISIGNUM(Y1).NE.ISIGNUM(Y2)) THEN 
     NNZX=1   
     X3=X1-Y1*(X2-X1)/(Y2-Y1)   
  ENDIF  
  
  
  !write(6,*) X1,Y1
  !write(6,*) X2,Y2
  !WRITE(6,*) NNZX,X3
  !WRITE(6,*)
  

       
    return
  end subroutine CalculNewAproxsimationEpsLagran




 ! ондопнцпюллю пюяверю медхюцнмюкэмшу лмнфхрекеи кюцпюмфю опълни пюявер медхюцнмюкэмшу лмнфхрекеи 
  subroutine CalculationEpsLagrangALFA(NfactorLagrang,NnumeroM,NumeroIS,M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran,Eesp,Rint)
    implicit none
    integer::NfactorLagrang,NnumeroM,M,IB,IS,N,NG,NU,NLST,IPOT 
   	real(8)::Z,H,EPS,ALFA,BET,TETMAX,RAL
    integer,dimension(:)::NN,L,Q,KK,IS1,IS2,NumeroIS
    real(8),dimension(:)::E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X 
    real(8),dimension(:,:)::RPOT,EepsLagran
    real(8),dimension(:,:,:)::Eesp,Rint
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IXKLS,IDPS,IDFGH,YDFGH,Nsum,ITIT,Nx,Lx,Qx,NcofLagER,ierr,IIYU,MAXCNVx,Nparam,INDEXEXIT
	! люйяхлюкэмне гмювемхъ вхякю жхйкнб опх онхяйх лмнфхрекеи кюцпюмфю
	integer,parameter::NMaxIteration=200
	real(8)::RGGJ,EcoffLag,Tx,Cx,TETMAXx
	! рнвмнярэ нпрнцнмюкхгхжхх тсмйжхи
	real(8),parameter::EpsInt=0.00000000001D0
	integer,allocatable,dimension(:,:)::INDEXKLUZ
    real(8),allocatable,dimension(:)::Pznl,Qznl,RfunNX,RfunNY,E3x,E4x,TETAx,R8x    
	real(8),allocatable,dimension(:,:)::Pz,Qz
    real(8),allocatable,dimension(:,:,:)::RIntort,EepsL
	real(8),dimension(6002)::R2x,R5x  
    

    !бшдекъел оюлърэ дкъ люяяхбнб 
   	allocate(RfunNX(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RfunNX" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(RfunNY(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RfunNY" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(INDEXKLUZ(IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "INDEXKLUZ" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(E3x(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "E3x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(E4x(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "E4x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif
	allocate(TETAx(IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "TETAx" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif  
    allocate(Pz(IS,N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Pz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(Qz(IS,N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Qz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(Pznl(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Pznl" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(Qznl(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "Qznl" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 

    allocate(RIntort(4,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "RIntort" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(EepsL(4,IS,IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "EepsL" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(R8x((N+2)*IS),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CalculationEpsLagrang'
      write(*,*) 'MEMORY ON THE FILE "R8x" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
     
    ! гюмскъел оепед пюявернл
    Pz=0.D0
	Qz=0.D0
	Pznl=0.D0
	Qznl=0.D0
    RIntort=0.D0
    EepsL=0.D0
    RfunNX=0.D0
	RfunNY=0.D0
	R8x=0.D0
    
	! нясыеярбкъел пюявер хмрецпюкнб нпрнцнмюкэмнярх х лмнфхрекеи кюцпюмфю дкъ мювюкэмнцн щрюою
    IF(M.EQ.2) THEN
	    DO IDFGH=1,IS
	       call MDRUM(IDFGH,2,N,RfunNX,R8)
	       DO YDFGH=1,IS
              IF(NN(IDFGH).NE.NN(YDFGH).AND.L(IDFGH).EQ.L(YDFGH)) THEN
                 call MDRUM(YDFGH,2,N,RfunNY,R8)
	             RGGJ=RIntOrtXX(RfunNX,RfunNY,N,H,RO1X) 
		         RIntort(1,IDFGH,YDFGH)=RGGJ
                 RIntort(2,IDFGH,YDFGH)=RGGJ
                 RIntort(3,IDFGH,YDFGH)=RGGJ
                 EepsL(1,IDFGH,YDFGH)=0.D0 
                 EepsL(2,IDFGH,YDFGH)=0.D0 
                 EepsL(3,IDFGH,YDFGH)=RGGJ
                 ! наыхи оюпюлерп
				 Rint(1,IDFGH,YDFGH)=RGGJ
                 Rint(2,IDFGH,YDFGH)=RGGJ
                 Rint(3,IDFGH,YDFGH)=RGGJ
	             Eesp(1,IDFGH,YDFGH)=0.D0
                 Eesp(2,IDFGH,YDFGH)=0.D0
				 Eesp(3,IDFGH,YDFGH)=RGGJ
			  ENDIF
	       ENDDO
	    ENDDO 
	   ELSE
        DO IDFGH=1,IS
	       call MDRUM(IDFGH,2,N,RfunNX,R8)
	       DO YDFGH=1,IS
              IF(NN(IDFGH).NE.NN(YDFGH).AND.L(IDFGH).EQ.L(YDFGH)) THEN
                 call MDRUM(YDFGH,2,N,RfunNY,R8)
	             RGGJ=RIntOrtXX(RfunNX,RfunNY,N,H,RO1X) 
		         RIntort(1,IDFGH,YDFGH)=Rint(1,IDFGH,YDFGH)
                 RIntort(2,IDFGH,YDFGH)=Rint(2,IDFGH,YDFGH)
                 RIntort(3,IDFGH,YDFGH)=RGGJ
                 EepsL(1,IDFGH,YDFGH)=Eesp(1,IDFGH,YDFGH)
                 EepsL(2,IDFGH,YDFGH)=Eesp(2,IDFGH,YDFGH)
                 EepsL(3,IDFGH,YDFGH)=Eesp(3,IDFGH,YDFGH)
                 WRITE(6,*) Eesp(1,IDFGH,YDFGH),Eesp(2,IDFGH,YDFGH),Eesp(3,IDFGH,YDFGH)
                 WRITE(6,*) Rint(1,IDFGH,YDFGH),Rint(2,IDFGH,YDFGH),RGGJ
                 WRITE(6,*)
 				 ! наыхи оюпюлерп
				 Eesp(1,IDFGH,YDFGH)=Eesp(2,IDFGH,YDFGH)
                 Eesp(2,IDFGH,YDFGH)=Eesp(3,IDFGH,YDFGH)
				 Rint(1,IDFGH,YDFGH)=Rint(2,IDFGH,YDFGH)
                 Rint(2,IDFGH,YDFGH)=RGGJ
                 
				

	          ENDIF
	       ENDDO
	    ENDDO 
	ENDIF  
	
	
     



	
	  
			
    ! тнплхпсел люяяхбш тсмйжхи Pz х Qz дкъ бяеу нанкнвей (кнйюкэмсч х мекнйюкэмсч вюярэ онремжхюкю)
	DO IXKLS=1,IS   
	   call MDRUM(IXKLS,2,N,RfunNX,R8)
       
	   ! пюявер опълни х налеммни вюярх онремжхюкю
       ! пюявер опълни вюярх
       call  Calcul_Vnl(M,N,IXKLS,Z,H,RAL,NG,NU,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,RPOT,P4,Pznl)
 
       ! пюявер налеммни вюярх
       call  Calcul_Xnl(M,M,N,IXKLS,IS,H,NG,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,RfunNX,Qznl)
	   
	   DO IDPS=1,N
          Pz(IXKLS,IDPS)=Pznl(IDPS)
	      Qz(IXKLS,IDPS)=Qznl(IDPS)
	   ENDDO        
    ENDDO

     
   

    
    ! жхйк дкъ онксвемхе лмнфхрекеи кюцпюмфю 
	EepsL=0.D0
    DO IDFGH=1,IS
	   call MDRUM(IDFGH,2,N,RfunNX,R8)
	   DO YDFGH=1,IS
          IF(NN(IDFGH).NE.NN(YDFGH).AND.L(IDFGH).EQ.L(YDFGH)) THEN
             IF(Q(IDFGH).NE.Q(YDFGH)) THEN
			    call MDRUM(YDFGH,2,N,RfunNY,R8)
	            RGGJ=CoffLagrang(IDFGH,YDFGH,RfunNX,RfunNY,Pz,Qz,N,H,R1,RO1X) 
		        EepsL(3,IDFGH,YDFGH)=FLOAT(Q(YDFGH)*Q(IDFGH))*RGGJ/FLOAT(Q(YDFGH)-Q(IDFGH)) 
             ENDIF                 
	      ENDIF
	   ENDDO
	ENDDO 
    



   !WRITE(6,*) 'EXIT',INDEXEXIT,Nsum
   ! онксвеммше лмнфхрекх кюцпюмфю
   DO IDFGH=1,IS
      DO YDFGH=1,IS 
         EepsLagran(IDFGH,YDFGH)=EepsL(3,IDFGH,YDFGH)
         Eesp(3,IDFGH,YDFGH)=EepsL(3,IDFGH,YDFGH)
      ENDDO
   ENDDO
     


  
    ! сдюкъел хг оюлърх люяяхбш	
	deallocate(INDEXKLUZ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "INDEXKLUZ" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(R8x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "R8x" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(TETAx,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "TETAx" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
	deallocate(E3x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "E3x" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(E4x,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "E4x" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
	deallocate(Pz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Pz" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(Qz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Qz" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
    deallocate(Pznl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Pznl" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(Qznl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "Qznl" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
	deallocate(RIntort,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RIntort" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(EepsL,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "EepsL" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
    deallocate(RfunNX,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RfunNX" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
  	deallocate(RfunNY,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'CalculationEpsLagrang'
       write(*,*) 'THE FILE "RfunNY" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif

    return
  end  subroutine CalculationEpsLagrangALFA


  ! ОНДОПНЦПЮЛЛЮ ОПЪЛНЦН ПЮЯВЕРЮ ЛМНФХРЕКЕИ КЮЦЮМФЮ

  real(8) function CoffLagrang(IX,IY,RX,RY,Vnl,Xnl,N,H,R1,RO1X)
     implicit none
     integer::IX,IY,N
	 real(8)::H
	 real(8),dimension(:)::RX,RY,R1,RO1X
     real(8),dimension(:,:)::Vnl,Xnl
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer::ierr,IIXD
	 real(8),allocatable,dimension(:)::RfunQ
    
	 !бшдекъел оюлърэ дкъ люяяхбнб
	 allocate(RfunQ(N),stat=ierr)
     if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'MEMORY ON THE FILE "RfunQ" IS NOT SELECTED'
       READ(*,*)
	   stop 
     endif 
     
     RfunQ=0.D0
	 DO IIXD=1,N
        RfunQ(IIXD)=2.D0*(Vnl(IX,IIXD)-Vnl(IY,IIXD))*RX(IIXD)*RY(IIXD)+Xnl(IX,IIXD)*RY(IIXD)-Xnl(IY,IIXD)*RX(IIXD)
        RfunQ(IIXD)=RfunQ(IIXD)/R1(IIXD)
	 ENDDO

	 CoffLagrang=SIMP(RfunQ,N,4,H,RO1X)
  	
	 ! сдюкемхе люяяхбнб хг оълърх 
	 deallocate(RfunQ,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'MAING'
        write(*,*) 'THE FILE "RfunQ" IS NOT REMOVED FROM MEMORY'
        READ(*,*)
	    stop 
     endif
  end function  CoffLagrang

 ! ондопнцпюллю пеьемхъ спюбмемхъ уюрпх-тнйю дкъ нанкнвйх хг йнмтхцспюжхх
 ! нохяюмхе ббндхлшу оюпюлерпнб
 ! M-мнлеп хрепюжхх
 ! IB-мнлеп нанкнвйх
 ! IS-вхякн нанкнвей
 ! N-вхякн рнвей (вхякн гмювемхи тсмйжхх)
 ! NG-вхякн налеммшу хмрецпюкнб
 ! NU-онкмне вхякн хмрецпюкнб
 ! NLST-оюпюлерп сйюгшбючыхи мю рхо ялеьхбюмхъ
 ! IPOT-оюпюлерп оевюрх пегскэрюрнб пюяверю
 ! NonEepsLag-оюпюлерп сйюгшбючыхи хяонкэгнбюрэ медхюцнмюкэмше лмнфхрекх кюцпюмфю б пюявере хкх мер
 ! Z-гюпъд ъдпю
 ! H-ьюц
 ! EPS-рнвмнярэ
 ! ALFA-оюпюлерп "мнбни" оепелеммни
 ! BET-оюпюлерп "мнбни" оепелеммни
 ! TETMAX-оюпюлерп сйюгшбючыхи мю бекхвхмс люйяхлюкэмнцн нрйкнмемхъ нр опедшдсыеи хрепюжхх    
 ! RAL-пюяярнъмхе дн кхцюмдю
 ! NN(IS)-люяяхб цкюбмшу йбюмрнбшу вхяек йнмтхцспюжхх
 ! L(IS)-люяяхб нпахрюкэмшу лнлемрнб йнмтхцспюжхх
 ! Q(IS)-люяяхб вхяек гюонкмемхъ йнмтхцспюжхх
 ! KK()-люяяхб гмювемхи лскэрхонкэмнярх   
 ! IS1()-люяяхбш гмювемхи мнлепнб нанкнвей яннрберярбсчыху KK()  
 ! IS2()-люяяхбш гмювемхи мнлепнб нанкнвей яннрберярбсчыху KK()  
 ! E3(IS)-люяяхб гмювемхи ндмнщкейрпнммшу щмепцхи опх дюммни хрепюжхх 
 ! E4(IS)-люяяхб гмювемхи ндмнщкейрпнммшу щмепцхи (опх дюммни хрепюжхх)
 ! TETA(IS)-люяяхб нрйкнмемхи 
 ! GAM()-люяяхб йнщттхжхемрнб йскнмнбяйнцн бгюхлндеиярбхъ
 ! R1(N)-люяяхб гмювемхи пюдхсяю
 ! P4(N)-люяяхб гмювемхи онремжхюкю йнмтхцспюжхх мю оепбнл ьюце (мскхбне опхакхфемхе)  
 ! R8(IS*(N+2))-люяяхб гмювемхи пюдхюкэмшу вюяреи бнкмнбшу тсмйжхи йнмтхцспюжхх
 ! RAB1(N)-бяонлнцюрекэмши люяяхб
 ! RAB2(N)-бяонлнцюрекэмши люяяхб
 ! RPOT(IS,N)-гмювемхе бмеьмецн онремжхюкю дкъ йюфдни нанкнвйх
 ! EepsLagran(IS,IS)-люяяхб медхюцнмюкэмшу лмнфхрекеи кюцпюмфю
 subroutine MAING(M,IB,IS,N,NG,NU,NLST,IPOT,NonEepsLag,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,EepsLagran)
    implicit none
    integer::M,IB,IS,N,NG,NU,NLST,IPOT,NonEepsLag      
   	real(8)::Z,H,EPS,ALFA,BET,TETMAX,RAL
    integer,dimension(:)::NN,L,Q,KK,IS1,IS2
    real(8),dimension(:)::E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X 
    real(8),dimension(:,:)::RPOT,EepsLagran
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::NI,LI,QI
	integer::ierr,IXXD,J,MAXCNV
	real(8)::EpsLagram,T,C 
    real(8),allocatable,dimension(:)::Pz,Qz,RfunQ
    real(8),dimension(6002)::R2,R5    
	 
    !бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RfunQ(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'MAING'
      write(*,*) 'MEMORY ON THE FILE "R2" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(Pz(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'MAING'
      write(*,*) 'MEMORY ON THE FILE "Pz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
	allocate(Qz(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'MAING'
      write(*,*) 'MEMORY ON THE FILE "Qz" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
  
    
	! гюмскъел оепед пюявернл
    Pz=0.D0
    Qz=0.D0
	RfunQ=0.D0

   
    
   
    call MDRUM(IB,2,N,R5,R8)
    ! гюохяшбюел опедьеярбсчысч хрепюжхч
	DO J=1,N+2
       R2(J)=R5(J)
    ENDDO
    	  
  
    NI=NN(IB)
    LI=L(IB)
    QI=Q(IB)

  
  ! WRITE(6,*) 'ETAP 2',IB 
   ! ондопнцпюллю пюяверю Pz х Qz -тсмйжхи менаундхлшу дкъ пеьемхъ дхттепемжхюкэмнцн спюбмемхъ
   call Calcul_Pz_Qz_nl(M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,Pz,Qz) 
     

  ! WRITE(6,*) 'ETAP 3',IB,M,NonEepsLag 
   ! опнбепъел менаундхлнярэ днаюбйх якюцюелшу я медхюцнмюкэмшлх лмнфхрекълх кюцпюмфю
   IF(M.NE.1.AND.NonEepsLag.NE.0) THEN
     ! WRITE(6,*) 'ETAP 4',IB 
	  DO IXXD=1,IS
         IF(NN(IXXD).NE.NN(IB).AND.L(IB).EQ.L(IXXD)) THEN
  	     	call MDRUM(IXXD,2,N,RfunQ,R8)
			EpsLagram=EepsLagran(IB,IXXD)/FLOAT(Q(IB))
			!WRITE(6,*) 'EPD',M,IXXD,EpsLagram
            DO J=1,N
               Qz(J+1)=Qz(J+1)+EpsLagram*RfunQ(J+2)/(DSQRT(RfunQ(2))*RO1X(J)**2)
            ENDDO 
		 ENDIF 
      ENDDO 
   ENDIF  
  
  !IF(IB.EQ.3) THEN 
  !   DO J=1,N
  !      WRITE(800+M,*) R1(J),Pz(J),Qz(J)
  ! 	 ENDDO 
  !ENDIF 
 

   
  ! пеьюел дхттепемжхюкэмне спюбмемхе лернднл бейрнпмни опнцнмйх мслепнбю  
  ! Pz+EPS=R6-люяяхб гмювемхи тсмйжхх  P (КНЙЮКЭМНИ ВЮЯРХ) Б СПЮБМЕМХХ   Y''=P*Y+q
  ! Qz-люяяхб гмювемхи тсмйжхх  q (МЕ КНЙЮКЭМНИ ВЮЯРХ) Б СПЮБМЕМХХ Y''=P*Y+q 
 ! WRITE(6,*) 'ETAP 5',IB 
   call PrimitiveDifferentialEquation(IB,NI,LI,N,M,NLST,MAXCNV,H,TETMAX,EPS,E3,E4,TETA,R2,Qz,R5,Pz,RO1X,T,C)      
 ! WRITE(6,*) 'ETAP 7',IB       
 
     

555 FORMAT(' N=',I4,'  Shell=',I2,'  E=',F10.4,' Teta=',D10.3,' Tener=',D10.3,' Imax=',I6)
556 FORMAT(' N=',I4,'    Shell=',I2,'    E=',F10.4,'    Teta=',D10.3)
      
    EepsLagran(IB,IB)=-R2(1)	
	IF(IPOT.EQ.3) write(6,555) M,IB,R2(1),T,C,MAXCNV
  
      
	write(*,556) M,IB,R2(1),T
    ! гюохяшбюел тсмйжхч
	call MDRUM(IB,7,N,R5,R8)
    


	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(RfunQ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "RfunNN" IS NOT REMOVED FROM MEMORY'
        READ(*,*)
	   stop 
    endif
  
	deallocate(Pz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "Pz" IS NOT REMOVED FROM MEMORY'
       READ(*,*)  
	   stop 
    endif
   
	deallocate(Qz,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'MAING'
       write(*,*) 'THE FILE "Qz" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
   
    return
  end subroutine MAING




 ! ондпнцпюллю тнплхпсер тсмйжхх Pz(ro) х Qz(ro) дкъ пеьемхъ дхттепемжхюкэмнцн спюбмемхъ дкъ дюммни нанкнвйх 
  subroutine Calcul_Pz_Qz_nl(M,IB,IS,N,NG,NU,NLST,IPOT,Z,H,EPS,ALFA,BET,TETMAX,RAL,NN,L,Q,KK,IS1,IS2,E3,E4,TETA,GAM,R1,P4,R8,RO1X,RO2X,RO3X,RPOT,Pz,Qz) 
    implicit none
    integer::M,IB,IS,N,NG,NU,NLST,IPOT 
   	real(8)::Z,H,EPS,ALFA,BET,TETMAX,RAL
    integer,dimension(:)::NN,L,Q,KK,IS1,IS2
    real(8),dimension(:)::E3,E4,TETA,GAM,R8,R1,P4,RO1X,RO2X,RO3X,Pz,Qz 
    real(8),dimension(:,:)::RPOT
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::NI,LI,QI
	integer::ierr,IISD,IXXD,J,MAXCNV
	real(8)::P,T,B,C,D,RO21,RO31,ROIF
    real(8),allocatable,dimension(:)::Vnl,Xnl,RfunQ
      
	 
    !бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RfunQ(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Pz_Qz_nl'
      write(*,*) 'MEMORY ON THE FILE "R2" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(Vnl(N),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Pz_Qz_nl'
      write(*,*) 'MEMORY ON THE FILE "Vnl" IS NOT SELECTED'
      READ(*,*)
	  stop 
    endif 
    allocate(Xnl(N),stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_nl'
       write(*,*) 'MEMORY ON THE FILE "Xnl" IS NOT SELECTED'
       READ(*,*)
	   stop 
    endif 
    
	! гюмскъел оепед пюявернл
    Pz=0.D0
    Qz=0.D0
	RfunQ=0.D0
	Vnl=0.D0
    Xnl=0.D0 

   
   
    
   
    call MDRUM(IB,2,N,RfunQ,R8)
    
    	  
  
    NI=NN(IB)
    LI=L(IB)
    QI=Q(IB)

   ! пюявер опълни х налеммни вюярх онремжхюкю
  
   ! пюявер опълни вюярх
   call  Calcul_Vnl(M,N,IB,Z,H,RAL,NG,NU,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,RPOT,P4,Vnl)
 
   ! пюявер налеммни вюярх
   call  Calcul_Xnl(100+M*10+IB,M,N,IB,IS,H,NG,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,RfunQ,Xnl)
        
 
    
   ! тнплхпсел кнйюкэмсч вюярэ ноепюрнпю тнйю  
   !T=BET*BET/4.D0
   !P=ALFA*BET
   D=FLOAT(LI*(LI+1))
   !DO J=1,N
   !   B=2.D0*R1(J)*Vnl(J)
   !   C=RAB1(J)*RAB1(J)
   !   Pz(J)=(D+B+(P*R1(J)+T)*C)*C
   !ENDDO
   DO J=1,N
      RO21=(RO2X(J)/RO1X(J))**2
	  RO31=RO3X(J)/RO1X(J)
	  ROIF=0.5D0*(RO31-1.5D0*RO21)/RO1X(J)**2
	  B=2.D0*Vnl(J)/(R1(J)*RO1X(J)*RO1X(J))
      C=D/(R1(J)*R1(J)*RO1X(J)*RO1X(J))
      Pz(J)=B+C+ROIF
   ENDDO 

  
   ! тнплхпсел мекнйюкэмсч (налеммсч) вюярэ онремжхюкю
   IF(M.NE.1) THEN
      !DO J=1,N
      !   Qz(J+1)=R1(J)*Xnl(J)*RAB1(J)**2
      !ENDDO
      DO J=1,N
         Qz(J+1)=Xnl(J)/(R1(J)*RO1X(J)*RO1X(J))
      ENDDO
   ENDIF  
   !100 FORMAT(4(1X,F20.15))
   !DO J=1,N    
   !   WRITE(100+M*10+IB,100) R1(J),Pz(J),Xnl(J)/(R1(J)*DSQRT(RO1X(J))),RfunQ(J+2)/DSQRT(RO1X(J))
   !ENDDO   


	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(RfunQ,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_nl'
       write(*,*) 'THE FILE "RfunNN" IS NOT REMOVED FROM MEMORY'
        READ(*,*)
	   stop 
    endif
    deallocate(Vnl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_nl'
       write(*,*) 'THE FILE "Vnl" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif
    deallocate(Xnl,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Pz_Qz_nl'
       write(*,*) 'THE FILE "Xnl" IS NOT REMOVED FROM MEMORY'
       	READ(*,*)  
	   stop 
    endif   
    return
  end subroutine Calcul_Pz_Qz_nl

     

  ! ондопнцпюллю пюяверю опълни (кнйюкэмни) вюярх онремжхюкю 
  ! нохяюмхе оюпюлерпнб ондопнцпюллш
  subroutine Calcul_Vnl(M,N,IB,Z,H,RAL,NG,NU,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,RPOT,P4,Vnl)
    implicit none  
    integer::M,N,IB,NG,NU
	integer,dimension(:)::L,Q,KK,IS1,IS2
	real(8)::Z,H,RAL
    real(8),dimension(:)::GAM,R1,RO1X,R8,P4,Vnl
	real(8),dimension(:,:)::RPOT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::LI,QI,IW,IA,J,ierr,MM,KM,NF,K,II
	real(8)::GAMA
    real(8),allocatable,dimension(:)::R5,R6
    

	IF(M.EQ.1) THEN
	   DO J=1,N
	      Vnl(J)=P4(J)
       ENDDO 
	   return
	ENDIF


    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R5(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Vnl'
      write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
      stop 
    endif 
    allocate(R6(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Vnl'
      write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
      stop 
    endif 
    
	! гюмскъел оепед пюявернл 
    R5=0.D0
    R6=0.D0

	LI=L(IB)
	QI=Q(IB)
	

	DO J=1,N
       Vnl(J)=-Z
    ENDDO
         
    MM=0
    KM=2*LI
    NF=NG+1
	
	! тнплхпсел щкейрпнммши онремжхюк (опълюъ вюярэ кнйюкэмюъ)
    DO K=MM,KM,2
       II=0
       DO J=1,N
          R6(J)=0.D0
       ENDDO
            
       DO IW=NF,NU
          IA=0
          IF(K.NE.KK(IW).OR.DABS(GAM(IW)).LE.1.D-8) THEN 
  	         CYCLE
          ENDIF
		  IF(IB.EQ.IS1(IW)) THEN 
		     IA=IS2(IW)
          ENDIF
          IF(IB.EQ.IS2(IW)) THEN 
		     IA=IS1(IW)
          ENDIF
		  IF(IA.EQ.0) THEN 
		     CYCLE
          ENDIF
			   
		  II=1
               
		  call MDRUM(IA,2,N,R5,R8)
               
		  
          IF(IB.EQ.IA) THEN 
		      GAMA=2.D0
             ELSE
              GAMA=1.D0
          ENDIF 
          
		  GAMA=GAM(IW)*GAMA/QI/R5(2)
               
		  DO J=1,N
             R6(J)=R6(J)+R5(J+2)**2*GAMA
          ENDDO	    
       ENDDO 
       
	   IF(II.EQ.0) THEN 
	       CYCLE
       ENDIF
	   
	   call POT(DBLE(K),N,H,R6,R1,RO1X)
       DO J=1,N
          Vnl(J)=Vnl(J)-R6(J)
       ENDDO
    ENDDO 

    ! днаюбкъел бмеьмхи онремжхюк
    IF(RAL.GT.1.D-3) THEN 
       DO J=1,N
          Vnl(J)=Vnl(J)+RPOT(IB,J)  
       ENDDO
    ENDIF
   
  




    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R5,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Vnl'
       write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R6,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Vnl'
       write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
   return
  end subroutine Calcul_Vnl



  ! ондопнцпюллю пюяверю налеммни (ме кнйюкэмни) вюярх онремжхюкю 
  ! нохяюмхе оюпюлерпнб ондопнцпюллш
  subroutine Calcul_Xnl(IIXX,M,N,IB,IS,H,NG,L,Q,KK,IS1,IS2,GAM,R1,RO1X,R8,R2,Xnl)
    implicit none  
    integer::IIXX,M,N,IB,NG,IS
	real(8)::H
	integer,dimension(:)::L,Q,KK,IS1,IS2
    real(8),dimension(:)::GAM,R1,RO1X,R8,R2,Xnl
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::LI,QI,IW,IA,J,ierr,II
	real(8)::GAMA,P,SSUM
    real(8),allocatable,dimension(:)::R5,R6

    100 FORMAT(4(1X,F20.15))
	! M-мнлеп хреппюжхх
	IF(M.EQ.1) THEN
	   Xnl=0.D0 
	   RETURN 
	ENDIF

    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R5(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Xnl'
      write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
      stop 
    endif 
    allocate(R6(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'Calcul_Xnl'
      write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
      stop 
    endif 

	LI=L(IB)
	QI=Q(IB)
    
	! гюмскъел оепед пюявернл 
    R5=0.D0
    R6=0.D0
    Xnl=0.D0

	! тнплхпсел налеммши щкейрпнммши онремжхюк (мекнйюкэмши)
    DO IW=1,NG
       IA=0
       IF(IB.EQ.IS1(IW)) THEN
	      IA=IS2(IW)
       ENDIF
       IF(IB.EQ.IS2(IW)) THEN
	      IA=IS1(IW)
       ENDIF
       IF(IA.EQ.0) THEN 
	      CYCLE 
       ENDIF
       call MDRUM(IA,2,N,R5,R8)
       DO J=1,N
          R6(J)=R5(J+2)*R2(J+2)
       ENDDO 
       call POT(DBLE(KK(IW)),N,H,R6,R1,RO1X)
	   GAMA=2.D0*GAM(IW)/QI/R5(2)
       DO J=1,N
          Xnl(J)=Xnl(J)+GAMA*R6(J)*R5(J+2)
	      !WRITE(IIXX,100) R1(J),-GAMA*R6(J)/R1(J),R2(J+2)*R5(J+2)/RO1X(J)**2
	   ENDDO 
	   !SSUM=0.D0
	   !DO J=1,N
	   !   SSUM=SSUM+R2(J+2)*R5(J+2)/RO1X(J)**2
       !ENDDO
       !WRITE(6,*) IIXX,SSUM*H/DSQRT(R5(2)*R2(2))

    ENDDO
    

    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R5,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Xnl'
       write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY'
       stop 
    endif
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R6,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'Calcul_Xnl'
       write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
    return
  end subroutine Calcul_Xnl



  ! ондопнцпюллю пеьемхъ дхттепемжхюкэмнцн спюбмемхъ (опнцнмйю+лернд мслепнбю)
  ! онксвемхе бнкмнбни тсмйжхх нанкнвйх   
  ! IB-мнлеп нанкнвйх б йнмтхцспюжхх
  ! NI-цкюбмне йбюмрнбне вхякн
  ! LI-нпахрюкэмши лнлемр
  ! N-вхякн рнвей
  ! M-мнлеп хрепюжхх
  ! NLST-йкчв сйюгшбючыхи мю опнжеяя опхлеьхбюмхъ опедшдсыеи йнмтхцспюжхх й дюммни
  ! MAXCNV-
  ! H-ьюц 
  ! TETMAX-люйяхлюкэмне нрйкюмемхе опедшдсыеи хрепюжхх нр дюммни 
  ! EPS-рнвмнярэ
  ! E3(IS)-люяяхб ндмнщкейрпнммшу щмепцхи (хрепюжхъ) 
  ! E4(IS)-люяяхб ндмнщкейрпнммшу щмепцхи (хрепюжхъ)
  ! TETA-люяяхб нрйкюмемхи йюфдни тсмйжхх 
  ! R4(N+2)-люяяхб гмювемхи Q(ro)-ТСМЙЖХХ Б УНДЪЫЕЕ Б СПЮБМЕМХЕ Б СПЮБМЕМХХ Y''=P*Y+Q-мекнйюкэмюъ вюярэ онремжхюкю 
  ! R3(N+2)-люяяхб гмювемхи P(ro)-ВЮЯРХ КНЦЮКЭМНИ ТСМЙЖХХ P БУНДЪЫЕИ Б СПЮБМЕМХЕ Y''=P*Y+Q-кнйюкэмюъ вюярэ онремжхюкю
  ! RAB2(N)-бяонлнцюрекэмши люяяхб дкъ пюяверю  
  ! T-оюпюлерпш сйюгшбючыхе мю йювеярбн янцкюянбюмхъ  
  ! C-оюпюлерпш сйюгшбючыхе мю йювеярбн янцкюянбюмхъ  
  ! R2(N+2)-люяяхб гмювемхи тсмйжхх дюммни нанкнвйх (ярюпни) оняке пюанрш ондопнцпюллш гюохяюмю мнбюъ тсмйжхъ дюммни хрепюжхх
  ! R5(N+2)-люяяхб гмювемхи тсмйжхх дюммни нанкнвйх (ярюпни) оняке пюанрш ондопнцпюллш гюохяшбюеряъ мнбюъ тсмйжхъ я опхлеяэч ярюпни тсмйжхх дкъ якедсчыеи хрепюжхх 
  subroutine PrimitiveDifferentialEquation(IB,NI,LI,N,M,NLST,MAXCNV,H,TETMAX,EPS,E3,E4,TETA,R2DDD,R4,R5,R3,RO1X,T,C)    
    implicit none
	integer::IB,NI,LI,N,M,NLST,MAXCNV
	real(8)::H,TETMAX,EPS
	real(8),dimension(:)::E3,E4,TETA,R2DDD,R4,R5,R3,RO1X
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::MM,J,NF,IA,IY,ierr,IIIXXX
    integer::M1,M2,IX,KSI,JJ,ISI,ISI7
	real(8)::A1,B1,A2,OMEGA,ET,GAMA,RO1,E1,E2,EZR
	real(8)::P,T,D,B,C,U,V,TAU1,TAU2
    real(8),allocatable,dimension(:)::R6
    
    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R6(N+2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'PrimitiveDifferentialEquation'
      write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
      stop 
    endif 

    ET=0.D0
	IY=1
   



    ! гюмскъел оепед пюявернл
    R6=0.D0
  !  IIIXXX=0
    ! щрюо 1. пЕЬЮЕЛ ДХТТЕПЕМЖХЮКЭМНЕ СПЮБМЕМХЕ БРНПНЦН ОНПЪДЙЮ (ОНКСВЮЕЛ БНКМНБСЧ ТСМЙЖХЧ)
101 P=H*H/12.D0
    D=LI+0.5D0
    M1=0
    M2=0
    IX=1
 !   IIIXXX=IIIXXX+1
 !   WRITE(*,*) IB,IIIXXX
    ! пюгахбюел хмрепбюк хмрецпхпнбюмхъ мю рпх свюярйю [0,M1],[M1,M2],[M2,N]
    DO J=1,N
       C=R3(J)+R2DDD(1)/RO1X(J)**2 !RAB2(J)
       R6(J+1)=C*P
       B=D*C
       D=C
       IF(B.GT.0) THEN
	      CYCLE 
       ENDIF
	   IF(IX.NE.2) THEN 
           M1=J
           IX=2
	      ELSE
           M2=J
       ENDIF
    ENDDO 

    
    IF(M1-M2) 14,13,19
 13 U=R6(2)
      
	DO J=1,N
       IF((R6(J+1)-U).LT.0.D0) THEN 
	      M1=J-2
          U=R6(J+1)
       ENDIF
    ENDDO 
      
	IF(U/P-5.0) 17,17,18
18  R2DDD(1)=R2DDD(1)*.8D0
    GOTO 101

17  M2=M1+4
    GOTO 14

19  R2DDD(1)=(R2DDD(1)+.1D0)*2.D0
    GOTO 101

14  R6(1)=0.D0
    R2DDD(2)=0.D0
    MM=M1-1
      
    DO  J=2,MM
        call AAB(J,R4,R6,A1,B1,A2,OMEGA,H)
        IF(J-2) 21,20,21
20      V=0.D0
22      U=V
        V=A2/(B1-A1*V)
        IF(DABS(V-U).GT.EPS) GOTO 22
        U=0.D0
        GOTO 23
21      B=1.D0/(B1-A1*V)
        V=A2*B
        U=(A1*U+OMEGA)*B
23      R6(J)=V
        R2DDD(J+1)=U
    ENDDO

 
    GAMA=1.D0
    RO1=U+V*GAMA
    KSI=0
    MM=M2-1

    DO J=M1,MM
       call AAB(J,R4,R6,A1,B1,A2,OMEGA,H)
       R2DDD(J+2)=GAMA
       V=RO1
       RO1=GAMA
       GAMA=(B1*RO1-A1*V-OMEGA)/A2
       IF(DSIGN(1.,RO1*GAMA).LT.0) THEN
  	      KSI=KSI+1
       ENDIF
    ENDDO

    ! нопедекъел вхякн сгкнб тсмйжхх
    IA=NI-LI-1

    IF(KSI.EQ.IA) THEN 
	   GOTO 26
    ENDIF
    IF(IY.NE.1) THEN 
	   GOTO 26
    ENDIF
    A1=.1D0*(KSI-IA)
102 R2DDD(1)=R2DDD(1)*(1.D0+A1/NI)
    GOTO 101

 26 MM=N-1
  
    DO JJ=1,MM
       J=N-JJ
       R6(J+3)=0.D0
       R2DDD(J+4)=0.D0
       IA=J
       IF(R6(J+1).LT.1.D0) THEN 
	      EXIT
       ENDIF
    ENDDO 



    DO JJ=M2,IA
       J=IA-JJ+M2
       call AAB(J,R4,R6,A1,B1,A2,OMEGA,H)
         
	   IF(IA-J) 30,31,30
 31       V=0.
 32       U=V
          V=A1/(B1-A2*V)
          IF(DABS(V-U).GT.EPS) GOTO 32
          U=0.D0
          GOTO 33
 30       B=1.D0/(B1-A2*V)
          V=A1*B
          U=(A2*U+OMEGA)*B
 33       R6(J+2)=V
          R2DDD(J+3)=U
    ENDDO




    A1=(GAMA-RO1*V-U)/H

    IF(DABS(A1).LT.DABS(EPS*GAMA).OR.DABS(ET-R2DDD(1)).LT.EPS) GOTO 103
    
	ET=R2DDD(1)
	IF(DSIGN(1.,A1)) 34,103,35
 35 TAU1=A1
    E1=R2DDD(1)
	GOTO(201,202,203),IY
 
 34 TAU2=A1
    E2=R2DDD(1)
    GOTO(201,202,203),IY

201 ISI=DSIGN(1.1,A1)
    IY=2
    GOTO 36

202 ISI7=DSIGN(1.1,A1)
    IF(ISI.NE.ISI7) GOTO 203
 36 A1=DSIGN(0.1,A1*(-1)**(NI-LI))
    GOTO 102

203 R2DDD(1)=(TAU2*E1-TAU1*E2)/(TAU2-TAU1)
    IY=3
    GOTO 101

103 MM=M1-1
    U=EPS*EPS
    DO 37 JJ=1,MM
       J=M1-JJ
       IF(DABS(R2DDD(J+3))-U) 91,91,92
 91    R2DDD(J+2)=0.D0
       GOTO 37
 92    R2DDD(J+2)=R2DDD(J+3)*R6(J)+R2DDD(J+1)
 37 CONTINUE
    
	MM=M2-1
    NF=N-1
      
	DO 38 J=MM,NF
       IF(DABS(R2DDD(J+2))-U) 93,93,94
 93    R2DDD(J+3)=0.D0
       GOTO 38
 94    R2DDD(J+3)=R2DDD(J+2)*R6(J+3)+R2DDD(J+4)
 38 CONTINUE

      
	! щрюо 2. яПЮБМХБЮЕЛ ОНКСВЕММСЧ БНКМНБСЧ ТСМЙЖХЧ Я ОНКСВЕММНИ МЮ ОПЕДШДСЫЕЛ ЩРЮОЕ ДКЪ МЮУНФДЕМХЪ ЙНЩТТХЖХЕМРЮ ЯЛЕЬХБЮМХЪ 
	! нопедекемхе йнщттхжхемрю опхлеьхбюмхъ C
	T=1.D0-R5(1)/R2DDD(1)
    IF(NLST.EQ.1) T=0.D0
    IF(M-3) 39,40,40
 39 C=0.D0
    GOTO 41
 40 IF((-1)**M) 42,42,43
 42 C=0.5D0
    GOTO 41

 43 EZR=R2DDD(1)-E4(IB)-R5(1)+E3(IB)
    IF(DABS(EZR).LT.1.D-8) EZR=1.D0
    C=(R2DDD(1)-E4(IB))/EZR
    IF(C.LT.0.D0.OR.C.GT.1.D0) C=0.5D0

 41 E3(IB)=R5(1)
    E4(IB)=R2DDD(1)
    IF(NLST.EQ.1) C=0.D0
	! гюохяэ ндмнщкейрпнммнцн гмювемхъ щмепцхх "мнбни" дюммни хреппюжхх
    R5(1)=C*R5(1)+(1.D0-C)*R2DDD(1)
    MAXCNV=0

    ! гюохяэ "мнбни" тсмйжхх дюммни нанкнвйх
    DO J=1,N
       B=R5(J+2)-R2DDD(J+2)
       IF(DABS(B).GT.DABS(T)) MAXCNV=J
       IF(DABS(B).GT.DABS(T)) T=B
       R5(J+2)=C*R5(J+2)+(1.D0-C)*R2DDD(J+2)
    ENDDO
	
	
	TETA(IB)=T

    ! пюявер гмювемхъ мнплш онксвеммни тсмйжхх  	
	DO J=1,N
       R6(J)=R5(J+2)*R5(J+2)
    ENDDO

	R5(2)=SIMP(R6,N,0,H,RO1X)


	! пюявер гмювемхъ мнплш онксвеммни тсмйжхх  	
	DO J=1,N
       R6(J)=R2DDD(J+2)*R2DDD(J+2)
    ENDDO

	R2DDD(2)=SIMP(R6,N,0,H,RO1X)
      
	TETMAX=TETMAX+DABS(T)
    
	

	! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R6,stat=ierr)
    if(ierr/=0) then
       write(*,*) 'PrimitiveDifferentialEquation'
       write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY'
	   stop 
    endif
    
	return
  end subroutine PrimitiveDifferentialEquation


  subroutine AAB(J,R4,R6,A1,B1,A2,OMEGA,H)
    implicit none
    integer:: J
    real(8)::A1,B1,A2,OMEGA,H
    real(8),dimension(:)::R4,R6 

    A1=1.D0-R6(J)
    B1=2.D0+10.D0*R6(J+1)
    A2=1.D0-R6(J+2)
    OMEGA=(R4(J+2)+10.D0*R4(J+1)+R4(J))*H*H/12.D0
    return
  end subroutine AAB
  
    
  integer function ISIGNUM(X)
    implicit none
	real(8)::X
	IF(X.EQ.0.D0) THEN
       ISIGNUM=0
   	ENDIF
	IF(X.GT.0.D0) THEN
	   ISIGNUM=1
	  ELSE
	   ISIGNUM=-1
	ENDIF    
	     
   
  end function  ISIGNUM

 
 ! оНДПНЦПЮЛЛЮ СМХРЮПМНЦН ОПЕНАПЮГНБЮМХЪ ТСМЙЖХИ АЮГХЯЮ
 subroutine UnitaryTransformation(IS,N,H,R8,RO1X,EepsLagran) 
   use diagonal,only:DIAG
   implicit none
   integer::IS,N
   real(8)::H
   real(8),dimension(:)::R8,RO1X
   real(8),dimension(:,:)::EepsLagran
   real(8),allocatable,dimension(:)::RValue,FunD,FunX 
   real(8),allocatable,dimension(:,:)::RSnorm,RVector,Rfun
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer::IXCV,JSDF,ISS,ierr
   real(8)::RcoffG
   !бшдекъел оюлърэ дкъ люяяхбнб
   allocate(RSnorm(IS,IS),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "RSnorm" IS NOT SELECTED'
	  stop 
   endif
   allocate(RValue(IS),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "RValue" IS NOT SELECTED'
	  stop 
   endif
   allocate(RVector(IS,IS),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "RVector" IS NOT SELECTED'
	  stop 
   endif
   allocate(Rfun(IS,N+2),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "Rfun" IS NOT SELECTED'
	  stop 
   endif
   allocate(FunD(N+2),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "FunD" IS NOT SELECTED'
	  stop 
   endif
   allocate(FunX(N),stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'MEMORY ON THE FILE "FunX" IS NOT SELECTED'
	  stop 
   endif
   RSnorm=0.D0
   
   DO IXCV=1,IS
      RSnorm(IXCV,IXCV)=1.D0
   ENDDO

   ! дхюцнмюкхгсел люрпхжс лмнфхрекеи кюцпюмфю
   call DIAG(IS,EepsLagran,RSnorm,RValue,RVector)
   

   WRITE(6,*) (RValue(IXCV),IXCV=1,IS)
   DO  JSDF=1,IS
       WRITE(6,*) (RVector(JSDF,IXCV),IXCV=1,IS)
   ENDDO
   ! гюохяшбюел тсмйжхх
   DO IXCV=1,IS
      call MDRUM(IXCV,2,N,FunD,R8)
      DO JSDF=1,N+2
         Rfun(IXCV,JSDF)=FunD(JSDF) 
	  ENDDO
   ENDDO

   ! тнплхпсел мнбше тсмйжхх
   DO IXCV=1,IS
      FunD=0.D0
      FunD(1)=RValue(IXCV)
      DO JSDF=1,IS
	     RcoffG=RVector(JSDF,IXCV)
         DO ISS=1,N
            FunD(ISS+2)=FunD(ISS+2)+RcoffG*Rfun(JSDF,ISS+2)
		 ENDDO
	  ENDDO
      ! пЮЯВЕР МНПЛХПНБНВМНЦН ЙНЩТТХЖХЕМРЮ
      DO ISS=1,N
         FunX(ISS)=FunD(ISS+2)*FunD(ISS+2)
      ENDDO

	  FunD(2)=SIMP(FunX,N,0,H,RO1X)
  
      call MDRUM(IXCV,7,N,FunD,R8)
   ENDDO
  



   



  ! сдюкемхе люяяхбнб хг оълърх 
  deallocate(FunD,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "FunD" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif 
   deallocate(FunX,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "FunX" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif 
  deallocate(Rfun,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "Rfun" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif 
   deallocate(RSnorm,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "RSnorm" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif
   deallocate(RValue,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "RValue" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif 
   deallocate(RVector,stat=ierr)
   if(ierr/=0) then
      write(*,*) 'UnitaryTransformation'
      write(*,*) 'THE FILE "RVector" IS NOT REMOVED FROM MEMORY' 
      stop 
   endif

   return
 end subroutine UnitaryTransformation
 
  real(8) function EMV2(LF,I,N,H,RAN,R8,R1,RAB1,RAB2,RO1X)
    implicit none
    integer::LF,I,N,J,I1,N0,L1,J1
    real(8)::H,RAN,P2,P12,P112
    real(8),dimension(:)::R8,R1,RAB1,RAB2,RO1X
 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::ierr
	real(8),allocatable,dimension(:)::R51,RP,G,G1,R5,R6
    
    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(R51(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "R51" IS NOT SELECTED'
	   stop 
	endif
    allocate(RP(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "RP" IS NOT SELECTED'
	   stop 
	endif
    allocate(G(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "G" IS NOT SELECTED'
	   stop 
	endif
	allocate(G1(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "G1" IS NOT SELECTED'
	   stop 
	endif
	allocate(R5(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
	   stop 
	endif
	allocate(R6(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'EMV2'
	   write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	   stop 
	endif
	! гюмскъел оепед пюявернл  
	R51=0.D0
    RP=0.D0
    G=0.D0
    G1=0.D0
    R5=0.D0
    R6=0.D0
		
      

     P2=0.D0
     P12=0.D0
     P112=0.D0
     J=N-2
    
     DO I1=1,N
        RP(I1)=R8(I+I1)*DSQRT(R1(I1)*RAB1(I1))
     ENDDO

	 DO I1=2,N
        IF(R1(I1).GT.RAN) THEN 
		   EXIT
        ENDIF
        N0=I1
     ENDDO



     R5(1)=0.D0
     N0=N0+2
     L1=LF+1
      
	 DO I1=2,N0
        G(I1)=RP(I1)/R1(I1)**L1
     ENDDO
	
	
	 N0=N0-2
     DO I1=2,N0
        G1(I1)=(-3.*G(I1)+4.*G(I1+1)-G(I1+2))/H/2./RAB2(I1)**0.5
        R5(I1)=G1(I1)*R1(I1)**L1+G(I1)*R1(I1)**LF*L1
     ENDDO
	
	 DO I1=N0,J
        R5(I1)=(RP(I1-2)-RP(I1+2)-8.*(RP(I1-1)-RP(I1+1)))/12./H/DSQRT(RAB2(I1))
	 ENDDO 

     R5(N-1)=0.D0
     R5(N)=0.D0
      
	 IF(LF.NE.0) THEN 
        DO I1=1,N
           R6(I1)=(R8(I+I1)/R1(I1)**2)**2*RAB2(I1)
        ENDDO
	    P2=2.*SIMP(R6,N,3,H,RO1X)
        DO I1=1,N
           R6(I1)=R5(I1)**2/R1(I1)*RAB1(I1)
        ENDDO
  	    P12=2.*SIMP(R6,N,3,H,RO1X)  
	 ENDIF  

     R51(1)=0.D0
     J1=N-1
     N0=N0-1
     
	 DO I1=2,N0
        R51(I1)=(-3.*G1(I1)+4.*G1(I1+1)-G1(I1+2))/H/2./DSQRT(RAB2(I1))*R1(I1)**L1+2.*G1(I1)*L1*R1(I1)**LF+G(I1)*L1*LF*R1(I1)**(LF-1)
     ENDDO
 	 
	 DO I1=N0,J1
        R51(I1)=(R5(I1-2)-R5(I1+2)-8.*R5(I1-1)+8.*R5(I1+1))/12./H/DSQRT(RAB2(I1))
     ENDDO

	 R51(N)=0.D0
      
	 DO I1=1,N
        R6(I1)=R51(I1)**2*R1(I1)*RAB1(I1)
     ENDDO

	 P112=2.*SIMP(R6,N,3,H,RO1X)
     EMV2=P112+2.*LF*L1*P12+LF*L1*(LF*L1-6)*P2
     

	 ! сдюкемхе люяяхбнб хг оълърх 
	 deallocate(R51,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "R51" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
     deallocate(RP,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "RP" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
     deallocate(G,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "G" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
     deallocate(G1,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "G1" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
     deallocate(R5,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
      deallocate(R6,stat=ierr)
     if(ierr/=0) then
        write(*,*) 'EMV2'
        write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY' 
 	    stop 
     endif 
	 return
   end function EMV2


   ! ондопнцпюллю пюяверю вхякю йскнмнбяйхху йнщттхжхемрнб 
   subroutine ENTR(NF,NG,NU,IS,L)
     implicit none
     integer(4)::NDP,NUF,NF,NG,NU,IS,I,J
     integer(4),dimension(:)::L 
      
     NF=0
     NG=0
     DO I=1,IS
        DO J=I,IS
           IF((I-J).EQ.0) THEN
               NF=NF+1+L(J)
              ELSE
               NG=NG+1+MIN0(L(I),L(J))
           ENDIF
        ENDDO
     ENDDO 
      
	 NU=NF+2*NG
      
	   
     return
   end subroutine ENTR

      
	
	
	
   subroutine KOEFRAL(NUF,NDP,N,IS,RAL,ACY,ISD,KD,DEL,R1,RPOT)
     implicit none
	 integer::NUF,NDP,N,IS
     real(8)::RAL,ACY
	 integer,dimension(:)::ISD,KD
	 real(8),dimension(:):: DEL,R1
     real(8),dimension(:,:)::RPOT
     !!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer::I,J,IW,K,M,ierr
     real(8)::D,DD
     real(4),allocatable,dimension(:)::RVINT
     real(8),allocatable,dimension(:)::R2	   
	 
     
      
	 IF(RAL.LT.1.D-3) THEN
	    return
     ENDIF
     !бшдекъел оюлърэ дкъ люяяхбнб
	 allocate(RVINT(N+2),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'KOEFRAL'
	    write(*,*) 'MEMORY ON THE FILE "RVINT" IS NOT SELECTED'
	    stop 
	 endif  
     allocate(R2(N+2),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'KOEFRAL'
	    write(*,*) 'MEMORY ON THE FILE "R2" IS NOT SELECTED'
	    stop 
	 endif  
    
     ! гюмскъел оепед пюявернл
     RVINT=0.D0
     R2=0.D0


     M=N*IS
	 
	 IF(NDP.GT.0) REWIND 2

     IF(NDP) 40,40,41
  41 IW=NDP-1
     IF(IW) 40,40,43
  43 DO J=1,IW
        READ (2)
     ENDDO
 40  DO I=1,IS
        ! гюмскъел люяяхб
 	    RVINT=0.D0
        R2=0.D0
        IF(NDP) 46,46,47
  47    READ(2)(RVINT(J),J=1,N)
        DO J=1,N
           R2(J)=RVINT(J)*R1(J)
        ENDDO
  46    CONTINUE

        DO 73 J=1,N
           IF(RAL-R1(J)) 74,75,75
  74       D=1
           GOTO 73
  75       D=R1(J)/RAL
  73       R2(J)=R2(J)-ACY*D
      
	       DO 78 IW=1,NUF
              IF(ISD(IW)-I) 78,77,78
  77             K=KD(IW)
                 DD=DEL(IW)
              DO J=1,N
                 IF((RAL-R1(J)).EQ.0.D0.AND.(RAL-R1(J)).LT.0.D0) THEN
		             D=(RAL/R1(J))**K
                     R2(J)=R2(J)-DD*D
		         ENDIF
		         IF((RAL-R1(J)).GT.0.D0) THEN 
	                 D=(R1(J)/RAL)**(K+1)
                     R2(J)=R2(J)-DD*D
                 ENDIF 
              ENDDO 
  78       CONTINUE
      
           DO J=1,N
              RPOT(I,J)=R2(J)
           ENDDO

      ENDDO

     ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(RVINT,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "RVINT" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
      deallocate(R2,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "R2" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
 
      
	return
  end subroutine KOEFRAL

      


  ! ОНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ ЙНЩТТХЖЕМРНБ ЙСКНМНБЯЙНЦН БГЮХЛНДЕИЯРБХЪ 
  subroutine KOEF(NU,NG,IS,L,Q,NST,IS1,IS2,KK,GAM,NFK,NGK,IXX,IYY,KX,DCOFF)
    implicit none
	integer::NU,NG,IS,NFK,NGK,IX,IY,I,K,J,MAX,IK,IIIX,IIIY,KKK
    real(8)::D
    integer,dimension(:)::L,Q,NST,IS1,IS2,KK 
	integer,dimension(:,:)::IXX,IYY,KX
	real(8),dimension(:)::GAM
    real(8),dimension(:,:)::DCOFF
   
    
	IX=NG

    IF(NG.EQ.0) IX=1
    IY=0

    DO 11 I=1,IS
       DO 11 J=I,IS
          MAX=MIN0(L(I),L(J))*2
          IK=0
          DO 12 K=IK,MAX,2
             IX=IX+1
             IS1(IX)=I
             IS2(IX)=J
             KK(IX)=K
             IF(K) 14,13,14
             13 D=Q(I)
             IF(I.EQ.J) D=(D-1.)*.5
             GAM(IX)=-D*Q(J)
             IF(I.NE.J.AND.NST(I)*NST(J).EQ.0) GOTO 15
             GOTO 12
             14  D=0.
             IF(I.EQ.J) D=FKK1(L(I),K,DBLE(Q(I)))
             GAM(IX)=D
          12  CONTINUE

          15  IF(I.EQ.J) GOTO 11
          IK=IABS(L(I)-L(J))
          MAX=L(I)+L(J)
          DO 16 K=IK,MAX,2
             IY=IY+1
             IS1(IY)=I
             IS2(IY)=J
             KK(IY)=K
          16  GAM(IY)=GKK1(L(I),L(J),K,DBLE(Q(I)),DBLE(Q(J)))
      11  CONTINUE
      NU=IX
      IF(NG.NE.0) GOTO 95
         NG=1
         IS1(1)=1
         IS2(1)=1
         KK(1)=0
         GAM(1)=0.
      IF(IX.NE.0) GOTO 95
         NU=2
         IS1(2)=1
         IS2(2)=1
         KK(2)=0
         GAM(2)=0.
  95  CONTINUE


	IF(NGK.NE.0) THEN
         DO I=1,NGK
	      IIIX=IXX(1,I)
            IIIY=IYY(1,I)
		  KKK=KX(1,I)
            DO J=1,NG
               IF(IIIX.EQ.IS1(J).AND.IIIY.EQ.IS2(J).AND.KKK.EQ.KK(J)) THEN
		          GAM(J)=DCOFF(1,I)
               ENDIF
            ENDDO
         ENDDO
    ENDIF

      
	IF(NFK.NE.0) THEN
        DO I=1,NFK
           IIIX=IXX(2,I)
           IIIY=IYY(2,I)
		 KKK=KX(2,I)
           IK=NG+1
           DO J=IK,NU
              IF(IIIX.EQ.IS1(J).AND.IIIY.EQ.IS2(J).AND.KKK.EQ.KK(J)) THEN
                 GAM(J)=DCOFF(2,I)
              ENDIF
           ENDDO
        ENDDO
    ENDIF


	return
  end subroutine KOEF
 


     
      ! ОНДОПНЦПЮЛЛШ ПЮЯВЕРЮ СЦКНБШУ ВЮЯРЕИ ОПЪЛНЦН Х НАЛЕММНЦН ХМРЕЦПЮКЮ
      real(8) function  FKK1(L,K,Q)
        implicit real(8)(F,C,Q)
	  integer:: L,K
      FKK1=CK(L,L,K)**2*Q*(Q-1.)/DBLE((4*L+2)*(4*L+1))
      return
      end function  FKK1


      real(8) function  GKK1(L1,L2,K,Q1,Q2)
       implicit real(8)(G,C,Q)
	 integer::L1,L2,K 
      GKK1=CK(L1,L2,K)**2*Q1*Q2/DBLE((4*L1+2)*(2*L2+1))
      return
      end function GKK1


      real(8) function CK(L1,L2,K)
       implicit real(8)(C,F)
	 integer::L1,L2,K,M
      M=(L1+L2+K)/2
      CK=DSQRT(DBLE((2*L1+1)*(2*L2+1))*FACT(L2-L1+K)*FACT(L1-L2+K)*FACT(L2+L1-K)/FACT(L1+L2+K+1))*FACT(M)/(FACT(M-L1)*FACT(M-L2)*FACT(M-K))
      return
      end function CK

      real(8) function FACT(N)
        implicit none
	  integer::I,N 
	  real(8)::F 
	
	  F=1.D0
        do I=1,N
           F=F*DBLE(I)
        enddo
        IF(N.LT.0) F=0.D0
        FACT=F
      return
      end function FACT

    


      ! ОНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ БЯОНЛЮЦЮРЕКЭМШУ ЛЮЯЯХБНБ ДКЪ ХМРЕЦППНБЮМХЪ
      subroutine VAR(NtipSetky,N,H,ALFA,BET,GAMMA,Ral,RO1,RO2,R,RAB1,RAB2,RO1X,RO2X,RO3X)
        implicit none
        integer::I1,I,N,NtipSetky
	    real(8)::H,ALFA,BET,GAMMA,RO1,RO2,RT,ROT,A,Ral
	    real(8),dimension(:)::R,RAB1,RAB2,RO1X,RO2X,RO3X
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		real(8)::R2,R3,RR,RR2,RDELTA,RDELTA2,RDELTA3,XXDF
          
		   
        ! юрнлмюъ яерйю
        IF(NtipSetky.EQ.1) THEN
          RT=(RO2-BET*DLOG(RO2))/ALFA
          DO I1=1,N
             I=N-I1+1
             ROT=RO1+FLOAT(I)*H
  3          A=(ALFA*RT+BET*DLOG(RT)-ROT)/(ALFA*RT+BET)
             RT=RT*(1.D0-A)
             IF(DABS(A)-1.D-12) 2,2,3  
  2          R(I)=RT
             RAB1(I)=1.D0/(ALFA*RT+BET)
             RAB2(I)=(RT/(ALFA*RT+BET))**2
		     RO1X(I)=(ALFA*RT+BET)/RT
             RO2X(I)=-BET/RT**2
             RO3X(I)=2.D0*BET/RT**3
          ENDDO

        ENDIF
        
		! лнкейскъпмюъ яерйю
        IF(NtipSetky.EQ.2) THEN
          RO1=RO1 !-2.D0
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
      end subroutine VAR 






    ! ондопнцпюллю нясыеярбкъер гюохяэ пегскэрюрнб б тюик 
    subroutine INTER(IZAP,IS,N,IOUT,IPOT,IFUNC,H,NN,L,Q,R1,RO1X,R8)
      implicit none
      integer::IZAP,IS,N,IOUT,IPOT,IFUNC
	  real(8)::H
	  integer,dimension(:)::NN,L,Q
	  real(8),dimension(:)::R1,RO1X,R8
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::I,J,IW,IEND,NU,IBGN,ND,K,QQ,IIIFD,ierr
      real(8)::D,P,T,HH,B,C
	  integer,dimension(IS)::IPRN
	  character,dimension(10)::LB
	  real,allocatable,dimension(:)::R5OUT
      real(8),allocatable,dimension(:)::R5,R6
      real(8),allocatable,dimension(:,:)::RIntOrt
    	   
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(R5(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'INTER'
	     write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
	     stop 
	  endif  
      allocate(R5OUT(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'INTER'
	     write(*,*) 'MEMORY ON THE FILE "R5OUT" IS NOT SELECTED'
	     stop 
	  endif  
      allocate(RIntOrt(IS,IS),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'INTER'
	     write(*,*) 'MEMORY ON THE FILE "RIntOrt" IS NOT SELECTED'
	     stop 
	  endif  


      DATA LB/'s','p','d','f','g','h','i','j','k','l'/
 
  55  FORMAT(/' Written from record  ',I5,' Continue from record',I5)
  1   FORMAT(/' Delta RO = ',F7.5,'         F - Function  of P-type ')
  2   FORMAT(/'      R   ',4('  F(',I2,A1,')',1X),10X,4('  F(',I2,A1,')',1X)/)
  3   FORMAT(/' Total and partial charges in sphere of definite radius')
  4   FORMAT(/' Total Q  ',4('  Q(',I2,A1,')',1X),10X,4('  Q(',I2,A1,')',1X)/)
  5   FORMAT(5F9.3,10X,4F9.3)
  98  FORMAT(2X,F15.7,10(1X,E15.7))
 110  FORMAT(/12X,'r',14X,10(I2,A1,13X)/)
 120  FORMAT(/2X,'Radial part Int Ort',6X,10(I2,A1,12X)/)
 123  FORMAT(15X,I2,A1,10(F15.7,1X))
      
	  ! пюявер хмрецпюкнб норнцнмюкэмнярх
	  DO I=1,IS
	     DO J=1,I
            DO IIIFD=1,N  
               R5(IIIFD)=R8(IIIFD+(N+2)*(I-1)+2)*R8(IIIFD+(N+2)*(J-1)+2)
			ENDDO
            RIntOrt(I,J)=SIMP(R5,N,0,H,RO1X)
		 ENDDO
	  ENDDO 
	  ! бшдювю хмрецпюкнб нпрнцнмюкэмнярх
	  WRITE(6,120)(NN(I),LB(L(I)+1),I=1,IS)
	  DO I=1,IS 
	     WRITE(6,123) NN(I),LB(L(I)+1),(RIntOrt(I,J),J=1,I) 
	  ENDDO 

 
	  
	  
	  DO I=1,IS
         IPRN(I)=I
      ENDDO

	  D=0.D0
      P=0.D0
      T=0.D0

      IBGN=IZAP-1
      IF(IBGN) 50,51,52
  52  REWIND 1
  
      DO I=1,IBGN
         READ (1)
      ENDDO

  51  IEND=IZAP+IS
      IW=N+2
      
	  ! гЮОХЯЭ ТСМЙЖХИ Б ТЮИК
  	  DO I=1,IS
         call MDRUM(I,2,N,R5,R8)
         DO J=1,IW
            R5OUT(J)=SNGL(R5(J))
         ENDDO
         WRITE(1)(R5OUT(J),J=1,IW)
      ENDDO

      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R5,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
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

 
      WRITE(6,55) IZAP,IEND
  50  CONTINUE

      ! бшдювю б тюик пюдхюкэмшу вюяреи бнкмнбшу тсмйжхи 
      WRITE(70,110)(NN(I),LB(L(I)+1),I=1,IS)
      DO J=1,N
         WRITE(70,98) R1(J),(R8(J+(N+2)*(I-1)+2)/DSQRT(RO1X(J)),I=1,IS)
      ENDDO

      
    



      ! ************  OUTPUT ON THE PRINTER  ************
      
	  IW=IOUT

      IF(IW.LE.0) RETURN
      IF(IPOT+IFUNC.EQ.0) RETURN

      ! *********** TRANSFORMATION OF FUNCTION INTO P-TYPE ********
      ! опенапюгнбюмхе бнкмнбшу тсмйжхи б P-рхо
	  DO I=1,IS
         K=(I-1)*(N+2)+2
         DO J=1,N
            R8(J+K)=R8(J+K)/DSQRT(RO1X(J))
         ENDDO
      ENDDO

	  IW=IOUT
      
	  ! ******** PRINT ****************************
	  IF(IFUNC.NE.0) THEN
	     HH=IW*H
         WRITE(6,1) HH
         NU=0
  10     ND=NU+1
         NU=NU+8
        IF(NU.GT.IS) NU=IS
        WRITE(6,2)(NN(IPRN(I)),LB(L(IPRN(I))+1),I=ND,NU)
        DO J=1,N,IW
           WRITE(6,5) R1(J),(R8(J+(N+2)*(I-1)+2),I=ND,NU)
        ENDDO
		IF(NU.LT.IS) GOTO 10  
	  ENDIF

      IF(IPOT.EQ.0) RETURN

      ! онкмши гюпъд х гюпъдш нанкнвей, пюяопедекемхе б ятепе нопедекеммнцн пюдхсяю 
      WRITE(6,3)
       
      allocate(R6(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'INTER'
	     write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	     stop 
	  endif  



	  R6=0.D0
      B=H/12.D0
 
      DO I=1,IS
         K=(I-1)*(N+2)+2
         QQ=Q(I)
         D=0.D0
         P=0.D0
         T=0.D0
         DO J=1,N
            C=R8(J+K)**2*QQ/RO1X(J)
            D=(8.D0*P-D+5.D0*C)*B
            T=T+D
            D=P
            P=C
            R8(J+K)=T
            R6(J)=R6(J)+T
         ENDDO
      ENDDO

      NU=0

  35  ND=NU+1
      NU=NU+8
      IF(NU.GT.IS) NU=IS
      WRITE(6,4)(NN(IPRN(I)),LB(L(IPRN(I))+1),I=ND,NU)
   
      DO J=1,N,IW
         WRITE(6,5) R6(J),(R8(J+(N+2)*(I-1)+2),I=ND,NU)
      ENDDO
      IF(NU.LT.IS) GOTO 35

      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R6,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'INTER'
         write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif

      return
    end subroutine INTER

  




      subroutine POT(AK,N,H,R6,R1,RO1X)
        implicit none
	    integer::I,I1,N
	    real(8)::AK,D,P,T,B,C,H,RRT
	    real(8),dimension(:)::R6
        real(8),dimension(:)::R1,RO1X

        D=0.D0
        P=0.D0
        T=0.D0
        B=H/12.D0

             
		do I=1,N
		   RRT=RO1X(I)*R1(I)
           C=R6(I)/(RO1X(I)*RO1X(I))
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
           C=(2.D0*AK+1.)*R6(I)/(R1(I)*RO1X(I))
           D=(8.D0*P-D-5.D0*C)*B
           T=T-D
           T=T/(1.D0+5.D0*(AK+1.)*B/RRT)
           D=P
           P=(AK+1.)*T/RRT-C
           R6(I)=T
        enddo

        return
      end subroutine POT


      subroutine MDRUM(M9,K,N,R5,R8)
        implicit none
        integer::M9,N,K,L9,MAX,I
	    real(8),dimension(:)::R8
        real(8),dimension(:)::R5

        L9=(M9-1)*(N+2)
        MAX=N+2
        IF(K-2) 1,2,1
  2        DO I=1,MAX
              R5(I)=R8(I+L9)
           ENDDO
           RETURN
  1        DO I=1,MAX
              R8(I+L9)=R5(I)
           ENDDO
	    RETURN
      end subroutine MDRUM


      subroutine  REWRIT(KEY,IS,N,R5,R8,Rfun)
        implicit none
        integer::I,IS,KEY,N,J
	    real(8),dimension(:)::R5,R8
        real(8),dimension(:,:)::Rfun



        DO I=1,IS
           IF(KEY.EQ.0) THEN
	          DO J=1,N+2 
	             R5(J)=Rfun(I,J)
	          ENDDO
	       ENDIF
	       IF(KEY.EQ.1) CALL MDRUM(I,2,N,R5,R8)
           IF(KEY.EQ.0) CALL MDRUM(I,7,N,R5,R8)
           IF(KEY.EQ.1) THEN
	          DO J=1,N+2 
	             Rfun(I,J)=R5(J)
	          ENDDO
	       ENDIF
        ENDDO
        IF(KEY.EQ.1) WRITE(6,20)
  20    FORMAT(' Frozen core ')

        return
      end subroutine  REWRIT


  ! ондопнцпюллю пюяверю йскнмнбяйху хмрецпюкнб
  ! нЯСЫЕЯРБКЪЕР РЮЙФЕ БШБНД ГМЮВЕМХИ ХМРЕЦПЮКНБ Б ТЮИК 
  real(8) function RCoulomb(K,I1,I2,I3,I4,IP,N,H,IW,NU,NG,NCN,NN,NST,R1,R8,RO1X,L)
    implicit none
    integer::K,I1,I2,I3,I4,IP,N,IW,NU,NG,NCN
    real(8)::H,CFS2
    integer,dimension(:)::NST,NN,L
    real(8),dimension(:)::R1,RO1X,R8
  	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::J1,J2,J3,J4,I,ierr
    integer,dimension(8)::KK9,II,JJ
    real(8),dimension(8)::RR
	real(8),allocatable,dimension(:)::R6
    character(1),dimension(10)::LB
    DATA LB/'s','p','d','f','g','h','i','j','k','l'/
      
	 
    
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(R6(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'RCoulomb'
	   write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	   stop 
	endif  

    ! гюмскъел оепед пюявернл
    R6=0.D0

      
	
    J1=(I1-1)*(N+2)+2
    J2=(I2-1)*(N+2)+2
    J3=(I3-1)*(N+2)+2
    J4=(I4-1)*(N+2)+2
    DO I=1,N
       R6(I)=R8(I+J1)*R8(I+J3)
    ENDDO
	  
    call POT(DBLE(K),N,H,R6,R1,RO1X)
      
	DO I=1,N
       R6(I)=R6(I)*R8(I+J2)*R8(I+J4)/R1(I)
    ENDDO

	RCoulomb=2*SIMP(R6,N,0,H,RO1X)

    IF(NST(I1)*NST(I2).NE.0.OR.((IW.EQ.NU.OR.IW.EQ.NG).AND.NCN.NE.0)) GOTO 3
  4 RETURN
  5 FORMAT('  F',I1,'(',2(I2,A1),')=',F10.7)
  6 FORMAT('  G',I1,'(',2(I2,A1),')=',F10.7)
  3 NCN=NCN+1
    KK9(NCN)=K
    II(NCN)=I1
    JJ(NCN)=I2
    RR(NCN)=RCoulomb
    IF(IP.EQ.0.AND.(NCN.EQ.5.OR.IW.EQ.NU)) GOTO 7
    IF(IP.EQ.1.AND.(NCN.EQ.5.OR.IW.EQ.NG)) GOTO 8
    RETURN
  7 DO I=1,NCN
       WRITE(6,5)KK9(I),NN(II(I)),LB(L(II(I))+1),NN(JJ(I)),LB(L(JJ(I))+1),RR(I)
    ENDDO
    NCN=0
    RETURN

  8 DO I=1,NCN 
       WRITE(6,6) KK9(I),NN(II(I)),LB(L(II(I))+1),NN(JJ(I)),LB(L(JJ(I))+1),RR(I)
    ENDDO
    NCN=0
    RETURN

  
    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(R6,stat=ierr)
    if(ierr/=0) then
         write(*,*) 'RCoulomb'
         write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif

      return
    end function RCoulomb










  ! оНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ ЙСКНМНБЯЙХУ ХМРЕЦПЮКНБ
  ! нЯСЫЕЯРБКЪЕР РЮЙФЕ БШБНД ГМЮВЕМХИ ХМРЕЦПЮКНБ Б ТЮИК 
  real(8) function RK(K,I1,I2,I3,I4,IP,ISWT,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)
    implicit none
    integer::K,I1,I2,I3,I4,IP,ISWT,N,IW,NU,NG,NCN
    real(8)::H,CFS2
    integer,dimension(:)::NST,NN,L
    real(8),dimension(:)::R1,RAB1,RAB2,RO1X,RO2X,RO3X,R8
  	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::J1,J2,J3,J4,I,ierr
    integer,dimension(8)::KK9,II,JJ
    real(8),dimension(8)::RR
	real(8),allocatable,dimension(:)::R6
    character(1),dimension(10)::LB
    DATA LB/'s','p','d','f','g','h','i','j','k','l'/
      
	 
    
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(R6(N+2),stat=ierr)
	if(ierr/=0) then
       write(*,*) 'RK'
	   write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	   stop 
	endif  

    ! гюмскъел оепед пюявернл
    R6=0.D0

      
	IF(ISWT.EQ.1) GOTO 12
      J1=(I1-1)*(N+2)+2
      J2=(I2-1)*(N+2)+2
      J3=(I3-1)*(N+2)+2
      J4=(I4-1)*(N+2)+2
      DO I=1,N
         R6(I)=R8(I+J1)*R8(I+J3)
      ENDDO
	  
      call POT(DBLE(K),N,H,R6,R1,RO1X)
      
	  DO I=1,N
         R6(I)=R6(I)*R8(I+J2)*R8(I+J4)/R1(I)
      ENDDO

	  RK=2*SIMP(R6,N,0,H,RO1X)

      IF(NST(I1)*NST(I2).NE.0.OR.((IW.EQ.NU.OR.IW.EQ.NG).AND.NCN.NE.0)) GOTO 3
  4   RETURN
  5   FORMAT('  F',I1,'(',2(I2,A1),')=',F10.7)
  6   FORMAT('  G',I1,'(',2(I2,A1),')=',F10.7)
  3   NCN=NCN+1
      KK9(NCN)=K
      II(NCN)=I1
      JJ(NCN)=I2
      RR(NCN)=RK
      IF(IP.EQ.0.AND.(NCN.EQ.5.OR.IW.EQ.NU)) GOTO 7
      IF(IP.EQ.1.AND.(NCN.EQ.5.OR.IW.EQ.NG)) GOTO 8
      RETURN
  7	  DO I=1,NCN
         WRITE(6,5)KK9(I),NN(II(I)),LB(L(II(I))+1),NN(JJ(I)),LB(L(JJ(I))+1),RR(I)
      ENDDO
	  NCN=0
      RETURN

  8   DO I=1,NCN 
         WRITE(6,6) KK9(I),NN(II(I)),LB(L(II(I))+1),NN(JJ(I)),LB(L(JJ(I))+1),RR(I)
      ENDDO
	  NCN=0
      RETURN




  12  J1=(I1-1)*(N+2)+2
      J2=(I2-1)*(N+2)+2
      J3=(I3-1)*(N+2)+2
      J4=(I4-1)*(N+2)+2
      DO I=1,N
      R6(I)=R8(I+J1)*R8(I+J2)*R8(I+J3)*R8(I+J4)/R1(I)*RAB1(I)
      ENDDO
      RK=SIMP(R6,N,0,H,RO1X)*CFS2*(2*K+1)/2


	  
	  ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R6,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'RK'
         write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif

      return
    end function RK




      real(8) function RMK0(L9,L1,KL,ISWK,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)
         implicit none
         integer::N,L9,L1,KL,ISWK
	     real(8)::H,CFS2
         real(8),dimension(:)::R1,RAB1,RAB2,RO1X
         real(8),dimension(:)::R8
         real(8),dimension(:)::R5,R6
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     real(8)::RIN
		 integer::I,J9
        
      
  10  FORMAT('+',T66,'    M0(',I2,I2,')=',E9.2)
  11  FORMAT('+',T88,'  M2(',I2,I2,')=',E9.2)
  12  FORMAT('+',T103,' M4(',I2,I2,')=',E9.2)


      R6(1)=0.D0
      J9=(L9+N)/(N+2)
      IF(ISWK.NE.2.AND.ISWK.NE.4) THEN
        DO I=1,N
           R5(I)=R8(L1+I)**2*RAB2(I)
        ENDDO
	  DO I=7,N
           R6(I)=SIMP(R5,I,3,H,RO1X)*R8(L9+I)**2/R1(I)*RAB1(I)**2
        ENDDO
	  DO I=2,6
           R6(I)=SIMP(R5,I,2,H,RO1X)*R8(L9+I)**2/R1(I)*RAB1(I)**2
        ENDDO
     	  RIN=SIMP(R6,N,3,H,RO1X)*CFS2/2.
        RMK0=RIN
        IF(KL.EQ.99) WRITE(6,10) J9,J9,RIN
	ENDIF

      IF(ISWK.EQ.2) THEN 
      
	   DO I=1,N
            R5(I)=R8(L1+I)**2*RAB2(I)*R1(I)**2
         ENDDO
	   DO I=7,N
            R6(I)=SIMP(R5,I,3,H,RO1X)*(R8(L9+I)*RAB1(I)/R1(I))**2/R1(I)
         ENDDO
	   DO I=2,6
            R6(I)=SIMP(R5,I,2,H,RO1X)*(R8(L9+I)*RAB1(I)/R1(I))**2/R1(I)
         ENDDO
	   RIN=SIMP(R6,N,3,H,RO1X)*CFS2/2.
         RMK0=RIN
      ENDIF

	IF(ISWK.EQ.4) THEN 
   
        DO I=1,N
           R5(I)=R8(L1+I)**2*RAB2(I)*R1(I)**4
        ENDDO
    	  DO I=7,N
          R6(I)=SIMP(R5,I,3,H,RO1X)/R1(I)**3*(R8(L9+I)*RAB1(I)/R1(I))**2
        ENDDO 
	  DO I=2,6
          R6(I)=SIMP(R5,I,2,H,RO1X)/R1(I)**3*(R8(L9+I)*RAB1(I)/R1(I))**2
        ENDDO
	  RMK0=SIMP(R6,N,3,H,RO1X)*CFS2/2.
      
	ENDIF
      
	end function RMK0


      real(8) function H31(L,R1,RAB1,R8)
             implicit none
	       integer::L
	       real(8),dimension(:)::R1,RAB1,R8
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             real(8)::B,B1,G,G1		    


      B=(R8(L+2)/R8(L+3))**2*RAB1(2)/RAB1(3)*R1(3)/R1(2)
      B1=(R8(L+3)/R8(L+4))**2*RAB1(3)/RAB1(4)*R1(4)/R1(3)
      G=DLOG(B)/(R1(3)-R1(2))
      G1=DLOG(B1)/(R1(4)-R1(3))
      B=(B+B1)*0.5
      G=(G+G1)*0.5
      H31=DEXP(G*R1(2))*R8(L+2)**2/R1(2)*RAB1(2)
     
      return
      end function H31



	  real(8) function RIntOrtXX(Fun1,Fun2,Npoints,H,RO1X)
        implicit none
		integer::Npoints
		real(8)::H
	    real(8),dimension(:)::Fun1,Fun2,RO1X
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		integer::ININSX,ierr
		real(8),allocatable,dimension(:)::RINT
    	   
	    !бшдекъел оюлърэ дкъ люяяхбнб
	    allocate(RINT(Npoints),stat=ierr)
	    if(ierr/=0) then
           write(*,*) 'RIntOrtXX'
	       write(*,*) 'MEMORY ON THE FILE "RINT" IS NOT SELECTED'
	       stop 
	    endif  
        
        DO ININSX=1,Npoints
           RINT(ININSX)=Fun1(ININSX+2)*Fun2(ININSX+2)/DSQRT(Fun1(2)*Fun2(2))
        ENDDO

        RIntOrtXX=SIMP(RINT,Npoints,4,H,RO1X)
        ! сдюкемхе люяяхбнб хг оълърх 
        deallocate(RINT,stat=ierr)
        if(ierr/=0) then
           write(*,*) 'RIntOrtXX'
           write(*,*) 'THE FILE "RINT" IS NOT REMOVED FROM MEMORY'
           stop 
        endif
      
        return
      end function RIntOrtXX 
	   

      ! ондопнцпюллю хмрецпхпнбюмхъ
      real(8) function SIMP(A,N9,ISIMP,H,RO1X)
        implicit none
	    integer::ISIMP,N9,I,K,I1
	    real(8)::SUM,H
	    real(8),dimension(:)::A,RO1X
	  
	
        SUM=0.D0

        IF(ISIMP.NE.1.AND.ISIMP.NE.2.AND.ISIMP.NE.3) THEN
           DO I=1,N9
              SUM=SUM+A(I)/RO1X(I)**2
           ENDDO
           SIMP=SUM*H  
        ENDIF

        IF(ISIMP.EQ.1) THEN
	       DO I=1,N9
              SUM=SUM+A(I)
           ENDDO
           SIMP=SUM*H
	    ENDIF

        IF(ISIMP.EQ.2) THEN
	       K=N9-1
           DO I=1,K
              SUM=SUM+A(I)+A(I+1)
           ENDDO
	       SIMP=SUM*H/2.D0
   	    ENDIF
      
        IF(ISIMP.EQ.3) THEN 
     	   I1=N9-3
           DO I=4,I1
              SUM=SUM+A(I)
           ENDDO
           SUM=SUM+(9.*(A(1)+A(N9))+28.*(A(2)+A(N9-1))+23.*(A(3)+A(N9-2)))/24.
           SIMP=SUM*H
        ENDIF      
     
        return
      end function SIMP


	! оНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ ОНКМНИ ЩМЕПЦХХ ЙНМТХЦСПЮЖХХ Я СВЕРНЛ ПЕКЪРХБХЯРЯЙХУ ОНОПЮБНЙ
    subroutine TOTAL(Z,N,H,NU,NG,IS,NCN,RAN,VAL,TETA,Q,L,KK,IS1,IS2,GAM,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X)
      implicit none
	  integer::N,NU,NG,IS,NCN
	  real(8)::Z,H,VAL,RAN
      integer,dimension(:)::Q,L,NN,NST,KK,IS1,IS2
      real(8),dimension(:)::TETA,GAM,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  character(1),dimension(10)::LB  
      integer::I,M9,QQ,IW,J,M1,ierr
	  real(8)::EYA,ELE,PE,TE,ESF,ENSF,RESF,RENSF,R,EKIN,VIR,CFS2,A1,B1
	  real(8)::DAR,RH2,RE,RKS,DAR1,RE1,DZ,C1,TER	   
	  real(8),parameter::CFS=0.0072973506D0 
	  real(8),allocatable,dimension(:)::R2,R4,R5,R6
    	   
	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(R2(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'TOTAL'
	     write(*,*) 'MEMORY ON THE FILE "R2" IS NOT SELECTED'
	     stop 
	  endif    
	  allocate(R4(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'TOTAL'
	     write(*,*) 'MEMORY ON THE FILE "R4" IS NOT SELECTED'
	     stop 
	  endif    	   
      allocate(R5(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'TOTAL'
	     write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
	     stop 
	  endif    
	  allocate(R6(N+2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'TOTAL'
	     write(*,*) 'MEMORY ON THE FILE "R6" IS NOT SELECTED'
	     stop 
	  endif    	  
      
	  DATA LB/'s','p','d','f','g','h','i','j','k','l'/
	
	
  99  FORMAT(/'  nl    Energy    Radius      Teta     Norm             Relat. mass,Ry    Darwin.      Total cor.      Spin-Orbit   ')
  7   FORMAT(I3,A1,2F10.4,2E11.3,6X,F13.6,2X,F13.6,3X,F13.6,3X,F12.6)
  9   FORMAT(47X,' Sum ',F13.6,2X,F13.6,7X,' Orbit -Orbit = ',F10.6)
      
	  WRITE(6,99)

      EYA=0.D0
      ELE=0.D0
      PE=0.D0
      TE=0.D0
      ESF=0.D0
      ENSF=0.D0
      RESF=0.D0
      RENSF=0.D0
      NCN=0
      DAR=0.D0
      RH2=0.D0
      RE=0.D0
      RKS=0.D0
      CFS2=CFS*CFS
      R4=0.D0
      R2=0.D0


      DO I=1,IS
         M9=(I-1)*(N+2)+2
         DO J=1,N
            R2(J)=R8(J+M9)**2*R1(J)
         ENDDO
	     VAL=SIMP(R2,N,0,H,RO1X)
         QQ=Q(I)
         TE=TE+QQ*R8(M9-1)
         DAR1=0.D0
         IF(L(I).EQ.0) DAR1=H31(M9,R1,RAB1,R8)*CFS2*Z/4
         DAR=DAR+DAR1*QQ
         RE1=-EMV2(L(I),M9,N,H,RAN,R8,R1,RAB1,RAB2,RO1X)*CFS2/8
         DZ=0.D0
      
	     IF(L(I).NE.0) THEN
            DO J=1,N
               R6(J)=R8(M9+J)**2/R1(J)*RAB1(J)**2
            ENDDO
            DZ=Z*CFS2*SIMP(R6,N,1,H,RO1X)
            DO J=1,IS
               M1=(J-1)*(N+2)+2
               IF(I.EQ.J) THEN 
	              DZ=DZ-RMK0(M9,M1,0,0,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)*(2*QQ-3)
               ENDIF
	           IF(I.NE.J) THEN 
	              DZ=DZ-RMK0(M9,M1,0,0,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)*2*Q(J)
               ENDIF
            ENDDO
         ENDIF
      
	     WRITE(6,7) NN(I),LB(L(I)+1),R8(M9-1),VAL,TETA(I),R8(M9),RE1,DAR1,RE1+DAR1,DZ

         RE=RE+RE1*QQ


         IF(L(I).EQ.0) GOTO 16
             A1=RMK0(M9,M9,1,0,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)
            IF(L(I).EQ.1) GOTO 13
               B1=RMK0(M9,M9,90,2,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)
            IF(L(I).EQ.2) GOTO 14
               C1=RMK0(M9,M9,90,4,N,H,CFS2,R1,R5,R6,R8,RAB1,RAB2,RO1X)
            IF(L(I).EQ.3) GOTO 15
  13           RH2=RH2+A1*QQ*(QQ-1)*0.8
               GOTO 16
  14           RH2=RH2+(A1+B1*2/7.)*QQ*(QQ-1)*4/3.
               GOTO 16
  15           RH2=RH2+8./13*QQ*(QQ-1)*(3*A1+B1+C1*5/11.)
  16     CONTINUE
      
         DO J=1,N
            R4(J)=R4(J)+R8(M9+J)**2*QQ
         ENDDO
	
      ENDDO




 
      WRITE(6,9) RE,DAR,RH2
      WRITE(6,7)

      M9=NG+1
  
      DO IW=M9,NU
         
	    IF(KK(IW).EQ.0) ESF=ESF-GAM(IW)*RK(0,IS1(IW),IS2(IW),IS1(IW),IS2(IW),0,0,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)
         
	    IF(KK(IW).NE.0) ENSF=ENSF-GAM(IW)*RK(KK(IW),IS1(IW),IS2(IW),IS1(IW),IS2(IW),0,0,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)
      RESF=RESF-RK(KK(IW),IS1(IW),IS2(IW),IS1(IW),IS2(IW),0,1,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)*GAM(IW)
      ENDDO
  
 
      NCN=0
 
      DO IW=1,NG
      ENSF=ENSF-GAM(IW)*RK(KK(IW),IS1(IW),IS2(IW),IS2(IW),IS1(IW),1,0,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)
      RENSF=RENSF-RK(KK(IW),IS1(IW),IS2(IW),IS2(IW),IS1(IW),1,1,N,H,CFS2,IW,NU,NG,NCN,NN,NST,R1,R8,RAB1,RAB2,RO1X,RO2X,RO3X,L)*GAM(IW)
      ENDDO
	 
      DO I=1,N
         R4(I)=R4(I)/R1(I)
      ENDDO
      EYA=-2.*Z*SIMP(R4,N,0,H,RO1X)


      R=ESF+ENSF
      RKS=RESF+RENSF
      PE=R+EYA
      TE=TE-R
      TER=TE+RE+DAR+RKS+RH2
      EKIN=TE-PE
      VIR=PE/EKIN


      ! бшбнд пегскэрюрнб
      
  6   FORMAT(/'  Total  energy      Nonrel. ',F10.3,5X,' Rel. ',F10.3/'  Kinetic energy             'F10.3,'    Virial ',F17.10)
 114  FORMAT(/'  Relat. correction of FK,GK ',F10.6)
      
      WRITE(6,114) RKS
      WRITE(6,6) TE,TER,EKIN,VIR
      
	  ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R2,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'TOTAL'
         write(*,*) 'THE FILE "R2" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
	  deallocate(R4,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'TOTAL'
         write(*,*) 'THE FILE "R4" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
	  deallocate(R5,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'TOTAL'
         write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
	  deallocate(R6,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'TOTAL'
         write(*,*) 'THE FILE "R6" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
      return
      end subroutine TOTAL



      ! оНДОПНЦПЮЛЛЮ ПЮЯВЕРЮ ОНКМНИ ЩМЕПЦХХ ЙНМТХЦСПЮЖХХ
      subroutine TOTALXF(Z,N,H,NU,NG,IS,Q,L,KK,IS1,IS2,GAM,NN,NST,R1,R8,RO1X)
        implicit none
	    integer::N,NU,NG,IS
	    real(8)::Z,H
        integer,dimension(:)::Q,L,NN,NST,KK,IS1,IS2
	    real(8),dimension(:)::GAM,R1,R8,RO1X 
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer::I,M9,QQ,IW,J,NCN,ierr
	    real(8)::EYA,ELE,PE,TE,ESF,ENSF,R,EKIN,VIR
	  	real(8),allocatable,dimension(:)::R4
    	   
	    !бшдекъел оюлърэ дкъ люяяхбнб
	    allocate(R4(N+2),stat=ierr)
	    if(ierr/=0) then
           write(*,*) 'TOTALXF'
	       write(*,*) 'MEMORY ON THE FILE "R4" IS NOT SELECTED'
	       stop 
	    endif  
   

        EYA=0.D0
        ELE=0.D0
        PE=0.D0
        TE=0.D0
        ESF=0.D0
        ENSF=0.D0
        NCN=0
        R4=0.D0
      
	 
	    DO I=1,IS
           M9=(I-1)*(N+2)+2
           QQ=Q(I)
           TE=TE+QQ*R8(M9-1)
           DO J=1,N
              R4(J)=R4(J)+R8(M9+J)**2*QQ
           ENDDO
        ENDDO

        M9=NG+1
        DO IW=M9,NU
           IF(KK(IW).EQ.0) THEN
	          ESF=ESF-GAM(IW)*RCoulomb(0,IS1(IW),IS2(IW),IS1(IW),IS2(IW),0,N,H,IW,NU,NG,NCN,NN,NST,R1,R8,RO1X,L)		                                    
		   ENDIF
	       IF(KK(IW).NE.0) THEN 
	         ENSF=ENSF-GAM(IW)*RCoulomb(KK(IW),IS1(IW),IS2(IW),IS1(IW),IS2(IW),0,N,H,IW,NU,NG,NCN,NN,NST,R1,R8,RO1X,L)
		   ENDIF
        ENDDO

        NCN=0
        DO IW=1,NG
           ENSF=ENSF-GAM(IW)*RCoulomb(KK(IW),IS1(IW),IS2(IW),IS2(IW),IS1(IW),1,N,H,IW,NU,NG,NCN,NN,NST,R1,R8,RO1X,L)                       
		ENDDO

        DO I=1,N
           R4(I)=R4(I)/R1(I)
        ENDDO

        EYA=-2.D0*Z*SIMP(R4,N,0,H,RO1X)

      
	    R=ESF+ENSF
       ! WRITE(6,*) 'EPS',TE
	   ! WRITE(6,*) 'Culomb',R
	   ! WRITE(6,*) 'Pot Z',EYA
		PE=R+EYA
	   !WRITE(6,*) 'POT',PE
	   !WRITE(6,*) 'KIN',TE-2.D0*R-EYA
	   !WRITE(6,*) 'COMPONET KIN',TE,-2.D0*R,-EYA
        TE=TE-R
        EKIN=TE-PE
        VIR=PE/EKIN


  100 FORMAT(2X,F14.6)
    6 FORMAT(/'  Total energy  ',T30,F14.6/'  Kinetic energy    ',T30,F14.6/'  Virial            ',T30,F17.10)
      
      WRITE(6,6) TE,EKIN,VIR
	 ! WRITE(6,*) DABS(TE)-DABS(EKIN)
	  WRITE(10,100) TE*0.5D0 ! б юрнлмшу едхмхжюу
	  
	  ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(R4,stat=ierr)
      if(ierr/=0) then
         write(*,*) 'TOTALXF'
         write(*,*) 'THE FILE "R4" IS NOT REMOVED FROM MEMORY' 
 	     stop 
      endif
     
     return
   end subroutine TOTALXF

    
      ! оНДОПНЦПЮЛЛЮ ОНКСВЮЕР МСКЕБНЕ ОПХАКХФМХЕ ДКЪ ПЕЬЕМХЪ ЯХЯРЕЛШ СПЮБМЕМХИ уЮПРПХ-тНЙЮ
     subroutine ZERO(Z,IS,NN,Q,L,NZ,NWR,N,RAL,ACY,VAL,R1,R8,P4,E3,Rfun)
      implicit none
	  integer::NZ,NWR,IS,N,K,I,I1
	  real(8)::Z,RAL,ACY,VAL,A,Z1,R
      integer,dimension(:)::NN,Q,L
	  real(8),dimension(:)::R8,P4,R1,E3
      real(8),dimension(:,:)::Rfun
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	  integer::ierr
	  real(8)::A1,B1
	  real(8),allocatable,dimension(:)::R5
    	   
	  
   
      
      ! оепегюохяэ тсмйжхи 
      IF(NZ.NE.1) THEN
         !бшдекъел оюлърэ дкъ люяяхбнб
	     allocate(R5(N+2),stat=ierr)
	     if(ierr/=0) then
            write(*,*) 'ZERO'
	        write(*,*) 'MEMORY ON THE FILE "R5" IS NOT SELECTED'
	        stop 
	     endif  

  	     call REWRIT(NWR,IS,N,R5,R8,Rfun)
         
		 ! сдюкемхе люяяхбнб хг оълърх 
         deallocate(R5,stat=ierr)
         if(ierr/=0) then
            write(*,*) 'ZERO'
            write(*,*) 'THE FILE "R5" IS NOT REMOVED FROM MEMORY' 
 	        stop 
         endif 
		 return
	  ENDIF
      
	  ! пЮЯВЕР МЮВЮКЭМНЦН БЮПХЮМРЮ НДМНЩКЕЙРПНММШУ ЩМЕПЦХИ
      K=(N+2)*IS
      A1=.5D0
      B1=.25D0
      
	  DO I=1,IS
         A=0
         DO I1=1,IS
            IF((NN(I)-NN(I1)).LT.0) THEN
               A=A+B1**(NN(I1)-NN(I))*Q(I1)           
            ENDIF
            IF((NN(I)-NN(I1)).EQ.0) THEN
               A=A+A1*Q(I1)           
            ENDIF
	        IF((NN(I)-NN(I1)).GT.0) THEN
               A=A+Q(I1)      
            ENDIF
         ENDDO
         A=Z-A-L(I)**2
         E3(I)=(A/NN(I))**2
         IF(E3(I).LT.1.OR.A.LT.0) THEN 
		    E3(I)=1+VAL
         ENDIF
	  ENDDO

 
      ! гюмскъел оепед пюявернл
      R8=0.D0
	  
      K=0

      DO I=1,IS
         R8(1+(I-1)*(N+2))=E3(I)
	     K=K+Q(I)
      ENDDO

	  K=IDINT(Z)-K+1+IDINT(ACY)
      IF(K.LT.1) THEN 
	     K=1
      ENDIF
      Z1=Z**(1./3.)
      DO I=1,N
          R=R1(I)*Z1
          P4(I)=-K-(Z-K)*DEXP(-0.2075*R)/(1.+1.19*R)-VAL/(1.+DABS(RAL-R1(I)))*R1(I)
      ENDDO
      
	  
       
	  return
     end subroutine ZERO








	 
	 
end module mAF
