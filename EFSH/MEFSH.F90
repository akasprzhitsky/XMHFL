    ! MODULE EXPANSION OF FUNCTION ON SPHERICAL HARMONICS VER 3.0 11.2019
	! VER 3.0 NEW  11,2019 цнд

	! SPHERICAL HARMONIC EXPANSION MODULE OF A FUNCTION VERSION 3.0
	

module mefsh
    implicit none	  
	      

	contains
    
! SUB-PROGRAM FOR DECOMPOSING THE ATOMIC WAVE FUNCTION INTO A SERIES IN SPHERICAL HARMONICS
! Nn-PRIMARY QUANTUM NUMBER
! L-ORBITAL MOMENT OF THE DEVELOPABLE FUNCTION
! ML-PROJECTION OF THE ORBITAL MOMENTUM OF THE DECOMPOSABLE FUNCTION
! Npoint-NUMBER OF DOTS OF THE ARGUMENT OF VALUES OF THE EXPANDABLE FUNCTION
! R (Npoint) - ARGUMENT VALUE ARRAY
! RFun (Npoint) - ARRAY OF FUNCTION VALUES
! RFunA (Npoint) - ARRAY OF APPROXIMATION VALUES OF RFun FUNCTION (FOR CHECKING APPROXIMATION)
! Lmin-ORBITAL MOMENT OF THE FIRST HARMONIC
! Lmax-ORBITAL MOMENT OF THE LAST HARMONIC
! A-DISTANCE FROM THE BEGINNING OF THE ATOMIC SYSTEM COORDINATES TO THE CENTER RELATIVE TO WHICH THE DECOMPOSITION IS CARRIED OUT
! SPECIAL REFERENCE TO THE TYPE OF DISPLACEMENT OF THE ATOMIC SYSTEM OF COORDINATES
! A> 0-offset in the positive direction of the OZ axis
! A <0-offset in the negative direction of the OZ axis
! THE DISTANCE SIGN INDICATES THE TYPE OF OFFSET
! NpointNew - NUMBER OF DOTS OF THE ARGUMENT OF VALUES OF THE ARGUMENT OF THE EXPANSION COEFFICIENT
! Rnew (NpointNew) - ARRAY OF VALUES OF THE ARGUMENT OF DECOMPOSITION COEFFICIENTS
! RcoffSH (Lmax-Lmin + 1, NpointNew) - ARRAY OF VALUES OF THE EXPANSION COEFFICIENTS OF THIS FUNCTION IN A SERIES IN SPHERICAL HARMONICS
! NRabParametrs () - ARRAY OF OPERATING PARAMETERS            
	subroutine EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS(Nn,L,ML,Npoint,R,RFun,RFunA,Lmin,Lmax,A,NpointNew,Rnew,RcoffSH,NRabParametrs) 
     implicit none
     integer::Nn,L,ML,Npoint,Lmin,Lmax,NpointNew
	 real(8)::A
	 integer,dimension(:)::NRabParametrs
	 real(8),dimension(:)::R,Rnew,RFun,RFunA
     real(8),dimension(:,:)::RcoffSH
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer::Ninterval,NpolAR,IS1,II,ILS,IZXC,ierr,NumbreIntTeta,IndexGarmon
	 integer,allocatable,dimension(:,:)::NumbreInter,NInterH,NInterK
     real(8),allocatable,dimension(:)::QL1,ZM,Acoff
     real(8),allocatable,dimension(:,:)::XlimF,ALFAPolin,AcoffPolinom,Zaa,RcoffSin,RcoffCos
	 real(8),allocatable,dimension(:,:,:,:,:,:)::XlimZFK,AcoffApro



	 ! бшдекъел оюлърэ онд люяяхбш юоопнйяхлюжхх
     ! ЛЮЙЯХЛЮКЭМНЕ ВХЯКН ХМРЕПБЮКНБ 200
     allocate(XlimF(2,200),stat=ierr)
	 if(ierr/=0) then
       write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
	   write(*,*)'MEMORY ON THE FILE "XlimF" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(ALFAPolin(3,200),stat=ierr)
	 if(ierr/=0) then
       write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
	   write(*,*)'MEMORY ON THE FILE "ALFAPolin" IS NOT SELECTED'
	   stop 
	 endif
     ! бшдекъел оюлърэ дкъ люяяхбнб +5-днонкмхрекэмне вхякн якюцюелшу
 	 allocate(AcoffPolinom(Nn-L+5+2,200),stat=ierr)
	 if(ierr/=0) then
       write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
	   write(*,*)'MEMORY ON THE FILE "AcoffPolinom" IS NOT SELECTED'
	   stop 
	 endif
     
	 allocate(QL1(L-IABS(ML)+1),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	    write(*,*) 'MEMORY ON THE FILE "QL1" IS NOT SELECTED'
	    stop 
	 endif
     allocate(ZM(IABS(ML)+1),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(*,*) 'MEMORY ON THE FILE "ZM" IS NOT SELECTED'
	   stop 
	 endif
     allocate(Zaa(L-IABS(ML)+1,L-IABS(ML)+1),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(*,*) 'MEMORY ON THE FILE "Zaa" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(RcoffSin(2,2*IABS(ML)+2),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(*,*) 'MEMORY ON THE FILE "RcoffSin" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(RcoffCos(L-IABS(ML)+1,L-IABS(ML)+1),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(*,*) 'MEMORY ON THE FILE "RcoffCos" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(Acoff(NpointNew),stat=ierr)
     if(ierr/=0) then
        write(*,*)'MEMORY ON THE FILE "Acoff" IS NOT SELECTED'
        stop 
     endif
	 ! гюмскъел оепед пюявернл
	 AcoffPolinom=0.D0
     XlimF=0.D0
     ALFAPolin=0.D0
	 !!!!!!!!!!!!!!!!!!!
	 QL1=0.D0
	 ZM=0.D0
	 Zaa=0.D0
	 RcoffSin=0.D0
	 RcoffCos=0.D0  
     Acoff=0.D0
     RCoffSH=0.D0

     ! тнплюрш бшдювх хмтнплюжхх
500  FORMAT(2X,'START',I4,1X,'Npoints= ',I5)


 

     ! щрюо 1. юопнйяхлхпсел пюдхюкэмсч вюярэ бнкмнбни тсмйжхх
	 call EFSH_APPROXIMATION_RADIAL_FUNCTION_ALFA(L,Npoint,R,RFun,10,Ninterval,NpolAR,XlimF,ALFAPolin,AcoffPolinom,RFunA) 
     
 


	 ! бшдекъел оюлърэ онд люяяхбш 
	 ! 15- люйяхлюкэмне вхякн хмрепбюкнб мю йнрнпше пюгахбюеряъ тсмйжхъ опх юоопнйяхлюжхх он сцкс 
	 ! вхякн хмрепбюкнб юоопнйяхлюжхх он сцкс нопедекъееряъ б ондопнцпюлле EFSH_PARAMETR_FUNCTION_RDFUN
	 allocate(XlimZFK(L-IABS(ML)+1,NpointNew,2,Ninterval,NpolAR+1,15),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	    write(*,*) 'MEMORY ON THE FILE "XlimZFK" IS NOT SELECTED'
	    stop 
	 endif
     allocate(AcoffApro(L-IABS(ML)+1,NpointNew,3,Ninterval,NpolAR+1,15),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	    write(*,*) 'MEMORY ON THE FILE "AcoffApro" IS NOT SELECTED'
	    stop 
	 endif
	 allocate(NumbreInter(L-IABS(ML)+1,NpointNew),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(*,*) 'MEMORY ON THE FILE "NumbreInter" IS NOT SELECTED'
	   stop 
	 endif
     allocate(NInterH(L-IABS(ML)+1,NpointNew),stat=ierr)
     if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(6,*) 'MEMORY ON THE FILE "NInterH" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(NInterK(L-IABS(ML)+1,NpointNew),stat=ierr)
     if(ierr/=0) then
       write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
	   write(6,*) 'MEMORY ON THE FILE "NInterK" IS NOT SELECTED'
	   stop 
	 endif
 
 
     ! гюмскъел оепед пюявернл
	 NumbreInter=0
     NInterH=0
     NInterK=0 
     XlimZFK=0.D0
     AcoffApro=0.D0
 
     
   
     ! щрюо 2.юоопнйяхлюжхъ пюгкнцюелни тсмйжхх он сцкс хмрецпхпнбюмхъ
     DO IS1=1,L-IABS(ML)+1 
	    DO II=1,NpointNew
           ! ОПНБНДХЛ ОНЯРПНЕМХЕ ЮООПНЙЯХЛЮЖХХ ХМРЕЦПХПСЕЛНИ ТСМЙЖХХ
	       call EFSH_PARAMETR_FUNCTION_ALFAZX(IS1,II,L,IABS(ML)+IS1-1,Rnew(II),A,Ninterval,NpolAR,XlimF,ALFAPolin,NumbreInter,NInterH,NInterK,NumbreIntTeta,XlimZFK,AcoffApro)
           WRITE(*,500) IS1,II
	    ENDDO
     ENDDO

     ! щрюо3.тнплхпсел бяонлнцюрекэмше люяхб дкъ пюгкнцюелни тсмйжхх 
     call EFSH_CALCULATION_COEFFICIENT_ALFAZX(L,ML,QL1,ZM,Zaa,RcoffSin,RcoffCos)
    
     ! жхйк он нпахрюкэмшл лнлемрюл пюгкнфемхъ
     IndexGarmon=0
	 DO ILS=Lmin,Lmax
	    IndexGarmon=IndexGarmon+1 
        ! ондопнцпюллю пюяверю йнщттхжхемрю пюгкнфемхъ пъдю 
        call EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA(L,ILS,ML,A,Ninterval,NpolAR,XlimF,AcoffPolinom,ALFAPolin,NpointNew,Rnew,ACoff,QL1,ZM,Zaa,RcoffSin,RcoffCos,NumbreInter,NInterH,NInterK,NumbreIntTeta,XlimZFK,AcoffApro)
	    ! гюохяшбюел онксвеммши йнщттхжхемр
	    DO IZXC=1,NpointNew
	       RCoffSH(IndexGarmon,IZXC)=ACoff(IZXC)
	    ENDDO  
	 ENDDO
    
	 ! сдюкъел люяяхбш хг оюлърх
     deallocate(AcoffPolinom,stat=ierr)
	 if(ierr/=0) then
        write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
		write(*,*) 'THE FILE "AcoffPolinom" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(XlimF,stat=ierr)
	 if(ierr/=0) then
        write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
		write(*,*) 'THE FILE "XlimF" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
   	 deallocate(ALFAPolin,stat=ierr)
	 if(ierr/=0) then
       write(*,*)'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS' 
	   write(*,*) 'THE FILE "ALFAPolin" IS NOT REMOVED FROM MEMORY'
	   stop 
	 endif
	 deallocate(QL1,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "QL1" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
     deallocate(ZM,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "ZM" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
     deallocate(Zaa,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "Zaa" IS NOT REMOVED FROM MEMORY'
	    stop 
     endif
	 deallocate(RcoffSin,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "RcoffSin" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(RcoffCos,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "RcoffCos" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(XlimZFK,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "XlimZFK" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(AcoffApro,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "AcoffApro" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(NumbreInter,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "NumbreInter" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(NInterH,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "NInterH" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
	 deallocate(NInterK,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
        write(*,*) 'THE FILE "NInterK" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
     deallocate(Acoff,stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS'
		write(*,*) 'THE FILE "Acoff" IS NOT REMOVED FROM MEMORY'
	    stop 
	 endif
         
     return
    end subroutine EFSH_FUNCTION_EXPANSION_SPHERICAL_HARMONICS 


! SUBPROGRAM FOR CALCULATING THE EXPANSION COEFFICIENT OF A FUNCTION IN A SERIES IN SPHERICAL HARMONICS
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL MOMENT OF FUNCTIONS DECOMPOSABLE IN A SERIES
! L2-ORBITAL MOMENT OF EXPANSION COEFFICIENT
! M-PROJECTION OF THE ORBITAL MOMENT
! Ninterval - NUMBER OF INTERVALS OF APPROXIMATED DECOMPOSABLE FUNCTION
! NpolA-degree of a polynomial with coefficients AcoffPolinom
! XlimF (2, Ninterval) - ARRAY OF APROXIMATION INTERVAL BOUNDARIES
! AcoffPolinom (NpolA + 1, Ninterval) - ARRAY OF COEFFICIENTS OF THE POLYNOM IN THE FUNCTION
! ALFAPolin (3, Ninterval) - ARRAY OF SECOND-ORDER POLYNOMA COEFFICIENTS ALFA
! Na-NUMBER OF DOTS OF DECOMPOSITION COEFFICIENT
! A-DISTANCE FROM THE BEGINNING OF THE ATOMIC SYSTEM COORDINATES TO THE CENTER RELATIVE TO WHICH THE DECOMPOSITION IS CARRIED OUT
! SPECIAL REFERENCE TO THE TYPE OF DISPLACEMENT OF THE ATOMIC SYSTEM OF COORDINATES
! A> 0-offset in the positive direction of the OZ axis
! A <0-offset in the negative direction of the OZ axis
! THE DISTANCE SIGN INDICATES THE TYPE OF OFFSET
! R (N) -ARRAY OF RADIUS VALUES
! RFUN (N) - ARRAY OF WAVE FUNCTION VALUES AT POINTS R (N)
! Ra (N) - ARRAY OF RADIUS VALUES IN THE SYSTEM OF COORDINATES REGARDING THE CENTER OF DECOMPOSITION
! ACoff (N) - VALUE OF THE EXPANSION COEFFICIENT AT POINTS Ra (N)
! QL1 (L1-IABS (M) +1) - ARRAY OF HARMONIC REPRESENTATION COEFFICIENTS OF THE DECOMPOSABLE FUNCTION
! ZM (IABS (M) +1) -BINOM COEFFICIENTS ARRAY
! Zaa (L1-IABS (M) + 1, L1-IABS (M) +1) -BINOMIAL COEFFICIENTS ARRAY
! RcoffSin (2,2 * IABS (M) +2) -ARRAY OF REPRESENTATION COEFFICIENTS sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
! RcoffCos (L1-IABS (M) + 1, L1-IABS (M) +1) -ARRAY OF REPRESENTATION COEFFICIENTS cos (x) ^ i = SUM (Rcoffp * cos (px)) FOR EACH DEGREE
! FROM 0 TO i
! NInt (L1-IABS (M) + 1, Na) - ARRAY OF INTERVAL NUMBERS FITTING IN THE INTEGRATION INTERVAL
! NInterH (L1-IABS (M) + 1, Na) -ARRAYS OF NUMBERS OF THE INITIAL APPROXIMATION INTERVAL INTO THE INTEGRATION INTERVAL
! NInterK (L1-IABS (M) + 1, Na) - ARRAY OF NUMBERS OF THE END INTERVAL OF APPROXIMATION OF THE FOLLOWED INTO THE INTEGRATION INTERVAL
! NumbreIntTeta-NUMBER OF INTERVALS OF APPROXIMATION OF FUNCTIONS BY ANGLE OF INTEGRATION
! XlimZFK (L1-IABS (M) + 1, Na, 2, NInt, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION INTERVALS BY THE INTEGRATION ANGLE
! AcoffApro (L1-IABS (M) + 1, Na, 3, NInt, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION COEFFICIENTS AT APPROXIMATION INTERVALS
    subroutine EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA(L1,L2,M,A,Ninterval,NpolA,XlimF,AcoffPolinom,ALFAPolin,Na,Ra,ACoff,QL1,ZM,Zaa,RcoffSin,RcoffCos,NumbreInter,NInterH,NInterK,NumbreIntTeta,XlimZFK,AcoffApro)
      implicit none

      integer::L1,L2,M,Na,Ninterval,NpolA,NumbreIntTeta
      real(8)::A
      integer,dimension(:,:)::NumbreInter,NInterH,NInterK
	  real(8),dimension(:)::Ra,ACoff,QL1,ZM
      real(8),dimension(:,:)::XlimF,AcoffPolinom,ALFAPolin
	  real(8),dimension(:,:)::Zaa,RcoffSin,RcoffCos
	  real(8),dimension(:,:,:,:,:,:)::XlimZFK,AcoffApro
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::I,K,IS1,IS2
      integer::ierr
	  real(8)::RcofA,SUM
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8),allocatable,dimension(:)::TL2,XXA,YYA
      real(8),allocatable,dimension(:,:)::RcoffAS
       
      


      ! опнбепъел ме опебшьюер кх гмювемхъ опнежхх лнлемрю мюд гмювемхе нпахрюкэмнцн   
	  ! лнлемрю йнщтхжхемрю пюгкнфемхъ
	  IF(IABS(M).GT.L2) THEN
         ACoff=0.d0
         return
	  ENDIF
  	
	  !бшдекъел оюлърэ дкъ люяяхбнб
      allocate(TL2(L2-IABS(M)+1),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
	     write(*,*) 'MEMORY ON THE FILE "TL2" IS NOT SELECTED'
	     stop 
	  endif
      allocate(RcoffAS(L1-IABS(M)+1,Na),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
	     write(*,*) 'MEMORY ON THE FILE "RcoffAS" IS NOT SELECTED'
	     stop 
	  endif
      allocate(XXA(L2-IABS(M)+1),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
	     write(*,*) 'MEMORY ON THE FILE "XXA" IS NOT SELECTED'
	     stop 
	  endif
      allocate(YYA(L2-IABS(M)+1),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
	     write(*,*) 'MEMORY ON THE FILE "YYA" IS NOT SELECTED'
	     stop 
	  endif
	  
      ! гюмскъел оепед пюявернл
      TL2=0.D0
      RcoffAS=0.D0
   	  XXA=0.D0
      YYA=0.D0
      ACoff=0.d0

     


	  ! щрюо 1. онксвюел йнщттхжемрш пюгкнфемхъ
	  ! онксвюел йнщттхжхемрш опедярюбкемхъ цюплнмхйх йнщттхжхемрю пюгкнфемхъ
      call EFSH_VALUE_TETA_FUNCTION(L2,IABS(M),L2-IABS(M)+1,XXA,YYA)
  	  call EFSH_COEFFICIENT_TETA_FUNCTION_APRO(IABS(M),L2-IABS(M)+1,XXA,YYA,TL2) 
      
	 
10000 FORMAT(2X,'N1= ',I3,' N2= ',I3,' Npoints=',I2,' RcoffA= ',E14.6,' Int= ',E14.6)
20000 FORMAT(2X,'Lharmonics=',I4,1X,'Npoints= ',I5)
      

      ! щрюо 2. пюявер йнщттхжхемрнб пюгкнфемхъ  
      DO IS1=1,L1-IABS(M)+1 
	     DO I=1,Na
            SUM=0.D0
	        ! бшъбкъел рхо ялеыемхъ
		    IF(A.GT.0) THEN
		       DO IS2=1,IS1
                  RcofA=(-1.D0)**(IS1+IS2-2)*Ra(I)**(IABS(M)+IS2-1)*QL1(IS1)*Zaa(IS1,IS2)*DABS(A)**(IS1-IS2)
                  SUM=SUM+RcofA*EFSH_COEFFICIENT_INTEGRAL_HARMONICS_ALFAZX(IS1,I,L1,L2,IABS(M),IS2-1,TL2,ZM,RcoffCos,RcoffSin,NumbreInter(IS1,I),NpolA,NInterH(IS1,I),NInterK(IS1,I),NumbreIntTeta,XlimZFK,AcoffApro,AcoffPolinom)
               ENDDO
              ELSE
               DO IS2=1,IS1
	              RcofA=Ra(I)**(IABS(M)+IS2-1)*DABS(A)**(IS1-IS2)*QL1(IS1)*Zaa(IS1,IS2)
                  SUM=SUM+RcofA*EFSH_COEFFICIENT_INTEGRAL_HARMONICS_ALFAZX(IS1,I,L1,L2,IABS(M),IS2-1,TL2,ZM,RcoffCos,RcoffSin,NumbreInter(IS1,I),NpolA,NInterH(IS1,I),NInterK(IS1,I),NumbreIntTeta,XlimZFK,AcoffApro,AcoffPolinom)
               ENDDO
               SUM=-SUM ! ябъгюмн я хглемхмхел мюопюбкемхъ 
		                ! хмрецпхпнбюмхъ (PI,0) дкъ яксвюъ ялеыемхе б нрпхжюрекэмнл мюопюбкемхх нях Z 
            ENDIF
	        WRITE(*,20000) L2,I
	        ! гюохяшбюел йнщттхжхемр
	        RcoffAS(IS1,I)=SUM 
	     ENDDO
      ENDDO

	  ! тнплхпсел люяяхб гмювемхи йнщттхжхемрю пюгкнфемхъ
	  DO I=1,Na
	     SUM=0.D0
	     DO IS1=1,L1-IABS(M)+1 
	        SUM=SUM+RcoffAS(IS1,I)
	     ENDDO 
	     ACoff(I)=SUM 
	  ENDDO   

  
	  ! сдюкемхе люяяхбнб хг оълърх 
	  deallocate(XXA,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
         write(*,*) 'THE FILE "XXA" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(YYA,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
         write(*,*) 'THE FILE "YYA" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(TL2,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
         write(*,*) 'THE FILE "TL2" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
      deallocate(RcoffAS,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA'
         write(*,*) 'THE FILE "RcoffAS" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
      
	  return
    end subroutine EFSH_COEFFICIENT_EXPANSION_SPHERICAL_HARMONICS_ALFA




! SUB-PROGRAM FOR CALCULATION OF AUXILIARY PARAMETERS
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL TORQUE OF THE EXPANDABLE FUNCTION
! M-PROJECTION OF THE MOMENT
! QL1 (L1-IABS (M) +1) - ARRAY OF HARMONIC REPRESENTATION COEFFICIENTS OF THE DECOMPOSABLE FUNCTION
! ZM (IABS (M) +1) -BINOM COEFFICIENTS ARRAY
! Zaa (L1-IABS (M) + 1, L1-IABS (M) +1) -BINOMIAL COEFFICIENTS ARRAY
! RcoffSin (2,2 * IABS (M) +2) - ARRAY OF REPRESENTATION COEFFICIENTS sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
! RcoffCos (L1-IABS (M) + 1, L1-IABS (M) +1) -ARRAY OF REPRESENTATION COEFFICIENTS cos (x) ^ i = SUM (Rcoffp * cos (px)) FOR EACH DEGREE
! FROM 0 TO i
    subroutine EFSH_CALCULATION_COEFFICIENT_ALFAZX(L1,M,QL1,ZM,Zaa,RcoffSin,RcoffCos)
     implicit none
     integer::L1,M
     real(8),dimension(:)::QL1,ZM
     real(8),dimension(:,:)::Zaa,RcoffSin,RcoffCos
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer::II,KK
     integer::ierr
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     real(8),allocatable,dimension(:)::Za
   	
	 ! бшдекъел оюлърэ онд люяяхб
	 allocate(Za(L1-IABS(M)+1),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_CALCULATION_COEFFICIENT_ALFAZX'
	    write(*,*) 'MEMORY ON THE FILE "Za" IS NOT SELECTED'
	    stop 
	 endif

     
     ! гюмскъел оепед пюявернл
	 QL1=0.D0
	 ZM=0.D0
	 Za=0.D0
     Zaa=0.D0
	 RcoffSin=0.D0
	 RcoffCos=0.D0

	 ! онксвюел йнщттхжхемрш опедярюбкемхъ цюплнмхйх пюгкнцюелни тсмйжхх 
     call EFSH_COEFFICIENT_TETA_FUNCTION_DIFFER(L1,IABS(M),QL1)
	 ! онксвюел йнщттхжемрш ахмнлю
	 call EFSH_BINOMIAL_COEFFICIENT(IABS(M),ZM)
	 ! онксвюел йнщтхжхемрш опедярюбкемхъ sin(x)^i=SUM(Rcoffp*cos(px))+SUM(Rcoffq*sin(qx))  
     call EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_SINN(2*IABS(M)+1,RcoffSin)
     ! онксвюел йнщттхжхемрш опедярюбкемхъ cos(x)^i=SUM(Rcoffp*cos(px))   
	 call EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COS(L1-IABS(M),RcoffCos)
     ! гюонкмъел люяяхб ахмнлхюкэмшу йнщттхжхемрнб 
     DO II=1,L1-IABS(M)+1
        call EFSH_BINOMIAL_COEFFICIENT(II-1,Za)
	    do KK=1,II
           Zaa(II,KK)=Za(KK)
	    enddo
	 ENDDO


     ! сдюкъел люяяхб хг оюлърх
	 deallocate(Za,stat=ierr)
	 if(ierr/=0) then
	    write(*,*) 'EFSH_CALCULATION_COEFFICIENT_ALFAZX'
        write(*,*) 'THE FILE "Za" IS NOT REMOVED FROM MEMORY'
	    stop 
     endif


 	 return
    end subroutine EFSH_CALCULATION_COEFFICIENT_ALFAZX



  ! SUB-PROGRAMS FOR REPLACING FUNCTIONS INCLUDED IN THE RADIAL PART OF THE WAVE FUNCTION
  ! POLYNOMA OF THE SECOND ORDER
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! IND1-FIRST ARRAY INDEX XlimZFK, AcoffApro
  ! IND2-SECOND ARRAY INDEX XlimZFK, AcoffApro
  ! L-ORBITAL MOMENT
  ! Nst-power of common factor r ^ Nst * sum (r ^ i * exp (-ALFA * r))
  ! RA-VALUE OF RADIUS IN THE NEW COORDINATE SYSTEM IN WHICH WE CARRY OUT THE CALCULATION
  ! A- OFFSET FROM THE ATOMIC SYSTEM OF COORDINATES
  ! INITIAL FUNCTION APPROXIMATION BY R
  ! Ninterval-number of intervals (function approximation by r)
  ! NpolA-POLYNOMA DEGREE WITH AcoffPolinom coefficients
  ! XlimF (2, Ninterval) - ARRAY OF APROXIMATION INTERVAL BOUNDARIES
  ! ALFA (3, Ninterval) - ARRAY OF SECOND-ORDER POLYNOMA COEFFICIENTS ALFA-EXPOTENTIAL DEGREE
  ! NInt-NUMBER OF INTERVALS FIT IN THE INTEGRATION INTERVAL
  ! NInterH-NUMBER OF THE STARTING INTERVAL OF APPROXIMATION OF THE INTERVAL IN THE INTEGRATION INTERVAL
  ! NInterK-NUMBER OF THE END INTERVAL OF APPROXIMATION OF THE INTERVAL INTO INTEGRATION INTERVAL
  ! NumbreIntTeta-NUMBER OF INTERVALS OF APPROXIMATION OF FUNCTIONS BY ANGLE OF INTEGRATION
  ! XlimZFK (IND1, IND2,2, NInt, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION INTERVALS BY THE INTEGRATION ANGLE
  ! AcoffApro (IND1, IND2,3, NInt, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION COEFFICIENTS ON APPROXIMATION INTERVALS
  subroutine EFSH_PARAMETR_FUNCTION_ALFAZX(IND1,IND2,L,Nst,RA,A,Ninterval,NpolA,XlimF,ALFA,NIntMas,NInterHMas,NInterKMas,NumbreIntTeta,XlimZFK,AcoffApro)
    implicit none
   
	integer::L,Nst,Ninterval,NumbreIntTeta
	integer::NpolA,IND1,IND2
	real(8)::RA,A
	integer,dimension(:,:)::NIntMas,NInterHMas,NInterKMas
    real(8),dimension(:,:)::XlimF,ALFA
    real(8),dimension(:,:,:,:,:,:)::XlimZFK,AcoffApro
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::IIDZ,IJDZ,ISUMDF,ISDFDA,ierr,NInt,NInterH,NInterK
    real(8)::RRF1,RRF2 
    real(8),allocatable,dimension(:,:)::RXlimAA,PXlimZF,APcoffZF



   


	! щрюо 1. нопедекъел мнлепю хмрепбюкнб б йнрнпше оноюдючр цпюмхжш хмрецпхпнбюмхъ
	DO IIDZ=1,Ninterval
	   ! нопедекъел мхфмчч цпюмхжс
	   RRF1=DABS(RA-DABS(A)) 
	   IF(RRF1.GE.XlimF(1,IIDZ).AND.RRF1.LE.XlimF(2,IIDZ)) THEN
            NInterH=IIDZ
	   ENDIF
	   ! нопедекъел бепумчч цпюмхжс
         RRF1=DABS(RA+DABS(A)) 
	   IF(RRF1.GE.XlimF(1,IIDZ).AND.RRF1.LE.XlimF(2,IIDZ)) THEN
            NInterK=IIDZ
	   ENDIF
	ENDDO  


	! нопедекъел вхякн хмрепбюкнб 
      NInt=NInterK-NInterH+1
	!WRITE(6,*) 'INTN',NInterH,NInterK,NInt
	!WRITE(6,*) XlimF(1,NInterH),XlimF(2,NInterH),DABS(RA-DABS(A)) 
	!WRITE(6,*) XlimF(1,NInterK),XlimF(2,NInterK),DABS(RA+DABS(A)) 

    !бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RXlimAA(2,NInt),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
	  write(*,*) 'MEMORY ON THE FILE "RXlimAA" IS NOT SELECTED'
	  stop 
	endif
	! люйяхлюкэмне вхякн хмрепбюкнб 100-вхякн хмрепбюкнб юоопнйяхюжхх он сцкс 
	! нопедекъеряъ б ондопнцпюлле EFSH_PARAMETR_FUNCTION_RDFUN 
	allocate(PXlimZF(2,100),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
	  write(*,*) 'MEMORY ON THE FILE "PXlimZF" IS NOT SELECTED'
	  stop 
	endif
	allocate(APcoffZF(3,100),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
	  write(*,*) 'MEMORY ON THE FILE "APcoffZF" IS NOT SELECTED'
	  stop 
	endif

    ! гюмскъел оепед хяонкэгнбюмхел
	PXlimZF=0.D0
    APcoffZF=0.D0


    
	
      

	! щрюо 2. тнплхпсел люяяхб опедекнб хмрепбюкю
    RXlimAA=0.D0
	! бшъямъел вхякн хмрепбюкнб
    IF(NInt.EQ.1) THEN
       IF(A.GT.0.D0) THEN
	      RXlimAA(1,NInt)=0.D0
          RXlimAA(2,NInt)=3.14159265358979D0
         ELSE
	      RXlimAA(1,NInt)=3.14159265358979D0
          RXlimAA(2,NInt)=0.D0
         ENDIF
	  ELSE
        ! оепбши хмрепбюк
	    IF(A.GT.0.D0) THEN
	       RXlimAA(1,1)=0.D0
          ELSE
	       RXlimAA(1,1)=3.14159265358979D0
        ENDIF
	    ! бшъямъел рхо ялеыемхъ 
	    IF(A.GT.0.D0) THEN
	       RRF1=(RA**2+DABS(A)**2-XlimF(2,NInterH)**2)/(2.D0*RA*DABS(A)) 
           RXlimAA(2,1)=DACOS(RRF1)
          ELSE
	       RRF1=(XlimF(2,NInterH)**2-(RA**2+DABS(A)**2))/(2.D0*RA*DABS(A)) 
           RXlimAA(2,1)=DACOS(RRF1)
	    ENDIF 
	    ! опнлефсрнвмше хмрепбюкш
        ! бшъямъел рхо ялеыемхъ 
	    IF(A.GT.0.D0) THEN 
	       DO IIDZ=2,NInt-1
              RRF1=(RA**2+DABS(A)**2-XlimF(1,NInterH+IIDZ-1)**2)/(2.D0*RA*DABS(A)) 
              RXlimAA(1,IIDZ)=DACOS(RRF1)
              RRF1=(RA**2+DABS(A)**2-XlimF(2,NInterH+IIDZ-1)**2)/(2.D0*RA*DABS(A)) 
              RXlimAA(2,IIDZ)=DACOS(RRF1)
	       ENDDO
	      ELSE
	  	   DO IIDZ=2,NInt-1
              RRF1=(XlimF(1,NInterH+IIDZ-1)**2-(RA**2+DABS(A)**2))/(2.D0*RA*DABS(A)) 
              RXlimAA(1,IIDZ)=DACOS(RRF1)
              RRF1=(XlimF(2,NInterH+IIDZ-1)**2-(RA**2+DABS(A)**2))/(2.D0*RA*DABS(A)) 
              RXlimAA(2,IIDZ)=DACOS(RRF1)
	       ENDDO
	    ENDIF 
	    ! онякедмхи хмрепбюк
	    ! бшъямъел рхо ялеыемхъ 
	    IF(A.GT.0.D0) THEN
	       RRF1=(RA**2+DABS(A)**2-XlimF(1,NInterK)**2)/(2.D0*RA*DABS(A)) 
           RXlimAA(1,NInt)=DACOS(RRF1)
           RXlimAA(2,NInt)=3.14159265358979D0
          ELSE
	       RRF1=(XlimF(1,NInterK)**2-(RA**2+DABS(A)**2))/(2.D0*RA*DABS(A)) 
           RXlimAA(1,NInt)=DACOS(RRF1)
           RXlimAA(2,NInt)=0.D0
	    ENDIF 
         
	ENDIF

   
	
	! щрюо 3. юоопнйяхлхпсел тсмйжхх
	! жхйк он хмрепбюкюл
	ISUMDF=0
	DO IIDZ=NInterH,NInterK
	   ISUMDF=ISUMDF+1
	   RRF1=RXlimAA(1,ISUMDF)
	   RRF2=RXlimAA(2,ISUMDF)
       ! жхйк он тсмйжхъл
	   DO IJDZ=1,NpolA+1
	      ! WRITE(6,*) 'SSSSST',L,Nst,IJDZ-1
          ! ондопнцпюллю юоопнйяхлюжхх тсмйжхх
          call EFSH_PARAMETR_FUNCTION_RDFUN(L-Nst+IJDZ-1,IIDZ,RA,A,ALFA,RRF1,RRF2,NumbreIntTeta,PXlimZF,APcoffZF) 
          ! гюохяшбюел пегскэрюрш юоопнйяхлюжхх 
          DO ISDFDA=1,NumbreIntTeta
	         ! опедекш
	         XlimZFK(IND1,IND2,1,ISUMDF,IJDZ,ISDFDA)=PXlimZF(1,ISDFDA)
             XlimZFK(IND1,IND2,2,ISUMDF,IJDZ,ISDFDA)=PXlimZF(2,ISDFDA)
	         !WRITE(6,*) 'PRED',XlimZFK(IND1,IND2,1,ISUMDF,IJDZ,ISDFDA),XlimZFK(IND1,IND2,2,ISUMDF,IJDZ,ISDFDA)
	         ! йнщттхжхемрш юоопнйяхлюжхх
             AcoffApro(IND1,IND2,1,ISUMDF,IJDZ,ISDFDA)=APcoffZF(1,ISDFDA)
	         AcoffApro(IND1,IND2,2,ISUMDF,IJDZ,ISDFDA)=APcoffZF(2,ISDFDA)
             AcoffApro(IND1,IND2,3,ISUMDF,IJDZ,ISDFDA)=APcoffZF(3,ISDFDA)
             !WRITE(6,*) AcoffApro(IND1,IND2,1,ISUMDF,IJDZ,ISDFDA),AcoffApro(IND1,IND2,2,ISUMDF,IJDZ,ISDFDA),AcoffApro(IND1,IND2,3,ISUMDF,IJDZ,ISDFDA)
		  ENDDO
	   ENDDO 
    ENDDO
      
	! гюохяшбюел дюммше 
	NIntMas(IND1,IND2)=NInt
	NInterHMas(IND1,IND2)=NInterH
	NInterKMas(IND1,IND2)=NInterK

   
	


	! сдюкемхе люяяхбнб хг оълърх  
	deallocate(RXlimAA,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
       write(*,*) 'THE FILE "RXlimAA" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif
    deallocate(PXlimZF,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
       write(*,*) 'THE FILE "PXlimZF" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif
	deallocate(APcoffZF,stat=ierr)
	if(ierr/=0) then
	  write(*,*) 'EFSH_PARAMETR_FUNCTION_ALFAZX'
      write(*,*) 'THE FILE "APcoffZF" IS NOT REMOVED FROM MEMORY'
	  stop 
	endif


    return  
  end subroutine EFSH_PARAMETR_FUNCTION_ALFAZX








! SUBPROGRAM FOR CALCULATING THE INTEGRAL OF THE DECOMPOSITION COEFFICIENT INCLUDED IN THE COMPOSITION
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! IND1-FIRST INDEX
! IND2-SECOND INDEX
! Lx-ORBITAL MOMENT OF THE EXPANDABLE FUNCTION
! L-ORBITAL MOMENT
! M-PROJECTION OF THE MOMENT
! NS1-PARAMETER OF INTEGRAL (INDICATOR)
! TL (L-IABS (M) +1) - COEFFICIENTS OF REPRESENTATION OF THETA FUNCTIONS IN THE FORM (sin (x)) ^ m * SUM (Ti * cos (i * x))
! ZZcoff (IABS (M) +1) - BINOMIAL COEFFICIENTS FOR DEGREE IABS (M)
! RcoffCos (NS + 1, NS1 + 1) - ARRAY OF REPRESENTATION COEFFICIENTS (COS (X)) ^ NS1 = SUM (Bp * cos (p * x))
! RcoffSin (2.2 * IABS (M) +2) - ARRAY OF REPRESENTATION COEFFICIENTS (SIN (X)) ^ (2 * IABS (M) +1) = SUM (Ap * sin (p * x)) + SUM (Bp * cos (p * x))
! NumbreInter-NUMBER OF INTERVALS INCLUDED IN THE INTERVAL OF INTEGRATION BY R
! NpolA-POLYNOMA DEGREE WITH AcoffPolinom coefficients
! NIntH-NUMBER OF START INTERVAL ON R
! NIntK-R END INTERVAL NUMBER
! NumbreIntTeta-NUMBER OF APPROXIMATION INTERVALS BY ANGLE
! XlimTeta (IND1, IND2,2, NumbreInter, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION INTERVALS BY THE INTEGRATION ANGLE
! AcoffAproTeta (IND1, IND2,3, NumbreInter, NpolA + 1, NumbreIntTeta) - ARRAY OF APPROXIMATION COEFFICIENTS AT APPROXIMATION INTERVALS
! AcoffPolinom (NpolA + 1, Ninterval) - ARRAY OF COEFFICIENTS OF THE POLYNOM IN THE FUNCTION
 	real(8) function EFSH_COEFFICIENT_INTEGRAL_HARMONICS_ALFAZX(IND1,IND2,Lx,L,M,NS1,TL,ZZcoff,RcoffCos,RcoffSin,NumbreInter,NpolA,NIntH,NIntK,NumbreIntTeta,XlimTeta,AcoffAproTeta,AcoffPolinom)
      implicit none
	
	  integer::L,Lx,M,NS1,NumbreInter,NIntH,NIntK,NumbreIntTeta,NpolA
	  integer::IND1,IND2
      real(8),dimension(:)::TL,ZZcoff
      real(8),dimension(:,:)::RcoffSin,RcoffCos,AcoffPolinom
      real(8),dimension(:,:,:,:,:,:)::XlimTeta,AcoffAproTeta
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer::ISD0,ISD1,ISD2,ISD3,ISD4,ISD01,ierr
      integer::NSP1,NSP2,NSP3,NSP4
      real(8)::SUMMM,Splus,Sminus,RFGH,XHlimZ,XKlimZ,RK,RH,AAASSS
	  real(8)::AAA,BBB,CCC,SUMAP,SUMBP,SUMCP,SUMAM,SUMBM,SUMCM,RA,RB,RC
    
      
      SUMMM=0.D0
	  NSP1=0
	  ! жхйк он мнлепюл хмрепбюкнб  
      DO ISD0=NIntH,NIntK
	     NSP1=NSP1+1
	     ! жхйк он мнлепюл тсмйжхи 
	     DO ISD1=1,NpolA+1
	        AAASSS=0.D0
            DO ISD01=1,NumbreIntTeta
               XHlimZ=XlimTeta(IND1,IND2,1,NSP1,ISD1,ISD01)
               XKlimZ=XlimTeta(IND1,IND2,2,NSP1,ISD1,ISD01)
               AAA=AcoffAproTeta(IND1,IND2,3,NSP1,ISD1,ISD01)
               BBB=AcoffAproTeta(IND1,IND2,2,NSP1,ISD1,ISD01)
               CCC=AcoffAproTeta(IND1,IND2,1,NSP1,ISD1,ISD01)
			   SUMAP=0.D0
			   SUMBP=0.D0
			   SUMCP=0.D0
			   SUMAM=0.D0
			   SUMBM=0.D0
			   SUMCM=0.D0
		       DO ISD2=1,L-IABS(M)+1   ! COS  
	              DO ISD3=1,NS1+1      ! COS 
	                 DO ISD4=2,2*IABS(M)+2 ! SIN
	                    ! опх ISD4=1- хмрецпюк гюмскъеряъ онщрнлс мювхмюел я ISD4=2
	                    RFGH=TL(ISD2)*RcoffCos(NS1+1,ISD3)*RcoffSin(1,ISD4)
                        ! write(6,*) 'SSSDF',TL(ISD2),RcoffCos(NS1+1,ISD3),RcoffSin(1,ISD4)
					    NSP2=ISD2-1
	                    NSP3=ISD3-1
                        NSP4=ISD4-1
                        ! WRITE(6,*) 'rcoff',RFGH
					    ! write(6,*) 'param',NSP2,NSP3,NSP4,XHlimZ,XKlimZ
					    RA=RFGH*EFSH_INT_SH_ALFA(2,NSP2,NSP3,NSP4,XHlimZ,XKlimZ)
                        RB=RFGH*EFSH_INT_SH_ALFA(1,NSP2,NSP3,NSP4,XHlimZ,XKlimZ)
                        RC=RFGH*EFSH_INT_SH_ALFA(0,NSP2,NSP3,NSP4,XHlimZ,XKlimZ)
	                    !WRITE(6,*) 'RR',RA,RB,RC
					    IF(RA.GE.0.D0) THEN
                           SUMAP=SUMAP+RA
	                      ELSE
                           SUMAM=SUMAM+RA
                        ENDIF
                        IF(RB.GE.0.D0) THEN
                           SUMBP=SUMBP+RB
	                      ELSE
                           SUMBM=SUMBM+RB
                        ENDIF
                        IF(RC.GE.0.D0) THEN
                           SUMCP=SUMCP+RC
	                      ELSE
                           SUMCM=SUMCM+RC
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

               AAASSS=AAASSS+AAA*(SUMAP+SUMAM)+BBB*(SUMBP+SUMBM)+CCC*(SUMCP+SUMCM)
			  ! WRITE(6,*) 'SSS',AAA,(SUMAP+SUMAM),BBB,(SUMBP+SUMBM),CCC,(SUMCP+SUMCM)
			ENDDO
		    SUMMM=SUMMM+AcoffPolinom(ISD1,ISD0)*AAASSS
	      ENDDO
      ENDDO

      EFSH_COEFFICIENT_INTEGRAL_HARMONICS_ALFAZX=SUMMM
        
     
      return
    end function







  ! FUNCTION CHANGE REGULATIONS R ^ NSS * EXP (-ALFA * R), R = SQRT (RA ^ 2 + A ^ 2 + (-) 2 * RA * A * COS (U))
  ! POLYNOMA OF THE SECOND ORDER
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! NSS DEGREE INDICATOR
  ! NALFA NUMBER OF APPROXIMATION INTERVAL
  ! RAD-FIRST FUNCTION PARAMETER
  ! AX-SECOND FUNCTION PARAMETER
  ! ALFA (3, NALFA) - EXPOTENTIAL FACTOR ARRAY
  ! RangleH-LOWER BOUNDARY OF APPROXIMATION INTERVAL
  ! RangleK-UPPER APPROXIMATE INTERVAL BORDER
  ! Ninterval-NUMBER OF INTERVALS INTO WHICH THE FULL INTERVAL IS BROKEN
  ! Xlimit (2, Ninterval) - ARRAY OF INTERVAL BOUNDARIES
  ! Xlimit (1, Ninterval) -BEGIN INTERVAL
  ! Xlimit (2, Ninterval) - END OF INTERVAL
  ! RcoffFUN (3, Ninterval) - ARRAY OF POLYNOMA COEFFICIENTS
  ! RcoffFUN (1, Ninterval) - COEFFICIENT CORRESPONDS TO ZERO POLYNOM DEGREE
  ! RcoffFUN (2, Ninterval) - COEFFICIENT CORRESPONDS TO THE FIRST DEGREE OF POLYNOMA
  ! RcoffFUN (3, Ninterval) - THE COEFFICIENT CORRESPONDS TO THE SECOND DEGREE OF POLYNOMA
  subroutine EFSH_PARAMETR_FUNCTION_RDFUN(NSS,NALFA,RAD,AX,ALFA,RangleH,RangleK,Ninterval,Xlimit,RcoffFUN) 
	implicit none
      
	integer::NSS,NALFA,Ninterval
	real(8)::RAD,AX,RangleH,RangleK
    real(8),dimension(:,:)::ALFA,Xlimit,RcoffFUN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer::IZXOP,ierr,NpointsA
    real(8)::Xpoint,Hdasd
	integer,allocatable,dimension(:)::NyzlovA
    real(8),allocatable,dimension(:)::XAproc,FunAproc,RcoffPolA
      


   
      !бшдекъел оюлърэ дкъ люяяхбнб
	allocate(XAproc(100),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
	  write(*,*) 'MEMORY ON THE FILE "XAproc" IS NOT SELECTED'
	  stop 
	endif
	allocate(FunAproc(100),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
	  write(*,*) 'MEMORY ON THE FILE "FunAproc" IS NOT SELECTED'
	  stop 
	endif
	allocate(NyzlovA(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
	  write(*,*) 'MEMORY ON THE FILE "NyzlovA" IS NOT SELECTED'
	  stop 
	endif
	allocate(RcoffPolA(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
	  write(*,*) 'MEMORY ON THE FILE "RcoffPolA" IS NOT SELECTED'
	  stop 
	endif

	! гюмскъел оепед пюявернл
	XAproc=0.D0
    FunAproc=0.D0
	NyzlovA=0
    RcoffPolA=0.D0

    ! вхякн рнвей опхундъыхуяъ мю онкмши хмрепбюк
	NpointsA=31  !101   !211   !1571  
    ! гюохяшбюел вхякн хмрепбюкнб
	Ninterval=(NpointsA-1)/2

    ! ьюц онксвемхъ гмювемхи тсмйжхх
    Hdasd=(RangleK-RangleH)/float(NpointsA-1)
   	! онксвюел мюанп гмювемхи тсмйжхх 
	!WRITE(6,*) 'NALFA',NALFA,ALFA(1,NALFA),ALFA(2,NALFA),ALFA(3,NALFA)
	DO IZXOP=1,NpointsA
       Xpoint=RangleH+Hdasd*float(IZXOP-1)
       IF(IZXOP.EQ.NpointsA) THEN 
	      Xpoint=RangleK
	   ENDIF
	   XAproc(IZXOP)=Xpoint
       FunAproc(IZXOP)=RDFEXP(NSS,NALFA,RAD,AX,ALFA,Xpoint)
	   !  WRITE(6,*) 'YFUN',XAproc(IZXOP), FunAproc(IZXOP)
	ENDDO

    ! онксвюел гмювемхъ йнщттхжхемрнб онкхмнлю
	DO IZXOP=1,Ninterval
       NyzlovA(1)=2*IZXOP-1
	   NyzlovA(2)=2*IZXOP
	   NyzlovA(3)=2*IZXOP+1
	   ! НОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ он лерндс йпюлепю
	   call EFSH_COEFFICIENT_POLINOM_KRAMERA2(NyzlovA,XAproc,FunAproc,RcoffPolA)
	   ! гюохяшбюел йнщттхжхемрш онкхмнлю
	   RcoffFUN(3,IZXOP)=RcoffPolA(1)
	   RcoffFUN(2,IZXOP)=RcoffPolA(2)
	   RcoffFUN(1,IZXOP)=RcoffPolA(3)
	   ! WRITE(6,*)  'PRED',XAproc(2*IZXOP-1),XAproc(2*IZXOP+1)
	   ! WRITE(6,*)  RcoffPolA(3),RcoffPolA(2),RcoffPolA(1)
	   ! гюохяшбюел опедекш хмрепбюкнб
	   Xlimit(1,IZXOP)=XAproc(2*IZXOP-1)
	   Xlimit(2,IZXOP)=XAproc(2*IZXOP+1)
    ENDDO

   
	! сдюкемхе люяяхбнб хг оълърх  
	deallocate(XAproc,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
       write(*,*) 'THE FILE "XAproc" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif
	deallocate(FunAproc,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
       write(*,*) 'THE FILE "FunAproc" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif
	deallocate(NyzlovA,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
       write(*,*) 'THE FILE "NyzlovA" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif
	deallocate(RcoffPolA,stat=ierr)
	if(ierr/=0) then
	   write(*,*) 'EFSH_PARAMETR_FUNCTION_RDFUN'
       write(*,*) 'THE FILE "RcoffPolA" IS NOT REMOVED FROM MEMORY'
	   stop 
	endif

      
	
	return  
  end subroutine EFSH_PARAMETR_FUNCTION_RDFUN


! SUB-PROGRAM FOR CALCULATING FUNCTION VALUES R ^ N1 * EXP (-ALFA * R), R = SQRT (RA ^ 2 + A ^ 2 + (-) 2 * RA * A * COS (U))
! N1-DEGREE OF RADIAL MULTIPLIER
! N2-INTERVAL NUMBER
! RA-VALUE OF RADIUS IN THE NEW COORDINATE SYSTEM (MOLECULAR COORDINATE SYSTEM)
! ALFA (3, N2) - THE ARRAY OF VALUES OF THE APPROXIMATION COEFFICIENTS BY THE SECOND-ORDER POLYNOMA
! XAS-ANGLE VALUE (RADIANS)
  real(8) function RDFEXP(N1,N2,RA,A,ALFA,XAS)
    implicit none
	integer::N1,N2
	real(8)::RA,A,XAS,RRR,ALFASSD
    real(8),dimension(:,:)::ALFA
	
	IF(A.GT.0.D0) THEN 
	   IF(DABS(RA**2+DABS(A)**2-2.D0*RA*DABS(A)*DCOS(XAS)).LT.10.D0**(-10)) THEN
           RRR=0.D0
		  ELSE
		   RRR=DSQRT(RA**2+DABS(A)**2-2.D0*RA*DABS(A)*DCOS(XAS))
	   ENDIF 
	  ELSE
       RRR=DSQRT(RA**2+DABS(A)**2+2.D0*RA*DABS(A)*DCOS(XAS))
	ENDIF
   
    ALFASSD=-(ALFA(1,N2)*RRR+ALFA(2,N2)*RRR*RRR+ALFA(3,N2)*RRR*RRR*RRR)
    RDFEXP=RRR**N1*DEXP(ALFASSD)
   	return
  end function





	




     

! INTEGRAL CALCULATION SUBPROGRAM
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! NpolAA-POLYNOMA DEGREE
! NSP2X- FREQUENCY OF HARMONIC FUNCTION
! NSP3X- HARMONIC FUNCTION FREQUENCY
! NSP4X- HARMONIC FUNCTION FREQUENCY
! XH-LOWER LIMIT OF INTEGRATION
! XK-UPPER INTEGRATION LIMIT
    real(8) function EFSH_INT_SH_ALFA(NpolAA,NSP2X,NSP3X,NSP4X,XH,XK)
      implicit none
      integer::NpolAA,NSP2X,NSP3X,NSP4X
      real(8)::XH,XK
 	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::NXZ1,NXZ2,NXZ3,NXZ4
      real(8)::SUM

	  NXZ1=NSP4X+NSP2X-NSP3X
      NXZ2=NSP2X+NSP3X-NSP4X
      NXZ3=NSP3X+NSP4X-NSP2X
      NXZ4=NSP3X+NSP4X+NSP2X

      SUM=Int_sin(NpolAA,NXZ1,XH,XK)
      SUM=SUM-Int_sin(NpolAA,NXZ2,XH,XK)
      SUM=SUM+Int_sin(NpolAA,NXZ3,XH,XK)
      SUM=SUM+Int_sin(NpolAA,NXZ4,XH,XK)
	 
	  EFSH_INT_SH_ALFA=0.25D0*SUM
      return
    end function



    ! ондопнцпюллю хмрецпхпнбюмхъ тсмйжхи бхдю X^Nvf*SIN(NSIN*X)
    ! Nvf-яреоемэ
	! NSIN-вюярнрю SIN
    ! XXH,XXK-цпюмхжш хмрепбюкю  
	real(8) function Int_sin(Nvf,NSIN,XXH,XXK)
      implicit none
	  integer::Nvf,NSIN
	  real(8)::XXH,XXK,RCDJK,RCDJH
	
      IF(NSIN.EQ.0) THEN
         Int_sin=0.D0
         return
	  ENDIF

	  IF(Nvf.EQ.0) THEN 
	     Int_sin=(DCOS(FLOAT(NSIN)*XXH)-DCOS(FLOAT(NSIN)*XXK))/FLOAT(NSIN)
      ENDIF

	  IF(Nvf.EQ.1) THEN 
         RCDJK=DSIN(FLOAT(NSIN)*XXK)-FLOAT(NSIN)*XXK*DCOS(FLOAT(NSIN)*XXK)
	     RCDJK=RCDJK/FLOAT(NSIN)**2

         RCDJH=DSIN(FLOAT(NSIN)*XXH)-FLOAT(NSIN)*XXH*DCOS(FLOAT(NSIN)*XXH)
         RCDJH=RCDJH/FLOAT(NSIN)**2
	     Int_sin=RCDJK-RCDJH
      ENDIF

	  IF(Nvf.EQ.2) THEN 
         RCDJK=2.D0*FLOAT(NSIN)*XXK*DSIN(FLOAT(NSIN)*XXK)-((FLOAT(NSIN)*XXK)**2-2.D0)*DCOS(FLOAT(NSIN)*XXK)
	     RCDJK=RCDJK/FLOAT(NSIN)**3

         RCDJH=2.D0*FLOAT(NSIN)*XXH*DSIN(FLOAT(NSIN)*XXH)-((FLOAT(NSIN)*XXH)**2-2.D0)*DCOS(FLOAT(NSIN)*XXH)
	     RCDJH=RCDJH/FLOAT(NSIN)**3
	     Int_sin=RCDJK-RCDJH
      ENDIF


	  return
    end function


     
      







! SUB-PROGRAM FOR OBTAINING THE VALUE OF THETA FUNCTION (RECURRENT METHOD)
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L-ORBITAL MOMENT FUNCTION
! M-PROJECTION OF THE ORBITAL MOMENTUM (THE PROGRAM CONSIDERS THE CASE (M> = 0)
! N-NUMBER OF POINTS WHICH WE FIND VALUES (IN THE INTERVAL (0, PI))
! X (N) -MASSIF OF ARGUMENT VALUES
! Y (N) -ARRAY OF FUNCTION VALUES
	  
    subroutine EFSH_VALUE_TETA_FUNCTION(L,M,N,X,Y)
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
	    RSDXF=DSQRT(EFSH_FACTORIAL(2*IABS(M)))/EFSH_FACTORIAL(IABS(M))
	    RSDXX=(DSIN(X(IAA1)*0.5D0)*DCOS(X(IAA1)*0.5D0))**IABS(M)
        RSDXX=RSDXX*RSDXF/DSQRT(2.D0)
	    Rcoff1=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+1))*RSDXX
        !  якюцюелне Teta M+1,M 
	    RSDXF=DSQRT(EFSH_FACTORIAL(2*IABS(M)+1))/EFSH_FACTORIAL(IABS(M))
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
    end subroutine EFSH_VALUE_TETA_FUNCTION
      



! SUBPROGRAM FOR CALCULATING THE COEFFICIENTS OF THE POLYNOMA OF THE 2ND ORDER
! BY CRAMERA METHOD
! DESCRIPTION OF SUBPROGRAM PARAMETERS
      
! Npoints (3) - ARRAY OF POINT NUMBERS
! R (NN) - ARGUMENT VALUE ARRAY
! FUN (NN) -ARRAY OF FUNCTION VALUES
! Acoff (3) -ARRAY OF POLYNOMA COEFFICIENTS
	subroutine EFSH_COEFFICIENT_POLINOM_KRAMERA2(Npoints,R,FUN,Acoff)
      implicit none
      
	  real(8)::X1R,X2R,X3R,Y1R,Y2R,Y3R,RDX1,RDX2,RDX3,RDX4
      integer,dimension(:)::Npoints
	  real(8),dimension(:)::R,FUN,Acoff


      ! мюундхл йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
      X1R=R(Npoints(1))
	  X2R=R(Npoints(2))
	  X3R=R(Npoints(3))
	  Y1R=FUN(Npoints(1))
	  Y2R=FUN(Npoints(2))
	  Y3R=FUN(Npoints(3))
	  
	  ! нопедекхрекэ яхярелш
	  RDX1=X1R**2*(X2R-X3R)+X2R**2*(X3R-X1R)+X3R**2*(X1R-X2R)
      ! нопедекхрекэ йнщттхжхемрю Acoff(1)
	  RDX2=Y1R*(X2R-X3R)+Y2R*(X3R-X1R)+Y3R*(X1R-X2R)
      ! нопедекхрекэ йнщттхжхемрю Acoff(2)
	  RDX3=Y1R*(X3R**2-X2R**2)+Y2R*(X1R**2-X3R**2)+Y3R*(X2R**2-X1R**2)
      ! нопедекхрекэ йнщттхжхемрю Acoff(3)
	  RDX4=X1R**2*(X2R*Y3R-Y2R*X3R)-X1R*(X2R**2*Y3R-Y2R*X3R**2)+Y1R*(X2R**2*X3R-X3R**2*X2R)
	
      ! гюохяшбюел йнщттхжхемрш
      Acoff(1)=RDX2/RDX1
      Acoff(2)=RDX3/RDX1
      Acoff(3)=RDX4/RDX1
         
	  return  
	end subroutine EFSH_COEFFICIENT_POLINOM_KRAMERA2




	
! SUBPROGRAM FOR CALCULATION OF N-ORDER POLYNOMA COEFFICIENTS
! BY CRAMERA METHOD
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! NL-DEGREE OF COEFFICIENT
! Np-POLYNOMA ORDER
! Npoints (N + 1) - ARRAY OF POINT NUMBERS
! R (NN) - ARGUMENT VALUE ARRAY
! FUN (NN) -ARRAY OF FUNCTION VALUES
! Acoff (N + 1) -ARRAY OF POLYNOMA COEFFICIENTS
! IN THE ARRAY OF COEFFICIENTS THE NUMBER OF THE ELEMENT CORRESPONDS THE NUMBER OF THE DEGREE BY THE RULE Ndegree = Nthe number of the element-1
! Filling occurs by decreasing the degree of the polynomial
	subroutine EFSH_COEFFICIENT_POLINOM_KRAMERA(NL,Np,Npoints,R,FUN,Acoff)
      implicit none

	integer::NL,Np,IJK,JJK,IJK1,IJK2,ierr 
	real(8)::RSysDet
      integer,dimension(:)::Npoints
	real(8),dimension(:)::R,FUN,Acoff
	real(8),allocatable,dimension(:,:)::Rmatix,RMatrixCoff

	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(Rmatix(Np+1,Np+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_POLINOM_KRAMERA'
	write(*,*) 'MEMORY ON THE FILE "Rmatix" IS NOT SELECTED'
	stop 
	endif
	allocate(RMatrixCoff(Np+1,Np+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_POLINOM_KRAMERA'
	write(*,*) 'MEMORY ON THE FILE "RMatrixCoff" IS NOT SELECTED'
	stop 
	endif


	! ГЮМСКЪЕЛ ОЕПЕД ХЯОНКЭГНБЮМХЕЛ
      Rmatix=0.d0
      RMatrixCoff=0.d0
	Acoff=0.d0

	! щрюо 1 ТНПЛХПСЕЛ ЛЮРПХЖС ЯХЯРЕЛШ КХМЕИМШУ СПЮБМЕМХИ
	do IJK=1,Np+1
         do JJK=Np+1,1,-1
            Rmatix(IJK,JJK)=R(Npoints(IJK))**(JJK-1+NL)
	   enddo
	enddo

      

      ! щрюо 2 МЮУНДХЛ НОПЕДЕКХРЕКЭ ЯХЯРЕЛШ
      RSysDet=EFSH_Determenant(Np+1,Rmatix)

	

      ! щрюо 3 МЮУНДХЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ
	do IJK=1,Np+1
         ! ТНПЛХПСЕЛ ЛЮРПХЖС ДКЪ ПЮЯВЕРЮ ЙНПМЪ ЯХЯРЕЛШ КХМЕИМШУ СПЮБМЕМХИ
         do IJK1=1,Np+1
            if(IJK1.EQ.IJK) then
		     do IJK2=1,Np+1 
	            RMatrixCoff(IJK2,IJK1)=FUN(Npoints(IJK2))
		     enddo  
		    else
			 do IJK2=1,Np+1 
	            RMatrixCoff(IJK2,IJK1)=Rmatix(IJK2,IJK1)
		     enddo  	  
	      endif
	   enddo
	  
	  ! IJK-ши йнпемэ яхярелш спюбмемхи нопедекъеряъ бшпюфемхел
          
	   Acoff(Np+2-IJK)=EFSH_Determenant(Np+1,RMatrixCoff)/RSysDet
	  ! WRITE(6,*) 'DD',Acoff(Np+2-IJK)
	enddo



       ! сдюкъел люяяхбш хг оюлърх
      deallocate(Rmatix,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_POLINOM_KRAMERA'
      write(*,*) 'THE FILE "Rmatix" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RMatrixCoff,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_POLINOM_KRAMERA'
      write(*,*) 'THE FILE "RMatrixCoff" IS NOT REMOVED FROM MEMORY'
	stop 
	endif

     
	return  
	end subroutine EFSH_COEFFICIENT_POLINOM_KRAMERA







	! ондопнцпюллю пюяверю нопедекхрекъ люрпхжш	

	real(8) function EFSH_Determenant(n,ra)
  	use msimsl,nouse=>fac
	implicit none
	integer::ierr,n,Lda,Ldfac,iiii,jjjj
	real(8)::det1,det2
	real(8),dimension(:,:)::ra
    real(8),allocatable,dimension(:,:)::a,fac
	integer,allocatable,dimension(:)::ipvt  
	
	! бшдекъел оълърэ онд люяяхбш
	allocate(a(n,n),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma MEMORY ON THE FILE "a" IS NOT SELECTED'
	stop 
	endif   
	allocate(fac(n,n),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma MEMORY ON THE FILE "fac" IS NOT SELECTED'
	stop 
	endif   
      allocate(ipvt(n),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma MEMORY ON THE FILE "ipvt" IS NOT SELECTED'
	stop 
	endif   

	Lda=n;Ldfac=n

	do iiii=1,n
	   do jjjj=1,n
	   a(iiii,jjjj)=ra(iiii,jjjj)
         enddo
      enddo


	!LU-пюгкнфемхе
	call DLufac2(n,a,fac,ipvt)
	call DLfdrg(n,fac,Ldfac,ipvt,det1,det2)
	EFSH_Determenant=det1*(10.D0)**det2 

	! сдюкемхе люяяхбнб хг оълърх 
	deallocate(a,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma THE FILE "a" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(fac,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma THE FILE "fac" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(ipvt,stat=ierr)
	if(ierr/=0) then
      write(*,*) 'cma THE FILE "ipvt" IS NOT REMOVED FROM MEMORY'
	stop 
	endif

	return    
      end function


	subroutine DLufac2(n,a,fac,ipvt)
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
	end subroutine DLufac2
	
	


! SUB-PROGRAM FOR CHECKING THE NORMALIZATION CONDITIONS OF TWO TETE FUNCTIONS
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L1-ORBITAL MOMENT OF FIRST FUNCTION
! L2-ORBITAL MOMENT OF THE SECOND FUNCTION
! M-PROJECTION OF THE MOMENT
! Q1 (L1-M + 1) - ARRAY OF FIRST FUNCTION COEFFICIENTS
! Q2 (L2-M + 1) - ARRAY OF FIRST FUNCTION COEFFICIENTS
      
	real(8) function EFSH_NORMIROVKA_TETA_FUNCTION(L1,L2,M,Q1,Q2)
      implicit none

	integer:: L1,L2,M,I1,I2,I3,I4,ierr
      real(8):: PROSS,PSUMPlus,PSUMMinus
      real(8),dimension(:)::Q1,Q2
      real(8),allocatable,dimension(:)::RcoffO
    
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(RcoffO(IABS(M)+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_NORMIROVKA_TETA_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RcoffO" IS NOT SELECTED'
	stop 
	endif


      ! гюмскъел оепед пюявернл
	RcoffO=0.D0
	
	! БШВХЯКЪЕЛ ЙНЩТТХЖХЕМРШ (ЙЮФДНИ ЯРПНЙХ)  
	call EFSH_BINOMIAL_COEFFICIENT(IABS(M),RcoffO)
	
      PSUMPlus=0.D0
	PSUMMinus=0.D0
      DO I1=1,L1-IABS(M)+1
         DO I2=1,L2-IABS(M)+1
            DO I3=1,IABS(M)+1
	         DO I4=1,IABS(M)+1
                  IF((1+(-1)**(I1+I2+I3+I4-4)).NE.0) THEN
	              IF(Q1(I1)*Q2(I2).NE.0.D0) THEN
				   IF(((-1.D0)**(I3-1)*Q1(I1)*Q2(I2)).GT.0.D0) THEN
      PROSS=Q1(I1)*Q2(I2)*RcoffO(I3)*RcoffO(I4)*2.D0
	PSUMPlus=PSUMPlus+(-1.D0)**(I3-1)*PROSS/FLOAT(I1+I2+I3+I4-3)	                
					 ELSE
      PROSS=Q1(I1)*Q2(I2)*RcoffO(I3)*RcoffO(I4)*2.D0
	PSUMMinus=PSUMMinus+(-1.D0)**(I3-1)*PROSS/FLOAT(I1+I2+I3+I4-3)	                 
	               ENDIF
	              ENDIF
				ENDIF  
               ENDDO
		  ENDDO
	   ENDDO 
      ENDDO 

	EFSH_NORMIROVKA_TETA_FUNCTION=PSUMPlus+PSUMMinus






	
	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(RcoffO,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_NORMIROVKA_TETA_FUNCTION'
      write(*,*) 'THE FILE "RcoffO" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
  
      return
      end function
	
	
	
		
      



! SUBPROGRAM FOR CALCULATING THE COEFFICIENTS OF THE THETA FUNCTION IN THE COMPOSITION OF A SPHERICAL FUNCTION
! (RECURRENT METHOD)
! coefficients of theta function for the case M> = 0
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L-VALUE OF ORBITAL MOMENT
! M-VALUE OF ORBITAL MOMENTUM PROJECTION
! Q (L-M + 1) - ARRAY OF INCLUSIONS IN THE STRUCTURE OF A SPHERICAL FUNCTION (sin (x)) ^ m * SUM (Qi * (cos (x)) ^ i)
	subroutine EFSH_COEFFICIENT_TETA_FUNCTION(L,M,Q)
	implicit none
	
	integer::I,K,L,M,ierr
      real(8)::PSUM,Rcoff1,Rcoff2,PSUMPlus,PSUMMinus
	real(8),dimension(:)::Q
      real(8),allocatable,dimension(:)::Na,Za
      real(8),allocatable,dimension(:,:)::Zaa

	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(Na(L-IABS(M)+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Na" IS NOT SELECTED'
	stop 
	endif
     	allocate(Za(L-IABS(M)+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Za" IS NOT SELECTED'
	stop 
	endif
      allocate(Zaa(L-IABS(M)+1,L-IABS(M)+1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Zaa" IS NOT SELECTED'
	stop 
	endif

      ! гюмскъел оепед пюявернл
      Na=0.D0
      Za=0.D0
	Zaa=0.D0
	      
	! нЯСЫЕЯРБКЪЕЛ ПЮЯВЕР ЙНЩТТХЖХЕМРНБ Na
	PSUM=1.D0
	DO I=L-IABS(M)+1,L+IABS(M)
         PSUM=PSUM*SQRT(float(I))
	ENDDO
	
	Rcoff1=(-1.D0)**IABS(M)*SQRT(float(2*L+1))*PSUM
      Rcoff2=EFSH_FACTORIAL(IABS(M))*2.d0**IABS(M)*DSQRT(2.D0)

	Na(1)=Rcoff1/Rcoff2

      DO I=2,L-IABS(M)+1 
	   K=I-2
	   Rcoff1=-float((L+IABS(M)+K+1)*(L-IABS(M)-K))    
         Rcoff2=float(2*(K+1)*(K+IABS(M)+1))
	   PSUM=Rcoff1/Rcoff2
	   Na(I)=Na(I-1)*PSUM
      ENDDO


      !WRITE(6,*) 'Na coff'
      !DO I=1,L-IABS(M)+1 
	!   write(6,*) Na(I) 
      !ENDDO
      
	!write(6,*)



      ! мюундхл йнщттхжхемрш
	DO I=1,L-IABS(M)+1
         ! БШВХЯКЪЕЛ ЙНЩТТХЖХЕМРШ (ЙЮФДНИ ЯРПНЙХ)  
	   call EFSH_BINOMIAL_COEFFICIENT(I-1,Za)
	   ! ГЮОНКМЪЕЛ ЛЮЯЯХБ АХМНЛХЮКЭМШУ ЙНЩТТХЖХЕМРНБ  
         do K=1,I
            Zaa(I,K)=Za(K)
	   enddo
	ENDDO
 
34545 FORMAT(1X,100(1X,F20.2))

      !write(6,*) 'Bionomial'
	!DO I=1,L-IABS(M)+1
      !   WRITE(6,34545)(Zaa(I,K),K=1,I) 
      !ENDDO
	!WRITE(6,*)
      
      ! мюундхл йнщттхжхемрш рерю тсмйжхх
	DO I=1,L-IABS(M)+1
         ! пюгдекэмне ясллхпнбюмхе ямхфюер онцпеьмнярэ б яксвюе гмюйноепелеммнцн ясллхпнбюмхъ
	   PSUMPlus=0.D0
	   PSUMMinus=0.D0
         DO K=I,L-IABS(M)+1
	      
            IF(Na(K).GT.0.D0) THEN
	  	     PSUMPlus=PSUMPlus+Na(K)*Zaa(K,I)
	        ELSE
               PSUMMinus=PSUMMinus+Na(K)*Zaa(K,I) 
	      ENDIF
	   ENDDO
	   PSUM=PSUMPlus+PSUMMinus
	   IF(DABS(PSUM).LE.(10.D0)**(-9)) PSUM=0.D0
	   Q(I)=PSUM*(-1.D0)**(I-1)
   	ENDDO


      
	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(Na,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
      write(*,*) 'THE FILE "Na" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Za,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
      write(*,*) 'THE FILE "Za" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Zaa,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION'
      write(*,*) 'THE FILE "Zaa" IS NOT REMOVED FROM MEMORY'
	stop 
      endif
	return
      end subroutine EFSH_COEFFICIENT_TETA_FUNCTION




! SUBPROGRAM FOR CALCULATING THE COEFFICIENTS OF THE THETA FUNCTION IN THE COMPOSITION OF A SPHERICAL FUNCTION
! (USING A DIFFERENTIAL REPRESENTATION)
! coefficients of theta function for the case M> = 0
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! L-VALUE OF ORBITAL MOMENT
! M-VALUE OF ORBITAL MOMENTUM PROJECTION
! Q (L-M + 1) - ARRAY OF INCLUSIONS IN THE STRUCTURE OF A SPHERICAL FUNCTION (sin (x)) ^ m * SUM (Qi * (cos (x)) ^ i)
	subroutine EFSH_COEFFICIENT_TETA_FUNCTION_DIFFER(L,M,Q)
	  implicit none
	
	  integer::I,K,L,M,ierr
      real(8)::PSUM,Rcoff1,Rcoff2
      real(8)::PSUMZZ
	  real(8),dimension(:)::Q
    

      ! гюмскъел оепед пюявернл
      Q=0.D0
   
	  PSUM=1.D0
	  DO I=L-IABS(M)+1,L+IABS(M)
         PSUM=PSUM*SQRT(float(I))
	  ENDDO
      	
	  Rcoff1=(-1.D0)**(IABS(M)+L)*SQRT(float(2*L+1))
      Rcoff2=2.d0**L*DSQRT(2.D0)*PSUM


      ! нопедекъел йнщттхжхемрш опедярюбкемхъ
      DO K=0,L
         IF((2*K-(L+IABS(M))).GE.0) THEN
            PSUM=1.D0
	        DO I=(2*K-(L+IABS(M)-1)),2*K
               PSUM=PSUM*float(I)
	        ENDDO
            PSUMZZ=PSUM*Rcoff1/(Rcoff2*EFSH_FACTORIAL(L-K)*EFSH_FACTORIAL(K))
            Q(2*K-(L+IABS(M))+1)=(-1.D0)**K*PSUMZZ
	     ENDIF
	  ENDDO
	

	  return
    end subroutine EFSH_COEFFICIENT_TETA_FUNCTION_DIFFER


! SUB-PROGRAM FOR OBTAINING THE SPECIES CONVERSION COEFFICIENTS
! cos (x) ^ i = SUM (Rcoffp * cos (px))
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! N-MAXIMUM POLYNOMA DEGREE
! Rcoff (N + 1, N + 1) -ARRAY OF COEFFICIENTS PRESENTATION cos (x) ^ i = SUM (Rcoffp * cos (px))
	subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COS(N,Rcoff)
	  implicit none

      integer::N,IYU,JYU
      real(8),dimension(:,:)::Rcoff

      ! гюмскъел оепед пюявернл
	  Rcoff=0.D0
      
	  IF(N.EQ.0) THEN
	     Rcoff(1,1)=1.d0
	     return
	  ENDIF
	
	
	  Rcoff(1,1)=1.d0
	  Rcoff(2,1)=0.d0 
	  Rcoff(2,2)=1.d0

	
	  ! жхйк он яреоемх янS(X)^IYU 
      DO IYU=1,N-1
         ! жхйк он йнщттхжхемрюл мнбнцн онкхмнлю (жхйк он яреоемъл)
	     Rcoff(IYU+2,2)=Rcoff(IYU+1,1) 
	     DO JYU=1,IYU  ! мювхмюеряъ я оепбни яреоемх р.й. "0"-СФЕ ГЮОХЯЮМЮ  
	        Rcoff(IYU+2,JYU+2)=Rcoff(IYU+2,JYU+2)+Rcoff(IYU+1,JYU+1)*0.5D0
	        Rcoff(IYU+2,JYU)=Rcoff(IYU+2,JYU)+Rcoff(IYU+1,JYU+1)*0.5D0   
	     ENDDO
	  ENDDO

  	  return
    end subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COS




! SUB-PROGRAM FOR OBTAINING THE SPECIES CONVERSION COEFFICIENTS
! cos (x) ^ i = SUM (Rcoffp * cos (px))
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! N-MAXIMUM POLYNOMA DEGREE
! RcoffN (N + 1) - COEFFICIENT ARRAY PRESENTATION cos (x) ^ i = SUM (Rcoffp * cos (px))
	subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COSN(N,RcoffN)
	  implicit none

      integer::N,IYU,JYU,ierr
      real(8),dimension(:)::RcoffN
      real(8),allocatable,dimension(:,:)::RSPS

	  !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(RSPS(N+1,N+1),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COSN'
	     write(*,*) 'MEMORY ON THE FILE "RSPS" IS NOT SELECTED'
	     stop 
	  endif

      ! гюмскъел оепед пюявернл
	  RSPS=0.D0
      RcoffN=0.D0
    
	  IF(N.EQ.0) THEN
	     RSPS(1,1)=1.d0
	    ELSE
	
	     RSPS(1,1)=1.d0
	     RSPS(2,1)=0.d0 
	     RSPS(2,2)=1.d0

	     ! жхйк он яреоемх янS(X)^IYU 
         DO IYU=1,N-1
            ! жхйк он йнщттхжхемрюл мнбнцн онкхмнлю (жхйк он яреоемъл)
	        RSPS(IYU+2,2)=RSPS(IYU+1,1) 
	        DO JYU=1,IYU  ! мювхмюеряъ я оепбни яреоемх р.й. "0"-СФЕ ГЮОХЯЮМЮ  
	           RSPS(IYU+2,JYU+2)=RSPS(IYU+2,JYU+2)+RSPS(IYU+1,JYU+1)*0.5D0
		       RSPS(IYU+2,JYU)=RSPS(IYU+2,JYU)+RSPS(IYU+1,JYU+1)*0.5D0   
	        ENDDO
	     ENDDO
	  ENDIF

      ! гюохяшбюел пегскэрюр пюяверю
	  DO IYU=1,N+1
         RcoffN(IYU)=RSPS(N+1,IYU)
      ENDDO

      ! сдюкемхе люяяхбнб хг оълърх   
      deallocate(RSPS,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COSN'
         write(*,*) 'THE FILE "RSPS" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif

      
	  return
    end subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_COSN







! SUBPROGRAM FOR OBTAINING THE POLYNOMY TRANSFORMATION COEFFICIENTS
! sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! N-MAXIMUM POLYNOMA DEGREE
! Rcoff (2, N + 1) COEFFICIENT ARRAY REPRESENTATION sin (x) ^ i = SUM (Rcoffp * cos (px)) + SUM (Rcoffq * sin (qx))
! Rcoff (1, N + 1) COEFFICIENTS AT SIN (PX)
! Rcoff (2, N + 1) COEFFICIENTS AT COS (PX)
	subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_SINN(N,Rcoff)
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
      end subroutine EFSH_COEFFICIENT_TRANSFORMATION_POLINOM_SINN







! SUB-PROGRAM FOR CALCULATING BINOMIAL COEFFICIENTS (EXPRESSION COEFFICIENTS (1 + X) ^ N)
! (RECURRENT METHOD)
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! N-MAXIMUM POLYNOMA DEGREE
! Z (N + 1) -ARRAY OF BINOMIAL COEFFICIENTS
      subroutine EFSH_BINOMIAL_COEFFICIENT(N,Z)
	implicit none

	integer::I,K,N
	real(8),dimension(:)::Z
      
	! ГЮМСКЕМХЕ ОЕПЕД ПЮЯВЕРНЛ
	Z=0.D0
      
	Z(1)=1.D0
	
	do I=2,N+1      
         K=I-1
         Z(I)=float(N-K+1)*Z(I-1)/float(K)
      enddo

      return
      end subroutine EFSH_BINOMIAL_COEFFICIENT


      ! ондопнцпюллю пюяверю тюйрнпхюкю
	! нохяюмхе оюпюлерпнб
      ! NOPI-жекне онкнфхрекэмне вхякн 
	real(8) function EFSH_FACTORIAL(NOPI)
      implicit none

	integer:: NOPI,IBKL
      real(8):: PROSS
      
	IF(NOPI.LT.0) THEN
        WRITE(*,*) 'ERROR FACTORIAL (N<0)'
        READ(*,*)
	  STOP
	ENDIF

      IF(NOPI.EQ.0) THEN
      EFSH_FACTORIAL=1.D0
      return
	ENDIF
      
	PROSS=1.D0
	do IBKL=1,NOPI
         PROSS=PROSS*float(IBKL)
	enddo

	EFSH_FACTORIAL=PROSS

   
      return
      end function
      
! READ AND WRITE RADIAL PARTS SUBROUTINE
! WAVE FUNCTIONS
! DESCRIPTION OF INPUT PARAMETERS OF THE SUBPROGRAM
! KL-PARAMETER (KL = 2-READ, KL = 1-WRITE)
! NF-FUNCTION NUMBER (SHELL NUMBER)
! R (M2Z) - ARRAY OF VALUES OF THE RADIAL PART OF THE WAVE FUNCTION
! RDM (IS * M2Z) - ARRAY OF VALUES OF RADIAL PARTS OF WAVE FUNCTIONS CONFIGURATION       
      subroutine  EFSH_RW(KL,NF,R,RDM,M2Z)
      implicit none
      integer::KL,NF,M2Z,I,J
      real(8),dimension(:)::R,RDM
  	
	    I=(NF-1)*M2Z
        do J=1,M2Z
           if(KL.EQ.2)  R(J)=RDM(I+J)    
           if(KL.EQ.1)  RDM(I+J)=R(J)   
        enddo
	    return
      end subroutine EFSH_RW


   

 

  ! ондопнцпюллю бшдювх гмювемхи юопнйяхлхпсчыеи тсмйжхх
  subroutine EFSH_APROCSIM_TETA(N,Ncoff,M,Rcoff,XX,YY) 
   implicit none
   integer::N,M,IASD,JASD,Ncoff
   real(8)::RHR,SUMCD
   real(8),dimension(:)::XX,YY,Rcoff


   ! гюмскъел оепед пюявернл
   XX=0.D0
   YY=0.D0
      

   ! тнплхпсел люяяхб гмювемхи юпцслемрю
   RHR=3.14159265358979D0/float(N+1)
   ! XX(1)=0.D0
   !XX(N)=3.14159265358979D0
   DO IASD=2,N+1
      XX(IASD-1)=RHR*float(IASD-1)
   ENDDO

   DO IASD=1,N
      SUMCD=0.D0
      DO JASD=0,Ncoff-1
         SUMCD=SUMCD+Rcoff(JASD+1)*DCOS(FLOAT(JASD)*XX(IASD))
      ENDDO
      YY(IASD)=DSIN(XX(IASD))**IABS(M)*SUMCD
   ENDDO
	
  
   return
  end subroutine EFSH_APROCSIM_TETA


! SUBPROGRAM FOR APPROXIMATING THETA FUNCTIONS BY FUNCTIONS OF THE FORM sin (x) ^ M * SUM (Qi * cos (i * x))
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! M-PROJECTION OF ORBITAL MOMENT
! N-NUMBER OF DOTS
! X (N) -MASSIF OF VALUES OF THE TETE FUNCTION ARGUMENT
! Y (N) -ARRAY OF VALUES OF THEETA FUNCTION
! Qcoff (N) -MASSIVE OF VALUES OF THE APPROXIMATION COEFFICIENTS
  subroutine EFSH_COEFFICIENT_TETA_FUNCTION_APRO(M,N,X,Y,Qcoff) 
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
   call EFSH_SYSTEM_LINEAR_EQUATIONS(N,AXXZ,Y,Qcoff)


      
	
   ! сдюкемхе люяяхбнб хг оълърх   
   deallocate(AXXZ,stat=ierr)
   if(ierr/=0) then
  	  write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION_APRO'
      write(*,*) 'THE FILE "AXXZ" IS NOT REMOVED FROM MEMORY'
	  stop 
   endif
   return
  end subroutine EFSH_COEFFICIENT_TETA_FUNCTION_APRO


! SUBPROGRAM FOR APPROXIMATING THETA FUNCTIONS BY FUNCTIONS OF THE FORM sin (x) ^ M * SUM (Qi * cos (x) ^ i)
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! M-PROJECTION OF ORBITAL MOMENT
! N-NUMBER OF DOTS
! X (N) -MASSIF OF VALUES OF THE TETE FUNCTION ARGUMENT
! Y (N) -ARRAY OF VALUES OF THEETA FUNCTION
! Qcoff (N) -MASSIVE OF VALUES OF THE APPROXIMATION COEFFICIENTS
	subroutine EFSH_COEFFICIENT_TETA_FUNCTION_APROTUS(M,N,X,Y,Qcoff) 
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
      write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION_APROTUS'
	write(*,*) 'MEMORY ON THE FILE "AXXZ" IS NOT SELECTED'
	stop 
	endif


      ! гюмскъел оепед пюявернл
      AXXZ=0.D0


	! гюонкмъел люяяхб
      DO JZCVX=1,N
	   RXCVB=X(JZCVX)
	   DO IZCVX=0,N-1
	AXXZ(JZCVX,IZCVX+1)=DSIN(RXCVB)**IABS(M)*DCOS(RXCVB)**IZCVX
         ENDDO
      ENDDO


	! ондопнцпюллю мюунфдемхъ йнпмеи яхярелш кхмеимшу спюбмемхи
	! лернд цюсяяю я опхлемемхел яуелш вюярхвмнцн бшанпю
      call EFSH_SYSTEM_LINEAR_EQUATIONS(N,AXXZ,Y,Qcoff)


      
	
	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(AXXZ,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_COEFFICIENT_TETA_FUNCTION_APROTUS'
      write(*,*) 'THE FILE "AXXZ" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      return
      end subroutine EFSH_COEFFICIENT_TETA_FUNCTION_APROTUS


! SUB-PROGRAM FOR APPROXIMATION OF THE RADIAL PART OF THE WAVE FUNCTION (DISCRETE WAVE FUNCTION)
! DESCRIPTION OF SUBPROGRAM PARAMETERS
! nn-principal quantum number
! l-orbital quantum number
! NpointFUN-number of points
! R (N) -array of radius values
! Rfun (N) -array of function values
! Nraz-parameter increasing by Nraz the initial number of intervals 2 (nn-l)
! Ninterval-number of intervals of approximation
! NpolAR-degree of a polynomial with ARCoffPolinom coefficients
! Xlim (2, Ninterval) -array of limits of approximation intervals
! Xlim (1, Ninterval) -lower limit
! Xlim (2, Ninterval) -upper limit
! A0 (3, Ninterval) -array of coefficients of the approximating function
! (A03 * x ^ 2 + A02 * x + A01) * x ^ (l + 1) * exp (- (ALFA3 * x ^ 2 + ALFA2 * x + ALFA3) * x) * (1 + sum (ak * x ^ k))
! ALFA (3, Ninterval) - the exponential factor of the approximating function depends on the interval
! ALFA (3, Ninterval) - CORRESPONDING TO DEGREE r ^ 2
! ALFA (2, Ninterval) - CORRESPONDING TO DEGREE r
! ALFA (1, Ninterval) - CORRESPONDING TO DEGREE r ^ 0
! ARCoffPolinom (2 + nn-l, Ninterval) -polomial multipliers (nn-l-1 multipliers) added 1 to make the array
! had at least one POLYNOMA PRODUCT RESULT (A03 * x ^ 2 + A02 * x + A01) * (1 + sum (ak * x ^ k))
! EPSfun-error in approximation (maximum deviation from the approximated function)
! EPSr-value of R-radius where maximum deviation is found
! RfunAro (N) -array of values ??of the approximated function (for comparison with the original)
	subroutine EFSH_APPROXIMATION_RADIAL_FUNCTION(nn,l,NpointFUN,R,Rfun,Nraz,Ninterval,NpolAR,Xlim,ALFA,ARCoffPolinom,EPSfun,EPSr,RfunAro) 
      implicit none
       
      integer::nn,l,Npoint,NpointFUN,Ninterval,Nraz,NpolAR
      real(8)::EPSfun,EPSr
      real(8),dimension(:)::R,Rfun,RfunAro
      real(8),dimension(:,:)::Xlim,ALFA,ARCoffPolinom
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer::IZXDC,IZBFG,IUTE,Npoint1,Npoint2,IMAXITER,ierr,NumbreZero
      integer::INDEX,NmaxSSD,NminSSD,NsredSSD,NtroisSSD,INTERIndex,IDDFG
      integer::INDEXXZ,IXDF,NpointSSD,INpoint1,INpoint2,INpoint3
      real(8)::RXDFG,EPSD,EPSQ,A,B,C,Xz1,Xz2,Xz3,RHAS,SUMDFG,F1,F2
	real(8)::A01,A02,ALFA12,RIntPLim1,RIntPLim2,Fsred,Xsred
 	integer,allocatable,dimension(:)::Nextremum,NpointZero,Nyzlov
      integer,allocatable,dimension(:)::NrazInt
      real(8),allocatable,dimension(:)::Xzero,RcoffPol,Ycoff,Rresh
      real(8),allocatable,dimension(:)::RSSD,ALFASSD,A0SSD,RRER,RRYRALFA
      real(8),allocatable,dimension(:)::RRYRA0,ACoffPolinom
      real(8),allocatable,dimension(:,:)::Amatrix,RlimTTR,A0 

	!щрюо 0 сярюмюбкхбюел рнвйс я йнрнпни тсмйжхъ гюмскъеряъ
      DO IZXDC=1,NpointFUN-1 
	  ! write(6,*) R(IZXDC),Rfun(IZXDC)
         IF(Rfun(IZXDC).EQ.0.D0.AND.Rfun(IZXDC+1).EQ.0.D0) THEN
            Npoint=IZXDC
		  EXIT 
	   ENDIF
	ENDDO
  
     	!бшдекъел оюлърэ дкъ люяяхбнб
      allocate(Xzero(nn-l),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Xzero" IS NOT SELECTED'
	stop 
	endif
	allocate(Nextremum(nn-l),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Nextremum" IS NOT SELECTED'
	stop 
	endif
      allocate(NpointZero(nn-l),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "NpointZero" IS NOT SELECTED'
	stop 
	endif
	allocate(ACoffPolinom(nn-l),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "ACoffPolinom" IS NOT SELECTED'
	stop 
	endif

      allocate(RcoffPol(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RcoffPol" IS NOT SELECTED'
	stop 
	endif
	allocate(Nyzlov(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Nyzlov" IS NOT SELECTED'
	stop 
	endif
	allocate(Amatrix(nn-l-1,nn-l-1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Amatrix" IS NOT SELECTED'
	stop 
	endif 
	allocate(Ycoff(nn-l-1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Ycoff" IS NOT SELECTED'
	stop 
	endif
	allocate(Rresh(nn-l-1),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Rresh" IS NOT SELECTED'
	stop 
	endif
      ! +1000-днонкмхрекэмне вхякн хмрепбюкнб(бнглнфмне вхякн хмрепбюкнб мю йнрнпше пюгахбюеряъ онякедмхи хмрепбюк) 
	allocate(RlimTTR(2,2*Nraz*(nn-l)+1000),stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RlimTTR" IS NOT SELECTED'
	stop 
	endif
      allocate(NrazInt(2*Nraz*(nn-l)+1000),stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "NrazInt" IS NOT SELECTED'
	stop 
	endif
	allocate(A0(3,2*Nraz*(nn-l)+1000),stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "A0" IS NOT SELECTED'
	stop 
	endif
  
      allocate(RSSD(Npoint),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RSSD" IS NOT SELECTED'
	stop 
	endif
	allocate(ALFASSD(Npoint),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "ALFASSD" IS NOT SELECTED'
	stop 
	endif
      allocate(A0SSD(Npoint),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "A0SSD" IS NOT SELECTED'
	stop 
	endif
      allocate(RRER(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RRER" IS NOT SELECTED'
	stop 
	endif
	allocate(RRYRALFA(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RRYRALFA" IS NOT SELECTED'
	stop 
	endif
	allocate(RRYRA0(3),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "RRYRA0" IS NOT SELECTED'
	stop 
	endif

	! гюмскъел оепед пюанрни
      Xzero=0.D0
      Nextremum=0
	NpointZero=0
      Nyzlov=0
      RcoffPol=0.D0
      Amatrix=0.D0
      Ycoff=0.D0
      Rresh=0.D0
	RlimTTR=0.D0
	A0=0.D0
	ALFA=0.D0
      Xlim=0.D0
      NrazInt=0
	RSSD=0.D0
	ALFASSD=0.D0
	A0SSD=0.D0
	RRER=0.D0
	RRYRALFA=0.D0
	RRYRA0=0.D0
	ACoffPolinom=0.D0
      RfunAro=0.D0

	
 
      ! щрюо 1. НОПЕДЕКЪЕЛ МСКХ ТСМЙЖХХ Pln(r)
	! ОПНБЕПЪЕЛ МЮКХВХЕ МСКЕИ С ТСМЙЖХХ
	! оепбши мскэ
	Xzero(1)=0.D0  ! яннрберярбсер гмюмскемхч тсмйжхх б мювюке хмрепбюкю
      NumbreZero=1  
	IF((nn-l-1).NE.0) THEN
        ! нясыеярбкъел онхяй мскеи
	  DO IZXDC=1,Npoint-2 ! ХЯЙКЧВЮЕЛ ХГ ПЮЯЯЛНРПЕМХЪ ЙПЮИМХЧ РНВЙС Npoint,
	                      ! Б ЩРНИ РНВЙХ ТСМЙЖХЪ АКХГЙЮ Й МСКЧ НДМЮЙН ЩРН МЕ ЯБЪГЮМН Я БНГМХЙМНБЕМХЕЛ ОЕПЕЯЕВЕМХЪ Я НЯЭЧ r   
	     ! бшъбкъел мскэ тсмйжхх
	IF(Rfun(IZXDC).GE.0.D0.AND.Rfun(IZXDC+1).LT.0.D0.OR.Rfun(IZXDC).LT.0.D0.AND.Rfun(IZXDC+1).GE.0.D0) THEN 
             ! мюидем мскэ тсмйжхх
	       NumbreZero=NumbreZero+1
             ! гюохяшбюел мнлеп рнвйх йнрнпюъ мюундхряъ бакхгх мскъ
	       NpointZero(NumbreZero)=IZXDC
	       ! мюундхл гмювемхе R опх йнрнпнл тсмйжхъ пюбмю мскч
		   Nyzlov(1)=IZXDC-1
		   Nyzlov(2)=IZXDC
		   Nyzlov(3)=IZXDC+1
	       ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,Rfun,RcoffPol)
	     ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
		   A=RcoffPol(1)
		   B=RcoffPol(2)
		   C=RcoffPol(3) 
             ! ОПХПЮБМХБЮЪ ОНКСВЕММШИ ОНКХМНЛ Й МСКЧ ЛШ РЕЛ ЯЮЛШЛ ОНКСВЮЕЛ ЙНПЕМЭ СПЮБМЕМХЪ (мскэ тсмйжхх)
	       Xz1=(-B+DSQRT(B**2-4.D0*A*C))/(2.D0*A)
		   Xz2=(-B-DSQRT(B**2-4.D0*A*C))/(2.D0*A)
	       ! БШЪЯМЪЕЛ ЙЮЙНИ ХГ ЩРХУ ГМЮВЕМХИ ЯННРБЕРЯРБСЕР ЙНПМЧ
		   INDEX=0
		   IF(Xz1.GT.R(IZXDC).AND.Xz1.LT.R(IZXDC+1)) THEN
		      Xzero(NumbreZero)=Xz1
	          INDEX=1
	       ENDIF 
             IF(Xz2.GT.R(IZXDC).AND.Xz2.LT.R(IZXDC+1)) THEN
		      Xzero(NumbreZero)=Xz2
	          INDEX=1
             ENDIF 
	       ! опнбепъел мюидем кх йнпемэ
	       IF(INDEX.EQ.0) THEN
                WRITE(*,*) 'INDEX=0, THE ROOT IS NOT FOUND'
	          READ(*,*)
	          STOP
             ENDIF
	     ENDIF
        ENDDO
	  ! опнбепъел яннрберярбхе онксвеммшу йнпмеи
	  IF((nn-l-1).NE.(NumbreZero-1)) THEN
          WRITE(*,*) 'DISCREPANCY OF NUMBER OF ROOTS'
		READ(*,*)
		STOP 
        ENDIF    

	ENDIF

      !WRITE(*,*) 'ZERO',NumbreZero
	!DO IZXDC=1,NumbreZero
      !   WRITE(*,*) Xzero(IZXDC) 
	!ENDDO
      !READ(*,*)



	 


	! щрюо 2. мюундхл щйярпелсл тсмйжхх (лхмхлслш х люйяхлслш)
	! бшъямъел вхякн щйярпелслнб
	IF((nn-l).EQ.1) THEN 
        ! ндхм щйярпелсл (люйяхлсл)
        NmaxSSD=1
        DO IZXDC=1,Npoint
         IF(Rfun(IZXDC).GT.Rfun(NmaxSSD)) NmaxSSD=IZXDC	
	  ENDDO
        ! гюохяшбюел мнлеп люйяхлюкэмни рнвйх
	  Nextremum(1)=NmaxSSD
      ELSE
        ! яксвюи меяйнкэйху щйярпелслнб
        ! оепбши хмрепбюк (дн оепбнцн мскъ) (люйяхлсл)
	  NmaxSSD=1
        DO IZXDC=2,NpointZero(2)
         IF(Rfun(IZXDC).GT.Rfun(NmaxSSD)) NmaxSSD=IZXDC	
	  ENDDO
        Nextremum(1)=NmaxSSD
        
	  ! хмрепбюкш опнлефсрнвмше 
	  DO IZXDC=2,NumbreZero-1
	     ! сярюмюбкхбюел лхмхлсл хкх люйяхлсл
	     INDEX=NpointZero(IZXDC)
	     IF((Rfun(INDEX+2)-Rfun(INDEX+1)).GT.0.D0) THEN
	        NmaxSSD=NpointZero(IZXDC)+1
              DO IZBFG=NpointZero(IZXDC)+1,NpointZero(IZXDC+1)
                 IF(Rfun(IZBFG).GT.Rfun(NmaxSSD)) NmaxSSD=IZBFG	
	        ENDDO
              Nextremum(IZXDC)=NmaxSSD
             ELSE
	        NminSSD=NpointZero(IZXDC)+1
              DO IZBFG=NpointZero(IZXDC)+1,NpointZero(IZXDC+1)
                 IF(Rfun(IZBFG).LT.Rfun(NminSSD)) NminSSD=IZBFG	
	        ENDDO
              Nextremum(IZXDC)=NminSSD
           ENDIF 
	  ENDDO
        
	  ! онякедмхи хмрепбюк
	  ! сярюмюбкхбюел лхмхлсл хкх люйяхлсл
	     INDEX=NpointZero(NumbreZero)
	     IF((Rfun(INDEX+2)-Rfun(INDEX+1)).GT.0.D0) THEN
	        NmaxSSD=NpointZero(NumbreZero)+1
              DO IZBFG=NpointZero(NumbreZero)+1,Npoint
                 IF(Rfun(IZBFG).GT.Rfun(NmaxSSD)) NmaxSSD=IZBFG	
	        ENDDO
              Nextremum(NumbreZero)=NmaxSSD
             ELSE
              NminSSD=NpointZero(NumbreZero)+1
              DO IZBFG=NpointZero(NumbreZero)+1,Npoint
                 IF(Rfun(IZBFG).LT.Rfun(NminSSD)) NminSSD=IZBFG	
	        ENDDO
              Nextremum(NumbreZero)=NminSSD
           ENDIF 
	ENDIF



      
      !WRITE(*,*) 'EXTREMUM',(nn-l)
	!DO IZXDC=1,(nn-l)
      !   WRITE(*,*) R(Nextremum(IZXDC)),Nextremum(IZXDC) 
	!ENDDO
      !READ(*,*)


	
 

	! щрюо 3. мюундхл йнщттхжхемрш онкхмнлю тсмйжхх
	! онкхмнл мскхбни яреоемх
	IF((nn-l-1).EQ.0) THEN
	  ACoffPolinom(1)=1.D0
	ENDIF
	! онкхмнл оепбни яреоемх
	IF((nn-l-1).EQ.1) THEN 
        ACoffPolinom(1)=1.D0
        ACoffPolinom(2)=-1.D0/Xzero(2)
      ENDIF
	! онкхмнл брнпни х анкэьеи яреоемх
      IF((nn-l-1).GE.2) THEN 
      
	DO IZXDC=1,nn-l-1
	   DO IZBFG=1,nn-l-1
	      Amatrix(IZXDC,IZBFG)=Xzero(IZXDC+1)**IZBFG
	   ENDDO  
	   Ycoff(IZXDC)=-1.D0  
      ENDDO
	! пеьюел яхярелс кхмеимшу спюбмемхи 
      call EFSH_SYSTEM_LINEAR_EQUATIONS(nn-l-1,Amatrix,Ycoff,Rresh)
	
	! гюохяшбюел йнщттхжхемрш онкхмнлю
	ACoffPolinom(1)=1.D0
	DO IZXDC=2,nn-l
	   ACoffPolinom(IZXDC)=Rresh(IZXDC-1)
	ENDDO 

      ENDIF 

      !WRITE(*,*) 'POLINIM',(nn-l)
	!DO IZXDC=1,(nn-l)
      !   WRITE(*,*) ACoffPolinom(IZXDC)
	!ENDDO
      !READ(*,*)




	
 
	! щрюо 4. тнплхпсел мюанп хмрепбюкнб юоопнйяхлюжхх
      ! тнплхпсел люяяхб хяундмшу хмрепбюкнб
      
	! тнплхпсел онякедмхи хмрепбюк(цкюбмше цпюмхжш хмрепбюкю)
      ! мхфмхи опедек
	RlimTTR(1,2*(nn-l))=R(Nextremum(nn-l))
	! бепумхи опедек
	RlimTTR(2,2*(nn-l))=R(Npoint) 
   	
	! тнплхпсел мевермше хмрепбюкш
	INDEX=0
	DO IZBFG=1,2*(nn-l)-1,2
	   ! явервхй мнлепнб цпюмхж
	   INDEX=INDEX+1
         ! мхфмъъ цпюмхжю
	   RlimTTR(1,IZBFG)=Xzero(INDEX)
         ! бепумъъ цпюмхжю
	   RlimTTR(2,IZBFG)=R(Nextremum(INDEX))
      ENDDO
      
	! тнплхпсел вермше хрепбюкш
	DO IZBFG=2,2*(nn-l)-1,2
	   ! мхфмъъ цпюмхжю
	   RlimTTR(1,IZBFG)=RlimTTR(2,IZBFG-1)
         ! бепумъъ цпюмхжю
	   RlimTTR(2,IZBFG)=RlimTTR(1,IZBFG+1)
      ENDDO
	
	! тнплхпсел онякедмхи хмрепбюк
	! мюундхл яйнкэйн пюг опедонякедмхи хмрепбюк сйкюдшбюеряъ б онякедмел
	Xz1=RlimTTR(2,2*(nn-l))-RlimTTR(1,2*(nn-l)) 
      Xz2=RlimTTR(2,2*(nn-l)-1)-RlimTTR(1,2*(nn-l)-1)
	! вхякн (опедонякедмху) хмрепбюкнб
	INTERIndex=IDINT(Xz1/Xz2)
	! нопедекъел вхякн хмрепбюкнб мю йнрнпше асдер пюгахр онякедмхи хмрепбюк
      INDEX=0
	IMAXITER=0
	DO IZBFG=1,INTERIndex
	   IF(IMAXITER+IZBFG.GT.INTERIndex) THEN
	      EXIT 
         ENDIF
	   INDEX=INDEX+1
	   IMAXITER=IMAXITER+IZBFG
      ENDDO

	! пюгахбюел онякедмхи хмрепбюк
      ! мювюкэмши хмрепбюк
      RlimTTR(1,2*(nn-l))=R(Nextremum(nn-l))
      RlimTTR(2,2*(nn-l))=RlimTTR(1,2*(nn-l))+Xz2 
      ! опнлефсрнвмше хмрепбюкш 
	DO IZBFG=2,INDEX-1
      RlimTTR(1,2*(nn-l)+IZBFG-1)=RlimTTR(2,2*(nn-l)+IZBFG-2)
      RlimTTR(2,2*(nn-l)+IZBFG-1)=RlimTTR(1,2*(nn-l))+FLOAT(IZBFG)*Xz2
	ENDDO
      ! онякедмхи хмрепбюк
      RlimTTR(1,2*(nn-l)+INDEX-1)=RlimTTR(2,2*(nn-l)+INDEX-2)
      RlimTTR(2,2*(nn-l)+INDEX-1)=R(Npoint) 



	! мюундхл люйяхлюкэмне вхякн Nraz дкъ йюфднцн хмрепбюкю мю йнрнпне
	! нм лнфер ашрэ пюгдекем
	! б хмрепбюке днкфмн яндепфюрэяъ ме лемэье рпеу рнвей
	IF(Nraz.NE.1) THEN 
        ! жхйк он хмрепбюкюл
	  DO IZBFG=1,2*(nn-l)+INDEX-1
	     ! жхйк он вхякс мю йнрнпне пюгахбюеряъ хмрепбюк
           DO IZXDC=Nraz,1,-1
              RHAS=(RlimTTR(2,IZBFG)-RlimTTR(1,IZBFG))/float(IZXDC) 
              ! пюяялюрпхбюел бнглнфмше хмрепбюкш
	        INDEXXZ=0 ! хмдейя сйюгшбючыхи мю рхо бшундю хг жхйкю
              DO INTERIndex=1,IZXDC
                 Xz1=RlimTTR(1,IZBFG)+RHAS*float(INTERIndex-1)
                 RIntPLim1=Xz1
                 RIntPLim2=Xz1+RHAS
	           IDDFG=0
	           DO IUTE=1,Npoint
                    ! нопедекъел лхмхлюкэмсч рнвйс
                    IF(R(IUTE).GE.RIntPLim1.AND.IDDFG.EQ.0) THEN
                    NminSSD=IUTE
	              IDDFG=1
		          ENDIF
	              ! нопедекъел люйяхлюкэмсч рнвйс
                    IF(R(IUTE).LE.RIntPLim2) THEN
                      NmaxSSD=IUTE
	              ENDIF
	           ENDDO
	           ! опнбепъел бшонкмемхе сякнбхъ (вхякн рнвей б хмрепбюке днкфмн ашрэ анкэье кхан пюбмн рпеу)
			   IF((NmaxSSD-NminSSD+1).LT.3) THEN
                   INDEXXZ=1
                   EXIT  
                 ENDIF
	        ENDDO
     	        ! опнбепъел сякнбхе бшундю хг жхйкю
              IF(INDEXXZ.EQ.0) THEN
                 NrazInt(IZBFG)=IZXDC
	           EXIT
	          ELSE
                 NrazInt(IZBFG)=IZXDC
              ENDIF
  		 ENDDO
        ENDDO
      ENDIF


      !WRITE(*,*) 'Nraz'
	!DO IZXDC=1,2*(nn-l)+INDEX-1
	!   WRITE(*,*) NrazInt(IZXDC)
	!ENDDO
      !READ(*,*)

      



      ! тнплхпсел йнмевмсч яерйс хмрепбюкнб 
      ! опнбепъел бн яйнкэйн пюг мсфмн сбекхвхрэ вхякн хмрепбюкнб 
      IF(Nraz.EQ.1) THEN
	  ! вхякн хмрепбюкнб 
        Ninterval=2*(nn-l)+INDEX-1
        DO IZBFG=1,2*(nn-l)+INDEX-1
	     Xlim(1,IZBFG)=RlimTTR(1,IZBFG)
           Xlim(2,IZBFG)=RlimTTR(2,IZBFG)
	  ENDDO
	ELSE
      
	Ninterval=0
	DO IZXDC=1,2*(nn-l)+INDEX-1
	   ! жхйк он сбекхвемхч вхякю хмрепбюкнб
         RHAS=(RlimTTR(2,IZXDC)-RlimTTR(1,IZXDC))/float(NrazInt(IZXDC))  
	   DO IZBFG=1,NrazInt(IZXDC)
            Ninterval=Ninterval+1
            Xz1=RlimTTR(1,IZXDC)+RHAS*float(IZBFG-1)
            Xlim(1,Ninterval)=Xz1
            Xlim(2,Ninterval)=Xz1+RHAS
	   ENDDO
      ENDDO

	ENDIF

      
      !WRITE(*,*) 'INTERVALS',Ninterval
	!DO IZXDC=1,Ninterval
      !   WRITE(*,*) Xlim(1,IZXDC),Xlim(2,IZXDC)
	!ENDDO
      !READ(*,*)
      !INDEXXZ,IXDF,NpointSSD

      
	! онксвюел гмювемхъ тсмйжхх ALFA and A0 б хгбеярмшу рнвйюу
      A0SSD=0.D0
      ALFASSD=0.D0
	NpointSSD=0
      ! ОНЯРПНЕМХЕ ТСМЙЖХХ ALFA,A0 НР R   RSSD,ALFASSD,A0SSD
	DO IZXDC=2,Npoint-1,2  ! йпюимхх рнвйх ме пюяялюрпхбюел 
	                       ! (хяйкчвюел пъд онякедму рнвей опхбндъыху й яхкэмшл хяйюфемъл)  
         NpointSSD=NpointSSD+1
	   ! опнбндхл юоопнйяхлюжхч 
	   ! гюохяшбюел гмювемхе пюдхсяю
	   RSSD(NpointSSD)=R(IZXDC)
	   ! ЩЙЯРПХЛЮКЭМЮЪ РНВЙЮ
		  SUMDFG=0.D0
		  DO IZBFG=1,nn-l
		     SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC+1)**(IZBFG-1)
	      ENDDO
		  F2=SUMDFG*R(IZXDC+1)**(l+1)
	      ! ЯПЕДМЪЪ РНВЙЮ
		  SUMDFG=0.D0
		  DO IZBFG=1,nn-l
		     SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC)**(IZBFG-1)
	      ENDDO
		  F1=SUMDFG*R(IZXDC)**(l+1)
	      ! НОПЕДЕКЪЕЛ ЮКЭТЮ Б ДЮММНЛ ЯКСВЮЕ
	      RXDFG=DLOG((Rfun(IZXDC)*F2)/(Rfun(IZXDC+1)*F1))
            ALFASSD(NpointSSD)=RXDFG/(R(IZXDC+1)-R(IZXDC))
	      ! нопедекъел A01 б рнвйе IZXDC
         	  SUMDFG=0.D0
		  DO IZBFG=1,nn-l
		     SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC)**(IZBFG-1)
	      ENDDO
		  F2=SUMDFG*R(IZXDC)**(l+1)
          A0SSD(NpointSSD)=Rfun(IZXDC)*DEXP(ALFASSD(NpointSSD)*R(IZXDC))/F2
                 
	   !WRITE(6,*) RSSD(NpointSSD),A0SSD(NpointSSD),ALFASSD(NpointSSD)
	ENDDO

	! мюундхл рнвйс янрберярбсчысч йнмжс (сярпюмъел бнгпюярюмхе мю йнмже хмрепбюкю)
      DO IZXDC=NpointSSD,2,-1
         IF(DABS(A0SSD(IZXDC)).LT.DABS(A0SSD(IZXDC-1))) THEN
            ! мюидемю йнмевмюъ рнвйю
		  NpointSSD=IZXDC
		  EXIT  
	   ENDIF
      ENDDO
      
      
	
	
	! щрюо 5 мюундхл йнщттхжхемрш гюбхяъыхе нр хмрепбюкю юоопнйяхлюжхх (ОНКХМНЛЮЛХ БРНПНЦН ОНПЪДЙЮ)
	! жхйк он хмрепбюкюл ALFA and A0 
	DO IZXDC=1,Ninterval
	   ! нопедекъел мнлепю рнвйх кефюыхх б дюммнл хмрепбюке 
         INDEX=0
 	   NsredSSD=1
	   Xz1=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))/4.D0
         Xz2=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))*0.5D0
         Xz3=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))*3.D0/4.D0
	   ! мюундхл рнвйх акхфюиьхх й дюммшл 
	   INpoint1=1
	   INpoint2=1
	   INpoint3=1
	   DO IZBFG=1,NpointSSD
            ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz1
            IF(DABS(RSSD(IZBFG)-Xz1).LE.DABS(RSSD(INpoint1)-Xz1)) THEN
              INpoint1=IZBFG
	      ENDIF
	      ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz2
            IF(DABS(RSSD(IZBFG)-Xz2).LE.DABS(RSSD(INpoint2)-Xz2)) THEN
              INpoint2=IZBFG
	      ENDIF
	      ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz3
            IF(DABS(RSSD(IZBFG)-Xz3).LE.DABS(RSSD(INpoint3)-Xz3)) THEN
              INpoint3=IZBFG
	      ENDIF
	   ENDDO

         ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint1.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint1-1
	      Nyzlov(2)=INpoint1
	      Nyzlov(3)=INpoint1+1
           ELSE
            Nyzlov(1)=INpoint1-2
	      Nyzlov(2)=INpoint1-1
	      Nyzlov(3)=INpoint1
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,ALFASSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz1 
	   RRER(1)=Xz1 
	   RRYRALFA(1)=A*Xz1**2+B*Xz1+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint1.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint1-1
	      Nyzlov(2)=INpoint1
	      Nyzlov(3)=INpoint1+1
           ELSE
            Nyzlov(1)=INpoint1-2
	      Nyzlov(2)=INpoint1-1
	      Nyzlov(3)=INpoint1
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,A0SSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz1
	   RRYRA0(1)=A*Xz1**2+B*Xz1+C
         
	  ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
        IF(INpoint2.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint2-1
	      Nyzlov(2)=INpoint2
	      Nyzlov(3)=INpoint2+1
           ELSE
            Nyzlov(1)=INpoint2-2
	      Nyzlov(2)=INpoint2-1
	      Nyzlov(3)=INpoint2
    	  ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,ALFASSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz2 
	   RRER(2)=Xz2 
	   RRYRALFA(2)=A*Xz2**2+B*Xz2+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint2.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint2-1
	      Nyzlov(2)=INpoint2
	      Nyzlov(3)=INpoint2+1
           ELSE
            Nyzlov(1)=INpoint2-2
	      Nyzlov(2)=INpoint2-1
	      Nyzlov(3)=INpoint2
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,A0SSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz2
	   RRYRA0(2)=A*Xz2**2+B*Xz2+C

        ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint3.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint3-1
	      Nyzlov(2)=INpoint3
	      Nyzlov(3)=INpoint3+1
           ELSE
            Nyzlov(1)=INpoint3-2
	      Nyzlov(2)=INpoint3-1
	      Nyzlov(3)=INpoint3
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,ALFASSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz3 
	   RRER(3)=Xz3 
	   RRYRALFA(3)=A*Xz3**2+B*Xz3+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint3.NE.NpointSSD) THEN
	      Nyzlov(1)=INpoint3-1
	      Nyzlov(2)=INpoint3
	      Nyzlov(3)=INpoint3+1
           ELSE
            Nyzlov(1)=INpoint3-2
	      Nyzlov(2)=INpoint3-1
	      Nyzlov(3)=INpoint3
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,A0SSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz3
	   RRYRA0(3)=A*Xz3**2+B*Xz3+C
         ! нопедекъел йнщттхжхемрш юопнйяхлюжхх ALFA  мю дюммнл хмрепбюке
	   Nyzlov(1)=1
	   Nyzlov(2)=2
	   Nyzlov(3)=3
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RRER,RRYRALFA,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   ALFA(3,IZXDC)=RcoffPol(1)
	   ALFA(2,IZXDC)=RcoffPol(2)
	   ALFA(1,IZXDC)=RcoffPol(3) 
         ! опнбепъел йнщттхжхемр опх люйяхлюкэмни яреоемх
         ! нрпхжюрекэмши йнщттхжхемр меопхелкел
	   IF(ALFA(3,IZXDC).LT.0.D0) THEN
	     ! нопедекъел йнщттхжхемрш юопнйяхлюжхх ALFA  мю дюммнл хмрепбюке
           Nyzlov(1)=1
	     Nyzlov(2)=3
	! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,1,Nyzlov,RRER,RRYRALFA,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	     ALFA(3,IZXDC)=0.D0
	     ALFA(2,IZXDC)=RcoffPol(1)
           ALFA(1,IZXDC)=RcoffPol(2) 
   	   ENDIF 

         ! нопедекъел йнщттхжхемрш юопнйяхлюжхх A0  мю дюммнл хмрепбюке
	   Nyzlov(1)=1
	   Nyzlov(2)=2
	   Nyzlov(3)=3
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RRER,RRYRA0,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A0(3,IZXDC)=RcoffPol(1)
	   A0(2,IZXDC)=RcoffPol(2)
	   A0(1,IZXDC)=RcoffPol(3)  
      ENDDO

      !WRITE(6,*) 'COFF',Ninterval 
      !DO IZXDC=1,Ninterval 
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC)
	!WRITE(6,*)  ALFA(1,IZXDC),ALFA(2,IZXDC),ALFA(3,IZXDC)
      !WRITE(6,*)  A0(3,IZXDC),A0(2,IZXDC),A0(1,IZXDC)
      !ENDDO




	! щрюо 6. онксвюел онкхмнл мю йюфднл хмрепбюке йюй опнхгбедемхе онкхмнлю A0 х ACoffPolinom  
      ARCoffPolinom=0.D0
	! жхйк он хмрепбюкюл
	DO IZXDC=1,Ninterval 
	   ! жхйк он йнщттхжхемрюл онкхмнлю A0
	   DO IZBFG=1,3
            !жхйк он йнщттхжхемрюл онкхмнлю ACoffPolinom 
		  DO INDEX=1,nn-l
      ARCoffPolinom(IZBFG+INDEX-1,IZXDC)=ARCoffPolinom(IZBFG+INDEX-1,IZXDC)+A0(IZBFG,IZXDC)*ACoffPolinom(INDEX)
	      ENDDO
         ENDDO
      ENDDO
      

      !WRITE(6,*) 'APOLINOM',(ACoffPolinom(INDEX),INDEX=1,nn-l)
	!DO IZXDC=1,Ninterval 
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC)
	!WRITE(6,*)  A0(1,IZXDC),A0(2,IZXDC),A0(3,IZXDC)
      !ENDDO
	!WRITE(6,*) 'REZ',2+nn-l
	!DO IZXDC=1,Ninterval 
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC)
	!WRITE(6,*)  (ARCoffPolinom(IZBFG,IZXDC),IZBFG=1,2+nn-l)
      !ENDDO







      EPSfun=0.D0
	EPSr=0.D0
      
	!WRITE(6,*)
	!WRITE(6,*)
      ! бшдюел пегскэрюр пюяверю 
	! нопедекъел онцпеьмнярэ юоопнйяхлюжхх (мюундхл люйяхлюкэмне нрйкнмемхе нр хяундмшу дюммшу)
      DO IZBFG=1,Npoint 
         DO IZXDC=1,Ninterval 
        IF(R(IZBFG).GE.Xlim(1,IZXDC).AND.R(IZBFG).LE.Xlim(2,IZXDC)) THEN
	      INpoint1=IZXDC
	  ENDIF
	   ENDDO 
         SUMDFG=0.D0
         DO IZXDC=1,2+nn-l
         SUMDFG=SUMDFG+ARCoffPolinom(IZXDC,INpoint1)*R(IZBFG)**(IZXDC-1)
	   ENDDO
         SUMDFG=SUMDFG*R(IZBFG)**(l+1)
	   ALFA12=0.D0
	   ALFA12=ALFA12-ALFA(3,INpoint1)*R(IZBFG)**3
         ALFA12=ALFA12-ALFA(2,INpoint1)*R(IZBFG)**2
         ALFA12=ALFA12-ALFA(1,INpoint1)*R(IZBFG)
	   F1=DEXP(ALFA12)*SUMDFG
	   RfunAro(IZBFG)=F1
	 !  WRITE(6,*) R(IZBFG),F1,Rfun(IZBFG)
	   ! БШЪЯМЪЕЛ МЕ ПЮБЕМЮ КХ ТСМЙЖХЪ МСКЧ Б ЩРНИ РНВЙХ 
         IF(Rfun(IZBFG).NE.0.D0) THEN
	      IF(DABS((F1-Rfun(IZBFG))/Rfun(IZBFG)).GE.EPSfun) THEN
              EPSfun=DABS((F1-Rfun(IZBFG))/Rfun(IZBFG))
			EPSr=R(IZBFG) 
	      ENDIF
         ENDIF 
      ENDDO

      ! ЯРЕОЕМЭ МНБНЦН ОНКХМНЛЮ ARCoffPolinom
      NpolAR=2+nn-l


      

     


      !WRITE(6,*) 'COFF',Ninterval
	!DO IZXDC=1,Ninterval
      !   WRITE(6,*) A0(IZXDC),ALFA(IZXDC)
	!   WRITE(6,*) Xlim(1,IZXDC),Xlim(2,IZXDC)
         !READ(*,*)
	!ENDDO
      !READ(*,*)
     

	
    

   

      ! сдюкемхе люяяхбнб хг оълърх 
	deallocate(ACoffPolinom,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "ACoffPolinom" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(A0,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "A0" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRER,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRER" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRYRALFA,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRYRALFA" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRYRA0,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRYRA0" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
      deallocate(RSSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RSSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(ALFASSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "ALFASSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(A0SSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "A0SSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(NrazInt,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "NrazInt" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RlimTTR,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RlimTTR" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(Rresh,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Rresh" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
      deallocate(Xzero,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Xzero" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(NpointZero,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "NpointZero" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(Nextremum,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Nextremum" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Nyzlov,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Nyzlov" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(RcoffPol,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RcoffPol" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(Amatrix,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Amatrix" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Ycoff,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Ycoff" IS NOT REMOVED FROM MEMORY'
	stop 
	endif


      return
      end subroutine EFSH_APPROXIMATION_RADIAL_FUNCTION 





  ! SUB-PROGRAM FOR APPROXIMATION OF THE RADIAL PART OF THE WAVE FUNCTION (DISCRETE WAVE FUNCTION)
  ! DESCRIPTION OF SUBPROGRAM PARAMETERS
  ! l-orbital quantum number
  ! NpointFUN-number of points
  ! R (N) -array of radius values
  ! Rfun (N) -array of function values
  ! Nraz-parameter increasing by Nraz the initial number of intervals 2 (nn-l)
  ! Ninterval-number of intervals of approximation
  ! NpolAR-degree of a polynomial with ARCoffPolinom coefficients
  ! Xlim (2, Ninterval) -array of limits of approximation intervals
  ! Xlim (1, Ninterval) -lower limit
  ! Xlim (2, Ninterval) -upper limit
  ! A0 (3, Ninterval) -array of coefficients of the approximating function
  ! (A03 * x ^ 2 + A02 * x + A01) * x ^ (l + 1) * exp (- (ALFA3 * x ^ 2 + ALFA2 * x + ALFA3) * x) * (1 + sum (ak * x ^ k))
  ! ALFA (3, Ninterval) - the exponential factor of the approximating function depends on the interval
  ! ALFA (3, Ninterval) - CORRESPONDING TO DEGREE r ^ 2
  ! ALFA (2, Ninterval) - CORRESPONDING TO DEGREE r
  ! ALFA (1, Ninterval) - CORRESPONDING TO DEGREE r ^ 0
  ! ARCoffPolinom (NpolAR + 1, Ninterval) -polomial factors RESULT OF PRODUCING POLYNOMS (A03 * x ^ 2 + A02 * x + A01) * (1 + sum (ak * x ^ k))
  ! RfunAro (N) -array of values ??of the approximated function (for comparison with the original)
	subroutine EFSH_APPROXIMATION_RADIAL_FUNCTION_ALFA(l,NpointFUN,R,Rfun,Nraz,Ninterval,NpolAR,Xlim,ALFA,ARCoffPolinom,RfunAro) 
     use dfimsl
	 implicit none
       
     integer::l,Npoint,NpointFUN,Ninterval,Nraz,NpolAR
     real(8),dimension(:)::R,Rfun,RfunAro
     real(8),dimension(:,:)::Xlim,ALFA,ARCoffPolinom
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IZXDC,IZBFG,IUTE,Npoint1,Npoint2,IMAXITER,ierr,NumbreZero,NminPoint
	 integer::NumbreExtrem,IMASIVFD,IIDDNpoin 
     integer::INDEX,NmaxSSD,NminSSD,NsredSSD,NtroisSSD,INTERIndex,IDDFG
     integer::INDEXXZ,IXDF,NpointSSD,INpoint1,INpoint2,INpoint3
     real(8)::RXDFG,EPSD,EPSQ,A,B,C,Xz1,Xz2,Xz3,RHAS,SUMDFG,F1,F2
	 real(8)::A01,A02,ALFA12,RIntPLim1,RIntPLim2,Fsred,Xsred,DELTA,Rk,Rb
     integer,allocatable,dimension(:)::Nyzlov
     integer,allocatable,dimension(:)::NrazInt
     real(8),allocatable,dimension(:)::Xzero,Xextremum,Xzeex,RcoffPol
  	 real(8),allocatable,dimension(:)::Ycoff
     real(8),allocatable,dimension(:)::RSSD,ALFASSD,A0SSD,RRER,RRYRALFA
     real(8),allocatable,dimension(:)::RRYRA0,ACoffPolinom,Rresh
 	 real(8),allocatable,dimension(:)::ALFArezF,A0rezF
     real(8),allocatable,dimension(:,:)::Amatrix,RlimTTR,A0 

     !DO IZXDC=1,NpointFUN 
     !   WRITE(6,*) R(IZXDC),Rfun(IZXDC)
     !ENDDO
      
  	 !щрюо 0 сярюмюбкхбюел рнвйс я йнрнпни тсмйжхъ гюмскъеряъ
	 Npoint=NpointFUN
	 DO IZXDC=NpointFUN,1,-1
        IF(Rfun(IZXDC).EQ.0.D0.AND.Rfun(IZXDC-1).NE.0.D0) THEN
           Npoint=IZXDC-1
	 	   EXIT 
	    ENDIF
	 ENDDO
     
	 !  сярюмюбкхбюел оепбсч рнвйс нркхвмсч нр мскъ я мювюкю хрепбюкю 
	 IF(Rfun(1).NE.0.D0) THEN
         NminPoint=1
	    ELSE
         DO IZXDC=1,NpointFUN
            IF(Rfun(IZXDC).EQ.0.D0.AND.Rfun(IZXDC+1).NE.0.D0) THEN
               NminPoint=IZXDC+1
		    ENDIF
         ENDDO
     ENDIF 
     
     

     !бшдекъел оюлърэ дкъ люяяхбнб
	 ! 1000-МСКЕИ ТСМЙЖХХ (люйяхлсл)
	 ! 1000-ЩЙЯРПЕЛСЛНБ   (люйяхлсл)
     allocate(Xzero(1000),stat=ierr)
	 if(ierr/=0) then
       write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	   write(*,*) 'MEMORY ON THE FILE "Xzero" IS NOT SELECTED'
	   stop 
	 endif
	 allocate(Xextremum(1000),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "Xextremum" IS NOT SELECTED'
	    stop 
	 endif
      

     allocate(RcoffPol(3),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "RcoffPol" IS NOT SELECTED'
	    stop 
	 endif
	 allocate(Nyzlov(3),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "Nyzlov" IS NOT SELECTED'
	    stop 
	 endif
     allocate(RSSD(Npoint),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "RSSD" IS NOT SELECTED'
	    stop 
	 endif
	 allocate(ALFASSD(Npoint),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "ALFASSD" IS NOT SELECTED'
	    stop 
	 endif
     allocate(A0SSD(Npoint),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "A0SSD" IS NOT SELECTED'
	    stop 
	 endif
     allocate(ALFArezF(Npoint),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "ALFArezF" IS NOT SELECTED'
	    stop 
	 endif
     allocate(A0rezF(Npoint),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "A0rezF" IS NOT SELECTED'
	    stop 
	 endif
     allocate(RRER(3),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "RRER" IS NOT SELECTED'
	    stop 
	 endif
	 allocate(RRYRALFA(3),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "RRYRALFA" IS NOT SELECTED'
	    stop 
	 endif
	 allocate(RRYRA0(3),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "RRYRA0" IS NOT SELECTED'
	    stop 
	 endif
      
	 ! гюмскъел оепед пюанрни
     Xzero=0.D0
     Xextremum=0
     Nyzlov=0
     RcoffPol=0.D0
     ALFA=0.D0
     Xlim=0.D0
     RSSD=0.D0
	 ALFASSD=0.D0
	 A0SSD=0.D0
	 ALFArezF=0.D0
	 A0rezF=0.D0
	 RRER=0.D0
	 RRYRALFA=0.D0
	 RRYRA0=0.D0
     RfunAro=0.D0

     
     
 
     ! щрюо 1. опнбндхл юмюкхг тсмйжхх сярюмюбкхбюел мскх х рнвйх щйярпелслю

	 call EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA(Npoint,R,RFUN,NumbreZero,Xzero,NumbreExtrem,Xextremum) 
       
     WRITE(6,*) 'The analysis of a radial part of function'
	 WRITE(6,*) 'NumbreZero= ',NumbreZero,' NumbreExtrem= ',NumbreExtrem
	 DO IZXDC=1,NumbreZero
        WRITE(6,*) 'ZERO',Xzero(IZXDC)
	 ENDDO
	 DO IZXDC=1,NumbreExtrem
        WRITE(6,*) 'EXTREMUM',Xextremum(IZXDC)
	 ENDDO

    

     !  бшдекъел оюлърэ онд люяяхбш
	 allocate(Xzeex(NumbreZero+NumbreExtrem+1),stat=ierr)
	 if(ierr/=0) then
        write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	    write(*,*) 'MEMORY ON THE FILE "Xzeex" IS NOT SELECTED'
	    stop 
	 endif
	 ! +1000-днонкмхрекэмне вхякн хмрепбюкнб(бнглнфмне вхякн хмрепбюкнб мю йнрнпше пюгахбюеряъ онякедмхи хмрепбюк) 
	 allocate(RlimTTR(2,Nraz*(NumbreZero+NumbreExtrem)+1000),stat=ierr) 
	 if(ierr/=0) then
       write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	   write(*,*) 'MEMORY ON THE FILE "RlimTTR" IS NOT SELECTED'
	   stop 
	 endif
     allocate(NrazInt(Nraz*(NumbreZero+NumbreExtrem)+1000),stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "NrazInt" IS NOT SELECTED'
	stop 
	endif
	allocate(A0(3,Nraz*(NumbreZero+NumbreExtrem)+1000),stat=ierr) 
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "A0" IS NOT SELECTED'
	stop 
	endif
	allocate(ACoffPolinom(NumbreZero),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "ACoffPolinom" IS NOT SELECTED'
	stop 
	endif
	allocate(Amatrix(NumbreZero,NumbreZero),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Amatrix" IS NOT SELECTED'
	stop 
	endif 
	allocate(Ycoff(NumbreZero),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Ycoff" IS NOT SELECTED'
	stop 
	endif
	allocate(Rresh(NumbreZero),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
	write(*,*) 'MEMORY ON THE FILE "Rresh" IS NOT SELECTED'
	stop 
	endif

      ! гюмскъел оепед пюявернл
	Xzeex=0.D0
	RlimTTR=0.D0
      NrazInt=0
	A0=0.D0
   	ACoffPolinom=0.D0
      Amatrix=0.D0
      Ycoff=0.D0
      Rresh=0.D0



    
    

	! щрюо 2. мюундхл йнщттхжхемрш онкхмнлю тсмйжхх (онкхмнл NumbreZero-1 яреоемх)
	! онкхмнл мскхбни яреоемх
	IF((NumbreZero-1).EQ.0) THEN
	  ACoffPolinom(1)=1.D0
	ENDIF
	! онкхмнл оепбни яреоемх
	IF((NumbreZero-1).EQ.1) THEN 
        ACoffPolinom(1)=1.D0
        ACoffPolinom(2)=-1.D0/Xzero(2)
      ENDIF
	! онкхмнл брнпни х анкэьеи яреоемх
      IF((NumbreZero-1).GE.2) THEN 
      
	DO IZXDC=1,NumbreZero-1
	   DO IZBFG=1,NumbreZero-1
	      Amatrix(IZXDC,IZBFG)=Xzero(IZXDC+1)**IZBFG
	   ENDDO  
	   Ycoff(IZXDC)=-1.D0  
      ENDDO
	! пеьюел яхярелс кхмеимшу спюбмемхи 
      call EFSH_SYSTEM_LINEAR_EQUATIONS(NumbreZero-1,Amatrix,Ycoff,Rresh)
	
	! гюохяшбюел йнщттхжхемрш онкхмнлю
	ACoffPolinom(1)=1.D0
	DO IZXDC=2,NumbreZero
	   ACoffPolinom(IZXDC)=Rresh(IZXDC-1)
	ENDDO 

      ENDIF 

      WRITE(6,*) 'POLINIM',NumbreZero
	DO IZXDC=1,NumbreZero
         WRITE(6,*) ACoffPolinom(IZXDC)
	ENDDO
      !READ(*,*)


  
   
	
	! щрюо 3. тнплхпсел мюанп хмрепбюкнб юоопнйяхлюжхх
      
	! тнплхпсел наыхи люяяхб опедекнб
	DO IZBFG=1,NumbreZero
         Xzeex(IZBFG)=Xzero(IZBFG)
      ENDDO
      DO IZBFG=1,NumbreExtrem
         Xzeex(IZBFG+NumbreZero)=Xextremum(IZBFG)
      ENDDO
         Xzeex(NumbreZero+NumbreExtrem+1)=R(Npoint) 

      !WRITE(6,*) 'Xalim',NumbreZero+NumbreExtrem+1
      !write(6,*) (Xzeex(IZBFG),IZBFG=1,NumbreZero+NumbreExtrem+1)


      ! янпрхпсел люяяхб он бнгпюярюмхч
	call dsvrgn(NumbreZero+NumbreExtrem+1,Xzeex,Xzeex)
      ! write(6,*) 'new'
	!write(6,*) (Xzeex(IZBFG),IZBFG=1,NumbreZero+NumbreExtrem+1)
	
	! тнплхпсел люяяхб хяундмшу хмрепбюкнб
      DO IZBFG=1,NumbreZero+NumbreExtrem
	   ! мхфмхи опедек
	   RlimTTR(1,IZBFG)=Xzeex(IZBFG)
	   ! бепумхи опедек
	   RlimTTR(2,IZBFG)=Xzeex(IZBFG+1)
	ENDDO

      !write(6,*) 'pred',NumbreZero+NumbreExtrem
      !DO IZBFG=1,NumbreZero+NumbreExtrem
      !   WRITE(6,*) RlimTTR(1,IZBFG),RlimTTR(2,IZBFG)
      !ENDDO
	
	! тнплхпсел онякедмхи хмрепбюк
	! мюундхл яйнкэйн пюг опедонякедмхи хмрепбюк сйкюдшбюеряъ б онякедмел
	  IMASIVFD=NumbreZero+NumbreExtrem
	  Xz1=RlimTTR(2,IMASIVFD)-RlimTTR(1,IMASIVFD) 
      Xz2=RlimTTR(2,IMASIVFD-1)-RlimTTR(1,IMASIVFD-1)
	! вхякн (опедонякедмху) хмрепбюкнб
	INTERIndex=IDINT(Xz1/Xz2)
	! нопедекъел вхякн хмрепбюкнб мю йнрнпше асдер пюгахр онякедмхи хмрепбюк
      INDEX=0
	IMAXITER=0
	DO IZBFG=1,INTERIndex
	   IF(IMAXITER+IZBFG.GT.INTERIndex) THEN
	      EXIT 
         ENDIF
	   INDEX=INDEX+1
	   IMAXITER=IMAXITER+IZBFG
      ENDDO

	! пюгахбюел онякедмхи хмрепбюк
      ! мювюкэмши хмрепбюк
      RlimTTR(2,IMASIVFD)=RlimTTR(1,IMASIVFD)+Xz2 
      ! опнлефсрнвмше хмрепбюкш 
	DO IZBFG=2,INDEX-1
      RlimTTR(1,IMASIVFD+IZBFG-1)=RlimTTR(2,IMASIVFD+IZBFG-2)
      RlimTTR(2,IMASIVFD+IZBFG-1)=RlimTTR(1,IMASIVFD)+FLOAT(IZBFG)*Xz2
	ENDDO
      ! онякедмхи хмрепбюк
      RlimTTR(1,IMASIVFD+INDEX-1)=RlimTTR(2,IMASIVFD+INDEX-2)
      RlimTTR(2,IMASIVFD+INDEX-1)=R(Npoint) 

      !write(6,*) 'pred',IMASIVFD+INDEX-1
      !DO IZBFG=1,IMASIVFD+INDEX-1
      !   WRITE(6,*) RlimTTR(1,IZBFG),RlimTTR(2,IZBFG)
      !ENDDO




	! мюундхл люйяхлюкэмне вхякн Nraz дкъ йюфднцн хмрепбюкю мю йнрнпне
	! нм лнфер ашрэ пюгдекем
	! б хмрепбюке днкфмн яндепфюрэяъ ме лемэье рпеу рнвей
	IF(Nraz.NE.1) THEN 
        ! жхйк он хмрепбюкюл
	  DO IZBFG=1,NumbreZero+NumbreExtrem+INDEX-1
	     ! жхйк он вхякс мю йнрнпне пюгахбюеряъ хмрепбюк
           DO IZXDC=Nraz,1,-1
              RHAS=(RlimTTR(2,IZBFG)-RlimTTR(1,IZBFG))/float(IZXDC) 
              ! пюяялюрпхбюел бнглнфмше хмрепбюкш
	        INDEXXZ=0 ! хмдейя сйюгшбючыхи мю рхо бшундю хг жхйкю
              DO INTERIndex=1,IZXDC
                 Xz1=RlimTTR(1,IZBFG)+RHAS*float(INTERIndex-1)
                 RIntPLim1=Xz1
                 RIntPLim2=Xz1+RHAS
	           IDDFG=0
	           DO IUTE=1,Npoint
                    ! нопедекъел лхмхлюкэмсч рнвйс
                    IF(R(IUTE).GE.RIntPLim1.AND.IDDFG.EQ.0) THEN
                    NminSSD=IUTE
	              IDDFG=1
		          ENDIF
	              ! нопедекъел люйяхлюкэмсч рнвйс
                    IF(R(IUTE).LE.RIntPLim2) THEN
                      NmaxSSD=IUTE
	              ENDIF
	           ENDDO
	           ! опнбепъел бшонкмемхе сякнбхъ (вхякн рнвей б хмрепбюке днкфмн ашрэ анкэье кхан пюбмн рпеу)
			   IF((NmaxSSD-NminSSD+1).LT.3) THEN
                   INDEXXZ=1
                   EXIT  
                 ENDIF
	        ENDDO
     	        ! опнбепъел сякнбхе бшундю хг жхйкю
              IF(INDEXXZ.EQ.0) THEN
                 NrazInt(IZBFG)=IZXDC
	           EXIT
	          ELSE
                 NrazInt(IZBFG)=IZXDC
              ENDIF
  		 ENDDO
        ENDDO
      ENDIF
       
      !DO IZXDC=1,NumbreZero+NumbreExtrem+INDEX-1
	!   NrazInt(IZXDC)=Nraz
	!ENDDO


      !WRITE(6,*) 'Nraz'
	!DO IZXDC=1,NumbreZero+NumbreExtrem+INDEX-1
	!   WRITE(6,*) NrazInt(IZXDC)
	!ENDDO
      !READ(*,*)

      



      ! тнплхпсел йнмевмсч яерйс хмрепбюкнб 
      ! опнбепъел бн яйнкэйн пюг мсфмн сбекхвхрэ вхякн хмрепбюкнб 
      IF(Nraz.EQ.1) THEN
	  ! вхякн хмрепбюкнб 
        Ninterval=NumbreZero+NumbreExtrem+INDEX-1
        DO IZBFG=1,Ninterval
	     Xlim(1,IZBFG)=RlimTTR(1,IZBFG)
           Xlim(2,IZBFG)=RlimTTR(2,IZBFG)
	  ENDDO
	ELSE
      
	Ninterval=0
	DO IZXDC=1,NumbreZero+NumbreExtrem+INDEX-1
	   ! жхйк он сбекхвемхч вхякю хмрепбюкнб
         RHAS=(RlimTTR(2,IZXDC)-RlimTTR(1,IZXDC))/float(NrazInt(IZXDC))  
	   DO IZBFG=1,NrazInt(IZXDC)
            Ninterval=Ninterval+1
            Xz1=RlimTTR(1,IZXDC)+RHAS*float(IZBFG-1)
            Xlim(1,Ninterval)=Xz1
            Xlim(2,Ninterval)=Xz1+RHAS
	   ENDDO
      ENDDO

	ENDIF

      
      !WRITE(6,*) 'INTERVALS',Ninterval
	!DO IZXDC=1,Ninterval
      !   WRITE(6,*) Xlim(1,IZXDC),Xlim(2,IZXDC)
	!ENDDO
      !READ(*,*)
      !INDEXXZ,IXDF,NpointSSD

      
	! онксвюел гмювемхъ тсмйжхх ALFA and A0 б хгбеярмшу рнвйюу
      A0SSD=0.D0
      ALFASSD=0.D0
	NpointSSD=0
	
      ! ОНЯРПНЕМХЕ ТСМЙЖХХ ALFA,A0 НР R   
	! онярпнемхе тсмйжхх ALFA НР R
	DO IZXDC=NminPoint,Npoint-2,2  ! пюяялюрпхбюел хмрепбюк б йнрнпнл тсмйжхъ нркхвмю нр мскъ  рн еярэ ( йпюимхх рнвйх ме пюяялюрпхбюел )
         NpointSSD=NpointSSD+1
	   ! опнбндхл юоопнйяхлюжхч 
	   ! гюохяшбюел гмювемхе пюдхсяю
	   RSSD(NpointSSD)=R(IZXDC)
	   ! ЩЙЯРПХЛЮКЭМЮЪ РНВЙЮ
		  SUMDFG=0.D0
		  DO IZBFG=1,NumbreZero
		     SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC+1)**(IZBFG-1)
	      ENDDO
		  F2=SUMDFG*R(IZXDC+1)**(l+1)
	      ! ЯПЕДМЪЪ РНВЙЮ
		  SUMDFG=0.D0
		  DO IZBFG=1,NumbreZero
		     SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC)**(IZBFG-1)
	      ENDDO
		  F1=SUMDFG*R(IZXDC)**(l+1)
	      ! НОПЕДЕКЪЕЛ ЮКЭТЮ Б ДЮММНЛ ЯКСВЮЕ
	      !WRITE(6,*) 'LOG',Rfun(IZXDC),F2,Rfun(IZXDC+1),F1
	      !READ(*,*)
		  RXDFG=DLOG((Rfun(IZXDC)*F2)/(Rfun(IZXDC+1)*F1))
          ALFASSD(NpointSSD)=RXDFG/(R(IZXDC+1)-R(IZXDC))
          !WRITE(6,*) 'OLD',NpointSSD,ALFASSD(NpointSSD),R(IZXDC) 
		  !write(6,*) 'old',R(IZXDC),ALFASSD(NpointSSD)
	       
	     
	!   WRITE(6,*) RSSD(NpointSSD),A0SSD(NpointSSD),ALFASSD(NpointSSD)
	ENDDO
    !write(*,*) 'ddd'
	!read(*,*)
    !WRITE(6,*)
    ! бшъбкъел нрйкнмемхъ
    DO IZXDC=2,NpointSSD-1
	   ! опнбепъел нрйкнмемхе
	   DELTA=100.D0*DABS((ALFASSD(IZXDC+1)-ALFASSD(IZXDC))/ALFASSD(IZXDC))
       ! еякх нрйкнмемхе анкке 1000% щрн гмювхр лш хлел цпсасч ньхайс
	   IF(DELTA.GT.1000.D0) THEN
          ! яцкюфхбюел цпсане нрйкнмемхе
          Rk=(ALFASSD(IZXDC)-ALFASSD(IZXDC-1))/(RSSD(IZXDC)-RSSD(IZXDC-1))
		  Rb=ALFASSD(IZXDC)-Rk*RSSD(IZXDC)
		  ALFASSD(IZXDC+1)=Rk*RSSD(IZXDC+1)+Rb 
	   ENDIF
  
    ENDDO
    
	! нясыеярбхл яцкюфхбюмхе онксвеммшу дюммшу
    !DO IZXDC=2,NpointSSD-1
    !   ALFASSD(IZXDC)=(ALFASSD(IZXDC-1)+ALFASSD(IZXDC)+ALFASSD(IZXDC+1))/3.D0
    !   write(6,*) 'NEW',RSSD(IZXDC),ALFASSD(IZXDC) 
	!ENDDO 
    !write(*,*) 'ddd'
	
	! реярхпсел онксвеммше тсмйжхч ALFA мю яксвюимше нрйкнмемхъ
	! йпхрепхел ъбкъеряъ мебнглнфмнярэ пегйни ялемш гмюйю опнхгбндмни 
	! нясыеярбкъел пюявер мю рпеу хмрепбюкю
	!DO IZXDC=2,NpointSSD
	   ! опнбепъел гмювемхе еярэ менцпюмхвеммн пюярсыее
	 !  IF(
	   
	  ! DELTA1=(ALFASSD(IZXDC-2)-ALFASSD(IZXDC-3))/(RSSD(IZXDC-2)-RSSD(IZXDC-3))   
      ! DELTA=(ALFASSD(IZXDC-1)-ALFASSD(IZXDC-2))/(RSSD(IZXDC-1)-RSSD(IZXDC-2))
	   !DELTA2=(ALFASSD(IZXDC-1)-ALFASSD(IZXDC-2))/(RSSD(IZXDC-1)-RSSD(IZXDC-2))   
	   ! опнбепъел врнаш 
    !ENDDO

	! онярпнемхе тсмйжхх A0 НР R   
	NpointSSD=0
	DO IZXDC=NminPoint,Npoint-2,2  ! пюяялюрпхбюел хмрепбюк б йнрнпнл тсмйжхъ нркхвмю нр мскъ  рн еярэ ( йпюимхх рнвйх ме пюяялюрпхбюел )
       NpointSSD=NpointSSD+1
	   ! нопедекъел A01 б рнвйе IZXDC
       SUMDFG=0.D0
	   DO IZBFG=1,NumbreZero
	      SUMDFG=SUMDFG+ACoffPolinom(IZBFG)*R(IZXDC)**(IZBFG-1)
	   ENDDO
	   F2=SUMDFG*R(IZXDC)**(l+1)
	   !WRITE(6,*) 'NEW',NpointSSD,ALFASSD(NpointSSD),R(IZXDC)
	   A0SSD(NpointSSD)=Rfun(IZXDC)*DEXP(ALFASSD(NpointSSD)*R(IZXDC))/F2            
    ENDDO

   
	! мюундхл рнвйс янрберярбсчысч йнмжс (сярпюмъел бнгпюярюмхе мю йнмже хмрепбюкю)
      DO IZXDC=NpointSSD,2,-1
         IF(DABS(A0SSD(IZXDC)).LT.DABS(A0SSD(IZXDC-1))) THEN
            ! мюидемю йнмевмюъ рнвйю
		  NpointSSD=IZXDC
		  EXIT  
	   ENDIF
      ENDDO

      ! ОПНБНДХЛ ЮОПНЙЯХЛЮЖХЧ ТСМЙЖХИ ALFA,A0  ДКЪ СБЕКХВЕМХЪ ЦСЯРНРШ РНВЕЙ
	ALFArezF=0.D0
	A0rezF=0.D0

      DO IZXDC=1,Npoint
         ! мюундхл рнвйс акхфмсчч й дюммни 
         IIDDNpoin=1 
	   DO IZBFG=1,NpointSSD
          F1=dabs(R(IZXDC)-RSSD(IZBFG))
	      F2=dabs(R(IZXDC)-RSSD(IIDDNpoin))
		  IF(F1.LT.F2) THEN 
	        IIDDNpoin=IZBFG
	      ENDIF 
	   ENDDO
	   ! пюяялюрпхбюел пюгмше бюпхюмрш юоопнйяхлюжхх
	   IF(IIDDNpoin.EQ.1) THEN
            ! тсмйжхъ ALFA
		  Nyzlov(1)=IIDDNpoin
	      Nyzlov(2)=IIDDNpoin+1
	      Nyzlov(3)=IIDDNpoin+2
         ENDIF
    	   IF(IIDDNpoin.EQ.NpointSSD) THEN
            ! тсмйжхъ ALFA
		  Nyzlov(1)=IIDDNpoin-2
	      Nyzlov(2)=IIDDNpoin-1
	      Nyzlov(3)=IIDDNpoin
         ENDIF
         IF(IIDDNpoin.NE.1.AND.IIDDNpoin.NE.NpointSSD) THEN
            ! тсмйжхъ ALFA
		  Nyzlov(1)=IIDDNpoin-1
	      Nyzlov(2)=IIDDNpoin
	      Nyzlov(3)=IIDDNpoin+1
         ENDIF
    	  
           ! тсмйжхъ ALFA
	     ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
        call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,ALFASSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3)
	   ALFArezF(IZXDC)=A*R(IZXDC)**2+B*R(IZXDC)+C
           ! тсмйжхъ A0
		 ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
        call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RSSD,A0SSD,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   A0rezF(IZXDC)=A*R(IZXDC)**2+B*R(IZXDC)+C  
       !  WRITE(6,*) R(IZXDC),A0rezF(IZXDC),ALFArezF(IZXDC)
      ENDDO
      
    
	
	! щрюо 4. мюундхл йнщттхжхемрш гюбхяъыхе нр хмрепбюкю юоопнйяхлюжхх (ОНКХМНЛЮЛХ БРНПНЦН ОНПЪДЙЮ)
	! жхйк он хмрепбюкюл ALFA and A0 
	DO IZXDC=1,Ninterval
	   ! нопедекъел мнлепю рнвйх кефюыхх б дюммнл хмрепбюке 
         INDEX=0
 	   NsredSSD=1
	   Xz1=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))/4.D0
         Xz2=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))*0.5D0
         Xz3=Xlim(1,IZXDC)+(Xlim(2,IZXDC)-Xlim(1,IZXDC))*3.D0/4.D0
	   ! мюундхл рнвйх акхфюиьхх й дюммшл 
	   INpoint1=1
	   INpoint2=1
	   INpoint3=1
	   DO IZBFG=1,Npoint
            ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz1
            IF(DABS(R(IZBFG)-Xz1).LE.DABS(R(INpoint1)-Xz1)) THEN
              INpoint1=IZBFG
	      ENDIF
	      ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz2
            IF(DABS(R(IZBFG)-Xz2).LE.DABS(R(INpoint2)-Xz2)) THEN
              INpoint2=IZBFG
	      ENDIF
	      ! нопедекъел рнвйс хмрепбюкю акхфмчч й Xz3
            IF(DABS(R(IZBFG)-Xz3).LE.DABS(R(INpoint3)-Xz3)) THEN
              INpoint3=IZBFG
	      ENDIF
	   ENDDO

         ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint1.NE.Npoint) THEN
	      Nyzlov(1)=INpoint1-1
	      Nyzlov(2)=INpoint1
	      Nyzlov(3)=INpoint1+1
           ELSE
            Nyzlov(1)=INpoint1-2
	      Nyzlov(2)=INpoint1-1
	      Nyzlov(3)=INpoint1
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,ALFArezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz1 
	   RRER(1)=Xz1 
	   RRYRALFA(1)=A*Xz1**2+B*Xz1+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint1.NE.Npoint) THEN
	      Nyzlov(1)=INpoint1-1
	      Nyzlov(2)=INpoint1
	      Nyzlov(3)=INpoint1+1
           ELSE
            Nyzlov(1)=INpoint1-2
	      Nyzlov(2)=INpoint1-1
	      Nyzlov(3)=INpoint1
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,A0rezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz1
	   RRYRA0(1)=A*Xz1**2+B*Xz1+C
         
	  ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
        IF(INpoint2.NE.Npoint) THEN
	      Nyzlov(1)=INpoint2-1
	      Nyzlov(2)=INpoint2
	      Nyzlov(3)=INpoint2+1
           ELSE
            Nyzlov(1)=INpoint2-2
	      Nyzlov(2)=INpoint2-1
	      Nyzlov(3)=INpoint2
    	  ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,ALFArezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz2 
	   RRER(2)=Xz2 
	   RRYRALFA(2)=A*Xz2**2+B*Xz2+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint2.NE.Npoint) THEN
	      Nyzlov(1)=INpoint2-1
	      Nyzlov(2)=INpoint2
	      Nyzlov(3)=INpoint2+1
           ELSE
            Nyzlov(1)=INpoint2-2
	      Nyzlov(2)=INpoint2-1
	      Nyzlov(3)=INpoint2
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,A0rezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz2
	   RRYRA0(2)=A*Xz2**2+B*Xz2+C

        ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint3.NE.Npoint) THEN
	      Nyzlov(1)=INpoint3-1
	      Nyzlov(2)=INpoint3
	      Nyzlov(3)=INpoint3+1
           ELSE
            Nyzlov(1)=INpoint3-2
	      Nyzlov(2)=INpoint3-1
	      Nyzlov(3)=INpoint3
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,ALFArezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz3 
	   RRER(3)=Xz3 
	   RRYRALFA(3)=A*Xz3**2+B*Xz3+C
        
	   ! нопедекъел гмювемхе тсмйжхх б щрху рнвйюу он япедярбюл юоопнйяхлюжхх
         IF(INpoint3.NE.Npoint) THEN
	      Nyzlov(1)=INpoint3-1
	      Nyzlov(2)=INpoint3
	      Nyzlov(3)=INpoint3+1
           ELSE
            Nyzlov(1)=INpoint3-2
	      Nyzlov(2)=INpoint3-1
	      Nyzlov(3)=INpoint3
    	   ENDIF
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,R,A0rezF,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3) 
	   ! онксвюел гмювемхе тсмйжхх б рнвйх Xz3
	   RRYRA0(3)=A*Xz3**2+B*Xz3+C
         ! нопедекъел йнщттхжхемрш юопнйяхлюжхх ALFA  мю дюммнл хмрепбюке
	   Nyzlov(1)=1
	   Nyzlov(2)=2
	   Nyzlov(3)=3
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RRER,RRYRALFA,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   ALFA(3,IZXDC)=RcoffPol(1)
	   ALFA(2,IZXDC)=RcoffPol(2)
	   ALFA(1,IZXDC)=RcoffPol(3) 
         ! опнбепъел йнщттхжхемр опх люйяхлюкэмни яреоемх
         ! нрпхжюрекэмши йнщттхжхемр меопхелкел
	!   IF(ALFA(3,IZXDC).LT.0.D0) THEN
	!     ! нопедекъел йнщттхжхемрш юопнйяхлюжхх ALFA  мю дюммнл хмрепбюке
      !     Nyzlov(1)=1
	!     Nyzlov(2)=3
	! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      !call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,1,Nyzlov,RRER,RRYRALFA,
      !*RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	!     ALFA(3,IZXDC)=0.D0
	!     ALFA(2,IZXDC)=RcoffPol(1)
      !     ALFA(1,IZXDC)=RcoffPol(2) 
   	!   ENDIF 

         ! нопедекъел йнщттхжхемрш юопнйяхлюжхх A0  мю дюммнл хмрепбюке
	   Nyzlov(1)=1
	   Nyzlov(2)=2
	   Nyzlov(3)=3
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA(0,2,Nyzlov,RRER,RRYRA0,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A0(3,IZXDC)=RcoffPol(1)
	   A0(2,IZXDC)=RcoffPol(2)
	   A0(1,IZXDC)=RcoffPol(3)  
      ENDDO

      !WRITE(6,*) 'COFF',Ninterval 
      !DO IZXDC=1,Ninterval 
	!WRITE(6,*)
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC),IZXDC
	!WRITE(6,*)  ALFA(1,IZXDC),ALFA(2,IZXDC),ALFA(3,IZXDC)
      !WRITE(6,*)  A0(1,IZXDC),A0(2,IZXDC),A0(3,IZXDC)
      !ENDDO



	! щрюо 5. онксвюел онкхмнл мю йюфднл хмрепбюке йюй опнхгбедемхе онкхмнлю A0 х ACoffPolinom  
      ARCoffPolinom=0.D0
	! жхйк он хмрепбюкюл
	DO IZXDC=1,Ninterval 
	   ! жхйк он йнщттхжхемрюл онкхмнлю A0
	   DO IZBFG=1,3
            !жхйк он йнщттхжхемрюл онкхмнлю ACoffPolinom 
		  DO INDEX=1,NumbreZero
      ARCoffPolinom(IZBFG+INDEX-1,IZXDC)=ARCoffPolinom(IZBFG+INDEX-1,IZXDC)+A0(IZBFG,IZXDC)*ACoffPolinom(INDEX)
	      ENDDO
         ENDDO
      ENDDO
      

      !WRITE(6,*) 'APOLINOM',(ACoffPolinom(INDEX),INDEX=1,NumbreZero)
	!DO IZXDC=1,Ninterval 
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC)
	!WRITE(6,*)  A0(1,IZXDC),A0(2,IZXDC),A0(3,IZXDC)
      !ENDDO
	!WRITE(6,*) 'REZ',2+NumbreZero
	!DO IZXDC=1,Ninterval 
      !WRITE(6,*)  Xlim(1,IZXDC),Xlim(2,IZXDC)
	!WRITE(6,*)  (ARCoffPolinom(IZBFG,IZXDC),IZBFG=1,2+NumbreZero)
      !ENDDO







  
      
	!WRITE(6,*)
	!WRITE(6,*)
      ! бшдюел пегскэрюр пюяверю 
      ! онксвюел юопнйяхлхпнбюммсч тсмйжхч (дкъ опнбепйх)
      DO IZBFG=1,Npoint 
         DO IZXDC=1,Ninterval 
        IF(R(IZBFG).GE.Xlim(1,IZXDC).AND.R(IZBFG).LE.Xlim(2,IZXDC)) THEN
	      INpoint1=IZXDC
	  ENDIF
	   ENDDO 
         SUMDFG=0.D0
         DO IZXDC=1,2+NumbreZero
         SUMDFG=SUMDFG+ARCoffPolinom(IZXDC,INpoint1)*R(IZBFG)**(IZXDC-1)
	   ENDDO
         SUMDFG=SUMDFG*R(IZBFG)**(l+1)
	   ALFA12=0.D0
	   ALFA12=ALFA12-ALFA(3,INpoint1)*R(IZBFG)**3
         ALFA12=ALFA12-ALFA(2,INpoint1)*R(IZBFG)**2
         ALFA12=ALFA12-ALFA(1,INpoint1)*R(IZBFG)
		 F1=DEXP(ALFA12)*SUMDFG
	     RfunAro(IZBFG)=F1
	     !WRITE(6,*) R(IZBFG),F1,Rfun(IZBFG)  
      ENDDO

      ! яреоемэ онкхмнлю я йнщттхжхемрюлх ARCoffPolinom
	NpolAR=1+NumbreZero


   

     
     

	
    

   

      ! сдюкемхе люяяхбнб хг оълърх 
      deallocate(ALFArezF,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "ALFArezF" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(A0rezF,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "A0rezF" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(Xzeex,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Xzeex" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(ACoffPolinom,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "ACoffPolinom" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(A0,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "A0" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRER,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRER" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRYRALFA,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRYRALFA" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RRYRA0,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RRYRA0" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
      deallocate(RSSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RSSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(ALFASSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "ALFASSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(A0SSD,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "A0SSD" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(NrazInt,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "NrazInt" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(RlimTTR,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RlimTTR" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
	deallocate(Rresh,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Rresh" IS NOT REMOVED FROM MEMORY'
	stop 
	endif 
      deallocate(Xzero,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Xzero" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
   	deallocate(Xextremum,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Xextremum" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Nyzlov,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Nyzlov" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(RcoffPol,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "RcoffPol" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(Amatrix,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Amatrix" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Ycoff,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_APPROXIMATION_RADIAL_FUNCTION'
      write(*,*) 'THE FILE "Ycoff" IS NOT REMOVED FROM MEMORY'
	stop 
	endif


      return
      end subroutine EFSH_APPROXIMATION_RADIAL_FUNCTION_ALFA 




     ! SUB-PROGRAM PERFORM ANALYSIS OF FUNCTION (ONE VARIABLE)
     ! SETS FUNCTION ZERO AND EXTREME POINTS
     ! DESCRIPTION OF SUBPROGRAM PARAMETERS
     ! Npoint-NUMBER OF DOTS FUNCTIONS
     ! R (Npoint) - ARGUMENT VALUE ARRAY
     ! RFUN (Npoint) -ARRAY OF FUNCTION VALUES
     ! NpointZERO-NUMBER OF ZERO FUNCTIONS
     ! Xzero (NpointZERO) - ARGUMENT OF VALUES OF A FUNCTION ARGUMENT IN WHICH IT VALUES TO ZERO
     ! NpointEXTREMA-NUMBER OF EXTREME POINTS
     ! Xextrema (NpointEXTREMA) - ARRAY OF VALUES OF EXTREME POINTS (ARGUMENT OF THE FUNCTION)
	subroutine EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA(Npoint,R,RFUN,NpointZERO,Xzero,NpointEXTREMA,Xextrema) 
      implicit none
	  integer::Npoint,NpointZERO,NpointEXTREMA
	  real(8),dimension(:)::R,RFUN,Xzero,Xextrema
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  integer::IZXDC,INDEX,ierr
	  real(8)::Xz1,Xz2,A,B,C
	  integer,allocatable,dimension(:)::Nyzlov
      real(8),allocatable,dimension(:)::RcoffPol,RFUNpro

      
      ! бшдекъел оюлърэ дкъ люяяхбнб
      allocate(RcoffPol(3),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
	     write(*,*) 'MEMORY ON THE FILE "RcoffPol" IS NOT SELECTED'
	     stop 
	  endif
  	  allocate(Nyzlov(3),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
	     write(*,*) 'MEMORY ON THE FILE "Nyzlov" IS NOT SELECTED'
	     stop 
	  endif
   	  allocate(RFUNpro(Npoint),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
	     write(*,*) 'MEMORY ON THE FILE "RFUNpro" IS NOT SELECTED'
	     stop 
	  endif

      ! гюмскъел оепед пюявернл 
	  Xzero=0.D0
	  Xextrema=0.D0
      
	  !WRITE(100,*) Npoint  
      !DO IZXDC=1,Npoint
      !   WRITE(100,*) R(IZXDC),RFUN(IZXDC)
	  !ENDDO
      


      ! оепбши мскэ
	  Xzero(1)=0.D0  ! яннрберярбсер гмюмскемхч тсмйжхх б мювюке хмрепбюкю
      NpointZERO=1   ! пюдхюкэмюъ вюярэ бнкмнбни тсмйжхх дхяйпермнцн яоейрпю 
	                 ! днкфмю гюмскъряъ б рнвйюу r=0 Х r=Infinity 

      ! нясыеярбкъел онхяй мскеи
	DO IZXDC=1,Npoint-2 ! ХЯЙКЧВЮЕЛ ХГ ПЮЯЯЛНРПЕМХЪ ЙПЮИМХЧ РНВЙС Npoint,
	                      ! Б ЩРНИ РНВЙХ ТСМЙЖХЪ АКХГЙЮ Й МСКЧ НДМЮЙН ЩРН МЕ ЯБЪГЮМН Я БНГМХЙМНБЕМХЕЛ ОЕПЕЯЕВЕМХЪ Я НЯЭЧ r   
	     ! бшъбкъел мскэ тсмйжхх
	IF(RFUN(IZXDC).GE.0.D0.AND.RFUN(IZXDC+1).LT.0.D0.OR.RFUN(IZXDC).LT.0.D0.AND.RFUN(IZXDC+1).GE.0.D0) THEN 
             ! мюидем мскэ тсмйжхх
	       NpointZERO=NpointZERO+1
             ! мюундхл гмювемхе R опх йнрнпнл тсмйжхъ пюбмю мскч
		   Nyzlov(1)=IZXDC-1
		   Nyzlov(2)=IZXDC
		   Nyzlov(3)=IZXDC+1
	       ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA2(Nyzlov,R,RFUN,RcoffPol)
	     ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
		   A=RcoffPol(1)
		   B=RcoffPol(2)
		   C=RcoffPol(3) 
             ! ОПХПЮБМХБЮЪ ОНКСВЕММШИ ОНКХМНЛ Й МСКЧ ЛШ РЕЛ ЯЮЛШЛ ОНКСВЮЕЛ ЙНПЕМЭ СПЮБМЕМХЪ (мскэ тсмйжхх)
	       Xz1=(-B+DSQRT(B**2-4.D0*A*C))/(2.D0*A)
		   Xz2=(-B-DSQRT(B**2-4.D0*A*C))/(2.D0*A)
	       ! БШЪЯМЪЕЛ ЙЮЙНИ ХГ ЩРХУ ГМЮВЕМХИ ЯННРБЕРЯРБСЕР ЙНПМЧ
		   INDEX=0
		   IF(Xz1.GT.R(IZXDC).AND.Xz1.LT.R(IZXDC+1)) THEN
		      Xzero(NpointZERO)=Xz1
	          INDEX=1
	       ENDIF 
             IF(Xz2.GT.R(IZXDC).AND.Xz2.LT.R(IZXDC+1)) THEN
		      Xzero(NpointZERO)=Xz2
	          INDEX=1
             ENDIF 
	       ! опнбепъел мюидем кх йнпемэ
	       IF(INDEX.EQ.0) THEN
                WRITE(*,*) 'INDEX=0, THE ROOT IS NOT FOUND'
	          READ(*,*)
	          STOP
             ENDIF
	     ENDIF
      ENDDO


      ! гюмскъел оепед пюанрни
      RFUNpro=0.D0
      ! ярпнхл оепбсч опнхгбндмсч тсмйжхх
	DO IZXDC=3,Npoint-2 ! хяйкчвюел хг пюялнрпемхъ йпюимхе рнвйх 
	                    ! б мху тсмйжхъ пюбмю мскч ( ме хяонкэгсел йпюимхх рнвйх б пюявере) 
         
         ! мюундхл гмювемхе R опх йнрнпнл тсмйжхъ пюбмю мскч
	   Nyzlov(1)=IZXDC-1
	   Nyzlov(2)=IZXDC
	   Nyzlov(3)=IZXDC+1
	   ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA2(Nyzlov,R,RFUN,RcoffPol)
	   ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
	   A=RcoffPol(1)
	   B=RcoffPol(2)
	   C=RcoffPol(3)
	   RFUNpro(IZXDC)=2.D0*A*R(IZXDC)+B
	!   WRITE(6,*)  R(IZXDC),RFUN(IZXDC),RFUNpro(IZXDC)
      ENDDO


       ! нясыеярбкъел онхяй мскеи опнхгбндмни ( рнвйх щйярпелслю) 
	NpointEXTREMA=0	
	DO IZXDC=3,Npoint-3     ! бшвепйхбюел йнмевмсч дкъ опнхгбндмни рнвйс 
	      ! бшъбкъел мскэ тсмйжхх
	IF(RFUNpro(IZXDC).GE.0.D0.AND.RFUNpro(IZXDC+1).LT.0.D0.OR.RFUNpro(IZXDC).LT.0.D0.AND.RFUNpro(IZXDC+1).GE.0.D0) THEN 
             ! мюидем мскэ тсмйжхх
	       NpointEXTREMA=NpointEXTREMA+1
             ! мюундхл гмювемхе R опх йнрнпнл тсмйжхъ пюбмю мскч
		   Nyzlov(1)=IZXDC-1
		   Nyzlov(2)=IZXDC
		   Nyzlov(3)=IZXDC+1
	       ! нОПЕДЕКЪЕЛ ЙНЩТТХЖХЕМРШ ОНКХМНЛЮ БРНПНЦН ОНПЪДЙЮ ( ОН ЛЕРНДС йПЮЛЕПЮ)
      call EFSH_COEFFICIENT_POLINOM_KRAMERA2(Nyzlov,R,RFUNpro,RcoffPol)
	     ! йнщттхжхемрш онкхмнлю брнпнцн онпъдйю
		   A=RcoffPol(1)
		   B=RcoffPol(2)
		   C=RcoffPol(3) 
             ! ОПХПЮБМХБЮЪ ОНКСВЕММШИ ОНКХМНЛ Й МСКЧ ЛШ РЕЛ ЯЮЛШЛ ОНКСВЮЕЛ ЙНПЕМЭ СПЮБМЕМХЪ (мскэ тсмйжхх)
	       Xz1=(-B+DSQRT(B**2-4.D0*A*C))/(2.D0*A)
		   Xz2=(-B-DSQRT(B**2-4.D0*A*C))/(2.D0*A)
	       ! БШЪЯМЪЕЛ ЙЮЙНИ ХГ ЩРХУ ГМЮВЕМХИ ЯННРБЕРЯРБСЕР ЙНПМЧ
		   INDEX=0
		   IF(Xz1.GT.R(IZXDC).AND.Xz1.LT.R(IZXDC+1)) THEN
		      Xextrema(NpointEXTREMA)=Xz1
	          INDEX=1
	       ENDIF 
             IF(Xz2.GT.R(IZXDC).AND.Xz2.LT.R(IZXDC+1)) THEN
		      Xextrema(NpointEXTREMA)=Xz2
	          INDEX=1
             ENDIF 
	       ! опнбепъел мюидем кх йнпемэ
	       IF(INDEX.EQ.0) THEN
                WRITE(*,*) 'INDEX=0, THE ROOT IS NOT FOUND'
	          READ(*,*)
	          STOP
             ENDIF
	     ENDIF
      ENDDO




	 





	deallocate(RcoffPol,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
      write(*,*) 'THE FILE "RcoffPol" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      deallocate(Nyzlov,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
      write(*,*) 'THE FILE "Nyzlov" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	deallocate(RFUNpro,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA'
      write(*,*) 'THE FILE "RFUNpro" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
      
	return
      end subroutine EFSH_ANALYSIS_FUNCTION_ZERO_AND_EXTREMA


    ! POLYNOMA ROOT CALCULATION SUBPROGRAM
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! Npol-MAXIMUM POLYNOMA DEGREE
    ! RcoffPol (Npol + 1) -ARRAY OF COEFFICIENTS
    ! Ns-NUMBER OF REAL ROOTS OF POLYNOMA
    ! X (Npol) -ROOT POLYNOMA
	subroutine EFSH_CALCULATION_ROOTS_POLYNOM(Npol,RcoffPol,Ns,X)
      use dfimsl
	implicit none
      integer::Npol,Ns,INDEXSUM,IYDS,ierr
      real(8),dimension(:)::RcoffPol,X
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	complex(8),allocatable,dimension(:)::Root
	
	!бшдекъел оюлърэ дкъ люяяхбнб
	allocate(Root(Npol),stat=ierr)
	if(ierr/=0) then
      write(*,*) 'EFSH_CALCULATION_ROOTS_POLYNOM'
	write(*,*) 'MEMORY ON THE FILE "Root" IS NOT SELECTED'
	stop 
	endif

      ! гюмскъел оепед пюявернл
	X=0.D0
      Root=0.D0

	! НЯСЫЕЯРБКЪЕЛ МЮУНФДЕМХЕ ЙНПМЕИ
	call dzplrc(Npol,RcoffPol,Root)

	! нопедекъел беыеярбеммши йнпмх
      INDEXSUM=0
      do IYDS=1,Npol
         ! ОПНБЕПЪЕЛ БЕЫЕЯРБЕММНЯРЭ ЙНПМЪ
	   if(dimag(Root(IYDS)).EQ.0.D0) then
	      INDEXSUM=INDEXSUM+1
	      X(INDEXSUM)=DREAL(Root(IYDS))
	   endif 
      enddo
      
	! гюохяшбюел вхякн беыеярбеммшу йнпмеи
      Ns=INDEXSUM


	! сдюкемхе люяяхбнб хг оълърх   
      deallocate(Root,stat=ierr)
	if(ierr/=0) then
	write(*,*) 'EFSH_CALCULATION_ROOTS_POLYNOM'
      write(*,*) 'THE FILE "Root" IS NOT REMOVED FROM MEMORY'
	stop 
	endif
	return
      end subroutine EFSH_CALCULATION_ROOTS_POLYNOM

    

   
    ! SUBPROGRAM FOR SOLVING A SYSTEM OF LINEAR EQUATIONS A * Xs = Y
    ! THE GAUSS METHOD USING A PARTIAL SELECTION SCHEME
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! N-NUMBER OF UNKNOWN
    ! A-MATRIX OF A SYSTEM OF LINEAR EQUATIONS
    ! Y-ARRAY FUNCTION VALUE
    ! Xs-ARRAY OF ROOT SYSTEM VALUES
      subroutine EFSH_SYSTEM_LINEAR_EQUATIONS(N,A,Y,Xs) 
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
	call DLufac2(N,A,fac,ipvt)
      ! пЕЬЮЕЛ КХМЕИМСЧ ЯХЯРЕЛС AXs=Y
	call useLU2(N,fac,ipvt,Y,Xs)  

	
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
      end subroutine EFSH_SYSTEM_LINEAR_EQUATIONS

      ! бяонлнцюрекэмюъ ондпнцпюллю
      subroutine useLU2(N,fac,ipvt,b,x)
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
      end subroutine useLU2


	
	
end module mefsh
