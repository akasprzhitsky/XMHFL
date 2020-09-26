      ! THE MODULE OF CALCULATION OF MOLECULAR ANGULAR VALUES VER 1.0 11.2019
	  ! VER 1.0 NEW  11,2019 √Œƒ

	! MODULE FOR CALCULATION OF MOLECULAR ANGULAR VALUES VERSION 1.0
	

module mcmav
  implicit none
	
      

	contains



    ! SUBPROGRAM FOR CALCULATING THE ANGULAR STRUCTURE OF THE INTERACTION OPERATOR
    ! NUCLEI OF A MOLECULE AND ELECTRONS (ELECTRON IN A CRYSTALLINE FIELD)
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! N-NUMBER OF NUCLEI
    ! L1-ORBITAL MOMENT OF THE FIRST ELECTRON
    ! ML1-PROJECTION OF THE ORBITAL MOMENTUM OF THE FIRST ELECTRON
    ! L2-ORBITAL MOMENT OF SECOND ELECTRON
    ! ML2-PROJECTION OF THE ORBITAL MOMENTUM OF THE SECOND ELECTRON
    ! AC (N, 2) - ARRAY OF CORNER COORDINATES OF NUCLEI (IN RADIANS)
    ! AC (N, 1) -ANGLE THETA (0 <X <= pi),
    ! AC (N, 2) -PI ANGLE (0 <Y <= 2 * pi)
    ! RcoffCF (N, Numbre, 2) - ARRAY OF ANGULAR COEFFICIENTS OF THE CRYSTALLINE FIELD
    ! RcoffCF (N, Numbre, 1) -REAL PART OF ANGULAR FACTOR (FIELD) 
    ! RcoffCF (N, Numbre, 2) -IME PART OF ANGULAR FACTOR (FIELD)
    ! EXPLANATION Numbre-NUMBER OF HARMONICS
	subroutine CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD(N,L1,ML1,L2,ML2,AC,RcoffCF)
	  use mc3js,only:C3JS_VAR,C3JS_Ckq,C3JS_CkqGG,C3JS_CkqISH
      implicit none
      
	  integer::N,L1,ML1,L2,ML2
      real(8),dimension(:,:)::AC
      real(8),dimension(:,:,:)::RcoffCF
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer,parameter::NmassivVar=800
	  integer ::Numbre,IX,IY,ierr
	  real(8):: FF,MatrixEleCkq,RAS
	  real(8),parameter::PInumbre=3.14159265358979D0
      real(8),allocatable,dimension(:)::Ykq,YSkq,F

      !¬€ƒ≈Àﬂ≈Ã œ¿Ãﬂ“‹ ƒÀﬂ Ã¿——»¬Œ¬
	  allocate(Ykq(2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "Ykq" IS NOT SELECTED'
	     stop 
	  endif
      allocate(YSkq(2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "Ykq" IS NOT SELECTED'
	     stop 
	  endif
      allocate(F(NmassivVar),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "F" IS NOT SELECTED'
	     stop 
	  endif


      ! «¿Õ”Àﬂ≈Ã œ≈–≈ƒ –¿—◊≈“ŒÃ
      RcoffCF=0.D0
	  Ykq=0.D0
	  YSkq=0.D0
      F=0.D0
      ! ¬—œŒÃŒ√¿“≈À‹Õ¿ﬂ œ–Œ÷≈ƒ”–¿ ƒÀﬂ –¿—◊≈“¿ Ã¿“–»◊Õ€’ ›À≈Ã≈Õ“Œ¬ —‘≈–»◊≈— Œ… √¿–ÃŒÕ» »
	  !call C3JS_VAR(NmassivVar,FF,F)


      Numbre=0
      ! ÷» À œŒ √¿–ÃŒÕ» 
      DO IX=IABS(L1-L2),L1+L2,2
	     Numbre=Numbre+1 
	     ! Œ—”Ÿ≈—“¬Àﬂ≈Ã –¿—◊≈“ «Õ¿◊≈Õ»… Ã¿“–»◊ÕŒ√Œ ›À≈Ã≈Õ“¿ —‘≈–»◊≈— Œ… √¿–ÃŒÕ» »
	     MatrixEleCkq=C3JS_CkqGG(L1,ML1,L2,ML2,IX) !C3JS_Ckq(L1,ML1,L2,ML2,IX,FF,F)
		 !write(6,*) 'ckq',MatrixEleCkq !,C3JS_CkqISH(L1,ML1,L2,ML2,IX) !,C3JS_Ckq(L1,ML1,L2,ML2,IX,FF,F)
	   	 !  Œ›‘‘»÷»≈Õ“
	     RAS=DSQRT(4.d0*PInumbre/float(2*IX+1))
         ! ÷» À œŒ ÕŒÃ≈–¿Ã ﬂƒ≈–
	     DO IY=1,N
            ! Œ—”Ÿ≈—“¬Àﬂ≈Ã –¿—◊≈“ «Õ¿◊≈Õ»… —‘≈–»◊≈— Œ… √¿–ÃŒÕ» »
            call CMAV_VALUE_SPHERICAL_FUNCTION(IX,ML1-ML2,AC(IY,1),AC(IY,2),Ykq,YSkq)
	        ! –≈¿À‹Õ¿ﬂ ◊¿—“‹  Œ›‘‘»÷»≈Õ“¿
	        RcoffCF(IY,Numbre,1)=MatrixEleCkq*RAS*YSkq(1)
	        ! ÃÕ»Ã¿ﬂ ◊¿—“‹  ¿›‘‘»÷»≈Õ“¿
            RcoffCF(IY,Numbre,2)=MatrixEleCkq*RAS*YSkq(2)
         ENDDO
      ENDDO


      ! ”ƒ¿À≈Õ»≈ Ã¿——»¬Œ¬ »« œﬂÃﬂ“»   
      deallocate(Ykq,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "Ykq" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(YSkq,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "YSkq" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(F,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "F" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif


      return
    end subroutine CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD










     ! SUBPROGRAM FOR CALCULATING THE NUMERICAL VALUE OF THE SPHERICAL HARMONY FOR K AND ALL PROJECTIONS
     ! DESCRIPTION OF SUBPROGRAM PARAMETERS
     ! SPHERICAL HARMONY K-NUMBER
     ! RTeta-angle in the theta spherical coordinate system (0 <RTeta <= pi)
     ! RFu-angle in spherical coordinates phi (0 <RFu <= 2 * pi)
     ! Ykq (2 * K + 1,2) -MASSIF OF VALUES OF THE SPHERICAL HARMONIC OF ALL PROJECTIONS FOR A GIVEN K
     ! Ykq (2 * K + 1,1) -REAL PART OF SPHERICAL HARMONICS
     ! Ykq (2 * K + 1,2) - IMAGINARY PART OF SPHERICAL HARMONICS
     ! YSkq (2 * K + 1,2) -MASSIF OF VALUES OF THE SPHERICAL HARMONIC OF COMPLEX-CONJUGATED ALL PROJECTIONS FOR A GIVEN K
	subroutine CMAV_SPHERICAL_FUNCTION(K,RTeta,RFu,Ykq,YSkq)
	  implicit none
      
	  integer::K
	  real(8)::RTeta,RFu
      real(8),dimension(:,:)::Ykq,YSkq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  real(8),parameter::PInum=3.14159265358979D0
	  integer::IZX,IZX1q
      real(8)::RTETAF
      
      ! «¿Õ”Àﬂ≈Ã œ≈–≈ƒ –¿—◊≈“ŒÃ
      Ykq=0.D0
      YSkq=0.D0

      IZX1q=-K
	  ! ÷» À œŒ œŒÀŒ∆»“≈À‹Õ€Ã œ–Œ≈ ÷»ﬂÃ √¿–ÃŒÕ» » ¬ Àﬁ◊¿ﬂ ÕŒÀ‹
      do IZX=K+1,2*K+1
         ! œ–Œ≈ ÷»ﬂ √¿–ÃŒÕ» »
	     IZX1q=IZX1q+(IZX-1)
         ! ¬€◊»—Àﬂ≈Ã «Õ¿◊≈Õ»ﬂ “≈“¿ ‘”Õ ÷»»
         RTETAF=CMAV_VALUE_TETA_FUNCTION_POINT(K,IZX1q,RTeta)
         ! «¿œ»—€¬¿≈Ã √¿–ÃŒÕ» ” — ÕŒÃ≈–ŒÃ K » œ–Œ≈ ÷»≈… IZX1q
	     ! –≈¿À‹Õ¿ﬂ ◊¿—“‹ 
	     Ykq(IZX,1)=DCOS(float(IZX1q)*RFu)*RTETAF/DSQRT(2.D0*PInum)
         ! ÃÕ»Ã¿ﬂ ◊¿—“‹
	     Ykq(IZX,2)=DSIN(float(IZX1q)*RFu)*RTETAF/DSQRT(2.D0*PInum)
         ! «¿œ»—€¬¿≈Ã √¿–ÃŒÕ» ” — ÕŒÃ≈–ŒÃ   » œ–Œ≈ ÷»≈…  -IZX1q
         ! –≈¿À‹Õ¿ﬂ ◊¿—“‹ 
	     Ykq(K+1-IZX1q,1)=(-1.D0)**(IZX1q)*Ykq(IZX,1)
         ! ÃÕ»Ã¿ﬂ ◊¿—“‹
	     Ykq(K+1-IZX1q,2)=(-1.D0)**(IZX1q+1)*Ykq(IZX,2)
	  enddo

      ! œŒÀ”◊¿≈Ã  ŒÃœÀ≈ —ÕŒ —Œœ–ﬂ∆≈ÕÕ€≈ √¿–ÃŒÕ» »
	  IZX1q=-K
	  do IZX=1,2*K+1
         ! œ–Œ≈ ÷»ﬂ √¿–ÃŒÕ» »
	     IZX1q=IZX1q+(IZX-1)
         ! –≈¿À‹Õ¿ﬂ ◊¿—“‹
	     YSkq(IZX,1)=DCOS(float(2*IZX1q)*RFu)*Ykq(IZX,1)
	     YSkq(IZX,1)=YSkq(IZX,1)+DSIN(float(2*IZX1q)*RFu)*Ykq(IZX,2)
	     ! ÃÕ»Ã¿ﬂ ◊¿—“‹
         YSkq(IZX,2)=DCOS(float(2*IZX1q)*RFu)*Ykq(IZX,2)
         YSkq(IZX,2)=YSkq(IZX,2)-DSIN(float(2*IZX1q)*RFu)*Ykq(IZX,1)
      enddo
	 
	  return
    end subroutine CMAV_SPHERICAL_FUNCTION


    ! SUB-PROGRAM FOR CALCULATING THE NUMERICAL VALUE OF THE SPHERICAL HARMONICS L, M
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! L-ORBITAL MOMENT OF SPHERICAL HARMONY
    ! M-PROJECTION OF THE ORBITAL MOMENT OF SPHERICAL HARMONIC
    ! RTeta-angle in the theta spherical coordinate system (0 <RTeta <= pi)
    ! RFu-angle in spherical coordinates phi (0 <RFu <= 2 * pi)
    ! Ykq (2) -ARRAY OF SPHERICAL HARMONIC VALUES FOR A GIVEN L AND PROJECTION M
    ! Ykq (1) -REAL PART OF SPHERICAL HARMONICS
    ! Ykq (2) -IME PART OF THE SPHERICAL HARMONIC
    ! YSkq (2) -ARRAY OF VALUES OF THE SPHERICAL HARMONICS OF THE INTEGRATED-CONNECTED FOR A GIVEN K AND PROJECTION M
    subroutine CMAV_VALUE_SPHERICAL_FUNCTION(L,M,RTeta,RFu,Ykq,YSkq)
	  implicit none
      
	  integer::L,M
	  real(8)::RTeta,RFu
      real(8),dimension(:)::Ykq,YSkq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	  real(8),parameter::PInum=3.14159265358979D0
      real(8)::RTETAF,REY,RIMY
       
      ! «¿Õ”Àﬂ≈Ã œ≈–≈ƒ –¿—◊≈“ŒÃ
      Ykq=0.D0
      YSkq=0.D0

      ! ¬€◊»—Àﬂ≈Ã «Õ¿◊≈Õ»ﬂ “≈“¿ ‘”Õ ÷»»
      RTETAF=CMAV_VALUE_TETA_FUNCTION_POINT(L,IABS(M),RTeta)
      REY=DCOS(float(IABS(M))*RFu)*RTETAF/DSQRT(2.D0*PInum)
      RIMY=DSIN(float(IABS(M))*RFu)*RTETAF/DSQRT(2.D0*PInum)
      ! ”—“¿Õ¿¬À»¬¿≈Ã «Õ¿  œ–Œ≈ ÷»»  
      IF(M.GE.0) THEN
         ! –≈¿À‹Õ¿ﬂ ◊¿—“‹ 
	     Ykq(1)=REY
         ! ÃÕ»Ã¿ﬂ ◊¿—“‹
	     Ykq(2)=RIMY
        ELSE
         ! –≈¿À‹Õ¿ﬂ ◊¿—“‹ 
	     Ykq(1)=(-1.D0)**IABS(M)*REY
         ! ÃÕ»Ã¿ﬂ ◊¿—“‹
	     Ykq(2)=(-1.D0)**(IABS(M)+1)*RIMY
      ENDIF
      

      ! œŒÀ”◊¿≈Ã  ŒÃœÀ≈ —ÕŒ —Œœ–ﬂ∆≈ÕÕ€≈ √¿–ÃŒÕ» »
      ! –≈¿À‹Õ¿ﬂ ◊¿—“‹
      YSkq(1)=DCOS(float(2*M)*RFu)*Ykq(1)
      YSkq(1)=YSkq(1)+DSIN(float(2*M)*RFu)*Ykq(2)
	  ! ÃÕ»Ã¿ﬂ ◊¿—“‹
      YSkq(2)=DCOS(float(2*M)*RFu)*Ykq(2)
      YSkq(2)=YSkq(2)-DSIN(float(2*M)*RFu)*Ykq(1)
      
      return
    end subroutine CMAV_VALUE_SPHERICAL_FUNCTION


 
    ! SUB-PROGRAM FOR OBTAINING THE VALUE OF THETA FUNCTION (RECURRENT METHOD)
    ! DESCRIPTION OF SUBPROGRAM PARAMETERS
    ! L-ORBITAL MOMENT FUNCTION
    ! M-PROJECTION OF THE ORBITAL MOMENTUM (THE PROGRAM CONSIDERS THE CASE (M> = 0)
    ! X-VALUE OF ARGUMENT (INTERVAL FROM (0, PI))
   real(8) function CMAV_VALUE_TETA_FUNCTION_POINT(L,M,X)
     implicit none
      
	 integer::L,M
     real(8)::X
     !!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IAA2
     real(8)::RHR,Rcoff1,Rcoff2,Rcoff3,RSDXF,RSDXX

     
     ! —À¿√¿≈ÃŒ≈ Teta M,M
	 RSDXF=DSQRT(CMAV_FACTORIAL(2*IABS(M)))/CMAV_FACTORIAL(IABS(M))
	 RSDXX=(DSIN(X*0.5D0)*DCOS(X*0.5D0))**IABS(M)
     RSDXX=RSDXX*RSDXF/DSQRT(2.D0)
	 Rcoff1=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+1))*RSDXX
     !  —À¿√¿≈ÃŒ≈ Teta M+1,M 
	 RSDXF=DSQRT(CMAV_FACTORIAL(2*IABS(M)+1))/CMAV_FACTORIAL(IABS(M))
	 RSDXX=(DSIN(X*0.5D0)*DCOS(X*0.5D0))**IABS(M)
     RSDXX=RSDXX*RSDXF*DCOS(X)/DSQRT(2.D0)
	 Rcoff2=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+3))*RSDXX
	 IF(IABS(M).EQ.L) Rcoff3=Rcoff1
     IF(IABS(M).EQ.L-1) Rcoff3=Rcoff2
     ! ÷» À œŒ ÃŒÃ≈Õ“¿Ã
	 DO IAA2=IABS(M),L-2
        RSDXF=SQRT(float(2*IAA2+3)*float(2*IAA2+5))
        RSDXF=RSDXF/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
          
	    RSDXX=SQRT(float(2*IAA2+5))
        RSDXX=RSDXX*SQRT(float(IAA2+1-IABS(M))*float(IAA2+1+IABS(M)))
        RSDXX=RSDXX/SQRT(float(2*IAA2+1))
        RSDXX=RSDXX/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
        ! –¿—◊≈“ —À≈ƒ”ﬁŸ≈√Œ ›À≈Ã≈Õ“¿
	    Rcoff3=RSDXF*DCOS(X)*Rcoff2-RSDXX*Rcoff1
	    ! œŒƒ√Œ“Œ¬ ¿   –¿—◊≈“” —À≈ƒ”ﬁŸ≈√Œ ›À≈Ã≈Õ“¿
        Rcoff1=Rcoff2
	    Rcoff2=Rcoff3
     ENDDO
      
	 ! «¿œ»—€¬¿≈Ã «Õ¿◊≈Õ»≈ ‘”Õ ÷»»
	 CMAV_VALUE_TETA_FUNCTION_POINT=Rcoff3 

     return
    end function CMAV_VALUE_TETA_FUNCTION_POINT


     
    ! œŒƒœ–Œ√–¿ÃÃ¿ –¿—◊≈“¿ ‘¿ “Œ–»¿À¿
	! Œœ»—¿Õ»≈ œ¿–¿Ã≈“–Œ¬
    ! NOPI-÷≈ÀŒ≈ œŒÀŒ∆»“≈À‹ÕŒ≈ ◊»—ÀŒ 
	real(8) function CMAV_FACTORIAL(NOPI)
      implicit none

	  integer:: NOPI,IBKL
      real(8):: PROSS
      
	  IF(NOPI.LT.0) THEN
         WRITE(*,*) 'ERROR FACTORIAL (N<0)'
         READ(*,*)
	     STOP
 	  ENDIF

      IF(NOPI.EQ.0) THEN
         CMAV_FACTORIAL=1.D0
         return
	  ENDIF
      
	  PROSS=1.D0
	  do IBKL=1,NOPI
         PROSS=PROSS*float(IBKL)
	  enddo

	  CMAV_FACTORIAL=PROSS

   
      return
    end function CMAV_FACTORIAL
	
end module mcmav
