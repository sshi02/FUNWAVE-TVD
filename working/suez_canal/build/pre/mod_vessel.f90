












!------------------------------------------------------------------------------------
!
!      FILE vessle.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!  1 is a module to model ship-wakes    
!
!  HISTORY :
!    10/07/2016  Fengyan Shi
!
!-------------------------------------------------------------------------------------


MODULE VESSEL_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,ETA0,Xco,Yco

  USE INPUT_READ
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
  IMPLICIT NONE
  SAVE

    INTEGER :: NumVessel,Kves,II,JJ
    INTEGER,DIMENSION(:),ALLOCATABLE :: VesselType
    CHARACTER(LEN=80),DIMENSION(:),ALLOCATABLE:: SourceType
    REAL(SP),DIMENSION(:),ALLOCATABLE :: Xvessel1,Yvessel1,Xvessel2,Yvessel2, &
                                       LengthVessel,WidthVessel, &
                                       AlphaVessel1,AlphaVessel2,BetaVessel,Pvessel, &
                                       TimeVessel1,TimeVessel2,ThetaVessel
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VesselPressureTotal,VesselPressureEach, &
                                       VesselPressureX,VesselPressureY
    LOGICAL :: OUT_VESSEL = .TRUE.
    REAL(SP),DIMENSION(:),ALLOCATABLE :: ResistanceX,ResistanceY,ResPosX,ResNegX, &
                                         ResPosY,ResNegY
    REAL(SP):: PLOT_INTV_VESSEL,PLOT_COUNT_VESSEL

    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VesselFluxGradient,VesselFluxGradientEach, &
                                           Gz
    REAL(SP),DIMENSION(:),ALLOCATABLE :: Uvel,Vvel

    REAL(SP) :: Xves,Yves,Lves,Wves,P_x,P_y,DetaX,DetaY

!   deep draft vessel

  
    LOGICAL :: MakeVesselDraft = .TRUE.


!INTERFACE READ_FOUR_TYPE_VALUES
!  Module Procedure VESSEL_INITIAL
!  Module Procedure VESSEL_FORCING
!END INTERFACE

CONTAINS
  
! READ 1

SUBROUTINE VESSEL_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL, INPUT_FILE_NAME
  USE GLOBAL,ONLY : iista,jjsta   !ykchoi Jan/23/2018
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::FILE_VESSEL=' '
  CHARACTER(LEN=80) :: VESSEL_FOLDER=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  CHARACTER(LEN=80)::VesselName = ' '
  INTEGER :: Ifile,ierr


! read vessel number and folder from input.txt
      FILE_NAME=INPUT_FILE_NAME

! vessel folder
      CALL READ_STRING(VESSEL_FOLDER,FILE_NAME,'VESSEL_FOLDER',ierr)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'VESSEL_FOLDER:', VESSEL_FOLDER

      CALL READ_INTEGER(NumVessel,FILE_NAME,'NumVessel',ierr)
      if (myid.eq.0) WRITE(3,'(A12,I3)') 'NumVessel = ',NumVessel

      CALL READ_LOGICAL(OUT_VESSEL,FILE_NAME,'OUT_VESSEL',ierr)

! deep draft vessel


! end deep draft vessel

      ALLOCATE (Xvessel1(NumVessel),Yvessel1(NumVessel),  &
                Xvessel2(NumVessel),Yvessel2(NumVessel),  &
                TimeVessel1(NumVessel),TimeVessel2(NumVessel), &
                LengthVessel(NumVessel),WidthVessel(NumVessel), &
                AlphaVessel1(NumVessel),AlphaVessel2(NumVessel), &
                BetaVessel(NumVessel),  &
                Pvessel(NumVessel),ThetaVessel(NumVessel),&
                ResistanceX(NumVessel),ResistanceY(NumVessel), &
                ResPosX(NumVessel),ResNegX(NumVessel), &
                ResPosY(NumVessel),ResNegY(NumVessel) )


      ALLOCATE (VesselFluxGradient(Mloc,Nloc),VesselFluxGradientEach(Mloc,Nloc), &
                Gz(Mloc,Nloc))
      ALLOCATE (Uvel(NumVessel),Vvel(NumVessel))
      ALLOCATE (SourceType(NumVessel),VesselType(NumVessel))
      
      Gz = ZERO
      Uvel = ZERO
      Vvel = ZERO

      ALLOCATE (VesselPressureTotal(Mloc,Nloc), VesselPressureEach(Mloc,Nloc),&
                VesselPressureX(Mloc,Nloc), &
                 VesselPressureY(Mloc,Nloc) )


! plot vessel intitial
     PLOT_COUNT_VESSEL = 0
     CALL READ_FLOAT(PLOT_INTV_VESSEL,FILE_NAME,'PLOT_INTV_VESSEL',ierr)
     IF(ierr==1)THEN
      if (myid.eq.0) WRITE(3,'(A50)')'PLOT_INTV_VESSEL not specified, use SMALL'
       PLOT_INTV_VESSEL = SMALL
     ENDIF
      
  DO Kves = 1, NumVessel

!  file name
    itmp1=mod(Kves/10000,10)
    itmp2=mod(Kves/1000,10)
    itmp3=mod(Kves/100,10)
    itmp4=mod(Kves/10,10)
    itmp5=mod(Kves,10)
    write(FILE_VESSEL(1:1),'(I1)')itmp1
    write(FILE_VESSEL(2:2),'(I1)')itmp2
    write(FILE_VESSEL(3:3),'(I1)')itmp3
    write(FILE_VESSEL(4:4),'(I1)')itmp4
    write(FILE_VESSEL(5:5),'(I1)')itmp5

    TMP_NAME = TRIM(VESSEL_FOLDER)//'vessel_'//TRIM(FILE_VESSEL)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), ' specified in ', TRIM(VESSEL_FOLDER), ' but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
  ENDIF

! open file
  Ifile=Kves+2000
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,'(A80)')  VesselName 
         READ(Ifile,*)  SourceType(Kves), VesselType(Kves)

         READ(Ifile,*)  ! length, width, etc
        IF(SourceType(Kves)(1:2)=='PR')THEN
         READ(Ifile,*)  LengthVessel(Kves), WidthVessel(Kves), &
                      AlphaVessel1(Kves),AlphaVessel2(Kves), &
                      BetaVessel(Kves),Pvessel(Kves)
!                for type 2 divid and volker 2017, alpha1 is cl,alpha2 is cb
!                beta is a

!   I remove the following protection, but will figure out how to avoid 
!   some crazy parameters
!          IF(VesselType(Kves)==2)THEN
!            IF(AlphaVessel1(Kves)<16.0_SP)AlphaVessel1(Kves)=16.0_SP
!            IF(AlphaVessel2(Kves)<2.0_SP)AlphaVessel2(Kves)=2.0_SP
!            IF(BetaVessel(Kves)<16.0_SP)BetaVessel(Kves) = 16.0_SP
!          ENDIF

         IF(VesselType(Kves)==1)THEN
            AlphaVessel1(Kves) = Max(SMALL, AlphaVessel1(Kves))
            BetaVessel(Kves) = Max(SMALL, BetaVessel(Kves))
            AlphaVessel1(Kves) = Min(1.0_SP, AlphaVessel1(Kves))
            BetaVessel(Kves) = Min(1.0_SP, BetaVessel(Kves))
         ENDIF

        ENDIF ! end pressure source

        IF(SourceType(Kves)(1:2)=='SL')THEN
         READ(Ifile,*)  LengthVessel(Kves), WidthVessel(Kves), &
                      AlphaVessel1(Kves),AlphaVessel2(Kves), &
                      BetaVessel(Kves),Pvessel(Kves)

        ENDIF
         READ(Ifile,*)  ! t, x, y
         READ(Ifile,*)  TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)

         TimeVessel1(Kves) = TimeVessel2(Kves)
         Xvessel1(Kves) = Xvessel2(Kves)
         Yvessel1(Kves) = Yvessel2(Kves)


   IF(MYID==0)THEN
   WRITE(3,*) '----- Vessel Name : ',  TRIM(VesselName)
   WRITE(3,*) 'Vessel Source Type: ',  TRIM(SourceType(Kves))
   WRITE(3,*) 'Vessel Type: ',  VesselType(Kves)
   WRITE(3,*) 'Vessel Length', LengthVessel(Kves)
   WRITE(3,*) 'Vessel Width', WidthVessel(Kves)
   WRITE(3,*) 'Vessel Alpha_1', AlphaVessel1(Kves)
   WRITE(3,*) 'Vessel Alpha_2', AlphaVessel2(Kves)
   WRITE(3,*) 'Vessel Beta', BetaVessel(Kves)
   WRITE(3,*) 'Vessel P', PVessel(Kves)
   WRITE(3,*) 'Initial Time, X, Y', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
   ENDIF

  ENDDO  ! end Kves

End SUBROUTINE VESSEL_INITIAL

SUBROUTINE VESSEL_FORCING
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO,DT,DEPTH
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J


  VesselFluxGradient = ZERO

  VesselPressureTotal = ZERO

  DO Kves = 1,NumVessel

    IF(TIME>TimeVessel1(Kves).AND.TIME>TimeVessel2(Kves)) THEN

         TimeVessel1(Kves)=TimeVessel2(Kves)
         Xvessel1(Kves) = Xvessel2(Kves)
         Yvessel1(Kves) = Yvessel2(Kves)

    Ifile = 2000 + Kves

! I add while to avoid dt is smaller than the time interval of vessel path 08/05/2019
    DO WHILE (TimeVessel2(Kves).LT.TIME+DT)
      READ(Ifile,*,END=120)  TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
    END DO

   IF(MYID==0)THEN
     WRITE(3,*)'Read Vessel # ', Kves
     WRITE(3,*)'T,X,Y = ', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
   ENDIF

    ThetaVessel(Kves) = ATAN2(Yvessel2(Kves)-Yvessel1(Kves),  &
                              Xvessel2(Kves)-Xvessel1(Kves))


! for panel source
    IF((TimeVessel2(Kves)-TimeVessel1(Kves))> ZERO)THEN
     Uvel(Kves) = (Xvessel2(Kves)-Xvessel1(Kves))/(TimeVessel2(Kves)-TimeVessel1(Kves))
     Vvel(Kves) = (Yvessel2(Kves)-Yvessel1(Kves))/(TimeVessel2(Kves)-TimeVessel1(Kves))
    ENDIF
! end panel source

    ENDIF ! end time > timevessel2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeVessel1(Kves))THEN
      IF(TimeVessel1(Kves).EQ.TimeVessel2(Kves))THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeVessel2(Kves)-TIME) &
            /MAX(SMALL, ABS(TimeVessel2(Kves)-TimeVessel1(Kves)))
      tmp1=1.0_SP - tmp2
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xves = Xvessel2(Kves)*tmp1 +Xvessel1(Kves)*tmp2
    Yves = Yvessel2(Kves)*tmp1 +Yvessel1(Kves)*tmp2


   IF(SourceType(Kves)(1:2)=='PR') THEN
   ! pressure source

     CALL PRESSURE_SOURCE (Xves,Yves)
    VesselPressureTotal = VesselPressureTotal+VesselPressureEach

   ENDIF

   IF(SourceType(Kves)(1:2)=='SL') THEN
   
     CALL SLENDER_BODY_SOURCE (Xves,Yves)
     VesselFluxGradient = VesselFluxGradient + VesselFluxGradientEach

   ENDIF

   IF(SourceType(Kves)(1:2)=='PA') THEN
   
     CALL GREEN_FUNCTION_SOURCE (Xves,Yves)
     VesselFluxGradient = VesselFluxGradient + VesselFluxGradientEach

   ENDIF   


120 CONTINUE  ! no more data for vessel Kves


  ENDDO  ! end Kves

! sourceX and sourceY


!   IF(SourceType(Kves)(1:2)=='PR') THEN   


    DO J=Jbeg,Jend
    DO I=Ibeg,Iend

!   I modified the term to negative 11/22/2016 fyshi

       VesselPressureX(I,J) = -Grav*H(I,J)*  &
               (VesselPressureTotal(I+1,J)-VesselPressureTotal(I-1,J))/2.0_SP  &
               /DX
       VesselPressureY(I,J) = -Grav*H(I,J)*  &
               (VesselPressureTotal(I,J+1)-VesselPressureTotal(I,J-1))/2.0_SP  &
               /DY

    ENDDO
    ENDDO

!  make initial draft 09/12/2018
    IF(MakeVesselDraft)THEN
       MakeVesselDraft = .FALSE.
        DO J=1,Nloc
        DO I=1,Mloc
          IF(ABS(VesselPressureTotal(I,J)).GT.SMALL)THEN      
            Eta(I,J) = - VesselPressureTotal(I,J)
            ETA0(I,J) = ETA(I,J)
          ENDIF
        ENDDO
        ENDDO
    ENDIF
!   ENDIF  ! end pressure source type 1 and 2

END SUBROUTINE VESSEL_FORCING

SUBROUTINE PRESSURE_SOURCE (Xves,Yves)
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO, H,U,V,  &
                    DX,DY,Ibeg,Iend,Jbeg,Jend,DEPTH
  IMPLICIT NONE
  INTEGER :: ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves
    REAL(SP) :: myvar

! rectangular
    VesselPressureEach = ZERO
    ResistanceX(Kves)=ZERO
    ResistanceY(Kves)=ZERO
    ResPosX(Kves) = ZERO
    ResNegX(Kves) = ZERO
    ResPosY(Kves) = ZERO
    ResNegY(Kves) = ZERO

    DO J=1,Nloc
    DO I=1,Mloc
      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))

   ! Ertekin et al. JFM 1986

    IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN


     IF(VesselType(Kves) == 1) THEN
      P_x = ZERO
      P_y = ZERO

      IF(Lves>0.5_SP*LengthVessel(Kves)*AlphaVessel1(Kves).AND. &
         Lves<0.5_SP*LengthVessel(Kves))THEN
          
          P_x = COS(PI*(Lves-0.5_SP*AlphaVessel1(Kves)*LengthVessel(Kves))  &
                   /((1.0_SP-AlphaVessel1(Kves))*LengthVessel(Kves)))**2
      ELSEIF(Lves<-0.5_SP*LengthVessel(Kves)*AlphaVessel2(Kves).AND. &
         Lves>-0.5_SP*LengthVessel(Kves))THEN
          P_x = COS(PI*(ABS(Lves)-0.5_SP*AlphaVessel2(Kves)*LengthVessel(Kves))  &
                   /((1.0_SP-AlphaVessel2(Kves))*LengthVessel(Kves)))**2

      ELSEIF(Lves<=0.5_SP*LengthVessel(Kves)*AlphaVessel1(Kves).AND. &
         Lves>=-0.5_SP*LengthVessel(Kves)*AlphaVessel2(Kves))THEN
          P_x = 1.0_SP
      ENDIF

      IF(ABS(Wves)>0.5_SP*WidthVessel(Kves)*BetaVessel(Kves).AND. &
         ABS(Wves)<0.5_SP*WidthVessel(Kves)) THEN
          P_y = COS(PI*(ABS(Wves)-0.5_SP*BetaVessel(Kves)*WidthVessel(Kves))  &
                  /((1.0_SP-BetaVessel(Kves))*WidthVessel(Kves)))**2
      ELSEIF(ABS(Wves)<=0.5_SP*WidthVessel(Kves)*BetaVessel(Kves)) THEN
          P_y = 1.0_SP
      ENDIF

         VesselPressureEach(I,J) = Pvessel(Kves)*P_x*P_y
     ENDIF  ! type 1

     IF(VesselType(Kves) == 2) THEN

        VesselPressureEach(I,J) = Pvessel(Kves)*  &
              (1.0_SP - AlphaVessel1(Kves)*(Lves/LengthVessel(Kves))**4)* &
              (1.0_SP - AlphaVessel2(Kves)*(Wves/WidthVessel(Kves))**2)* &
              EXP(-BetaVessel(Kves)*(Wves/WidthVessel(Kves))**2)

     ENDIF ! end type 2     

    ENDIF  ! end inside ship rectangule
     
    ENDDO
    ENDDO  ! end grid

! calculate resistance. In the initial code, the resistance was calculated
! inside the last loop. It turns out that resistance was re-counted 
! because of including ghost cells (thanks to Jeff Harris)    

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend

         DetaX = (ETA(I+1,J)-ETA(I-1,J))/2.0_SP/DX
         DetaY = (ETA(I,J+1)-ETA(I,J-1))/2.0_SP/DY

         IF(DetaX>=0)THEN
           ResPosX(Kves) = ResPosX(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaX*DY
         ELSE
           ResNegX(Kves) = ResNegX(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaX*DY           
         ENDIF

         IF(DetaY>=0)THEN
           ResPosY(Kves) = ResPosY(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaY*DX
         ELSE
           ResNegY(Kves) = ResNegY(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaY*DX          
         ENDIF

         ResistanceX(Kves) = ResPosX(Kves) + ResNegX(Kves)
         ResistanceY(Kves) = ResPosY(Kves) + ResNegY(Kves)
      
    ENDDO
    ENDDO

  ! pressure source

     call MPI_ALLREDUCE(ResistanceX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResistanceX(Kves) = myvar
     call MPI_ALLREDUCE(ResistanceY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResistanceY(Kves) = myvar
     call MPI_ALLREDUCE(ResPosX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResPosX(Kves) = myvar
     call MPI_ALLREDUCE(ResPosY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResPosY(Kves) = myvar
     call MPI_ALLREDUCE(ResNegX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResNegX(Kves) = myvar
     call MPI_ALLREDUCE(ResNegY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResNegY(Kves) = myvar


END SUBROUTINE PRESSURE_SOURCE

SUBROUTINE GREEN_FUNCTION_SOURCE (Xves,Yves)
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO, H,U,V,  &
                    DX,DY,Ibeg,Iend,Jbeg,Jend,DEPTH
  IMPLICIT NONE
  INTEGER :: ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves


  VesselFluxGradientEach = ZERO

  

! assume single value of submerged vessel body
! we only record Gz at z_prime

    DO J=1,Nloc
    DO I=1,Mloc
      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))

      
      IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN
         VesselFluxGradientEach(I,J) = Pvessel(Kves)  &
                  *SIN(2.0*PI*Lves/(LengthVessel(Kves))) &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2   
!   print*,i,j,VesselFluxGradientEach(I,J), SIN(2.0*PI*Lves/(LengthVessel(Kves))), &
!          COS(PI*Wves/(WidthVessel(Kves)))**2
      ENDIF
    ENDDO
    ENDDO

  !  end realistic vessel body, otherwise slender

END SUBROUTINE GREEN_FUNCTION_SOURCE


SUBROUTINE SLENDER_BODY_SOURCE (Xves,Yves)
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO, H,U,V,  &
                    DX,DY,Ibeg,Iend,Jbeg,Jend,DEPTH
  IMPLICIT NONE
  INTEGER :: ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves
  REAL(SP) :: x_right_1, x_right_2, x_left_1, x_left_2


  VesselFluxGradientEach = ZERO

! initial parabolic function by Tanimoto et al. 2000
  
    DO J=1,Nloc
    DO I=1,Mloc
      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))

      
      IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN

        IF(VesselType(kves)==1)THEN
         VesselFluxGradientEach(I,J) = Pvessel(Kves)  &
                  *SIN(2.0*PI*Lves/(LengthVessel(Kves))) &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2   
        ENDIF  ! type 1

        IF(VesselType(kves)==2)THEN
         
         x_right_1=0.5_SP*LengthVessel(Kves)*AlphaVessel1(Kves)
         x_right_2=0.5_SP*LengthVessel(Kves)
         x_left_1=-0.5_SP*LengthVessel(Kves)*AlphaVessel2(Kves)
         x_left_2=-0.5_SP*LengthVessel(Kves)

         IF(Lves>x_right_1.AND. &
            Lves<x_right_2)THEN
          
          VesselFluxGradientEach(I,J) = Pvessel(Kves)  &
                  *SIN(PI*(Lves-x_right_1)/(x_right_2-x_right_1)) &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2 

         ELSEIF(Lves<x_left_1.AND. &
            Lves>x_left_2)THEN
          VesselFluxGradientEach(I,J) = Pvessel(Kves)  &
                  *SIN(PI*(Lves-x_left_1)/(x_left_1-x_left_2)) &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2 

         ELSEIF(Lves<=x_right_1.AND. &
            Lves>=x_left_1)THEN
          VesselFluxGradientEach(I,J) = ZERO
         ENDIF

        ENDIF  ! type 2

      ENDIF
    ENDDO
    ENDDO


!open(99,file='tmp.txt')
!do j=1,Nloc
! write(99,108)(VesselFluxGradientEach(i,j),i=1,Mloc)
!enddo
!close(99)
!108 format(5000E12.3)
!stop

END SUBROUTINE SLENDER_BODY_SOURCE

END MODULE VESSEL_MODULE

! end vessel
