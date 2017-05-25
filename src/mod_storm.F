!------------------------------------------------------------------------------------
!
!      FILE mod_storm.F
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
!  STORM is a module to model wind and pressure effects    
!
!  HISTORY :
!    05/16/2017  Fengyan Shi
!
!-------------------------------------------------------------------------------------

# if defined (STORM)

MODULE STORM_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,WindU2D,WindV2D
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    REAL(SP) :: Xstorm1,Ystorm1,Xstorm2,Ystorm2, &
                Pn_storm1, Pc_storm1,A_storm1,B_storm1, &
                Pn_storm2, Pc_storm2,A_storm2,B_storm2, &
                TimeStorm1,TimeStorm2,ThetaStorm
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: StormPressureTotal, &
                                           StormPressureX,StormPressureY
    REAL(SP),DIMENSION(:),ALLOCATABLE ::  Xco,Yco
    LOGICAL :: OUT_STORM = .TRUE.
    REAL(SP):: PLOT_INTV_STORM,PLOT_COUNT_STORM

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif    

CONTAINS
  
! READ STORM

SUBROUTINE STORM_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,SMALL,LARGE
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::STORM_FILE =' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: StormName

! read storm from input.txt
      FILE_NAME='input.txt'

! storm file
      CALL READ_STRING(STORM_FILE,FILE_NAME,'STORM_FILE',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# else
      WRITE(3,'(A15,A50)')'STORM_FILE:', STORM_FILE
# endif

      CALL READ_LOGICAL(OUT_STORM,FILE_NAME,'OUT_STORM',ierr)

      ALLOCATE (Xco(Mloc),Yco(Nloc))

      ALLOCATE (StormPressureTotal(Mloc,Nloc),&
                StormPressureX(Mloc,Nloc), &
                StormPressureY(Mloc,Nloc))

! plot intitial
     PLOT_COUNT_STORM = 0
     CALL READ_FLOAT(PLOT_INTV_STORM,FILE_NAME,'PLOT_INTV_STORM',ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A50)')'PLOT_INTV_STORM not specified, use LARGE'
# else
      WRITE(3,'(A50)')'PLOT_INTV_STORM not specified, use LARGE'
# endif
       PLOT_INTV_STORM = LARGE
     ENDIF
      
! Xco, and Yco

# if defined (PARALLEL)
     Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO

    TMP_NAME = TRIM(STORM_FILE)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), 'CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=300
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  StormName  !  name
         READ(Ifile,*)  ! t,x,y, Pn, Pc, A and B
         READ(Ifile,*)  TimeStorm2,Xstorm2,Ystorm2,  &
                        Pn_storm2,Pc_storm2, A_storm2, B_storm2

         TimeStorm1 = TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2
         Pn_storm1 = Pn_storm2
         Pc_storm1 = Pc_storm2
         A_storm1 = A_storm2
         B_storm1 = B_storm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   ENDIF
# else
   WRITE(*,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(*,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   WRITE(3,*) 'Storm Name: ',  TRIM(StormName)
   WRITE(3,*) 'Initial Time, X, Y', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
# endif

End SUBROUTINE STORM_INITIAL

SUBROUTINE STORM_FORCING
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Xstorm,Ystorm,Px,Py,DetaX,DetaY
  REAL(SP) :: Pn_storm, Pc_storm,A_storm,B_storm,Pw,Vw,Rdis,AngleLocal

  StormPressureTotal = ZERO

    IF(TIME>TimeStorm1.AND.TIME>TimeStorm2) THEN

         TimeStorm1=TimeStorm2
         Xstorm1 = Xstorm2
         Ystorm1 = Ystorm2

    Ifile = 300

    READ(Ifile,*,END=120)  TimeStorm2,Xstorm2,Ystorm2,  &
                        Pn_storm2,Pc_storm2, A_storm2, B_storm2

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   ENDIF
# else
   WRITE(*,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(*,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
   WRITE(3,*)'T,X,Y = ', TimeStorm2,Xstorm2,Ystorm2
   WRITE(3,*) 'Storm Pn,Pc,A, B: ',Pn_storm2,Pc_storm2, A_storm2, B_storm2
# endif

    ThetaStorm = ATAN2(Ystorm2-Ystorm1,  &
                              Xstorm2-Xstorm1)


    ENDIF ! end time > timestorm2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeStorm1)THEN
      IF(TimeStorm1.EQ.TimeStorm2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeStorm2-TIME) &
            /MAX(SMALL, ABS(TimeStorm2-TimeStorm1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xstorm = Xstorm2*tmp1 +Xstorm1*tmp2
    Ystorm = Ystorm2*tmp1 +Ystorm1*tmp2
    Pn_storm = Pn_storm2*tmp1 + Pn_storm1*tmp2
    Pc_storm = Pc_storm2*tmp1 + Pc_storm1*tmp2
    A_storm = A_storm2*tmp1 + A_storm1*tmp2
    B_storm = B_storm2*tmp1 + B_storm1*tmp2

120 CONTINUE  ! no more data for vessel Kves

! sourceX and sourceY


    DO J=1,Nloc
    DO I=1,Mloc
       Rdis=SQRT((Xco(I)-Xstorm)**2+(Yco(J)-Ystorm)**2)/1000.0_SP
       Rdis=MAX(SMALL,Rdis)
          ! Rdis is in km, Pw is in mb
       Pw=Pc_storm +(Pn_storm-Pc_storm)*EXP(-A_storm/Rdis**B_storm)

       StormPressureTotal(I,J) = Pw/100.0_SP  ! convert from cm to meter

       Vw=SQRT(A_storm*B_storm*100.0_SP*ABS(Pn_storm-Pc_storm)  &
           *EXP(-A_storm/Rdis**B_storm) / Rho_air/Rdis**B_storm)

       AngleLocal = ATAN2((Xco(I)-Xstorm),(Yco(J)-Ystorm))

       WindU2D(I,J)=-Vw*COS(AngleLocal)
       WindV2D(I,J)=Vw*SIN(AngleLocal)

    ENDDO
    ENDDO

#ifndef CARTESIAN
#ifdef STORM
#error "Please compile without -DSTORM when -DCARTESIAN is not set."
#endif
#endif

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
       StormPressureX(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I+1,J)-StormPressureTotal(I-1,J))/2.0_SP  &
               /DX
       StormPressureY(I,J) = -Grav*H(I,J)*  &
               (StormPressureTotal(I,J+1)-StormPressureTotal(I,J-1))/2.0_SP  &
               /DY
    ENDDO
    ENDDO

END SUBROUTINE STORM_FORCING

END MODULE STORM_MODULE

# endif 
! end storm
