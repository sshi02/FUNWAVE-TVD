!------------------------------------------------------------------------------------
!
!      FILE mod_sediment.F
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
!  SEDIMENT is a module to model sediment transport and morpho changes    
!
!  HISTORY :
!    03/14/2016 Babak Tehranirad
!    06/09/2017 Fengyan Shi
!    02/05/2018 put in the sediment module 
!-------------------------------------------------------------------------------------

# if defined (SEDIMENT)

MODULE SEDIMENT_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,U,V,P,Q,MinDepth,MASK,DT
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    CHARACTER(LEN=80) Mask_s_FILE
    CHARACTER(LEN=80) Sed_Scheme

    LOGICAL :: Bed_Change=.FALSE.
    LOGICAL :: Avalanche=.FALSE.
    LOGICAL :: IN_Mask_s=.FALSE.

    INTEGER :: Counter_s=0, Counter_ava=0
    INTEGER :: Morph_step,c_dum

    REAL(SP) :: D50,Sdensity,n_porosity,DT1,Counter,WS,Kappa1,Kappa2
    
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Zb,P_ave,D_ave,Mask_s,D,Pickup,Delta_c
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::CH,CHX,CHY,CHH,CHHxR,CHHxL,CHXxL,CHXxR
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::CHHyR,CHHyL,CHYyL,CHYyR,ZbOld,dZb
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::FsL,FsR,GsL,GsR,Fs,Gs,CH0,Zs,ava

    REAL(SP) :: Shields_cr,A_sed,Tau_cr,D_sed,htt,k_s,u_c,ustar_c,ec_max
    REAL(SP) :: beta_d,reduction,tau_xy,c_a,c_b,ec_max_w,ec_max_c
    REAL(SP) :: u_c1,u_c2,u_c3,u_c4,ustar_c1,ustar_c2,ustar_c3,ustar_c4
    REAL(SP) :: k1,k2,k3,k4,F1,F2,F3,F4,tan_phi,Dstar,viscosity

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif    


CONTAINS
  
SUBROUTINE SEDIMENT_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL, &
                    iista,jjsta   !ykchoi Jan/23/2018
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr

! sediment parameters 
      viscosity=0.000001_SP


! read vessel number and folder from input.txt
      FILE_NAME='input.txt'

! input is based on Babaks original but make default values
      CALL GET_STRING_VAL(Sed_Scheme,FILE_NAME,'Sed_Scheme',ierr)

      IF(ierr==1)THEN
        Sed_Scheme = 'Upwinding'
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
         WRITE(3,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
      endif
# else
         WRITE(*,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
         WRITE(3,'(A50)')'You dont specify Sed_Scheme, use Upwinding.'
# endif
      ENDIF

      CALL GET_Float_VAL(D50,FILE_NAME,'D50',ierr)

      IF(ierr==1)THEN
        D50 = 0.0005_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40)')'You dont specify D50, use default: 0.5mm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 0.5mm'
      endif
# else
         WRITE(*,'(A40)')'You dont specify D50, use default: 0.5mm'
         WRITE(3,'(A40)')'You dont specify D50, use default: 0.5mm'
# endif
      ENDIF 

      CALL GET_Float_VAL(Sdensity,FILE_NAME,'Sdensity',ierr)

      IF(ierr==1)THEN
        Sdensity = 2.68_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Sdensity, use default: 2.68'
         WRITE(3,'(A60)')'You dont specify Sdensity, use default: 2.68'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Sdensity, use default: 2.68'
         WRITE(3,'(A60)')'You dont specify Sdensity, use default: 2.68'
# endif
      ENDIF 

      CALL GET_LOGICAL_VAL(BED_CHANGE,FILE_NAME,'Bed_Change',ierr)

      IF(ierr == 1)THEN
       BED_CHANGE = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'Bed_Change not defined, Default: False'
       WRITE(*,'(A40)')'Bed_Change not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'Bed_Change not defined, Default: False'
      WRITE(*,'(A40)')'Bed_Change not defined, Default: False'
# endif 
      ENDIF

      CALL GET_Float_VAL(n_porosity,FILE_NAME,'n_porosity',ierr)

      IF(ierr==1)THEN
        n_porosity = 0.47_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify n_porosity, use default: 0.47'
         WRITE(3,'(A60)')'You dont specify n_porosity, use default: 0.47'
      endif
# else
         WRITE(*,'(A60)')'You dont specify n_porosity, use default: 0.47'
         WRITE(3,'(A60)')'You dont specify n_porosity, use default: 0.47'
# endif
      ENDIF 

      CALL GET_Float_VAL(WS,FILE_NAME,'WS',ierr)

      IF(ierr==1)THEN
        WS = 0.125_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify WS, use default: 0.125'
         WRITE(3,'(A60)')'You dont specify WS, use default: 0.125'
      endif
# else
         WRITE(*,'(A60)')'You dont specify WS, use default: 0.125'
         WRITE(3,'(A60)')'You dont specify WS, use default: 0.125'
# endif
      ENDIF

      CALL GET_Float_VAL(Shields_cr,FILE_NAME,'Shields_cr',ierr)

      IF(ierr==1)THEN
        Shields_cr = 0.055_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Shields_cr, use default: 0.055'
         WRITE(3,'(A60)')'You dont specify Shields_cr, use default: 0.055'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Shields_cr, use default: 0.055'
         WRITE(3,'(A60)')'You dont specify Shields_cr, use default: 0.055'
# endif
      ENDIF

      CALL GET_INTEGER_VAL(Morph_step,FILE_NAME,'Morph_step',ierr)

      IF(ierr==1)THEN
        Morph_step = 25
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Morph_step, use default: 25'
         WRITE(3,'(A60)')'You dont specify Morph_step, use default: 25'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Morph_step, use default: 25'
         WRITE(3,'(A60)')'You dont specify Morph_step, use default: 25'
# endif
      ENDIF

      CALL GET_LOGICAL_VAL(Avalanche,FILE_NAME,'Avalanche',ierr)
      IF(ierr == 1)THEN
       Avalanche = .FALSE. 
# if defined (PARALLEL) 
      if (myid.eq.0)then
       WRITE(3,'(A40)')'Avalanche not defined, Default: False'
       WRITE(*,'(A40)')'Avalanche not defined, Default: False'
      endif
# else
      WRITE(3,'(A40)')'Avalanche not defined, Default: False'
      WRITE(*,'(A40)')'Avalanche not defined, Default: False'
# endif 
      ENDIF

      CALL GET_Float_VAL(tan_phi,FILE_NAME,'Tan_phi',ierr)

      IF(ierr==1)THEN
        Tan_phi = 0.7_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Tan_phi, use default: 0.7'
         WRITE(3,'(A60)')'You dont specify Tan_phi, use default: 0.7'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Tan_phi, use default: 0.7'
         WRITE(3,'(A60)')'You dont specify Tan_phi, use default: 0.7'
# endif
      ENDIF

      CALL GET_LOGICAL_VAL(IN_Mask_s,FILE_NAME,'Hard_bottom',ierr)

      IF(ierr==1)THEN
        IN_Mask_s = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A60)')'You dont specify Hard_bottom, use default: False'
         WRITE(3,'(A60)')'You dont specify Hard_bottom, use default: False'
      endif
# else
         WRITE(*,'(A60)')'You dont specify Hard_bottom, use default: False'
         WRITE(3,'(A60)')'You dont specify Hard_bottom, use default: False'
# endif
      ENDIF

     IF(IN_Mask_s)THEN
      CALL GET_STRING_VAL(Mask_s_File,FILE_NAME,'Hard_bottom_file',ierr) 

      IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
         WRITE(3,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
      endif
       call MPI_FINALIZE ( ier )
# else
         WRITE(*,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
         WRITE(3,'(A80)')'You use Hard_bottom, Hard_bottom_file NOT FOUND, STOP'
# endif
        STOP
      ENDIF ! end ierr

     ENDIF ! end IN_Mask_s


      CALL GET_Float_VAL(Kappa1,FILE_NAME,'Kappa1',ierr)

      IF(ierr==1)THEN
        Kappa1 = 0.3333_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
         WRITE(3,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
      endif
# else
         WRITE(*,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
         WRITE(3,'(A80)')'You dont specify Kappa1, use default 4th-order: 0.3333'
# endif
      ENDIF

      CALL GET_Float_VAL(Kappa2,FILE_NAME,'Kappa2',ierr)
      IF(ierr==1)THEN
        Kappa2 = 1.0_SP
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
         WRITE(3,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
      endif
# else
         WRITE(*,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
         WRITE(3,'(A80)')'You dont specify Kappa2, use default 4th-order: 1.0'
# endif
      ENDIF

! ALLOCATION

        ALLOCATE(Zb(Mloc,Nloc))
        ALLOCATE(CHH(Mloc,Nloc))
        ALLOCATE(CHHxR(Mloc1,Nloc),CHHxL(Mloc1,Nloc),  &
                 CHHyR(Mloc,Nloc1),CHHyL(Mloc,Nloc1))
        ALLOCATE(CHY(Mloc,Nloc),CH(Mloc,Nloc), &
                 CHXxL(Mloc1,Nloc),CHXxR(Mloc1,Nloc))
        ALLOCATE(CHX(Mloc,Nloc),CHYyL(Mloc,Nloc1),CHYyR(Mloc,Nloc1)) 
        ALLOCATE(Mask_s(Mloc,Nloc))
        ALLOCATE(D(Mloc,Nloc),Pickup(Mloc,Nloc))
        ALLOCATE(Delta_c(Mloc,Nloc),D_ave(Mloc,Nloc),P_ave(Mloc,Nloc))
        ALLOCATE(FsL(Mloc1,Nloc),FsR(Mloc1,Nloc), &
                 GsL(Mloc,Nloc1),GsR(Mloc,Nloc1))
        ALLOCATE(Gs(Mloc,Nloc1),Fs(Mloc1,Nloc),CH0(Mloc,Nloc))
        ALLOCATE(ZbOld(Mloc,Nloc),dZb(Mloc,Nloc), &
                 Zs(Mloc,Nloc),ava(Mloc,Nloc))

        ava=ZERO
        Zs=ZERO		
        ZbOld=ZERO
        dZb=ZERO
        CH0=ZERO
        CHXxR=ZERO
        CHXxL=ZERO
        CH=ZERO
        CHX=ZERO
        CHXxR=ZERO
        CHXxL=ZERO
        CHYyR=ZERO
        CHYyL=ZERO   
        CHY=ZERO
        CHH=ZERO
        CHHxR=ZERO
        CHHxL=ZERO
        CHHyR=ZERO
        CHHyL=ZERO
        FsL=ZERO
        FsR=ZERO
        GsL=ZERO
        GsR=ZERO
        Fs=ZERO
        Gs=ZERO

        D=ZERO
        Pickup=ZERO
        D_ave=ZERO
        P_ave=ZERO
        Delta_c=ZERO
        Mask_s=ZERO
        Zb=ZERO		

End SUBROUTINE SEDIMENT_INITIAL


