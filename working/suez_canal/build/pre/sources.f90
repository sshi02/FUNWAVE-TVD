












!------------------------------------------------------------------------------------
!
!      FILE sources.F
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
!    SourceTerms is subroutine for all source terms
!
!    HISTORY: 
!       05/01/2010 Fengyan Shi
!       09/26/2013 Babak Tehranirad, added 2D Cd
!       08/18/2015 YoungKwang Choi, modified viscosity breaking
!       02/08/2016 Fengyan Shi, corrected wavemaker corresponding to 
!                               conservative form of momentum equations
!
! --------------------------------------------------
SUBROUTINE SourceTerms
     USE GLOBAL
     USE VESSEL_MODULE

     IMPLICIT NONE
     REAL,DIMENSION(Mloc,Nloc) :: nu_vis
     LOGICAL :: PQ_scheme = .FALSE.
     REAL(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
     INTEGER :: kd,kf


     DXg=DX
     DYg=DY


! wavemaker fyshi 02/08/2016

          DO J=1,Nloc
          DO I=1,Mloc

![---ykchoi Jan/23/2018
!            xmk=(I-Ibeg)*DXg+npx*(Mloc-2*Nghost)*DXg
!            ymk=(J-Jbeg)*DYg+npy*(Nloc-2*Nghost)*DYg
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
!---ykchoi Jan/23/2018]
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

       IF(WAVEMAKER(1:6)=='WK_REG')THEN          
          WaveMaker_Mass(I,J)=TANH(PI/(Time_ramp*Tperiod)*TIME)*D_gen &
                 *EXP(-Beta_gen*(xmk-Xc_WK)**2)&
                 *SIN(rlamda*(ymk-ZERO)-2.0_SP*PI/Tperiod*TIME)    
       ENDIF

!  added TMA_1D and JON_1D 02/24/2017

        IF(WaveMaker(1:6)=='WK_IRR'.OR.WaveMaker(1:6)=='TMA_1D'  &
           .OR.WaveMaker(1:6)=='JON_1D'.OR.WaveMaker(1:6)=='JON_2D')THEN
          WK_Source=ZERO
          DO kf=1,Nfreq
           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
                       *COS(OMGN_IR(KF)*TIME) &
                       +Sm(I,J,kf)*SIN(OMGN_IR(KF)*TIME))
          ENDDO
          WaveMaker_Mass(I,J)=WK_Source

       ENDIF

       IF(WAVEMAKER(1:7)=='WK_TIME')THEN
           WK_Source=ZERO
           DO kf=1,NumWaveComp
             WK_Source=WK_Source &
               +TANH(PI/(Time_ramp*PeakPeriod)*TIME)*D_genS(kf) &
                 *EXP(-Beta_genS(kf)*(xmk-Xc_WK)**2)&
                 *COS(2.0_SP*PI/WAVE_COMP(kf,1)*TIME-WAVE_COMP(kf,3)) 
           ENDDO
          WaveMaker_Mass(I,J)=WK_Source
       ENDIF

       IF(WAVEMAKER(1:9)=='WK_DATA2D')THEN

        ! make the efficient scheme using Cm and Sm 06/07/2016 fyshi
          WK_Source=ZERO

          DO kf=1,NumFreq
           WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm(I,J,kf) &
                       *COS(OMGN2D(KF)*TIME) &
                       +Sm(I,J,kf)*SIN(OMGN2D(KF)*TIME))
          ENDDO
          WaveMaker_Mass(I,J)=WK_Source
       ENDIF

      ENDIF ! end maker domain

         ENDDO
         ENDDO

! ----

     nu_vis=ZERO	
	! ykchoi(08.18.2015)
	! new variable :: WAVEMAKER_VIS
        ! should include ghost cells 02/28/2017 fyshi

      !IF(VISCOSITY_BREAKING)THEN
      IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS)THEN

        nu_vis = nu_break


!       DO J=Jbeg,Jend
!       DO I=Ibeg,Iend
!         nu_vis(I,J)=nu_break(I,J)
!       ENDDO
!       ENDDO

      ENDIF


     IF(DIFFUSION_SPONGE)THEN
!    should include ghost cells 02/28/2017, fyshi

       nu_vis = nu_vis + nu_sponge

!       DO J=Jbeg,Jend
!       DO I=Ibeg,Iend
!         nu_vis(I,J)=nu_vis(I,J)+nu_sponge(I,J)
!       ENDDO
!       ENDDO          

     ENDIF

107  format(500f12.6)

! depth gradient term
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend


! second order, move the second term to left-hand side
       SourceX(I,J)=GRAV*(Eta(I,J))/DX*(Depthx(I+1,J)-Depthx(I,J))*MASK(I,J) &
 ! friction
                   -Cd(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &

                       ! dispersion
                       ! (h+eta)(u*\nabla V4 + V4 * \nabla u - v1pp-v2-v3)
                   + Gamma1*MASK9(I,J)*Max(H(I,J),MinDepthFrc)*(         & 
                     U(I,J)*0.5_SP*(U4(I+1,J)-U4(I-1,J))/DX+V(I,J)*0.5_SP*(U4(I,J+1)-U4(I,J-1))/DY &
                     +U4(I,J)*0.5_SP*(U(I+1,J)-U(I-1,J))/DX+V4(I,J)*0.5_SP*(U(I,J+1)-U(I,J-1))/DY  &
                     -Gamma2*MASK9(I,J)*(U1pp(I,J)+U2(I,J)+U3(I,J)) &
                     )    &
                        ! Ht(-V4+V1p) = div(M)*(U4-U1p)
                    +Gamma1*MASK9(I,J)*((P(I+1,J)-P(I,J))/DX+(Q(I,J+1)-Q(I,J))/DY) &
                      *(U4(I,J)-U1p(I,J)) &
! wavemaker
                   +WaveMaker_Mass(I,J)*U(I,J)

       IF(FRICTION_SPONGE) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceX(I,J) = SourceX(I,J) &
                 - CD_4_SPONGE(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
                   *Depth(I,J)
       ENDIF


       IF(WaveMakerCurrentBalance)THEN
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
            SourceX(I,J) = SourceX(I,J) &
                  -WaveMakerCd*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
         ENDIF
       ENDIF ! current balance

       IF(BREAKWATER) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceX(I,J) = SourceX(I,J) &
                 - CD_breakwater(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
                   *Depth(I,J)
       ENDIF
          
       SourceY(I,J)=GRAV*(Eta(I,J))/DY*(Depthy(I,J+1)-Depthy(I,J))*MASK(I,J) &
                          ! friction
                   -Cd(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
		          ! dispersion
                          ! (h+eta)(u*\nabla V4 + V4 * \nabla u -v1pp-v2-v3)
                   + Gamma1*MASK9(I,J)*Max(H(I,J),MinDepthFrc)*(         & 
                     U(I,J)*0.5_SP*(V4(I+1,J)-V4(I-1,J))/DX+V(I,J)*0.5_SP*(V4(I,J+1)-V4(I,J-1))/DY &
                     +U4(I,J)*0.5_SP*(V(I+1,J)-V(I-1,J))/DX+V4(I,J)*0.5_SP*(V(I,J+1)-V(I,J-1))/DY  &
                     -Gamma2*MASK9(I,J)*(V1pp(I,J)+V2(I,J)+V3(I,J)) &
                     )    &
                          ! Ht(-V4+V1p) = div(Q)*(V4-V1p)
                    +Gamma1*MASK9(I,J)*((P(I+1,J)-P(I,J))/DX+(Q(I,J+1)-Q(I,J))/DY) &
                      *(V4(I,J)-V1p(I,J))  &
! wavemaker
                   +WaveMaker_Mass(I,J)*V(I,J)
       IF(FRICTION_SPONGE) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceY(I,J) = SourceY(I,J) &
                   -CD_4_SPONGE(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
                    *Depth(I,J)             
       ENDIF


       IF(WaveMakerCurrentBalance)THEN
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
            SourceY(I,J) = SourceY(I,J) &
                  -WaveMakerCd*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
         ENDIF
       ENDIF ! current balance


       IF(BREAKWATER) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceY(I,J) = SourceY(I,J) &
                   -CD_BREAKWATER(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
                    *Depth(I,J)             
       ENDIF




! eddy viscosity breaking
   ! ykchoi(08.18.2015)
   ! new variable :: WAVEMAKER_VIS
   !IF(VISCOSITY_BREAKING.OR.DIFFUSION_SPONGE)THEN
   IF(VISCOSITY_BREAKING.OR.DIFFUSION_SPONGE.OR.WAVEMAKER_VIS)THEN

     IF(PQ_scheme)THEN

!      it turns out P and Q are not exchanged at processor interface
!      it affects edges, make PQ_scheme=false

       SourceX(I,J) = SourceX(I,J) + 0.5_SP/DX*( &
                       nu_vis(I+1,J)* &
                       1.0_SP/DX*(P(I+2,J)-P(I+1,J)) &
                      -nu_vis(I-1,J)* &
                       1.0_SP/DX*(P(I,J)-P(I-1,J)) ) &
!
                                   + 1.0_SP/DY*( &
                       0.5_SP*(nu_vis(I,J+1)+nu_vis(I,J))* &                 
                       0.5_SP/DY*(P(I,J+1)+P(I+1,J+1)-P(I,J)-P(I+1,J)) &
                      -0.5_SP*(nu_vis(I,J-1)+nu_vis(I,J))* &
                       0.5_SP/DY*(P(I,J)+P(I+1,J)-P(I,J-1)-P(I+1,J-1)) )

       SourceY(I,J) = SourceY(I,J) + 0.5_SP/DY*( &
                       nu_vis(I,J+1)* &
                       1.0_SP/DY*(Q(I,J+2)-Q(I,J+1)) &
                      -nu_vis(I,J-1)* &
                       1.0_SP/DY*(Q(I,J)-Q(I,J-1)) ) &
!
                                   + 1.0_SP/DX*( &
                       0.5_SP*(nu_vis(I+1,J)+nu_vis(I,J))* &                 
                       0.5_SP/DX*(Q(I+1,J)+Q(I+1,J+1)-Q(I,J)-Q(I,J+1)) &
                      -0.5_SP*(nu_vis(I-1,J)+nu_vis(I,J))* &
                       0.5_SP/DX*(Q(I,J)+Q(I,J+1)-Q(I-1,J)-Q(I-1,J+1)) )

     ELSE
       SourceX(I,J) = SourceX(I,J) + 0.5_SP/DX*( &
                      (nu_vis(I+1,J)+nu_vis(I,J)) &
                      *1.0_SP/DX*(HU(I+1,J)-HU(I,J)) &
                     -(nu_vis(I-1,J)+nu_vis(I,J)) &
                      *1.0_SP/DX*(HU(I,J)-HU(I-1,J)) ) &
                                   + 0.5_SP/DY*( &
                      (nu_vis(I,J+1)+nu_vis(I,J)) &
                      *1.0_SP/DY*(HU(I,J+1)-HU(I,J)) &
                     -(nu_vis(I,J-1)+nu_vis(I,J)) &
                      *1.0_SP/DY*(HU(I,J)-HU(I,J-1)) )


                     

       SourceY(I,J) = SourceY(I,J) + 0.5_SP/DX*( &
                      (nu_vis(I+1,J)+nu_vis(I,J)) &
                      *1.0_SP/DX*(HV(I+1,J)-HV(I,J)) &
                     -(nu_vis(I-1,J)+nu_vis(I,J)) &
                      *1.0_SP/DX*(HV(I,J)-HV(I-1,J)) ) &
                                   + 0.5_SP/DY*( &
                      (nu_vis(I,J+1)+nu_vis(I,J)) &
                      *1.0_SP/DY*(HV(I,J+1)-HV(I,J)) &
                     -(nu_vis(I,J-1)+nu_vis(I,J)) &
                      *1.0_SP/DY*(HV(I,J)-HV(I,J-1)) )

     ENDIF ! end pq_scheme


    ENDIF  ! end eddy viscosity breaking

     ENDDO
     ENDDO

      SourceX = SourceX + VesselPressureX
      SourceY = SourceY + VesselPressureY


END SUBROUTINE SourceTerms

