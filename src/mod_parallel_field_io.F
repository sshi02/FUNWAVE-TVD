!------------------------------------------------------------------------------------
!
!      FILE mod_parallel_field_io.F
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
!    PARALLEL_FIELD_IO is a module that contains subroutines for outputting
!                      field data in various data format
!
!    HISTORY:
!      12/06/2017 Michael-Angelo Y.-H. Lam
!-------------------------------------------------------------------------------------

MODULE PARALLEL_FIELD_IO
USE GLOBAL

CONTAINS

SUBROUTINE PutFileASCII(FILE,PHI)
    
     IMPLICIT NONE
     CHARACTER(LEN=80), INTENT(IN) :: FILE
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI

![-------ykchoi (06/May/2017)
      !INTEGER :: l
      ! could be max. procs
      !INTEGER,DIMENSION(NumberProcessor) :: npxs,npys
	!REAL(SP),DIMENSION(NumberProcessor) :: xx
     INTEGER :: irank, lenx, leny, lenxy, ireq
     INTEGER, ALLOCATABLE :: iistas(:), iiends(:), jjstas(:), jjends(:)
     INTEGER :: istatus(mpi_status_size)

     REAL(SP), ALLOCATABLE :: LocalPHI(:,:)
     REAL(SP), ALLOCATABLE :: xx(:,:)
     
      !REAL(SP),DIMENSION(MGlob+2*Nghost,NGlob+2*Nghost) :: PHIGLOB   
     REAL(SP),DIMENSION( MGlob, NGlob ) :: PHIGLOB
!-------ykchoi (06/May/2017)]
     LOGICAL :: FirstCallPutFile = .TRUE.
     SAVE  FirstCallPutFile
     
! first time call 
     IF(FirstCallPutFile)THEN
        FirstCallPutFile = .FALSE.
! format length
        write(FORMAT_LEN(1:1),'(A1)') '('
        write(FORMAT_LEN(2:8),'(I7)') Mglob
        write(FORMAT_LEN(9:13),'(A5)') 'E16.6'
        write(FORMAT_LEN(14:14),'(A1)') ')'
     ENDIF

![-------ykchoi (06/May/2017)
      !call MPI_Gather(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,&
      !    0,MPI_COMM_WORLD,ier)
      !call MPI_Gather(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,&
      !    0,MPI_COMM_WORLD,ier)
     ALLOCATE( LocalPHI(Mloc-2*Nghost, Nloc-2*Nghost) )
     LocalPHI=PHI(Ibeg:Iend,Jbeg:Jend)

     allocate( iistas(nprocs), iiends(nprocs), jjstas(nprocs), jjends(nprocs) ) 

     call MPI_Gather( iista, 1, MPI_INTEGER, iistas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( iiend, 1, MPI_INTEGER, iiends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( jjsta, 1, MPI_INTEGER, jjstas, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
     call MPI_Gather( jjend, 1, MPI_INTEGER, jjends, 1, MPI_INTEGER, &
                      0, MPI_COMM_WORLD, ier )
      
     if( myid == 0 ) then
	 PHIGLOB( iista:iiend, jjsta:jjend ) = LocalPHI
     endif

     do irank=1, px*py-1  	  
	  
	  if( myid == 0 ) then
	    lenx = iiends(irank+1) - iistas(irank+1) + 1 
	    leny = jjends(irank+1) - jjstas(irank+1) + 1
	    lenxy = lenx*leny
	    allocate( xx(lenx, leny) )
	  
	    call mpi_irecv( xx, lenxy, mpi_sp, irank, 1, mpi_comm_world, ireq, ier )
	    call mpi_wait( ireq, istatus, ier )
	    
		PHIGLOB( iistas(irank+1):iiends(irank+1), jjstas(irank+1):jjends(irank+1) ) = xx
		deallocate( xx )

        elseif( myid == irank ) then
	    
		lenxy = ( iiend-iista+1 )*( jjend-jjsta+1 )
	    call mpi_isend( LocalPHI, lenxy, mpi_sp, 0, 1, mpi_comm_world, ireq, ier )
	    call mpi_wait( ireq, istatus, ier )

        endif	  

     enddo

     deallocate(LocalPHI, iistas, iiends, jjstas, jjends)

     if (myid.eq.0) then
        OPEN(1,FILE=TRIM(FILE))
	  DO J=1,NGlob,OUTPUT_RES
	     WRITE(1,FORMAT_LEN)(real(PHIGLOB(I,J)),I=1,MGlob,OUTPUT_RES)
	  ENDDO
        !DO J=Nghost+1,NGlob+NGhost,OUTPUT_RES
        !   WRITE(1,FORMAT_LEN)(real(PHIGLOB(I,J)),I=Nghost+1,MGlob+Nghost,OUTPUT_RES)
        !ENDDO
!100  FORMAT(5000E16.6)
!100   FORMAT(FORMAT_LEN)
        CLOSE(1)
     endif
!-------ykchoi (06/May/2017)]

END SUBROUTINE PutfileASCII

END MODULE PARALLEL_FIELD_IO
