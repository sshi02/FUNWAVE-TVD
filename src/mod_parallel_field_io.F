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
!      10/10/2019 Zhouteng Ye, add HDF5 support
!-------------------------------------------------------------------------------------

# if defined (PARALLEL)

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

! Code derived from William D. Gropp lecture slides 33 
! http://wgropp.cs.illinois.edu/courses/cs598-s16/lectures/lecture33.pdf
SUBROUTINE PutFileBinary(FILE_NAME,PHI)
    
     IMPLICIT NONE
     CHARACTER(LEN=80), INTENT(IN) :: FILE_NAME
     REAL(SP),DIMENSION(Mloc,Nloc), INTENT(IN) :: PHI
     REAL(SP),DIMENSION(Mloc,Nloc) :: PHITEST

     INTEGER ::  file_handler, ierr , istatus  , mem_type_loc , file_type_glob
     INTEGER(8) :: offset
     INTEGER , DIMENSION(ndims) ::  sizes_loc , subsizes_loc , start_index_loc
     INTEGER , DIMENSION(ndims) ::  sizes_glob , subsizes_glob , start_index_glob

     offset = 0


     ! Defines each nodes subdomain 'PHI' and the surrounding ghost points
     ! Used to extract data minus ghost points
     sizes_loc(1) = Mloc 
     sizes_loc(2) = Nloc

     subsizes_loc(1) =  sizes_loc(1) - 2*Nghost
     subsizes_loc(2) =  sizes_loc(2) - 2*Nghost

     start_index_loc(1) = Nghost  ! 0 start index even in Fortran
     start_index_loc(2) = Nghost  ! 0 start index even in Fortran

     ! Defines the mapping between subdomains (minus ghost points)
     ! to the global domain
     sizes_glob(1) = Mglob
     sizes_glob(2) = Nglob

     subsizes_glob(1) =   subsizes_loc(1) 
     subsizes_glob(2) =   subsizes_loc(2)

     start_index_glob(1) = iista - 1 ! 0 start index even in Fortran 
     start_index_glob(2) = jjsta - 1 ! 0 start index even in Fortran


     CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes_loc , subsizes_loc , &
                                   start_index_loc , MPI_ORDER_FORTRAN, mpi_sp , mem_type_loc, ierr)
     CALL MPI_TYPE_COMMIT(mem_type_loc, ierr)


     CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes_glob , subsizes_glob , &
                                   start_index_glob , MPI_ORDER_FORTRAN, mpi_sp , file_type_glob , ierr)
     CALL MPI_TYPE_COMMIT(file_type_glob, ierr)


     CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FILE_NAME , & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, file_handler , ierr) 

     CALL MPI_FILE_SET_VIEW( file_handler , offset , mpi_sp , file_type_glob , 'native' , MPI_INFO_NULL, ierr)


     CALL MPI_FILE_WRITE_ALL(file_handler, PHI , 1, mem_type_loc , MPI_STATUS_IGNORE, ierr)

     CALL MPI_FILE_CLOSE(file_handler, ierr )


END SUBROUTINE PutfileBinary

#if defined (HDF5IO)
Subroutine H5Init(h5file)
  use hdf5
  use global
  Implicit None
  INTEGER(HID_T) :: plist_id      ! Property list identifier 
  Type(hdf5file) :: h5file

  Integer :: h5error

  h5_total_dims(1) = Mglob
  h5_total_dims(2) = Nglob
  h5_block_dims(1) = Mloc - 2 * Nghost
  h5_block_dims(2) = Nloc - 2 * Nghost


  h5_offset(1) = npx*int(mglob/px) + min(npx, mod( mglob, px ))
  h5_offset(2) = npy*int(nglob/py) + min(npy, mod( nglob, py ))


  if (myid .eq. 0 ) Call system('rm -f'//trim(h5file%filename))

  Call MPI_barrier(MPI_COMM_WORLD, h5error)

  Call h5open_f(h5error)
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5error)
  CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5error)

  ! Create the file collectively.
  CALL h5fcreate_f(h5file%filename, H5F_ACC_TRUNC_F, h5file%file_id, h5error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, h5error)


  IF(OUT_DEPTH.OR.BREAKWATER)THEN
    h5group_dep%groupname   = 'dep'
    Call HDF5CreateGroup(h5file, h5group_dep)
    h5group_breakwater%groupname   = 'breakwater'
    Call HDF5CreateGroup(h5file, h5group_breakwater)
  ENDIF

  IF(OUT_ETA)THEN
    h5group_eta%groupname = 'eta'
    Call HDF5CreateGroup(h5file, h5group_eta)
  end IF

  IF(OUT_Hmax)THEN
    h5group_hmax%groupname = 'Hmax'
    Call HDF5CreateGroup(h5file, h5group_hmax)
  ENDIF

  IF(OUT_Hmin)THEN
    h5group_hmin%groupname = 'Hmin'
    Call HDF5CreateGroup(h5file, h5group_hmin)
  ENDIF

  IF(OUT_Umax)THEN
    h5group_Umax%groupname = 'Umax'
    call HDF5CreateGroup(h5file, h5group_umax)
  ENDIF

  IF(OUT_MFmax)THEN
    h5group_MFmax%groupname = 'MFmax'
    call HDF5CreateGroup(h5file, h5group_MFmax)
  ENDIF

  IF(OUT_VORmax)THEN
    h5group_VORmax%groupname = 'VORmax'
    call HDF5CreateGroup(h5file, h5group_VORmax)
  ENDIF

  IF(OUT_MASK)THEN
    h5group_mask%groupname = 'mask'
    call HDF5CreateGroup(h5file, h5group_mask)
  ENDIF

  IF(OUT_MASK9)THEN
    h5group_mask9%groupname = 'mask9'
    call HDF5CreateGroup(h5file, h5group_mask9)
  ENDIF

  IF(OUT_U)THEN
    h5group_u%groupname   = 'U'
    Call HDF5CreateGroup(h5file, h5group_u)
  EndIf
  IF(OUT_V)THEN
    h5group_v%groupname   = 'V'
    Call HDF5CreateGroup(h5file, h5group_v)
  end IF
  IF(OUT_P)THEN
    h5group_P%groupname   = 'P'
    Call HDF5CreateGroup(h5file, h5group_P)
  end IF
  IF(OUT_Q)THEN
    h5group_Q%groupname   = 'Q'
    Call HDF5CreateGroup(h5file, h5group_Q)
  end IF

  
  IF(OUT_AGE)THEN
    IF(SHOW_BREAKING)THEN
      h5group_age%groupname   = 'age'
      call HDF5CreateGroup(h5file, h5group_age)
    ENDIF
  ENDIF

  IF(OUT_ROLLER)THEN
    h5group_roller%groupname   = 'roller'
    call HDF5CreateGroup(h5file, h5group_roller)
  ENDIF

  IF(OUT_UNDERTOW)THEN
    h5group_undertow_u%groupname   = 'undertow_u'
    h5group_undertow_v%groupname   = 'undertow_v'
    call HDF5CreateGroup(h5file, h5group_undertow_u)
    call HDF5CreateGroup(h5file, h5group_undertow_v)
  ENDIF

  IF(VISCOSITY_BREAKING)THEN
    IF(OUT_NU)THEN
      h5group_nubrk%groupname   = 'nu_break'
      call HDF5CreateGroup(h5file, h5group_nubrk)
    ENDIF
  ENDIF

# if defined (METEO)
   IF(OUT_METEO)THEN

    IF(WindHollandModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      h5group_Ustorm%groupname   = 'U_storm'
      h5group_Vstorm%groupname   = 'V_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
      call HDF5CreateGroup(h5file, h5group_Ustorm)
      call HDF5CreateGroup(h5file, h5group_Vstorm)
    ENDIF

    IF(MeteoGausian)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
    ENDIF

    IF(SlideModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CreateGroup(h5file, h5group_Pstorm)
    ENDIF

  ENDIF
# endif
     
# if defined (VESSEL)
     IF(OUT_VESSEL)THEN
# if defined (VESSEL_PANEL_SOURCE)
       h5group_Fves%groupname   = 'Fves'
       call HDF5CreateGroup(h5file, h5group_Fves)
# else
       h5group_Pves%groupname   = 'Pves'
       call HDF5CreateGroup(h5file, h5group_Pves)
# endif
   ! end panel source

     ENDIF
# endif

! sediment
# if defined (SEDIMENT)
     h5group_c%groupname   = 'c'
     h5group_pick%groupname   = 'pick'
     h5group_dego%groupname   = 'depo'
     h5group_Dchgs%groupname   = 'Dchgs'
     h5group_DchgB%groupname   = 'DchgB'
     h5group_BedFx%groupname   = 'Bed_Fx'
     h5group_BedFy%groupname   = 'Bed_Fy'
     call HDF5CreateGroup(h5file, h5group_c)
     call HDF5CreateGroup(h5file, h5group_pick)
     call HDF5CreateGroup(h5file, h5group_depo)
     call HDF5CreateGroup(h5file, h5group_Dchgs)
     call HDF5CreateGroup(h5file, h5group_DchgB)
     call HDF5CreateGroup(h5file, h5group_BedFx)
     call HDF5CreateGroup(h5file, h5group_BedFy)
# endif


End Subroutine H5Init

Subroutine H5Finalize(h5file)
  use hdf5
  use global
  Implicit None
  Integer :: error
  Type(hdf5file) :: h5file

  Call HDF5_WRITE_PARAMETERS(h5file)

  IF(OUT_DEPTH.OR.BREAKWATER)THEN
    h5group_dep%groupname   = 'dep'
    Call HDF5CloseGroup( h5group_dep)
    h5group_breakwater%groupname   = 'breakwater'
    Call HDF5CloseGroup( h5group_breakwater)
  ENDIF

  IF(OUT_ETA)THEN
    h5group_eta%groupname = 'eta'
    Call HDF5CloseGroup( h5group_eta)
  end IF

  IF(OUT_Hmax)THEN
    h5group_hmax%groupname = 'Hmax'
    Call HDF5CloseGroup( h5group_hmax)
  ENDIF

  IF(OUT_Hmin)THEN
    h5group_hmin%groupname = 'Hmin'
    Call HDF5CloseGroup( h5group_hmin)
  ENDIF

  IF(OUT_Umax)THEN
    call HDF5CloseGroup( h5group_umax)
  ENDIF

  IF(OUT_MFmax)THEN
    call HDF5CloseGroup( h5group_MFmax)
  ENDIF

  IF(OUT_VORmax)THEN
    call HDF5CloseGroup( h5group_VORmax)
  ENDIF

  IF(OUT_MASK)THEN
    call HDF5CloseGroup( h5group_mask)
  ENDIF

  IF(OUT_MASK9)THEN
    call HDF5CloseGroup( h5group_mask9)
  ENDIF

  IF(OUT_U)THEN
    h5group_u%groupname   = 'U'
    Call HDF5CloseGroup( h5group_u)
  EndIf
  IF(OUT_V)THEN
    h5group_v%groupname   = 'V'
    Call HDF5CloseGroup( h5group_v)
  end IF
  IF(OUT_P)THEN
    h5group_P%groupname   = 'P'
    Call HDF5CloseGroup( h5group_P)
  end IF
  IF(OUT_Q)THEN
    h5group_Q%groupname   = 'Q'
    Call HDF5CloseGroup( h5group_Q)
  end IF

  
  IF(OUT_AGE)THEN
    IF(SHOW_BREAKING)THEN
      h5group_age%groupname   = 'age'
      call HDF5CloseGroup( h5group_age)
    ENDIF
  ENDIF

  IF(OUT_ROLLER)THEN
    h5group_roller%groupname   = 'roller'
    call HDF5CloseGroup( h5group_roller)
  ENDIF

  IF(OUT_UNDERTOW)THEN
    h5group_undertow_u%groupname   = 'undertow_u'
    h5group_undertow_v%groupname   = 'undertow_v'
    call HDF5CloseGroup( h5group_undertow_u)
    call HDF5CloseGroup( h5group_undertow_v)
  ENDIF

  IF(VISCOSITY_BREAKING)THEN
    IF(OUT_NU)THEN
      h5group_nubrk%groupname   = 'nu_break'
      call HDF5CloseGroup( h5group_nubrk)
    ENDIF
  ENDIF

# if defined (METEO)
   IF(OUT_METEO)THEN

    IF(WindHollandModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      h5group_Ustorm%groupname   = 'U_storm'
      h5group_Vstorm%groupname   = 'V_storm'
      call HDF5CloseGroup( h5group_Pstorm)
      call HDF5CloseGroup( h5group_Ustorm)
      call HDF5CloseGroup( h5group_Vstorm)
    ENDIF

    IF(MeteoGausian)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CloseGroup( h5group_Pstorm)
    ENDIF

    IF(SlideModel)THEN
      h5group_Pstorm%groupname   = 'P_storm'
      call HDF5CloseGroup( h5group_Pstorm)
    ENDIF

  ENDIF
# endif
     
# if defined (VESSEL)
     IF(OUT_VESSEL)THEN
# if defined (VESSEL_PANEL_SOURCE)
       h5group_Fves%groupname   = 'Fves'
       call HDF5CloseGroup( h5group_Fves)
# else
       h5group_Pves%groupname   = 'Pves'
       call HDF5CloseGroup( h5group_Pves)
# endif
   ! end panel source

     ENDIF
# endif

! sediment
# if defined (SEDIMENT)
     h5group_c%groupname   = 'c'
     h5group_pick%groupname   = 'pick'
     h5group_dego%groupname   = 'depo'
     h5group_Dchgs%groupname   = 'Dchgs'
     h5group_DchgB%groupname   = 'DchgB'
     h5group_BedFx%groupname   = 'Bed_Fx'
     h5group_BedFy%groupname   = 'Bed_Fy'
     call HDF5CloseGroup( h5group_c)
     call HDF5CloseGroup( h5group_pick)
     call HDF5CloseGroup( h5group_depo)
     call HDF5CloseGroup( h5group_Dchgs)
     call HDF5CloseGroup( h5group_DchgB)
     call HDF5CloseGroup( h5group_BedFx)
     call HDF5CloseGroup( h5group_BedFy)
# endif



  Call h5fclose_f(h5file%file_id, error)

End Subroutine H5Finalize

Subroutine HDF5CreateGroup(h5file, h5group)
  use global
  Implicit None

  Type(hdf5file) :: h5file
  Type(hdf5group) :: h5group

  Integer :: h5error

  CALL h5gcreate_f(h5file%file_id, h5group%groupname, h5group%group_id, h5error)

End Subroutine HDF5CreateGroup

Subroutine HDF5CloseGroup(h5group)
  use global
  Implicit None

  Type(hdf5group) :: h5group

  Integer :: h5error

  CALL h5gclose_f(h5group%group_id, h5error)

End Subroutine HDF5CloseGroup

SUBROUTINE HDF5PutFile(h5group, data, data_name)
  USE GLOBAL
  IMPLICIT NONE

  Type(hdf5group) :: h5group
  Real(SP), Intent(In) :: data(Mloc,Nloc)
  Character(80) :: data_name

  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier 
  INTEGER(HID_T) :: dset_id       ! Data set identifier 

  Integer :: mblock, nblock

  integer :: rank = 2

  Integer :: h5error


  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(rank, h5_total_dims, filespace, h5error)

  ! Create the dataset with default properties.
# if defined DOUBLE_PRECISION
  CALL h5dcreate_f(h5group%group_id, data_name, H5T_NATIVE_DOUBLE, filespace, &
      dset_id, h5error)
#else
  CALL h5dcreate_f(h5group%group_id, data_name, H5T_NATIVE_REAL, filespace, &
      dset_id, h5error)
#endif
  CALL h5sclose_f(filespace, h5error)

  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(rank, h5_block_dims, memspace, h5error)


  ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, h5error)
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, h5_offset, h5_block_dims, h5error)

    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5error)


# if defined DOUBLE_PRECISION
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data(ibeg:iend,jbeg:jend), h5_total_dims ,error, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)
# else
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data(ibeg:iend,jbeg:jend), h5_total_dims ,h5error, &
          file_space_id = filespace, mem_space_id = memspace, xfer_prp=plist_id)
# endif


    ! Close dataspaces.
    CALL h5sclose_f(filespace, h5error)
    CALL h5sclose_f(memspace, h5error)
    
    ! Close the dataset.
    CALL h5dclose_f(dset_id, h5error)

    Call h5pclose_f(plist_id, h5error)


end SUBROUTINE HDF5PutFile

SUBROUTINE HDF5_WRITE_PARAMETERS(h5file)
  USE GLOBAL
  Implicit None
  Type(hdf5file) :: h5file
  Character(80) :: tmp_name 


  ! write the total grid information
  tmp_name = 'nx'
  Call add_Parameter_integer(h5file%file_id, tmp_name, mglob, 1)
  tmp_name = 'ny'
  Call add_Parameter_integer(h5file%file_id, tmp_name, nglob, 1)
  ! write the grid resolution
  tmp_name = 'dx'
  Call add_Parameter_float(h5file%file_id, tmp_name, dx, 1)
  tmp_name = 'dy'
  Call add_Parameter_float(h5file%file_id, tmp_name, dy, 1)
  ! write the count of data
  tmp_name = 'out_interval'
  Call add_Parameter_float(h5file%file_id, tmp_name, PLOT_INTV,1)
  tmp_name = 'out_count'
  Call add_Parameter_integer(h5file%file_id, tmp_name, icount, 1)

END SUBROUTINE HDF5_WRITE_PARAMETERS

Subroutine add_Parameter_character(id, p_name, p_value)
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Character(80) :: p_name
  Character(80) :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  Integer(Size_T) :: attrlen=80
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  CALL h5tset_size_f(atype_id, attrlen, error)
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_character

Subroutine add_Parameter_Integer(id, p_name, p_value, dim)
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Integer :: dim
  Character(80) :: p_name
  integer :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  adims = dim
  data_dims = dim

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_Integer

Subroutine add_Parameter_Float(id, p_name, p_value, dim)
  use param
  use HDF5
  Implicit None
  Integer(HID_T) :: id
  Integer :: dim
  Character(80) :: p_name
  Real(SP) :: p_value
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims = (/1/)
  INTEGER     ::   arank = 1                      ! Attribure rank

  Integer :: error

  adims = dim
  data_dims = dim

  CALL h5screate_simple_f(arank, adims, aspace_id, error)
# if defined (DOUBLE_PRECISION)
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
# else
  CALL h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
# endif
  CALL h5acreate_f(id, p_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, p_value, data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5tclose_f(atype_id, error)
  CALL h5sclose_f(aspace_id, error)

End Subroutine add_Parameter_Float




#endif

END MODULE PARALLEL_FIELD_IO


# endif
 ! end parallel
