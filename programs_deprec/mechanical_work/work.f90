program work

    use omp_lib
    use netcdf
    use libinterp, only: check
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    ! integer, parameter :: chunk_size = 384                                ! chunk size along y axis
    ! integer, parameter :: chunk_size = 24                                ! chunk size along y axis
    integer, parameter :: chunk_size = 1                                ! chunk size along y axis


    double precision, dimension(:,:,:,:), allocatable :: dz_buffer         ! layer geometric depth
    double precision, dimension(:,:,:,:), allocatable :: temp_buffer       ! temperature
    double precision, dimension(:,:,:,:), allocatable :: omega_buffer      ! omega
    double precision, dimension(:,:,:,:), allocatable :: qv_buffer         ! qv
    double precision, dimension(:,:,:,:), allocatable :: omt_buffer        ! omT
    double precision, dimension(:,:,:,:), allocatable :: omqv_buffer       ! omqv
    double precision, dimension(:,:,:), allocatable   :: work_out_buffer   ! mechanical work (output buffer)
    double precision, dimension(:,:,:), allocatable   :: work_th_out_buffer   ! thermal part of mechanical work (output buffer)
    double precision, dimension(:,:,:), allocatable   :: work_qv_out_buffer   ! vapor part of mechanical work (output buffer)
    double precision                                  :: geometric_thickness


    integer :: nx, ny, nplev, nt                                           ! note: nplev is full levels
    integer :: t, p, y, x, yc, nchunks, ystart
    integer :: tmp_unit, ios
    integer :: ncid_dz, ncid_temp, ncid_omega, ncid_qv, ncid_omt, ncid_omqv, ncid_work_out, ncid_work_th_out, ncid_work_qv_out
    integer :: varid_dz, varid_temp, varid_omega, varid_qv, varid_omt, varid_omqv, varid_work_out, varid_work_th_out, varid_work_qv_out

    character(len=255) :: path_dz, path_temp, path_omega, path_qv, path_omt, path_omqv, path_work_out, path_work_th_out, path_work_qv_out

    namelist /config/ path_dz, path_temp, path_omega, path_qv, path_omt, path_omqv, path_work_out, path_work_th_out, path_work_qv_out, nx, ny, nplev, nt



    ! get configuration namelist from command line and read it afterwards
    call get_command_argument(1, filenml)
    open(newunit=tmp_unit, file=trim(adjustl(filenml)), status='old', iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to open configuration namelist, iomsg='//trim(msg)
        stop
    endif
    read(unit=tmp_unit, nml=config, iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to assign namelist objects, iomsg='//trim(msg)
        stop
    endif
    close(tmp_unit)


    ! allocate buffers
    allocate(dz_buffer(nx, chunk_size, nplev, 1))
    allocate(temp_buffer(nx, chunk_size, nplev, 1))
    allocate(omega_buffer(nx, chunk_size, nplev, 1))
    allocate(qv_buffer(nx, chunk_size, nplev, 1))
    allocate(omt_buffer(nx, chunk_size, nplev, 1))
    allocate(omqv_buffer(nx, chunk_size, nplev, 1))
    allocate(work_out_buffer(nx, chunk_size, 1))
    allocate(work_th_out_buffer(nx, chunk_size, 1))
    allocate(work_qv_out_buffer(nx, chunk_size, 1))
    

    ! get variable IDs for layer geometric depth
    call check(nf90_open(trim(adjustl(path_dz)), nf90_write, ncid_dz))
    call check(nf90_inq_varid(ncid_dz, "DZ", varid_dz))

    ! get variable IDs for temperature
    call check(nf90_open(trim(adjustl(path_temp)), nf90_write, ncid_temp))
    call check(nf90_inq_varid(ncid_temp, "temp_coarse", varid_temp))

    ! get variable IDs for omega
    call check(nf90_open(trim(adjustl(path_omega)), nf90_write, ncid_omega))
    call check(nf90_inq_varid(ncid_omega, "ptend_coarse", varid_omega))

    ! get variable IDs for qv
    call check(nf90_open(trim(adjustl(path_qv)), nf90_write, ncid_qv))
    call check(nf90_inq_varid(ncid_qv, "sphum_coarse", varid_qv))

    ! get variable IDs for omt
    call check(nf90_open(trim(adjustl(path_omt)), nf90_write, ncid_omt))
    call check(nf90_inq_varid(ncid_omt, "omT_coarse", varid_omt))

    ! get variable IDs for omqv
    call check(nf90_open(trim(adjustl(path_omqv)), nf90_write, ncid_omqv))
    call check(nf90_inq_varid(ncid_omqv, "omqv_coarse", varid_omqv))

    ! get variable IDs for work_out (output variable)
    call check(nf90_open(trim(adjustl(path_work_out)), nf90_write, ncid_work_out))
    call check(nf90_inq_varid(ncid_work_out, "work", varid_work_out))

    ! get variable IDs for work_th_out (output variable)
    call check(nf90_open(trim(adjustl(path_work_th_out)), nf90_write, ncid_work_th_out))
    call check(nf90_inq_varid(ncid_work_th_out, "work_th", varid_work_th_out))

    ! get variable IDs for work_qv_out (output variable)
    call check(nf90_open(trim(adjustl(path_work_qv_out)), nf90_write, ncid_work_qv_out))
    call check(nf90_inq_varid(ncid_work_qv_out, "work_qv", varid_work_qv_out))

    nchunks = ny/chunk_size

    !$omp parallel private(t,yc,y,x,p,ystart)
    do t=1,nt
        !$omp master
        print *,'timestep',t,'of',nt
        !$omp end master
        do yc=1,nchunks
            !$omp master
            ystart = (yc-1)*chunk_size+1
            call check(nf90_get_var(ncid_dz, varid_dz, dz_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_temp, varid_temp, temp_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omega, varid_omega, omega_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qv, varid_qv, qv_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omt, varid_omt, omt_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqv, varid_omqv, omqv_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            !$omp end master
            !$omp barrier
            !$omp do
            do y=1,chunk_size
                do x=1,nx
                    work_out_buffer(x,y,1) = 0
                    work_th_out_buffer(x,y,1) = 0
                    work_qv_out_buffer(x,y,1) = 0
                    do p=1,nplev
                        geometric_thickness = -dz_buffer(x,y,p,1)
                        work_out_buffer(x,y,1) = work_out_buffer(x,y,1) - (omt_buffer(x,y,p,1) / temp_buffer(x,y,p,1) + 1.61 * (omqv_buffer(x,y,p,1) + omega_buffer(x,y,p,1) * qv_buffer(x,y,p,1))) * geometric_thickness
                        work_th_out_buffer(x,y,1) = work_th_out_buffer(x,y,1) - omt_buffer(x,y,p,1) / temp_buffer(x,y,p,1) * geometric_thickness
                        work_qv_out_buffer(x,y,1) = work_qv_out_buffer(x,y,1) - (1.61 * (omqv_buffer(x,y,p,1) + omega_buffer(x,y,p,1) * qv_buffer(x,y,p,1))) * geometric_thickness
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp barrier
            !$omp master
            call check(nf90_put_var(ncid_work_out, varid_work_out, work_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_work_th_out, varid_work_th_out, work_th_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_work_qv_out, varid_work_qv_out, work_qv_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            !$omp end master 
        enddo
    enddo
    !$omp end parallel

    call check(nf90_close(ncid_dz))
    call check(nf90_close(ncid_temp))
    call check(nf90_close(ncid_omega))
    call check(nf90_close(ncid_qv))
    call check(nf90_close(ncid_omt))
    call check(nf90_close(ncid_omqv))
    call check(nf90_close(ncid_work_out))
    call check(nf90_close(ncid_work_th_out))
    call check(nf90_close(ncid_work_qv_out))

end program

