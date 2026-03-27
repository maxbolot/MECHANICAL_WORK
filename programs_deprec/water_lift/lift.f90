program lift

    use omp_lib
    use netcdf
    use libinterp, only: check
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    ! integer, parameter :: chunk_size = 384                                ! chunk size along y axis
    integer, parameter :: chunk_size = 24                                ! chunk size along y axis
    ! integer, parameter :: chunk_size = 1                                ! chunk size along y axis


    double precision, dimension(:,:,:,:), allocatable :: dz_buffer         ! layer geometric depth
    double precision, dimension(:,:,:,:), allocatable :: omega_buffer      ! omega
    double precision, dimension(:,:,:,:), allocatable :: qv_buffer      ! qv
    double precision, dimension(:,:,:,:), allocatable :: qw_buffer      ! qw
    double precision, dimension(:,:,:,:), allocatable :: qr_buffer      ! qr
    double precision, dimension(:,:,:,:), allocatable :: qi_buffer      ! qi
    double precision, dimension(:,:,:,:), allocatable :: qs_buffer      ! qs
    double precision, dimension(:,:,:,:), allocatable :: qg_buffer      ! qg
    double precision, dimension(:,:,:,:), allocatable :: omqv_buffer       ! omqv
    double precision, dimension(:,:,:,:), allocatable :: omqw_buffer       ! omqw
    double precision, dimension(:,:,:,:), allocatable :: omqr_buffer       ! omqr
    double precision, dimension(:,:,:,:), allocatable :: omqi_buffer       ! omqi
    double precision, dimension(:,:,:,:), allocatable :: omqs_buffer       ! omqs
    double precision, dimension(:,:,:,:), allocatable :: omqg_buffer       ! omqg
    double precision, dimension(:,:,:), allocatable   :: lift_out_buffer   ! mechanical work to lift water (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qv_out_buffer   ! mechanical work to lift vapor (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qw_out_buffer   ! mechanical work to lift cloud liquid water (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qr_out_buffer   ! mechanical work to lift rain (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qi_out_buffer   ! mechanical work to cloud ice (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qs_out_buffer   ! mechanical work to snow (output buffer)
    double precision, dimension(:,:,:), allocatable   :: lift_qg_out_buffer   ! mechanical work to graupel (output buffer)
    double precision                                  :: geometric_thickness


    integer :: nx, ny, nplev, nt                                           ! note: nplev is full levels
    integer :: t, p, y, x, yc, nchunks, ystart
    integer :: tmp_unit, ios
    integer :: ncid_dz, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg
    integer :: ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg
    integer :: ncid_lift_out, ncid_lift_qv_out, ncid_lift_qw_out, ncid_lift_qr_out, ncid_lift_qi_out, ncid_lift_qs_out, ncid_lift_qg_out

    integer :: varid_dz, varid_omega, varid_qv, varid_qw, varid_qr, varid_qi, varid_qs, varid_qg
    integer :: varid_omqv, varid_omqw, varid_omqr, varid_omqi, varid_omqs, varid_omqg 
    integer :: varid_lift_out, varid_lift_qv_out, varid_lift_qw_out, varid_lift_qr_out, varid_lift_qi_out, varid_lift_qs_out, varid_lift_qg_out

    character(len=255) :: path_dz, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg
    character(len=255) :: path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg
    character(len=255) :: path_lift_out, path_lift_qv_out, path_lift_qw_out, path_lift_qr_out, path_lift_qi_out, path_lift_qs_out, path_lift_qg_out

    namelist /config/ path_dz, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg, &
        & path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, &
        & path_lift_out, path_lift_qv_out, path_lift_qw_out, path_lift_qr_out, path_lift_qi_out, path_lift_qs_out, path_lift_qg_out, &
        & nx, ny, nplev, nt



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
    allocate(omega_buffer(nx, chunk_size, nplev, 1))
    allocate(qv_buffer(nx, chunk_size, nplev, 1))
    allocate(qw_buffer(nx, chunk_size, nplev, 1))
    allocate(qr_buffer(nx, chunk_size, nplev, 1))
    allocate(qi_buffer(nx, chunk_size, nplev, 1))
    allocate(qs_buffer(nx, chunk_size, nplev, 1))
    allocate(qg_buffer(nx, chunk_size, nplev, 1))
    allocate(omqv_buffer(nx, chunk_size, nplev, 1))
    allocate(omqw_buffer(nx, chunk_size, nplev, 1))
    allocate(omqr_buffer(nx, chunk_size, nplev, 1))
    allocate(omqi_buffer(nx, chunk_size, nplev, 1))
    allocate(omqs_buffer(nx, chunk_size, nplev, 1))
    allocate(omqg_buffer(nx, chunk_size, nplev, 1))
    allocate(lift_out_buffer(nx, chunk_size, 1))
    allocate(lift_qv_out_buffer(nx, chunk_size, 1))
    allocate(lift_qw_out_buffer(nx, chunk_size, 1))
    allocate(lift_qr_out_buffer(nx, chunk_size, 1))
    allocate(lift_qi_out_buffer(nx, chunk_size, 1))
    allocate(lift_qs_out_buffer(nx, chunk_size, 1))
    allocate(lift_qg_out_buffer(nx, chunk_size, 1))
    

    ! get variable IDs for layer geometric depth
    call check(nf90_open(trim(adjustl(path_dz)), nf90_write, ncid_dz))
    call check(nf90_inq_varid(ncid_dz, "DZ", varid_dz))

    ! get variable IDs for omega
    call check(nf90_open(trim(adjustl(path_omega)), nf90_write, ncid_omega))
    call check(nf90_inq_varid(ncid_omega, "ptend_coarse", varid_omega))

    ! get variable IDs for qv
    call check(nf90_open(trim(adjustl(path_qv)), nf90_write, ncid_qv))
    call check(nf90_inq_varid(ncid_qv, "sphum_coarse", varid_qv))

    ! get variable IDs for qw
    call check(nf90_open(trim(adjustl(path_qw)), nf90_write, ncid_qw))
    call check(nf90_inq_varid(ncid_qw, "liq_wat_coarse", varid_qw))

    ! get variable IDs for qr
    call check(nf90_open(trim(adjustl(path_qr)), nf90_write, ncid_qr))
    call check(nf90_inq_varid(ncid_qr, "rainwat_coarse", varid_qr))

    ! get variable IDs for qi
    call check(nf90_open(trim(adjustl(path_qi)), nf90_write, ncid_qi))
    call check(nf90_inq_varid(ncid_qi, "ice_wat_coarse", varid_qi))

    ! get variable IDs for qs
    call check(nf90_open(trim(adjustl(path_qs)), nf90_write, ncid_qs))
    call check(nf90_inq_varid(ncid_qs, "snowwat_coarse", varid_qs))

    ! get variable IDs for qg
    call check(nf90_open(trim(adjustl(path_qg)), nf90_write, ncid_qg))
    call check(nf90_inq_varid(ncid_qg, "graupel_coarse", varid_qg))

    ! get variable IDs for omqv
    call check(nf90_open(trim(adjustl(path_omqv)), nf90_write, ncid_omqv))
    call check(nf90_inq_varid(ncid_omqv, "omqv_coarse", varid_omqv))

    ! get variable IDs for omqw
    call check(nf90_open(trim(adjustl(path_omqw)), nf90_write, ncid_omqw))
    call check(nf90_inq_varid(ncid_omqw, "omql_coarse", varid_omqw))

    ! get variable IDs for omqr
    call check(nf90_open(trim(adjustl(path_omqr)), nf90_write, ncid_omqr))
    call check(nf90_inq_varid(ncid_omqr, "omqr_coarse", varid_omqr))

    ! get variable IDs for omqi
    call check(nf90_open(trim(adjustl(path_omqi)), nf90_write, ncid_omqi))
    call check(nf90_inq_varid(ncid_omqi, "omqi_coarse", varid_omqi))

    ! get variable IDs for omqs
    call check(nf90_open(trim(adjustl(path_omqs)), nf90_write, ncid_omqs))
    call check(nf90_inq_varid(ncid_omqs, "omqs_coarse", varid_omqs))

    ! get variable IDs for omqg
    call check(nf90_open(trim(adjustl(path_omqg)), nf90_write, ncid_omqg))
    call check(nf90_inq_varid(ncid_omqg, "omqg_coarse", varid_omqg))

    ! get variable IDs for lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_out)), nf90_write, ncid_lift_out))
    call check(nf90_inq_varid(ncid_lift_out, "water_lift", varid_lift_out))

    ! get variable IDs for vapor_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qv_out)), nf90_write, ncid_lift_qv_out))
    call check(nf90_inq_varid(ncid_lift_qv_out, "vapor_lift", varid_lift_qv_out))

    ! get variable IDs for liqwat_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qw_out)), nf90_write, ncid_lift_qw_out))
    call check(nf90_inq_varid(ncid_lift_qw_out, "liqwat_lift", varid_lift_qw_out))

    ! get variable IDs for rain_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qr_out)), nf90_write, ncid_lift_qr_out))
    call check(nf90_inq_varid(ncid_lift_qr_out, "rain_lift", varid_lift_qr_out))

    ! get variable IDs for icewat_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qi_out)), nf90_write, ncid_lift_qi_out))
    call check(nf90_inq_varid(ncid_lift_qi_out, "icewat_lift", varid_lift_qi_out))

    ! get variable IDs for snow_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qs_out)), nf90_write, ncid_lift_qs_out))
    call check(nf90_inq_varid(ncid_lift_qs_out, "snow_lift", varid_lift_qs_out))

    ! get variable IDs for graupel_lift (output variable)
    call check(nf90_open(trim(adjustl(path_lift_qg_out)), nf90_write, ncid_lift_qg_out))
    call check(nf90_inq_varid(ncid_lift_qg_out, "graupel_lift", varid_lift_qg_out))



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
            call check(nf90_get_var(ncid_omega, varid_omega, omega_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qv, varid_qv, qv_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qw, varid_qw, qw_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qr, varid_qr, qr_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qi, varid_qi, qi_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qs, varid_qs, qs_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qg, varid_qg, qg_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqv, varid_omqv, omqv_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqw, varid_omqw, omqw_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqr, varid_omqr, omqr_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqi, varid_omqi, omqi_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqs, varid_omqs, omqs_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqg, varid_omqg, omqg_buffer, (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            !$omp end master
            !$omp barrier
            !$omp do
            do y=1,chunk_size
                do x=1,nx
                    lift_out_buffer(x,y,1) = 0
                    lift_qv_out_buffer(x,y,1) = 0
                    lift_qw_out_buffer(x,y,1) = 0
                    lift_qr_out_buffer(x,y,1) = 0
                    lift_qi_out_buffer(x,y,1) = 0
                    lift_qs_out_buffer(x,y,1) = 0
                    lift_qg_out_buffer(x,y,1) = 0
                    do p=1,nplev
                        geometric_thickness = -dz_buffer(x,y,p,1)
                        lift_out_buffer(x,y,1) = lift_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * (qv_buffer(x,y,p,1) + qw_buffer(x,y,p,1) + qr_buffer(x,y,p,1) + qi_buffer(x,y,p,1) + qs_buffer(x,y,p,1) + qg_buffer(x,y,p,1)) &
                            & + omqv_buffer(x,y,p,1) + omqw_buffer(x,y,p,1) + omqr_buffer(x,y,p,1) + omqi_buffer(x,y,p,1) + omqs_buffer(x,y,p,1) + omqg_buffer(x,y,p,1) ) * geometric_thickness
                        lift_qv_out_buffer(x,y,1) = lift_qv_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qv_buffer(x,y,p,1) + omqv_buffer(x,y,p,1)) * geometric_thickness
                        lift_qw_out_buffer(x,y,1) = lift_qw_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qw_buffer(x,y,p,1) + omqw_buffer(x,y,p,1)) * geometric_thickness
                        lift_qr_out_buffer(x,y,1) = lift_qr_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qr_buffer(x,y,p,1) + omqr_buffer(x,y,p,1)) * geometric_thickness
                        lift_qi_out_buffer(x,y,1) = lift_qi_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qi_buffer(x,y,p,1) + omqi_buffer(x,y,p,1)) * geometric_thickness
                        lift_qs_out_buffer(x,y,1) = lift_qs_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qs_buffer(x,y,p,1) + omqs_buffer(x,y,p,1)) * geometric_thickness
                        lift_qg_out_buffer(x,y,1) = lift_qg_out_buffer(x,y,1) - (omega_buffer(x,y,p,1) * qg_buffer(x,y,p,1) + omqg_buffer(x,y,p,1)) * geometric_thickness
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp barrier
            !$omp master
            call check(nf90_put_var(ncid_lift_out, varid_lift_out, lift_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qv_out, varid_lift_qv_out, lift_qv_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qw_out, varid_lift_qw_out, lift_qw_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qr_out, varid_lift_qr_out, lift_qr_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qi_out, varid_lift_qi_out, lift_qi_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qs_out, varid_lift_qs_out, lift_qs_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            call check(nf90_put_var(ncid_lift_qg_out, varid_lift_qg_out, lift_qg_out_buffer, (/1,ystart,t/), (/nx,chunk_size,1/)))
            !$omp end master 
        enddo
    enddo
    !$omp end parallel

    call check(nf90_close(ncid_dz))
    call check(nf90_close(ncid_omega))
    call check(nf90_close(ncid_qv))
    call check(nf90_close(ncid_qw))
    call check(nf90_close(ncid_qr))
    call check(nf90_close(ncid_qi))
    call check(nf90_close(ncid_qs))
    call check(nf90_close(ncid_qg))
    call check(nf90_close(ncid_omqv))
    call check(nf90_close(ncid_omqw))
    call check(nf90_close(ncid_omqr))
    call check(nf90_close(ncid_omqi))
    call check(nf90_close(ncid_omqs))
    call check(nf90_close(ncid_omqg))
    call check(nf90_close(ncid_lift_out))
    call check(nf90_close(ncid_lift_qv_out))
    call check(nf90_close(ncid_lift_qw_out))
    call check(nf90_close(ncid_lift_qr_out))
    call check(nf90_close(ncid_lift_qi_out))
    call check(nf90_close(ncid_lift_qs_out))
    call check(nf90_close(ncid_lift_qg_out))

end program

