program compute_work_omp_pipeline

    use omp_lib
    use netcdf
    use libinterp
    use libnc
    use libhist
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    integer, parameter :: chunk_size = 24  ! chunk size along y axis
    integer, parameter :: nbuf = 2          ! number of pipeline buffers

    double precision, dimension(:,:,:,:,:), allocatable :: dz_buffer       ! layer geometric depth
    double precision, dimension(:,:,:,:,:), allocatable :: temp_buffer     ! temperature
    double precision, dimension(:,:,:,:,:), allocatable :: omega_buffer    ! omega
    double precision, dimension(:,:,:,:,:), allocatable :: qv_buffer       ! qv
    double precision, dimension(:,:,:,:,:), allocatable :: qw_buffer       ! qw
    double precision, dimension(:,:,:,:,:), allocatable :: qr_buffer       ! qr
    double precision, dimension(:,:,:,:,:), allocatable :: qi_buffer       ! qi
    double precision, dimension(:,:,:,:,:), allocatable :: qs_buffer       ! qs
    double precision, dimension(:,:,:,:,:), allocatable :: qg_buffer       ! qg
    double precision, dimension(:,:,:,:,:), allocatable :: omt_buffer      ! omT
    double precision, dimension(:,:,:,:,:), allocatable :: omqv_buffer     ! omqv
    double precision, dimension(:,:,:,:,:), allocatable :: omqw_buffer     ! omqw
    double precision, dimension(:,:,:,:,:), allocatable :: omqr_buffer     ! omqr
    double precision, dimension(:,:,:,:,:), allocatable :: omqi_buffer     ! omqi
    double precision, dimension(:,:,:,:,:), allocatable :: omqs_buffer     ! omqs
    double precision, dimension(:,:,:,:,:), allocatable :: omqg_buffer     ! omqg
    double precision, dimension(:,:,:,:), allocatable :: work_out_buffer   ! mechanical work (output buffer)
    double precision, dimension(:,:,:,:), allocatable :: lift_out_buffer   ! mechanical work to lift water (output buffer)
      
    double precision :: geometric_thickness      ! geometric thickness of layer 
    double precision :: work_acc, lift_acc
    integer :: nx, ny, nplev, nt                 ! note: nplev is full levels
    integer :: t, p, y, x, yc, nchunks, ystart
    integer :: ibuf, istate
    integer :: buf_state(nbuf), buf_ystart(nbuf)
    logical :: active
    integer :: tmp_unit, ios
    integer :: ncid_dz, ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg
    integer :: ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg
    integer :: ncid_work_out, ncid_lift_out
    integer :: varid_dz, varid_temp, varid_omega, varid_qv, varid_qw, varid_qr, varid_qi, varid_qs, varid_qg
    integer :: varid_omt, varid_omqv, varid_omqw, varid_omqr, varid_omqi, varid_omqs, varid_omqg
    integer :: varid_work_out, varid_lift_out
    character(len=255) :: path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg
    character(len=255) :: path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg
    character(len=255) :: path_work_out, path_lift_out

    namelist /config/ path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg, &
                     path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, &
                     path_work_out, path_lift_out, nx, ny, nplev, nt

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
    allocate(dz_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(temp_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omega_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qv_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qw_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qr_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qi_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qs_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(qg_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omt_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqv_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqw_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqr_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqi_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqs_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(omqg_buffer(nx, chunk_size, nplev, 1, nbuf))
    allocate(work_out_buffer(nx, chunk_size, 1, nbuf))
    allocate(lift_out_buffer(nx, chunk_size, 1, nbuf))

    ! populate file IDs
    call check(nf90_open(trim(adjustl(path_dz)), nf90_write, ncid_dz))
    call check(nf90_open(trim(adjustl(path_temp)), nf90_write, ncid_temp))
    call check(nf90_open(trim(adjustl(path_omega)), nf90_write, ncid_omega))
    call check(nf90_open(trim(adjustl(path_qv)), nf90_write, ncid_qv))
    call check(nf90_open(trim(adjustl(path_qw)), nf90_write, ncid_qw))
    call check(nf90_open(trim(adjustl(path_qr)), nf90_write, ncid_qr))
    call check(nf90_open(trim(adjustl(path_qi)), nf90_write, ncid_qi))
    call check(nf90_open(trim(adjustl(path_qs)), nf90_write, ncid_qs))
    call check(nf90_open(trim(adjustl(path_qg)), nf90_write, ncid_qg))
    call check(nf90_open(trim(adjustl(path_omt)), nf90_write, ncid_omt))
    call check(nf90_open(trim(adjustl(path_omqv)), nf90_write, ncid_omqv))
    call check(nf90_open(trim(adjustl(path_omqw)), nf90_write, ncid_omqw))
    call check(nf90_open(trim(adjustl(path_omqr)), nf90_write, ncid_omqr))
    call check(nf90_open(trim(adjustl(path_omqi)), nf90_write, ncid_omqi))
    call check(nf90_open(trim(adjustl(path_omqs)), nf90_write, ncid_omqs))
    call check(nf90_open(trim(adjustl(path_omqg)), nf90_write, ncid_omqg))
    call check(nf90_open(trim(adjustl(path_work_out)), nf90_write, ncid_work_out))
    call check(nf90_open(trim(adjustl(path_lift_out)), nf90_write, ncid_lift_out))

    ! populate variable IDs
    call check(nf90_inq_varid(ncid_dz, "DZ", varid_dz))
    call check(nf90_inq_varid(ncid_temp, "temp_coarse", varid_temp))
    call check(nf90_inq_varid(ncid_omega, "ptend_coarse", varid_omega))
    call check(nf90_inq_varid(ncid_qv, "sphum_coarse", varid_qv))
    call check(nf90_inq_varid(ncid_qw, "liq_wat_coarse", varid_qw))
    call check(nf90_inq_varid(ncid_qr, "rain_coarse", varid_qr))
    call check(nf90_inq_varid(ncid_qi, "ice_coarse", varid_qi))
    call check(nf90_inq_varid(ncid_qs, "snow_coarse", varid_qs))
    call check(nf90_inq_varid(ncid_qg, "graupel_coarse", varid_qg))
    call check(nf90_inq_varid(ncid_omt, "omT_coarse", varid_omt))
    call check(nf90_inq_varid(ncid_omqv, "omqv_coarse", varid_omqv))
    call check(nf90_inq_varid(ncid_omqw, "omqw_coarse", varid_omqw))
    call check(nf90_inq_varid(ncid_omqr, "omqr_coarse", varid_omqr))
    call check(nf90_inq_varid(ncid_omqi, "omqi_coarse", varid_omqi))
    call check(nf90_inq_varid(ncid_omqs, "omqs_coarse", varid_omqs))
    call check(nf90_inq_varid(ncid_omqg, "omqg_coarse", varid_omqg))
    call check(nf90_inq_varid(ncid_work_out, "work", varid_work_out))
    call check(nf90_inq_varid(ncid_lift_out, "lift", varid_lift_out))

    nchunks = ny/chunk_size

    !$omp parallel default(shared) private(t,yc,ibuf,istate,ystart,y,x,p,geometric_thickness,work_acc,lift_acc,active)
    !$omp master
    do t = 1, nt
        print *, 'timestep', t, 'of', nt

        buf_state = 0
        buf_ystart = 0

        do yc = 1, nchunks
            ibuf = mod(yc - 1, nbuf) + 1

            do
                !$omp atomic read
                istate = buf_state(ibuf)
                if (istate /= 1) exit
                !$omp taskyield
            end do

            if (istate == 2) then
                call check(nf90_put_var(ncid_work_out, varid_work_out, work_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                call check(nf90_put_var(ncid_lift_out, varid_lift_out, lift_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                !$omp atomic write
                buf_state(ibuf) = 0
            end if

            ystart = (yc - 1) * chunk_size + 1
            buf_ystart(ibuf) = ystart

            call check(nf90_get_var(ncid_dz, varid_dz, dz_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_temp, varid_temp, temp_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omega, varid_omega, omega_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qv, varid_qv, qv_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qw, varid_qw, qw_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qr, varid_qr, qr_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qi, varid_qi, qi_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qs, varid_qs, qs_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_qg, varid_qg, qg_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omt, varid_omt, omt_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqv, varid_omqv, omqv_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqw, varid_omqw, omqw_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqr, varid_omqr, omqr_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqi, varid_omqi, omqi_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqs, varid_omqs, omqs_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
            call check(nf90_get_var(ncid_omqg, varid_omqg, omqg_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))

            !$omp atomic write
            buf_state(ibuf) = 1

            !$omp task firstprivate(ibuf) private(y,x,p,geometric_thickness,work_acc,lift_acc)
            do y = 1, chunk_size
                do x = 1, nx
                    work_acc = 0.0d0
                    lift_acc = 0.0d0
                    do p = 1, nplev
                        geometric_thickness = -dz_buffer(x,y,p,1,ibuf)
                        work_acc = work_acc - &
                            (omt_buffer(x,y,p,1,ibuf) / temp_buffer(x,y,p,1,ibuf) + 1.61d0 * (omqv_buffer(x,y,p,1,ibuf) + omega_buffer(x,y,p,1,ibuf) * qv_buffer(x,y,p,1,ibuf))) * geometric_thickness
                        lift_acc = lift_acc - &
                            (omega_buffer(x,y,p,1,ibuf) * (qv_buffer(x,y,p,1,ibuf) + qw_buffer(x,y,p,1,ibuf) + qr_buffer(x,y,p,1,ibuf) + qi_buffer(x,y,p,1,ibuf) + qs_buffer(x,y,p,1,ibuf) + qg_buffer(x,y,p,1,ibuf)) &
                            + omqv_buffer(x,y,p,1,ibuf) + omqw_buffer(x,y,p,1,ibuf) + omqr_buffer(x,y,p,1,ibuf) + omqi_buffer(x,y,p,1,ibuf) + omqs_buffer(x,y,p,1,ibuf) + omqg_buffer(x,y,p,1,ibuf)) * geometric_thickness
                    end do
                    work_out_buffer(x,y,1,ibuf) = work_acc
                    lift_out_buffer(x,y,1,ibuf) = lift_acc
                end do
            end do
            !$omp atomic write
            buf_state(ibuf) = 2
            !$omp end task
        end do

        do
            active = .false.
            do ibuf = 1, nbuf
                !$omp atomic read
                istate = buf_state(ibuf)
                if (istate == 1) then
                    active = .true.
                else if (istate == 2) then
                    call check(nf90_put_var(ncid_work_out, varid_work_out, work_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                    call check(nf90_put_var(ncid_lift_out, varid_lift_out, lift_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                    !$omp atomic write
                    buf_state(ibuf) = 0
                    active = .true.
                end if
            end do
            if (.not. active) exit
            !$omp taskyield
        end do

        !$omp taskwait
    end do
    !$omp end master
    !$omp end parallel

    call check(nf90_close(ncid_dz))
    call check(nf90_close(ncid_temp))
    call check(nf90_close(ncid_omega))
    call check(nf90_close(ncid_qv))
    call check(nf90_close(ncid_qw))
    call check(nf90_close(ncid_qr))
    call check(nf90_close(ncid_qi))
    call check(nf90_close(ncid_qs))
    call check(nf90_close(ncid_qg))
    call check(nf90_close(ncid_omt))
    call check(nf90_close(ncid_omqv))
    call check(nf90_close(ncid_omqw))
    call check(nf90_close(ncid_omqr))
    call check(nf90_close(ncid_omqi))
    call check(nf90_close(ncid_omqs))
    call check(nf90_close(ncid_omqg))
    call check(nf90_close(ncid_work_out))
    call check(nf90_close(ncid_lift_out))

end program compute_work_omp_pipeline
