program compute_work

    use omp_lib
    use netcdf
    use interp
    use nc
    use hist
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    integer, parameter :: chunk_size = 24  ! chunk size along y axis
    integer, parameter :: nbuf = 2          ! number of pipeline buffers
    integer, parameter :: lat_south = -30, lat_north = 30, lon_west = 0, lon_east = 360  ! boundaries of domain for histogramming
    integer, parameter :: npr_edges = 400
    integer, parameter :: nwork_edges = 5001

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
    double precision, dimension(:), allocatable :: lon                     ! longitude
    double precision, dimension(:), allocatable :: lat                     ! latitude
    double precision, dimension(:), allocatable :: time_vals               ! time coordinate values
    double precision, dimension(:,:), allocatable :: cell_area             ! grid-cell area (nx x ny)
    double precision, dimension(:), allocatable :: pr_edges                ! precipitation histogram edges
    double precision, dimension(:), allocatable :: work_edges              ! work histogram edges
    double precision, dimension(:,:,:,:), allocatable :: pr_buffer         ! precipitation rate
    double precision, dimension(:,:,:,:), allocatable :: work_out_buffer   ! mechanical work (output buffer)
    double precision, dimension(:,:,:,:), allocatable :: lift_out_buffer   ! mechanical work to lift water (output buffer)
    double precision, dimension(:), allocatable :: hist_area_out           ! histogram output (area-weighted)
    double precision, dimension(:,:), allocatable :: hist2d_work_out       ! 2D histogram output (area-weighted)
    double precision, dimension(:,:), allocatable :: hist2d_lift_out       ! 2D histogram output (area-weighted)
    double precision, dimension(:), allocatable :: hist_area_chunk         ! per-chunk area histogram
    double precision, dimension(:,:), allocatable :: hist2d_work_chunk     ! per-chunk work histogram
    double precision, dimension(:,:), allocatable :: hist2d_lift_chunk     ! per-chunk lift histogram

      
    double precision :: geometric_thickness      ! geometric thickness of layer 
    double precision :: work_acc, lift_acc
    double precision :: pr_log_step
    integer :: nx, ny, nplev, nt                 ! note: nplev is full levels
    integer :: t, p, y, x, yc, nchunks, ystart, iedge
    ! Global clipped-domain index bounds used only for histogram accumulation.
    integer :: lat_clip_start, lat_clip_end, lon_clip_start, lon_clip_end
    ! Overlap between clipped latitude range and current chunk.
    integer :: hist_y_start, hist_y_end, y_local_start, y_local_end
    integer :: ibuf, istate
    integer :: buf_state(nbuf), buf_ystart(nbuf)
    logical :: active
    logical :: found_lat_start, found_lat_end, found_lon_start, found_lon_end
    integer :: tmp_unit, ios
    integer :: ncid_dz, ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg
    integer :: ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr
    integer :: ncid_work_out
    integer :: ncid_hist_out
    integer :: dimid_out_lon, dimid_out_lat, dimid_out_time
    integer :: tmp_dimids(4)
    integer :: ncstatus
    integer :: varid_dz, varid_temp, varid_omega, varid_qv, varid_qw, varid_qr, varid_qi, varid_qs, varid_qg
    integer :: varid_omt, varid_omqv, varid_omqw, varid_omqr, varid_omqi, varid_omqs, varid_omqg, varid_pr
    integer :: varid_lon, varid_lat, varid_time
    integer :: varid_out_lon, varid_out_lat, varid_out_time
    integer :: varid_work_out, varid_lift_out
    integer :: varid_hist_area, varid_hist2d_work, varid_hist2d_lift, varid_pr_edges, varid_work_edges
    integer :: dimid_nbin_pr, dimid_nbin_work, dimid_nedges_pr, dimid_nedges_work
    character(len=255) :: path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg
    character(len=255) :: path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, path_pr
    character(len=255) :: path_work_out, path_hist_out
    character(len=255) :: time_units, time_calendar

    namelist /config/ path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg, &
                     path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, & 
                     path_pr, &
                     path_work_out, path_hist_out

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

    ! open path_dz and read grid dimensions
    call check(nf90_open(trim(adjustl(path_dz)), nf90_nowrite, ncid_dz))
    call check(nf90_inq_varid(ncid_dz, "DZ", varid_dz))
    call check(nf90_inquire_variable(ncid_dz, varid_dz, dimids=tmp_dimids))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(1), len=nx))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(2), len=ny))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(3), len=nplev))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(4), len=nt))

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
    allocate(lon(nx))
    allocate(lat(ny))
    allocate(time_vals(nt))
    allocate(pr_edges(npr_edges))
    allocate(work_edges(nwork_edges))
    allocate(pr_buffer(nx, chunk_size, 1, nbuf))
    allocate(work_out_buffer(nx, chunk_size, 1, nbuf))
    allocate(lift_out_buffer(nx, chunk_size, 1, nbuf))
    allocate(hist_area_out(npr_edges - 1))
    allocate(hist2d_work_out(npr_edges - 1, nwork_edges - 1))
    allocate(hist2d_lift_out(npr_edges - 1, nwork_edges - 1))

    hist_area_out = 0.0d0
    hist2d_work_out = 0.0d0
    hist2d_lift_out = 0.0d0

    ! read coordinates and compute cell areas
    call check(nf90_inq_varid(ncid_dz, "grid_xt_coarse", varid_lon))
    call check(nf90_inq_varid(ncid_dz, "grid_yt_coarse", varid_lat))
    call check(nf90_get_var(ncid_dz, varid_lon, lon))
    call check(nf90_get_var(ncid_dz, varid_lat, lat))
    time_units = "hours since 1900-01-01 00:00:00"
    time_calendar = "standard"
    do t = 1, nt
        time_vals(t) = dble(t - 1)
    end do
    ncstatus = nf90_inq_varid(ncid_dz, "time", varid_time)
    if (ncstatus == nf90_noerr) then
        call check(nf90_get_var(ncid_dz, varid_time, time_vals))
        ncstatus = nf90_get_att(ncid_dz, varid_time, "units", time_units)
        if (ncstatus /= nf90_noerr) time_units = "hours since 1900-01-01 00:00:00"
        ncstatus = nf90_get_att(ncid_dz, varid_time, "calendar", time_calendar)
        if (ncstatus /= nf90_noerr) time_calendar = "standard"
    end if
    call compute_cell_areas(lon, lat, cell_area)

    ! Resolve longitude clipping indices on the model grid.
    found_lon_start = .false.
    found_lon_end = .false.
    do x = 1, nx
        if (.not. found_lon_start .and. lon(x) >= dble(lon_west)) then
            lon_clip_start = x
            found_lon_start = .true.
        end if
        if (lon(x) <= dble(lon_east)) then
            lon_clip_end = x
            found_lon_end = .true.
        end if
    end do

    if (.not. found_lon_start) lon_clip_start = 1
    if (.not. found_lon_end) lon_clip_end = nx

    ! Resolve latitude clipping indices on the model grid.
    found_lat_start = .false.
    found_lat_end = .false.
    do y = 1, ny
        if (.not. found_lat_start .and. lat(y) >= dble(lat_south)) then
            lat_clip_start = y
            found_lat_start = .true.
        end if
        if (lat(y) <= dble(lat_north)) then
            lat_clip_end = y
            found_lat_end = .true.
        end if
    end do

    if (.not. found_lat_start) lat_clip_start = 1
    if (.not. found_lat_end) lat_clip_end = ny

    ! open remaining input files
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
    call check(nf90_open(trim(adjustl(path_pr)), nf90_write, ncid_pr))

    ! create work output file
    call check(nf90_create(trim(adjustl(path_work_out)), nf90_clobber, ncid_work_out))
    call check(nf90_put_att(ncid_work_out, nf90_global, "Conventions", "CF-1.8"))
    call check(nf90_def_dim(ncid_work_out, "lon", nx, dimid_out_lon))
    call check(nf90_def_dim(ncid_work_out, "lat", ny, dimid_out_lat))
    call check(nf90_def_dim(ncid_work_out, "time", nt, dimid_out_time))
    call check(nf90_def_var(ncid_work_out, "lon", nf90_double, (/dimid_out_lon/), varid_out_lon))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, "long_name", "longitude"))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, "standard_name", "longitude"))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, "units", "degrees_east"))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, "axis", "X"))
    call check(nf90_def_var(ncid_work_out, "lat", nf90_double, (/dimid_out_lat/), varid_out_lat))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, "long_name", "latitude"))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, "standard_name", "latitude"))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, "units", "degrees_north"))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, "axis", "Y"))
    call check(nf90_def_var(ncid_work_out, "time", nf90_double, (/dimid_out_time/), varid_out_time))
    call check(nf90_put_att(ncid_work_out, varid_out_time, "long_name", "time"))
    call check(nf90_put_att(ncid_work_out, varid_out_time, "standard_name", "time"))
    call check(nf90_put_att(ncid_work_out, varid_out_time, "units", trim(time_units)))
    call check(nf90_put_att(ncid_work_out, varid_out_time, "calendar", trim(time_calendar)))
    call check(nf90_put_att(ncid_work_out, varid_out_time, "axis", "T"))
    call check(nf90_def_var(ncid_work_out, "work", nf90_double, (/dimid_out_lon, dimid_out_lat, dimid_out_time/), varid_work_out))
    call check(nf90_put_att(ncid_work_out, varid_work_out, "long_name", "work"))
    call check(nf90_put_att(ncid_work_out, varid_work_out, "units", "Watts per square meter"))
    call check(nf90_put_att(ncid_work_out, varid_work_out, "coordinates", "lon lat"))
    call check(nf90_def_var(ncid_work_out, "lift", nf90_double, (/dimid_out_lon, dimid_out_lat, dimid_out_time/), varid_lift_out))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, "long_name", "lift work"))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, "units", "Watts per square meter"))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, "coordinates", "lon lat"))
    call check(nf90_enddef(ncid_work_out))
    call check(nf90_put_var(ncid_work_out, varid_out_lon, lon))
    call check(nf90_put_var(ncid_work_out, varid_out_lat, lat))
    call check(nf90_put_var(ncid_work_out, varid_out_time, time_vals))

    ! populate variable IDs for input files
    call check(nf90_inq_varid(ncid_temp, "temp_coarse", varid_temp))
    call check(nf90_inq_varid(ncid_omega, "ptend_coarse", varid_omega))
    call check(nf90_inq_varid(ncid_qv, "sphum_coarse", varid_qv))
    call check(nf90_inq_varid(ncid_qw, "liq_wat_coarse", varid_qw))
    call check(nf90_inq_varid(ncid_qr, "rainwat_coarse", varid_qr))
    call check(nf90_inq_varid(ncid_qi, "ice_wat_coarse", varid_qi))
    call check(nf90_inq_varid(ncid_qs, "snowwat_coarse", varid_qs))
    call check(nf90_inq_varid(ncid_qg, "graupel_coarse", varid_qg))
    call check(nf90_inq_varid(ncid_omt, "omT_coarse", varid_omt))
    call check(nf90_inq_varid(ncid_omqv, "omqv_coarse", varid_omqv))
    call check(nf90_inq_varid(ncid_omqw, "omql_coarse", varid_omqw))
    call check(nf90_inq_varid(ncid_omqr, "omqr_coarse", varid_omqr))
    call check(nf90_inq_varid(ncid_omqi, "omqi_coarse", varid_omqi))
    call check(nf90_inq_varid(ncid_omqs, "omqs_coarse", varid_omqs))
    call check(nf90_inq_varid(ncid_omqg, "omqg_coarse", varid_omqg))
    call check(nf90_inq_varid(ncid_pr, "PRATEsfc_coarse", varid_pr))

    ! initialize precipitation histogram edges (logarithmically spaced)
    pr_log_step = 5.0d0 / dble(npr_edges - 1)
    do iedge = 1, npr_edges
        pr_edges(iedge) = 10.0d0 ** (-6.0d0 + dble(iedge - 1) * pr_log_step)
    end do

    ! initialize work histogram edges (linearly spaced)
    do iedge = 1, nwork_edges
        work_edges(iedge) = -2500.0d0 + dble(iedge - 1)
    end do

    ! determine number of y-chunks to process based on chunk size and total y dimension; assume ny is divisible by chunk_size for simplicity
    nchunks = ny/chunk_size

    ! Buffer state convention:
    !   0 = free slot
    !   1 = compute task is running on this buffer
    !   2 = compute finished; chunk is ready to be written to output NetCDF
    ! Pipeline over y-chunks: load, compute work/lift, accumulate clipped histograms, then flush outputs.
    !$omp parallel default(shared) private(t,yc,ibuf,istate,ystart,y,x,p,geometric_thickness,work_acc,lift_acc,active,hist_y_start,hist_y_end,y_local_start,y_local_end)
    !$omp master
    do t = 1, nt
        print *, 'timestep', t, 'of', nt

        ! Reset buffer metadata for this timestep.
        buf_state = 0
        buf_ystart = 0

        ! Producer loop: load one y-chunk into a reusable pipeline buffer.
        do yc = 1, nchunks
            ibuf = mod(yc - 1, nbuf) + 1

            ! If this buffer is still computing (state 1), keep yielding until it is reusable.
            do
                !$omp atomic read
                istate = buf_state(ibuf)
                if (istate /= 1) exit
                !$omp taskyield
            end do

            ! Opportunistic drain: if a previous chunk in this buffer is ready, write it now.
            if (istate == 2) then
                call check(nf90_put_var(ncid_work_out, varid_work_out, work_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                call check(nf90_put_var(ncid_work_out, varid_lift_out, lift_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                !$omp atomic write
                buf_state(ibuf) = 0
            end if

            ystart = (yc - 1) * chunk_size + 1
            buf_ystart(ibuf) = ystart

            print *, 'loading y-chunk', yc, 'into pipeline'

            ! Fill this buffer with all required input fields for one chunk.
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
            call check(nf90_get_var(ncid_pr, varid_pr, pr_buffer(:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,1,1/)))

            ! Mark as in-flight before launching the compute task.
            !$omp atomic write
            buf_state(ibuf) = 1

            ! Per-chunk compute task: derive work/lift columns and update histogram chunks.
            !$omp task firstprivate(ibuf,ystart) private(y,x,p,geometric_thickness,work_acc,lift_acc,hist_area_chunk,hist2d_work_chunk,hist2d_lift_chunk,hist_y_start,hist_y_end,y_local_start,y_local_end)
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

            ! Convert clipped global y-bounds to local chunk indices before histogramming.
            hist_y_start = max(ystart, lat_clip_start)
            hist_y_end = min(ystart + chunk_size - 1, lat_clip_end)
            if (hist_y_start <= hist_y_end .and. lon_clip_start <= lon_clip_end) then
                y_local_start = hist_y_start - ystart + 1
                y_local_end = hist_y_end - ystart + 1

                call bin_sumarea(pr_buffer(lon_clip_start:lon_clip_end,y_local_start:y_local_end,:,ibuf), &
                                 cell_area(lon_clip_start:lon_clip_end,hist_y_start:hist_y_end), pr_edges, hist_area_chunk)
                call joint_bin_sumarea(pr_buffer(lon_clip_start:lon_clip_end,y_local_start:y_local_end,:,ibuf), &
                                       work_out_buffer(lon_clip_start:lon_clip_end,y_local_start:y_local_end,:,ibuf), &
                                       cell_area(lon_clip_start:lon_clip_end,hist_y_start:hist_y_end), pr_edges, work_edges, hist2d_work_chunk)
                call joint_bin_sumarea(pr_buffer(lon_clip_start:lon_clip_end,y_local_start:y_local_end,:,ibuf), &
                                       lift_out_buffer(lon_clip_start:lon_clip_end,y_local_start:y_local_end,:,ibuf), &
                                       cell_area(lon_clip_start:lon_clip_end,hist_y_start:hist_y_end), pr_edges, work_edges, hist2d_lift_chunk)
            end if

            ! Some chunks may not intersect the clipped domain; only accumulate allocated chunk histograms.
            !$omp critical(histogram_accumulation)
            if (allocated(hist_area_chunk)) hist_area_out = hist_area_out + hist_area_chunk
            if (allocated(hist2d_work_chunk)) hist2d_work_out = hist2d_work_out + hist2d_work_chunk
            if (allocated(hist2d_lift_chunk)) hist2d_lift_out = hist2d_lift_out + hist2d_lift_chunk
            !$omp end critical(histogram_accumulation)

            if (allocated(hist_area_chunk)) deallocate(hist_area_chunk)
            if (allocated(hist2d_work_chunk)) deallocate(hist2d_work_chunk)
            if (allocated(hist2d_lift_chunk)) deallocate(hist2d_lift_chunk)

            ! Hand off completed chunk to writer side of the pipeline.
            !$omp atomic write
            buf_state(ibuf) = 2
            !$omp end task
        end do

        ! Final drain for this timestep: flush all buffers that complete after producer loop ends.
        do
            active = .false.
            do ibuf = 1, nbuf
                !$omp atomic read
                istate = buf_state(ibuf)
                if (istate == 1) then
                    active = .true.
                else if (istate == 2) then
                    call check(nf90_put_var(ncid_work_out, varid_work_out, work_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                    call check(nf90_put_var(ncid_work_out, varid_lift_out, lift_out_buffer(:,:,:,ibuf), (/1,buf_ystart(ibuf),t/), (/nx,chunk_size,1/)))
                    !$omp atomic write
                    buf_state(ibuf) = 0
                    active = .true.
                end if
            end do
            if (.not. active) exit
            !$omp taskyield
        end do

        ! Ensure all tasks for this timestep are complete before advancing t.
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

    ! write histogram outputs to NetCDF
    call check(nf90_create(trim(adjustl(path_hist_out)), nf90_clobber, ncid_hist_out))
    call check(nf90_def_dim(ncid_hist_out, "nbin_pr", npr_edges - 1, dimid_nbin_pr))
    call check(nf90_def_dim(ncid_hist_out, "nbin_work", nwork_edges - 1, dimid_nbin_work))
    call check(nf90_def_dim(ncid_hist_out, "nedges_pr", npr_edges, dimid_nedges_pr))
    call check(nf90_def_dim(ncid_hist_out, "nedges_work", nwork_edges, dimid_nedges_work))
    call check(nf90_def_var(ncid_hist_out, "hist_area", nf90_double, (/dimid_nbin_pr/), varid_hist_area))
    call check(nf90_def_var(ncid_hist_out, "hist2d_work", nf90_double, (/dimid_nbin_pr, dimid_nbin_work/), varid_hist2d_work))
    call check(nf90_def_var(ncid_hist_out, "hist2d_lift", nf90_double, (/dimid_nbin_pr, dimid_nbin_work/), varid_hist2d_lift))
    call check(nf90_def_var(ncid_hist_out, "pr_edges", nf90_double, (/dimid_nedges_pr/), varid_pr_edges))
    call check(nf90_def_var(ncid_hist_out, "work_edges", nf90_double, (/dimid_nedges_work/), varid_work_edges))
    call check(nf90_put_att(ncid_hist_out, nf90_global, "clip_lon_west", lon(lon_clip_start)))
    call check(nf90_put_att(ncid_hist_out, nf90_global, "clip_lon_east", lon(lon_clip_end)))
    call check(nf90_put_att(ncid_hist_out, nf90_global, "clip_lat_south", lat(lat_clip_start)))
    call check(nf90_put_att(ncid_hist_out, nf90_global, "clip_lat_north", lat(lat_clip_end)))
    call check(nf90_enddef(ncid_hist_out))
    call check(nf90_put_var(ncid_hist_out, varid_hist_area, hist_area_out))
    call check(nf90_put_var(ncid_hist_out, varid_hist2d_work, hist2d_work_out))
    call check(nf90_put_var(ncid_hist_out, varid_hist2d_lift, hist2d_lift_out))
    call check(nf90_put_var(ncid_hist_out, varid_pr_edges, pr_edges))
    call check(nf90_put_var(ncid_hist_out, varid_work_edges, work_edges))
    call check(nf90_close(ncid_hist_out))

end program compute_work
