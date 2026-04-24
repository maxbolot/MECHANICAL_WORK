! async fork of compute_work_async.f90 with precipitation-threshold masking
!
! Key structural behavior is unchanged from the async implementation:
! NetCDF reads are launched in !$omp task blocks (serialized inside
! !$omp critical(netcdf_io)) so I/O and compute can overlap.
!
! Additional behavior in this fork:
! - Thresholds are read from an ASCII percentile file.
! - work/lift are computed once per grid point and then masked for each
!   percentile threshold using PRATEsfc_coarse >= threshold.
! - Daily event counts are tracked and work/lift are written as
!   conditional means over threshold-exceedance events.
! - Outputs include a leading percentile dimension.
!
! Buffer state convention:
!   0 = free
!   2 = compute finished for current chunk
!   3 = in use (being read OR computed)
!
! Task dependency token array fetch_ready(nbuf) is used only as an OMP
! dependency address; its integer value is never read.

program compute_work_async_prate_threshold

    use omp_lib
    use netcdf
    use nc
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    integer, parameter :: chunk_size = 144  ! chunk size along y axis
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

    double precision, dimension(:), allocatable :: lon                     ! longitude
    double precision, dimension(:), allocatable :: lat                     ! latitude
    double precision, dimension(:), allocatable :: time_vals               ! time coordinate values
    double precision, dimension(:), allocatable :: day_time_vals           ! representative time value for each aggregated day
    double precision, dimension(:), allocatable :: prate_percentiles       ! percentile levels loaded from ASCII
    double precision, dimension(:), allocatable :: prate_thresholds        ! precipitation thresholds loaded from ASCII

    double precision, dimension(:,:,:,:), allocatable :: pr_buffer         ! precipitation rate
    double precision, dimension(:,:,:,:,:), allocatable :: work_out_buffer ! thresholded work output buffer
    double precision, dimension(:,:,:,:,:), allocatable :: lift_out_buffer ! thresholded lift output buffer
    double precision, dimension(:,:,:,:,:), allocatable :: count_out_buffer ! thresholded count of timesteps exceeding threshold (for conditional mean calculation)
    double precision, dimension(:,:,:,:), allocatable :: work_accum_full   ! full-grid time accumulation of work by percentile
    double precision, dimension(:,:,:,:), allocatable :: lift_accum_full   ! full-grid time accumulation of lift by percentile
    double precision, dimension(:,:,:,:), allocatable :: count_accum_full  ! full-grid time accumulation of threshold exceedance count by percentile

    double precision :: geometric_thickness      ! geometric thickness of layer
    double precision :: work_acc, lift_acc       ! per-gridpoint unmasked work/lift

    integer :: nx, ny, nplev, nt                 ! note: nplev is full levels
    integer :: ndays
    integer :: npercentiles                      ! number of loaded percentile thresholds
    integer :: t, p, y, x, yc, nchunks, ystart, iday
    integer :: ip                                ! percentile index
    integer :: ibuf, istate
    integer :: buf_state(nbuf), buf_ystart(nbuf)
    integer :: fetch_ready(nbuf)                 ! task dependency token (address only, never read)
    integer, dimension(:), allocatable :: day_index_of_t
    integer, dimension(:), allocatable :: day_first_idx
    integer :: tmp_unit, ios

    integer :: ncid_dz, ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg
    integer :: ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr
    integer :: ncid_work_out

    integer :: dimid_out_pct, dimid_out_lon, dimid_out_lat, dimid_out_time
    integer :: tmp_dimids(4)
    integer :: ncstatus

    integer :: varid_dz, varid_temp, varid_omega, varid_qv, varid_qw, varid_qr, varid_qi, varid_qs, varid_qg
    integer :: varid_omt, varid_omqv, varid_omqw, varid_omqr, varid_omqi, varid_omqs, varid_omqg, varid_pr
    integer :: varid_lon, varid_lat, varid_time
    integer :: varid_out_pct, varid_out_lon, varid_out_lat, varid_out_time
    integer :: varid_work_out, varid_lift_out, varid_count_out

    character(len=255) :: path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg
    character(len=255) :: path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, path_pr
    character(len=255) :: path_work_out, path_thresholds
    character(len=255) :: time_units, time_calendar
    character(len=64) :: attr_name               ! global attribute key for per-percentile threshold

    namelist /config/ path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg, &
                     path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, &
                     path_pr, path_work_out, path_thresholds

    ! get configuration namelist from command line and read it afterwards
    call get_command_argument(1, filenml)
    open(newunit=tmp_unit, file=trim(adjustl(filenml)), status='old', iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to open configuration namelist, iomsg='//trim(msg)
        stop 1
    end if

    read(unit=tmp_unit, nml=config, iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to assign namelist objects, iomsg='//trim(msg)
        stop 1
    end if
    close(tmp_unit)

    ! load percentile thresholds (percentile, threshold_prate) from ASCII output
    call read_thresholds_ascii(trim(adjustl(path_thresholds)), prate_percentiles, prate_thresholds)
    npercentiles = size(prate_percentiles)
    if (npercentiles <= 0) then
        error stop 'No precipitation thresholds loaded from threshold file.'
    end if

    ! open path_dz and read grid dimensions
    call check(nf90_open(trim(adjustl(path_dz)), nf90_nowrite, ncid_dz))
    call check(nf90_inq_varid(ncid_dz, 'DZ', varid_dz))
    call check(nf90_inquire_variable(ncid_dz, varid_dz, dimids=tmp_dimids))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(1), len=nx))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(2), len=ny))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(3), len=nplev))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(4), len=nt))

    ! allocate async pipeline buffers
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

    allocate(pr_buffer(nx, chunk_size, 1, nbuf))
    allocate(work_out_buffer(npercentiles, nx, chunk_size, 1, nbuf))
    allocate(lift_out_buffer(npercentiles, nx, chunk_size, 1, nbuf))
    allocate(count_out_buffer(npercentiles, nx, chunk_size, 1, nbuf))

    allocate(work_accum_full(npercentiles, nx, ny, nt))
    allocate(lift_accum_full(npercentiles, nx, ny, nt))
    allocate(count_accum_full(npercentiles, nx, ny, nt))
    work_accum_full = 0.0d0
    lift_accum_full = 0.0d0
    count_accum_full = 0.0d0

    call check(nf90_inq_varid(ncid_dz, 'grid_xt_coarse', varid_lon))
    call check(nf90_inq_varid(ncid_dz, 'grid_yt_coarse', varid_lat))
    call check(nf90_get_var(ncid_dz, varid_lon, lon))
    call check(nf90_get_var(ncid_dz, varid_lat, lat))

    call check(nf90_open(trim(adjustl(path_temp)), nf90_nowrite, ncid_temp))
    time_units = 'hours since 1900-01-01 00:00:00'
    time_calendar = 'standard'
    do t = 1, nt
        time_vals(t) = dble(t - 1)
    end do

    ncstatus = nf90_inq_varid(ncid_temp, 'time', varid_time)
    if (ncstatus == nf90_noerr) then
        call check(nf90_get_var(ncid_temp, varid_time, time_vals))
        ncstatus = nf90_get_att(ncid_temp, varid_time, 'units', time_units)
        if (ncstatus /= nf90_noerr) time_units = 'hours since 1900-01-01 00:00:00'
        ncstatus = nf90_get_att(ncid_temp, varid_time, 'calendar', time_calendar)
        if (ncstatus /= nf90_noerr) time_calendar = 'standard'
    end if
    ! Convert native model timesteps into contiguous day bins so outputs are
    ! aggregated by day (not by raw input timestep index).
    call build_day_groups(time_vals, time_units, day_index_of_t, day_first_idx, day_time_vals, ndays)

    call check(nf90_open(trim(adjustl(path_omega)), nf90_nowrite, ncid_omega))
    call check(nf90_open(trim(adjustl(path_qv)), nf90_nowrite, ncid_qv))
    call check(nf90_open(trim(adjustl(path_qw)), nf90_nowrite, ncid_qw))
    call check(nf90_open(trim(adjustl(path_qr)), nf90_nowrite, ncid_qr))
    call check(nf90_open(trim(adjustl(path_qi)), nf90_nowrite, ncid_qi))
    call check(nf90_open(trim(adjustl(path_qs)), nf90_nowrite, ncid_qs))
    call check(nf90_open(trim(adjustl(path_qg)), nf90_nowrite, ncid_qg))
    call check(nf90_open(trim(adjustl(path_omt)), nf90_nowrite, ncid_omt))
    call check(nf90_open(trim(adjustl(path_omqv)), nf90_nowrite, ncid_omqv))
    call check(nf90_open(trim(adjustl(path_omqw)), nf90_nowrite, ncid_omqw))
    call check(nf90_open(trim(adjustl(path_omqr)), nf90_nowrite, ncid_omqr))
    call check(nf90_open(trim(adjustl(path_omqi)), nf90_nowrite, ncid_omqi))
    call check(nf90_open(trim(adjustl(path_omqs)), nf90_nowrite, ncid_omqs))
    call check(nf90_open(trim(adjustl(path_omqg)), nf90_nowrite, ncid_omqg))
    call check(nf90_open(trim(adjustl(path_pr)), nf90_nowrite, ncid_pr))

    ! create thresholded work output file
    call check(nf90_create(trim(adjustl(path_work_out)), nf90_clobber, ncid_work_out))
    call check(nf90_put_att(ncid_work_out, nf90_global, 'Conventions', 'CF-1.8'))
    call check(nf90_put_att(ncid_work_out, nf90_global, 'threshold_file', trim(adjustl(path_thresholds))))
    call check(nf90_put_att(ncid_work_out, nf90_global, 'daily_stat', 'mixed_by_variable'))
    call check(nf90_put_att(ncid_work_out, nf90_global, 'daily_aggregation', 'mean over native timesteps within each day'))

    call check(nf90_def_dim(ncid_work_out, 'percentile', npercentiles, dimid_out_pct))
    call check(nf90_def_dim(ncid_work_out, 'lon', nx, dimid_out_lon))
    call check(nf90_def_dim(ncid_work_out, 'lat', ny, dimid_out_lat))
    ! Output time axis is daily, with one representative timestamp per day.
    call check(nf90_def_dim(ncid_work_out, 'time', ndays, dimid_out_time))

    call check(nf90_def_var(ncid_work_out, 'percentile', nf90_double, (/dimid_out_pct/), varid_out_pct))
    call check(nf90_put_att(ncid_work_out, varid_out_pct, 'long_name', 'precipitation percentile threshold'))
    call check(nf90_put_att(ncid_work_out, varid_out_pct, 'units', '1'))

    call check(nf90_def_var(ncid_work_out, 'lon', nf90_double, (/dimid_out_lon/), varid_out_lon))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, 'long_name', 'longitude'))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, 'standard_name', 'longitude'))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, 'units', 'degrees_east'))
    call check(nf90_put_att(ncid_work_out, varid_out_lon, 'axis', 'X'))

    call check(nf90_def_var(ncid_work_out, 'lat', nf90_double, (/dimid_out_lat/), varid_out_lat))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, 'long_name', 'latitude'))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, 'standard_name', 'latitude'))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, 'units', 'degrees_north'))
    call check(nf90_put_att(ncid_work_out, varid_out_lat, 'axis', 'Y'))

    call check(nf90_def_var(ncid_work_out, 'time', nf90_double, (/dimid_out_time/), varid_out_time))
    call check(nf90_put_att(ncid_work_out, varid_out_time, 'long_name', 'time'))
    call check(nf90_put_att(ncid_work_out, varid_out_time, 'standard_name', 'time'))
    call check(nf90_put_att(ncid_work_out, varid_out_time, 'units', trim(time_units)))
    call check(nf90_put_att(ncid_work_out, varid_out_time, 'calendar', trim(time_calendar)))
    call check(nf90_put_att(ncid_work_out, varid_out_time, 'axis', 'T'))

    call check(nf90_def_var(ncid_work_out, 'work', nf90_double, (/dimid_out_pct, dimid_out_lon, dimid_out_lat, dimid_out_time/), varid_work_out))
    call check(nf90_put_att(ncid_work_out, varid_work_out, 'long_name', 'work filtered by precipitation threshold'))
    call check(nf90_put_att(ncid_work_out, varid_work_out, 'units', 'Watts per square meter'))
    call check(nf90_put_att(ncid_work_out, varid_work_out, 'coordinates', 'percentile lon lat'))
    call check(nf90_put_att(ncid_work_out, varid_work_out, 'time_stat', 'daily_conditional_mean'))
    call check(nf90_put_att(ncid_work_out, varid_work_out, 'time_aggregation', 'mean over threshold exceedance events within each day'))
    call check(nf90_put_att(ncid_work_out, varid_work_out, '_FillValue', -9999.0d0))

    call check(nf90_def_var(ncid_work_out, 'lift', nf90_double, (/dimid_out_pct, dimid_out_lon, dimid_out_lat, dimid_out_time/), varid_lift_out))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, 'long_name', 'lift work filtered by precipitation threshold'))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, 'units', 'Watts per square meter'))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, 'coordinates', 'percentile lon lat'))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, 'time_stat', 'daily_conditional_mean'))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, 'time_aggregation', 'mean over threshold exceedance events within each day'))
    call check(nf90_put_att(ncid_work_out, varid_lift_out, '_FillValue', -9999.0d0))

    call check(nf90_def_var(ncid_work_out, 'event_count', nf90_double, (/dimid_out_pct, dimid_out_lon, dimid_out_lat, dimid_out_time/), varid_count_out))
    call check(nf90_put_att(ncid_work_out, varid_count_out, 'long_name', 'number of events exceeding precipitation threshold'))
    call check(nf90_put_att(ncid_work_out, varid_count_out, 'units', 'count'))
    call check(nf90_put_att(ncid_work_out, varid_count_out, 'coordinates', 'percentile lon lat'))
    call check(nf90_put_att(ncid_work_out, varid_count_out, 'time_stat', 'daily_count'))
    call check(nf90_put_att(ncid_work_out, varid_count_out, 'time_aggregation', 'sum of threshold exceedance events over native timesteps within each day'))

    ! Store threshold values in global attributes for traceability
    ! (for example: p99_threshold, p9995_threshold, ...).
    do ip = 1, npercentiles
        attr_name = percentile_attr_name(prate_percentiles(ip))
        call check(nf90_put_att(ncid_work_out, nf90_global, trim(attr_name), prate_thresholds(ip)))
    end do

    call check(nf90_enddef(ncid_work_out))

    call check(nf90_put_var(ncid_work_out, varid_out_pct, prate_percentiles))
    call check(nf90_put_var(ncid_work_out, varid_out_lon, lon))
    call check(nf90_put_var(ncid_work_out, varid_out_lat, lat))
    call check(nf90_put_var(ncid_work_out, varid_out_time, day_time_vals))

    call check(nf90_inq_varid(ncid_temp, 'temp_coarse', varid_temp))
    call check(nf90_inq_varid(ncid_omega, 'ptend_coarse', varid_omega))
    call check(nf90_inq_varid(ncid_qv, 'sphum_coarse', varid_qv))
    call check(nf90_inq_varid(ncid_qw, 'liq_wat_coarse', varid_qw))
    call check(nf90_inq_varid(ncid_qr, 'rainwat_coarse', varid_qr))
    call check(nf90_inq_varid(ncid_qi, 'ice_wat_coarse', varid_qi))
    call check(nf90_inq_varid(ncid_qs, 'snowwat_coarse', varid_qs))
    call check(nf90_inq_varid(ncid_qg, 'graupel_coarse', varid_qg))
    call check(nf90_inq_varid(ncid_omt, 'omT_coarse', varid_omt))
    call check(nf90_inq_varid(ncid_omqv, 'omqv_coarse', varid_omqv))
    call check(nf90_inq_varid(ncid_omqw, 'omql_coarse', varid_omqw))
    call check(nf90_inq_varid(ncid_omqr, 'omqr_coarse', varid_omqr))
    call check(nf90_inq_varid(ncid_omqi, 'omqi_coarse', varid_omqi))
    call check(nf90_inq_varid(ncid_omqs, 'omqs_coarse', varid_omqs))
    call check(nf90_inq_varid(ncid_omqg, 'omqg_coarse', varid_omqg))
    call check(nf90_inq_varid(ncid_pr, 'PRATEsfc_coarse', varid_pr))

    ! determine number of y-chunks
    nchunks = ny / chunk_size

    ! -----------------------------------------------------------------------
    ! Async pipeline over y-chunks.
    !
    ! For each chunk the master launches two dependent tasks:
    !   1. read_task    : fills one buffer from NetCDF under a critical
    !                     section so file reads remain serialized.
    !   2. compute_task : waits on read completion, computes work/lift, then
    !                     applies precipitation thresholds for each percentile.
    !
    ! The anti-dependency on fetch_ready(ibuf) prevents buffer reuse races.
    ! -----------------------------------------------------------------------

    !$omp parallel default(shared) private(t,yc,ibuf,istate,ystart,y,x,p,ip,geometric_thickness,work_acc,lift_acc)
    !$omp master

    do t = 1, nt
        print *, 'timestep', t, 'of', nt
        ! Map this timestep to its daily accumulation slot.
        iday = day_index_of_t(t)

        buf_state = 0
        buf_ystart = 0
        fetch_ready = 0

        do yc = 1, nchunks
            ibuf = mod(yc - 1, nbuf) + 1
            ystart = (yc - 1) * chunk_size + 1

            ! Wait until this buffer slot is fully free.
            do
                !$omp atomic read
                istate = buf_state(ibuf)
                if (istate == 0) exit
                if (istate == 2) then
                    !$omp atomic write
                    buf_state(ibuf) = 0
                    exit
                end if
                !$omp taskyield
            end do

            buf_ystart(ibuf) = ystart
            !$omp atomic write
            buf_state(ibuf) = 3

            ! Read task: fills buffers from NetCDF under critical I/O section.
            !$omp task depend(out: fetch_ready(ibuf)) firstprivate(ibuf, ystart, t)
                !$omp critical(netcdf_io)
                    call check(nf90_get_var(ncid_dz,    varid_dz,    dz_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_temp,  varid_temp,  temp_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omega, varid_omega, omega_buffer(:,:,:,:,ibuf), (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qv,    varid_qv,    qv_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qw,    varid_qw,    qw_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qr,    varid_qr,    qr_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qi,    varid_qi,    qi_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qs,    varid_qs,    qs_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_qg,    varid_qg,    qg_buffer(:,:,:,:,ibuf),    (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omt,   varid_omt,   omt_buffer(:,:,:,:,ibuf),   (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqv,  varid_omqv,  omqv_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqw,  varid_omqw,  omqw_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqr,  varid_omqr,  omqr_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqi,  varid_omqi,  omqi_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqs,  varid_omqs,  omqs_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_omqg,  varid_omqg,  omqg_buffer(:,:,:,:,ibuf),  (/1,ystart,1,t/), (/nx,chunk_size,nplev,1/)))
                    call check(nf90_get_var(ncid_pr,    varid_pr,    pr_buffer(:,:,:,ibuf),       (/1,ystart,1,t/), (/nx,chunk_size,1,1/)))
                !$omp end critical(netcdf_io)
            !$omp end task

            ! Compute task: waits for read completion and applies thresholds.
            !$omp task depend(in: fetch_ready(ibuf)) firstprivate(ibuf, ystart, t, iday) &
            !$omp&    private(y,x,p,ip,geometric_thickness,work_acc,lift_acc)

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

                    ! Threshold masking per percentile rank.
                    do ip = 1, npercentiles
                        if (pr_buffer(x,y,1,ibuf) >= prate_thresholds(ip)) then
                            work_out_buffer(ip,x,y,1,ibuf) = work_acc
                            lift_out_buffer(ip,x,y,1,ibuf) = lift_acc
                            count_out_buffer(ip,x,y,1,ibuf) = 1.0d0
                        else
                            work_out_buffer(ip,x,y,1,ibuf) = 0.0d0
                            lift_out_buffer(ip,x,y,1,ibuf) = 0.0d0
                            count_out_buffer(ip,x,y,1,ibuf) = 0.0d0
                        end if
                    end do
                end do
            end do

            ! Accumulate into daily bins. Multiple timesteps in the same day
            ! are summed into the same iday slice.
            work_accum_full(:,:,ystart:ystart+chunk_size-1,iday) = work_accum_full(:,:,ystart:ystart+chunk_size-1,iday) + work_out_buffer(:,:,:,1,ibuf)
            lift_accum_full(:,:,ystart:ystart+chunk_size-1,iday) = lift_accum_full(:,:,ystart:ystart+chunk_size-1,iday) + lift_out_buffer(:,:,:,1,ibuf)
            count_accum_full(:,:,ystart:ystart+chunk_size-1,iday) = count_accum_full(:,:,ystart:ystart+chunk_size-1,iday) + count_out_buffer(:,:,:,1,ibuf)

            !$omp atomic write
            buf_state(ibuf) = 2

            !$omp end task
        end do

        !$omp taskwait
        do ibuf = 1, nbuf
            buf_state(ibuf) = 0
        end do
    end do

    !$omp end master
    !$omp end parallel

    ! Compute conditional means; keep count at 0.0 for no-event cases.
    do iday = 1, ndays
        do y = 1, ny
            do x = 1, nx
                do ip = 1, npercentiles
                    if (count_accum_full(ip,x,y,iday) > 0.0d0) then
                        work_accum_full(ip,x,y,iday) = work_accum_full(ip,x,y,iday) / count_accum_full(ip,x,y,iday)
                        lift_accum_full(ip,x,y,iday) = lift_accum_full(ip,x,y,iday) / count_accum_full(ip,x,y,iday)
                    else
                        work_accum_full(ip,x,y,iday) = -9999.0d0
                        lift_accum_full(ip,x,y,iday) = -9999.0d0
                    end if
                end do
            end do
        end do
    end do

    call check(nf90_put_var(ncid_work_out, varid_work_out, work_accum_full(:,:,:,1:ndays), start=(/1,1,1,1/), count=(/npercentiles,nx,ny,ndays/)))
    call check(nf90_put_var(ncid_work_out, varid_lift_out, lift_accum_full(:,:,:,1:ndays), start=(/1,1,1,1/), count=(/npercentiles,nx,ny,ndays/)))
    call check(nf90_put_var(ncid_work_out, varid_count_out, count_accum_full(:,:,:,1:ndays), start=(/1,1,1,1/), count=(/npercentiles,nx,ny,ndays/)))

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
    call check(nf90_close(ncid_pr))
    call check(nf90_close(ncid_work_out))

contains

    subroutine read_thresholds_ascii(path_in, percentiles, thresholds)
        character(len=*), intent(in) :: path_in
        real(8), allocatable, intent(out) :: percentiles(:)
        real(8), allocatable, intent(out) :: thresholds(:)

        integer :: unit_in, io, n, i
        character(len=256) :: line
        real(8) :: pval, tval

        n = 0
        open(newunit=unit_in, file=trim(path_in), status='old', action='read', iostat=io, iomsg=msg)
        if (io /= 0) then
            write(error_unit,*) 'Failed to open threshold file, iomsg='//trim(msg)
            stop 1
        end if

        do
            read(unit_in, '(A)', iostat=io) line
            if (io /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            read(line, *, iostat=io) pval, tval
            if (io /= 0) cycle
            n = n + 1
        end do

        if (n <= 0) then
            close(unit_in)
            error stop 'Threshold file contains no valid data rows.'
        end if

        rewind(unit_in)
        allocate(percentiles(n), thresholds(n))
        i = 0
        do
            read(unit_in, '(A)', iostat=io) line
            if (io /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            read(line, *, iostat=io) pval, tval
            if (io /= 0) cycle
            i = i + 1
            percentiles(i) = pval
            thresholds(i) = tval
        end do

        close(unit_in)
    end subroutine read_thresholds_ascii

    function percentile_attr_name(p) result(name)
        real(8), intent(in) :: p
        character(len=64) :: name
        character(len=32) :: pstr
        character(len=32) :: packed
        integer :: i, j, last
        real(8) :: pct

        pct = 100.0d0 * p
        write(pstr, '(F10.4)') pct
        pstr = adjustl(pstr)

        last = len_trim(pstr)
        do while (last > 1 .and. pstr(last:last) == '0')
            last = last - 1
        end do
        if (pstr(last:last) == '.') last = last - 1

        packed = ''
        j = 0
        do i = 1, last
            if (pstr(i:i) /= '.') then
                j = j + 1
                packed(j:j) = pstr(i:i)
            end if
        end do

        name = 'p'//trim(packed)//'_threshold'
    end function percentile_attr_name

    ! Build day grouping metadata from raw time values.
    ! day_of_t(i) gives the daily index for timestep i,
    ! and day_times(k) stores the first timestamp encountered for day k.
    subroutine build_day_groups(time_in, units, day_of_t, day_first, day_times, ndays)
        real(8), intent(in) :: time_in(:)
        character(len=*), intent(in) :: units
        integer, allocatable, intent(out) :: day_of_t(:), day_first(:)
        real(8), allocatable, intent(out) :: day_times(:)
        integer, intent(out) :: ndays

        integer :: i
        integer, allocatable :: day_key(:)
        real(8) :: scale_to_days

        ! Convert native time units to days before integer day bucketing.
        scale_to_days = units_to_days(units)
        allocate(day_key(size(time_in)))
        do i = 1, size(time_in)
            day_key(i) = floor(time_in(i) * scale_to_days + 1.0d-9)
        end do

        ndays = 1
        do i = 2, size(day_key)
            if (day_key(i) /= day_key(i-1)) ndays = ndays + 1
        end do

        allocate(day_of_t(size(time_in)))
        allocate(day_first(ndays))
        allocate(day_times(ndays))

        ndays = 1
        day_of_t(1) = 1
        day_first(1) = 1
        day_times(1) = time_in(1)
        do i = 2, size(day_key)
            if (day_key(i) /= day_key(i-1)) then
                ndays = ndays + 1
                day_first(ndays) = i
                day_times(ndays) = time_in(i)
            end if
            day_of_t(i) = ndays
        end do
    end subroutine build_day_groups

    ! Return multiplicative factor that converts values in `units` to days.
    ! Unknown prefixes default to 1 (already in days).
    real(8) function units_to_days(units)
        character(len=*), intent(in) :: units
        character(len=255) :: lower
        integer :: i, n

        lower = ''
        n = len_trim(units)
        do i = 1, n
            if (units(i:i) >= 'A' .and. units(i:i) <= 'Z') then
                lower(i:i) = achar(iachar(units(i:i)) + 32)
            else
                lower(i:i) = units(i:i)
            end if
        end do

        units_to_days = 1.0d0
        if (index(lower, 'hour') == 1) then
            units_to_days = 1.0d0 / 24.0d0
        else if (index(lower, 'minute') == 1) then
            units_to_days = 1.0d0 / 1440.0d0
        else if (index(lower, 'second') == 1) then
            units_to_days = 1.0d0 / 86400.0d0
        end if
    end function units_to_days

end program compute_work_async_prate_threshold
