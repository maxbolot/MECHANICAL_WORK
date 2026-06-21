! histogram-only async fork of compute_work_async.f90
!
! Outputs monthly histogram products only:
!   hist_area(time, nbin_pr, nlat)
!   hist_work(time, nbin_pr, nlat)
!   hist_lift(time, nbin_pr, nlat)
!   hist2d_work(time, nbin_work, nbin_pr, nlat)
!   hist2d_lift(time, nbin_work, nbin_pr, nlat)
!
! Latitude bands are configurable via namelist.  Precipitation bins reuse the
! piecewise-log definition from compute_prate_thresholds_by_lat_band.f90.
! Aggregation grouping (daily or monthly) is based on the NetCDF time
! coordinate and supports standard, gregorian, and julian calendars.

program compute_work_async_histograms_by_lat_band

    use omp_lib
    use netcdf
    use interp
    use nc
    use hist
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg
    character(len=255) :: path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg
    character(len=255) :: path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, path_pr
    character(len=255) :: path_hist_out
    character(len=255) :: date_list_file
    character(len=255) :: target_period_key
    character(len=255) :: time_units, time_calendar
    character(len=255) :: source_dir_current
    character(len=255) :: source_date_current
    character(len=255) :: self_date, self_root
    character(len=512) :: incident_log_path
    character(len=512) :: skip_detail

    integer, parameter :: chunk_size = 144
    integer, parameter :: nbuf = 2
    integer, parameter :: npr_edges = 1201
    integer, parameter :: nwork_edges = 5001
    integer, parameter :: max_lat_band_bounds = 181
    integer, parameter :: max_segments_per_band = 8

    double precision, parameter :: fill_value = -9999.0d0

    double precision, dimension(:,:,:,:,:), allocatable :: dz_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: temp_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omega_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qv_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qw_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qr_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qi_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qs_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: qg_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omt_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqv_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqw_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqr_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqi_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqs_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: omqg_buffer
    double precision, dimension(:,:,:,:), allocatable :: pr_buffer
    double precision, dimension(:,:,:,:), allocatable :: work_buffer
    double precision, dimension(:,:,:,:), allocatable :: lift_buffer

    double precision, dimension(:), allocatable :: lon
    double precision, dimension(:), allocatable :: lat
    double precision, dimension(:,:), allocatable :: cell_area
    double precision, dimension(:), allocatable :: time_vals
    double precision, dimension(:), allocatable :: month_time_vals
    integer, dimension(:), allocatable :: month_index_of_t
    integer, dimension(:), allocatable :: month_first_idx
    integer, dimension(:), allocatable :: month_step_count
    integer :: nmonths

    integer :: nx, ny, nplev, nt
    integer :: nt_src
    integer :: nt_read_src
    integer :: nlat_bands
    integer :: lat_band_bounds_count
    logical :: use_custom_lat_band_bounds
    logical :: use_custom_lat_band_segments
    double precision :: lat_band_bounds(max_lat_band_bounds)
    integer :: lat_band_segment_count(max_lat_band_bounds - 1)
    double precision :: lat_band_segment_south(max_segments_per_band, max_lat_band_bounds - 1)
    double precision :: lat_band_segment_north(max_segments_per_band, max_lat_band_bounds - 1)
    double precision :: lat_band_south(max_lat_band_bounds - 1)
    double precision :: lat_band_north(max_lat_band_bounds - 1)
    integer :: lat_band_start(max_segments_per_band, max_lat_band_bounds - 1)
    integer :: lat_band_end(max_segments_per_band, max_lat_band_bounds - 1)

    integer :: t, p, x, y, yc, ibuf, istate, iday, ilat, iseg
    integer :: nchunks, ystart
    integer :: nx_chunk, ny_chunk, y_local_start, y_local_end
    integer :: hist_y_start, hist_y_end
    integer :: buf_state(nbuf), buf_ystart(nbuf)
    integer :: fetch_ready(nbuf)
    logical :: found_lat_start, found_lat_end, found_lon_start, found_lon_end
    integer :: tmp_unit, ios
    integer :: ncid_dz, ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg
    integer :: ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr
    integer :: ncid_hist_out
    integer :: tmp_dimids(4)
    integer :: ncstatus
    integer :: varid_dz, varid_temp, varid_omega, varid_qv, varid_qw, varid_qr, varid_qi, varid_qs, varid_qg
    integer :: varid_omt, varid_omqv, varid_omqw, varid_omqr, varid_omqi, varid_omqs, varid_omqg, varid_pr
    integer :: nt_omega_src, nt_qv_src, nt_qw_src, nt_qr_src, nt_qi_src, nt_qs_src, nt_qg_src
    integer :: nt_omt_src, nt_omqv_src, nt_omqw_src, nt_omqr_src, nt_omqi_src, nt_omqs_src, nt_omqg_src
    integer :: nx_var, ny_var, nz_var, nt_var
    integer :: varid_lon, varid_lat, varid_time
    integer :: varid_out_lon, varid_out_lat, varid_out_time, varid_lat_band
    integer :: varid_lat_band_south, varid_lat_band_north
    integer :: varid_lat_band_segment_count
    integer :: varid_lat_band_segment_south, varid_lat_band_segment_north
    integer :: varid_hist_area, varid_hist_work, varid_hist_lift
    integer :: varid_hist2d_work, varid_hist2d_lift
    integer :: dimid_out_lon, dimid_out_lat, dimid_out_time, dimid_out_pr, dimid_out_work
    integer :: dimid_out_band
    integer :: dimid_nbin_pr, dimid_nbin_work, dimid_nedges_pr, dimid_nedges_work
    integer :: dimid_hist_time
    integer :: dimid_hist_nlat
    integer :: dimid_hist_nsegment
    integer :: iedge
    double precision :: pr_edges(npr_edges)
    double precision :: work_edges(nwork_edges)
    double precision :: geometric_thickness, work_acc, lift_acc
    double precision :: month_norm
    double precision :: lat_south, lat_north
    character(len=64) :: attr_name
    character(len=32) :: lat_mode
    character(len=5) :: continue_mode_attr
    character(len=32) :: aggregation_mode
    character(len=255) :: lower_units, ref_string
    integer :: ref_year, ref_month, ref_day
    integer :: ref_hour, ref_minute, ref_second
    double precision :: scale_to_days
    integer :: month_key_prev, month_key_curr
    integer :: ymd_year, ymd_month, ymd_day
    integer :: band_n
    integer :: isrc, n_sources, n_targets
    integer :: list_unit, line_count
    integer :: incident_log_unit
    integer :: sources_processed, sources_skipped
    integer :: got_nx, got_ny, got_nz, got_nt
    logical :: use_period_filter
    logical :: continue_on_source_error
    logical :: first_time_set
    logical :: matches_period
    logical :: band_has_data

    character(len=255), allocatable :: source_dates(:)
    character(len=255), allocatable :: source_roots(:)

    double precision, dimension(:,:,:), allocatable :: hist_area_out
    double precision, dimension(:,:,:), allocatable :: hist_work_out
    double precision, dimension(:,:,:), allocatable :: hist_lift_out
    double precision, dimension(:,:,:,:), allocatable :: hist2d_work_out
    double precision, dimension(:,:,:,:), allocatable :: hist2d_lift_out

    double precision, dimension(:), allocatable :: hist_area_chunk
    double precision, dimension(:), allocatable :: hist_work_chunk
    double precision, dimension(:), allocatable :: hist_lift_chunk
    double precision, dimension(:,:), allocatable :: hist2d_work_chunk
    double precision, dimension(:,:), allocatable :: hist2d_lift_chunk

    namelist /config/ path_dz, path_temp, path_omega, path_qv, path_qw, path_qr, path_qi, path_qs, path_qg, &
                     path_omt, path_omqv, path_omqw, path_omqr, path_omqi, path_omqs, path_omqg, path_pr, &
                     path_hist_out, nlat_bands, lat_band_bounds_count, lat_band_bounds, use_custom_lat_band_bounds, aggregation_mode, &
                     date_list_file, target_period_key, use_custom_lat_band_segments, lat_band_segment_count, &
                     lat_band_segment_south, lat_band_segment_north, continue_on_source_error

    path_dz = ''
    path_temp = ''
    path_omega = ''
    path_qv = ''
    path_qw = ''
    path_qr = ''
    path_qi = ''
    path_qs = ''
    path_qg = ''
    path_omt = ''
    path_omqv = ''
    path_omqw = ''
    path_omqr = ''
    path_omqi = ''
    path_omqs = ''
    path_omqg = ''
    path_pr = ''
    path_hist_out = ''
    date_list_file = ''
    target_period_key = ''
    nlat_bands = 18
    lat_band_bounds_count = 0
    use_custom_lat_band_bounds = .false.
    use_custom_lat_band_segments = .false.
    continue_on_source_error = .true.
    aggregation_mode = 'monthly'
    lat_band_bounds = 1.0d99
    lat_band_segment_count = 0
    lat_band_segment_south = 1.0d99
    lat_band_segment_north = 1.0d99

    ! Read runtime configuration from the namelist file passed on argv(1).
    call get_command_argument(1, filenml)
    open(newunit=tmp_unit, file=trim(adjustl(filenml)), status='old', iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to open configuration namelist, iomsg='//trim(msg)
        stop 1
    end if
    read(unit=tmp_unit, nml=config, iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to read configuration namelist, iomsg='//trim(msg)
        stop 1
    end if
    close(tmp_unit)

    if (len_trim(path_hist_out) == 0) then
        write(error_unit,*) 'path_hist_out is required.'
        stop 1
    end if

    use_period_filter = (len_trim(target_period_key) > 0)

    if (len_trim(date_list_file) > 0) then
        open(newunit=list_unit, file=trim(adjustl(date_list_file)), status='old', iostat=ios, iomsg=msg)
        if (ios /= 0) then
            write(error_unit,*) 'Failed to open date_list_file, iomsg='//trim(msg)
            stop 1
        end if

        line_count = 0
        do
            read(list_unit, '(A)', iostat=ios) msg
            if (ios /= 0) exit
            if (len_trim(msg) == 0) cycle
            line_count = line_count + 1
        end do
        rewind(list_unit)

        if (line_count < 1) then
            write(error_unit,*) 'date_list_file is empty: ', trim(date_list_file)
            stop 1
        end if

        allocate(source_dates(line_count))
        allocate(source_roots(line_count))
        n_sources = 0
        do
            read(list_unit, '(A)', iostat=ios) msg
            if (ios /= 0) exit
            if (len_trim(msg) == 0) cycle
            n_sources = n_sources + 1
            call split_date_root_line(trim(msg), source_dates(n_sources), source_roots(n_sources))
        end do
        close(list_unit)

        if (.not. use_period_filter) then
            write(error_unit,*) 'date_list_file mode requires target_period_key for one-pass period aggregation.'
            stop 1
        end if

        if (len_trim(path_dz) == 0 .or. len_trim(path_temp) == 0 .or. len_trim(path_omega) == 0 .or. len_trim(path_qv) == 0 .or. &
            len_trim(path_qw) == 0 .or. len_trim(path_qr) == 0 .or. len_trim(path_qi) == 0 .or. len_trim(path_qs) == 0 .or. &
            len_trim(path_qg) == 0 .or. len_trim(path_omt) == 0 .or. len_trim(path_omqv) == 0 .or. len_trim(path_omqw) == 0 .or. &
            len_trim(path_omqr) == 0 .or. len_trim(path_omqi) == 0 .or. len_trim(path_omqs) == 0 .or. len_trim(path_omqg) == 0 .or. &
            len_trim(path_pr) == 0) then
            source_dir_current = trim(source_roots(1))//'/'//trim(source_dates(1))
            path_dz = trim(source_dir_current)//'/DZ_C3072_1440x720.fre.nc'
            path_temp = trim(source_dir_current)//'/temp_coarse_C3072_1440x720.fre.nc'
            path_omega = trim(source_dir_current)//'/ptend_coarse_C3072_1440x720.fre.nc'
            path_qv = trim(source_dir_current)//'/sphum_coarse_C3072_1440x720.fre.nc'
            path_qw = trim(source_dir_current)//'/liq_wat_coarse_C3072_1440x720.fre.nc'
            path_qr = trim(source_dir_current)//'/rainwat_coarse_C3072_1440x720.fre.nc'
            path_qi = trim(source_dir_current)//'/ice_wat_coarse_C3072_1440x720.fre.nc'
            path_qs = trim(source_dir_current)//'/snowwat_coarse_C3072_1440x720.fre.nc'
            path_qg = trim(source_dir_current)//'/graupel_coarse_C3072_1440x720.fre.nc'
            path_omt = trim(source_dir_current)//'/omT_coarse_C3072_1440x720.fre.nc'
            path_omqv = trim(source_dir_current)//'/omqv_coarse_C3072_1440x720.fre.nc'
            path_omqw = trim(source_dir_current)//'/omql_coarse_C3072_1440x720.fre.nc'
            path_omqr = trim(source_dir_current)//'/omqr_coarse_C3072_1440x720.fre.nc'
            path_omqi = trim(source_dir_current)//'/omqi_coarse_C3072_1440x720.fre.nc'
            path_omqs = trim(source_dir_current)//'/omqs_coarse_C3072_1440x720.fre.nc'
            path_omqg = trim(source_dir_current)//'/omqg_coarse_C3072_1440x720.fre.nc'
            path_pr = trim(source_dir_current)//'/PRATEsfc_coarse_C3072_1440x720.fre.nc'
        end if
    else
        n_sources = 1
    end if

    write(*,'(A)') '=== compute_work_async_histograms_by_lat_band startup ==='
    write(*,'(A)') 'config_file='//trim(filenml)
    write(*,'(A)') 'aggregation_mode='//trim(to_lower(aggregation_mode))
    if (use_period_filter) then
        write(*,'(A)') 'target_period_key='//trim(target_period_key)
    else
        write(*,'(A)') 'target_period_key=<none>'
    end if
    if (len_trim(date_list_file) > 0) then
        write(*,'(A)') 'date_list_file='//trim(date_list_file)
    else
        write(*,'(A)') 'date_list_file=<none>'
    end if
    write(*,'(A,I0)') 'source_count=', n_sources
    write(*,'(A)') 'path_hist_out='//trim(path_hist_out)
    write(*,'(A,I0)') 'nlat_bands=', nlat_bands
    write(*,'(A)') '========================================================='

    if (nlat_bands <= 0) then
        write(error_unit,*) 'nlat_bands must be positive.'
        stop 1
    end if
    if (nlat_bands > max_lat_band_bounds - 1) then
        write(error_unit,*) 'nlat_bands exceeds supported maximum.'
        stop 1
    end if

    call check(nf90_open(trim(adjustl(path_dz)), nf90_nowrite, ncid_dz))
    call check(nf90_inq_varid(ncid_dz, 'DZ', varid_dz))
    call check(nf90_inquire_variable(ncid_dz, varid_dz, dimids=tmp_dimids))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(1), len=nx))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(2), len=ny))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(3), len=nplev))
    call check(nf90_inquire_dimension(ncid_dz, tmp_dimids(4), len=nt))

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
    allocate(pr_buffer(nx, chunk_size, 1, nbuf))
    allocate(work_buffer(nx, chunk_size, 1, nbuf))
    allocate(lift_buffer(nx, chunk_size, 1, nbuf))
    allocate(lon(nx))
    allocate(lat(ny))
    allocate(time_vals(nt))

    if (use_custom_lat_band_bounds .and. use_custom_lat_band_segments) then
        write(error_unit,*) 'Choose either use_custom_lat_band_bounds or use_custom_lat_band_segments, not both.'
        stop 1
    end if

    if (use_custom_lat_band_segments) then
        do ilat = 1, nlat_bands
            if (lat_band_segment_count(ilat) <= 0 .or. lat_band_segment_count(ilat) > max_segments_per_band) then
                write(error_unit,*) 'lat_band_segment_count out of range for band ', ilat
                stop 1
            end if
            do iseg = 1, lat_band_segment_count(ilat)
                if (lat_band_segment_south(iseg, ilat) >= lat_band_segment_north(iseg, ilat)) then
                    write(error_unit,*) 'lat band segment has south >= north for band/segment ', ilat, iseg
                    stop 1
                end if
                if (lat_band_segment_south(iseg, ilat) < -90.0d0 .or. lat_band_segment_north(iseg, ilat) > 90.0d0) then
                    write(error_unit,*) 'lat band segment bounds must be inside [-90,90] for band/segment ', ilat, iseg
                    stop 1
                end if
            end do
        end do
    else if (use_custom_lat_band_bounds) then
        if (lat_band_bounds_count /= nlat_bands + 1) then
            write(error_unit,*) 'lat_band_bounds_count must equal nlat_bands + 1 when custom bounds are enabled.'
            stop 1
        end if
    end if

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
    lower_units = to_lower(trim(time_units))
    call parse_time_units(trim(time_units), ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second, scale_to_days)
    if (index(lower_units, 'standard') == 0 .and. index(lower_units, 'gregorian') == 0 .and. index(lower_units, 'hour') == 0 .and. index(lower_units, 'day') == 0 .and. index(lower_units, 'minute') == 0 .and. index(lower_units, 'second') == 0) then
        ! Units are still acceptable if parse_time_units succeeded, but keep the
        ! message explicit for unsupported calendars.
    end if
    if (.not. is_supported_calendar(time_calendar)) then
        write(error_unit,*) 'Unsupported calendar: ', trim(time_calendar)
        stop 1
    end if

    if (use_period_filter) then
        nmonths = 1
        allocate(month_index_of_t(nt))
        allocate(month_first_idx(1))
        allocate(month_time_vals(1))
        allocate(month_step_count(1))
        month_index_of_t = 0
        month_first_idx(1) = 1
        month_time_vals(1) = 0.0d0
        month_step_count(1) = 0
    else
        select case (trim(to_lower(aggregation_mode)))
        case ('monthly')
            ! Map each native timestep to a month index and track per-month sample count.
            call build_month_groups(time_vals, time_units, time_calendar, month_index_of_t, month_first_idx, month_time_vals, month_step_count, nmonths)
        case ('daily')
            ! Map each native timestep to a day index and track per-day sample count.
            call build_day_groups(time_vals, time_units, time_calendar, month_index_of_t, month_first_idx, month_time_vals, month_step_count, nmonths)
        case default
            write(error_unit,*) 'Unsupported aggregation_mode: ', trim(aggregation_mode)
            stop 1
        end select
    end if
    if (nmonths <= 0) then
        write(error_unit,*) 'No month groups were created.'
        stop 1
    end if
    call check(nf90_close(ncid_temp))

    allocate(hist_area_out(nmonths, npr_edges - 1, nlat_bands))
    allocate(hist_work_out(nmonths, npr_edges - 1, nlat_bands))
    allocate(hist_lift_out(nmonths, npr_edges - 1, nlat_bands))
    allocate(hist2d_work_out(nmonths, nwork_edges - 1, npr_edges - 1, nlat_bands))
    allocate(hist2d_lift_out(nmonths, nwork_edges - 1, npr_edges - 1, nlat_bands))
    hist_area_out = 0.0d0
    hist_work_out = 0.0d0
    hist_lift_out = 0.0d0
    hist2d_work_out = 0.0d0
    hist2d_lift_out = 0.0d0

    call compute_cell_areas(lon, lat, cell_area)

    call resolve_lat_band_segments(nlat_bands, lat_band_bounds_count, use_custom_lat_band_bounds, lat_band_bounds, &
                                   use_custom_lat_band_segments, lat_band_segment_count, lat_band_segment_south, lat_band_segment_north, &
                                   lat_band_south, lat_band_north)
    do ilat = 1, nlat_bands
        do iseg = 1, lat_band_segment_count(ilat)
            call resolve_lat_band_indices(lat, lat_band_segment_south(iseg, ilat), lat_band_segment_north(iseg, ilat), lat_band_start(iseg, ilat), lat_band_end(iseg, ilat))
        end do
    end do

    first_time_set = .false.
    sources_processed = 0
    sources_skipped = 0
    incident_log_unit = -1
    incident_log_path = trim(path_hist_out)//'.skipped_sources.tsv'
    open(newunit=incident_log_unit, file=trim(incident_log_path), status='replace', action='write', iostat=ios, iomsg=msg)
    if (ios == 0) then
        write(incident_log_unit,'(A)') 'source_date'//char(9)//'stage'//char(9)//'variable'//char(9)//'detail'
    else
        incident_log_unit = -1
        write(error_unit,*) 'Warning: failed to open incident log file, iomsg='//trim(msg)
    end if

    source_loop: do isrc = 1, n_sources
        source_date_current = 'single_source'
        ncid_temp = -1
        ncid_omega = -1
        ncid_qv = -1
        ncid_qw = -1
        ncid_qr = -1
        ncid_qi = -1
        ncid_qs = -1
        ncid_qg = -1
        ncid_omt = -1
        ncid_omqv = -1
        ncid_omqw = -1
        ncid_omqr = -1
        ncid_omqi = -1
        ncid_omqs = -1
        ncid_omqg = -1
        ncid_pr = -1
        if (len_trim(date_list_file) > 0) then
            source_dir_current = trim(source_roots(isrc))//'/'//trim(source_dates(isrc))
            source_date_current = trim(source_dates(isrc))
            path_dz = trim(source_dir_current)//'/DZ_C3072_1440x720.fre.nc'
            path_temp = trim(source_dir_current)//'/temp_coarse_C3072_1440x720.fre.nc'
            path_omega = trim(source_dir_current)//'/ptend_coarse_C3072_1440x720.fre.nc'
            path_qv = trim(source_dir_current)//'/sphum_coarse_C3072_1440x720.fre.nc'
            path_qw = trim(source_dir_current)//'/liq_wat_coarse_C3072_1440x720.fre.nc'
            path_qr = trim(source_dir_current)//'/rainwat_coarse_C3072_1440x720.fre.nc'
            path_qi = trim(source_dir_current)//'/ice_wat_coarse_C3072_1440x720.fre.nc'
            path_qs = trim(source_dir_current)//'/snowwat_coarse_C3072_1440x720.fre.nc'
            path_qg = trim(source_dir_current)//'/graupel_coarse_C3072_1440x720.fre.nc'
            path_omt = trim(source_dir_current)//'/omT_coarse_C3072_1440x720.fre.nc'
            path_omqv = trim(source_dir_current)//'/omqv_coarse_C3072_1440x720.fre.nc'
            path_omqw = trim(source_dir_current)//'/omql_coarse_C3072_1440x720.fre.nc'
            path_omqr = trim(source_dir_current)//'/omqr_coarse_C3072_1440x720.fre.nc'
            path_omqi = trim(source_dir_current)//'/omqi_coarse_C3072_1440x720.fre.nc'
            path_omqs = trim(source_dir_current)//'/omqs_coarse_C3072_1440x720.fre.nc'
            path_omqg = trim(source_dir_current)//'/omqg_coarse_C3072_1440x720.fre.nc'
            path_pr = trim(source_dir_current)//'/PRATEsfc_coarse_C3072_1440x720.fre.nc'
            write(*,'(A,I0,A,A)') 'source ', isrc, ': ', trim(source_dates(isrc))
        end if

        ! In date-list mode, reopen DZ for each source date so the time
        ! dimension follows the current source cadence.
        if (ncid_dz > 0) call check(nf90_close(ncid_dz))
        ncid_dz = -1
        ncstatus = nf90_open(trim(adjustl(path_dz)), nf90_nowrite, ncid_dz)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'DZ', trim(nf90_strerror(ncstatus)))
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        call check(nf90_inq_varid(ncid_dz, 'DZ', varid_dz))

        ncstatus = nf90_open(trim(adjustl(path_temp)), nf90_nowrite, ncid_temp)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'temp_coarse', trim(nf90_strerror(ncstatus)))
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        call check(nf90_inq_varid(ncid_temp, 'temp_coarse', varid_temp))
        call check(nf90_inquire_variable(ncid_temp, varid_temp, dimids=tmp_dimids))
        call check(nf90_inquire_dimension(ncid_temp, tmp_dimids(4), len=nt_src))

        ! Date-list mode may mix sources with different record counts. Use
        ! source-local nt to avoid out-of-bounds NetCDF reads.
        if (.not. use_period_filter .and. nt_src /= nt) then
            write(error_unit,*) 'In non-period mode, all sources must have the same nt. expected/got=', nt, nt_src
            stop 1
        end if
        if (nt_src > size(time_vals)) then
            deallocate(time_vals)
            allocate(time_vals(nt_src))
        end if

        time_units = 'hours since 1900-01-01 00:00:00'
        time_calendar = 'standard'
        do t = 1, nt_src
            time_vals(t) = dble(t - 1)
        end do
        ncstatus = nf90_inq_varid(ncid_temp, 'time', varid_time)
        if (ncstatus == nf90_noerr) then
            call check(nf90_get_var(ncid_temp, varid_time, time_vals(1:nt_src)))
            ncstatus = nf90_get_att(ncid_temp, varid_time, 'units', time_units)
            if (ncstatus /= nf90_noerr) time_units = 'hours since 1900-01-01 00:00:00'
            ncstatus = nf90_get_att(ncid_temp, varid_time, 'calendar', time_calendar)
            if (ncstatus /= nf90_noerr) time_calendar = 'standard'
        end if

        ncstatus = nf90_open(trim(adjustl(path_omega)), nf90_nowrite, ncid_omega)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'ptend_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qv)), nf90_nowrite, ncid_qv)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'sphum_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qw)), nf90_nowrite, ncid_qw)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'liq_wat_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qr)), nf90_nowrite, ncid_qr)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'rainwat_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qi)), nf90_nowrite, ncid_qi)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'ice_wat_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qs)), nf90_nowrite, ncid_qs)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'snowwat_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_qg)), nf90_nowrite, ncid_qg)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'graupel_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omt)), nf90_nowrite, ncid_omt)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omT_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqv)), nf90_nowrite, ncid_omqv)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omqv_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqw)), nf90_nowrite, ncid_omqw)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omql_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqr)), nf90_nowrite, ncid_omqr)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omqr_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqi)), nf90_nowrite, ncid_omqi)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omqi_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqs)), nf90_nowrite, ncid_omqs)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omqs_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_omqg)), nf90_nowrite, ncid_omqg)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'omqg_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        ncstatus = nf90_open(trim(adjustl(path_pr)), nf90_nowrite, ncid_pr)
        if (ncstatus /= nf90_noerr) then
            call log_skip_incident(incident_log_unit, source_date_current, 'open', 'PRATEsfc_coarse', trim(nf90_strerror(ncstatus)))
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if

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

        if (.not. check_var_shape(ncid_dz, varid_dz, 'DZ', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'DZ', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_temp, varid_temp, 'temp_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'temp_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omega, varid_omega, 'ptend_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'ptend_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qv, varid_qv, 'sphum_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'sphum_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qw, varid_qw, 'liq_wat_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'liq_wat_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qr, varid_qr, 'rainwat_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'rainwat_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qi, varid_qi, 'ice_wat_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'ice_wat_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qs, varid_qs, 'snowwat_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'snowwat_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_qg, varid_qg, 'graupel_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'graupel_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omt, varid_omt, 'omT_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omT_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqv, varid_omqv, 'omqv_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omqv_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqw, varid_omqw, 'omql_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omql_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqr, varid_omqr, 'omqr_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omqr_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqi, varid_omqi, 'omqi_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omqi_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqs, varid_omqs, 'omqs_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omqs_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_omqg, varid_omqg, 'omqg_coarse', nx, ny, nplev, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, nplev, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'omqg_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. check_var_shape(ncid_pr, varid_pr, 'PRATEsfc_coarse', nx, ny, 1, nt_src, source_date_current, got_nx, got_ny, got_nz, got_nt)) then
            write(skip_detail,'(A,4(I0,1X),A,4(I0,1X))') 'expected=', nx, ny, 1, nt_src, ' got=', got_nx, got_ny, got_nz, got_nt
            call log_skip_incident(incident_log_unit, source_date_current, 'shape', 'PRATEsfc_coarse', skip_detail)
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if

        call get_var_time_len(ncid_omega, varid_omega, nt_omega_src)
        call get_var_time_len(ncid_qv, varid_qv, nt_qv_src)
        call get_var_time_len(ncid_qw, varid_qw, nt_qw_src)
        call get_var_time_len(ncid_qr, varid_qr, nt_qr_src)
        call get_var_time_len(ncid_qi, varid_qi, nt_qi_src)
        call get_var_time_len(ncid_qs, varid_qs, nt_qs_src)
        call get_var_time_len(ncid_qg, varid_qg, nt_qg_src)
        call get_var_time_len(ncid_omt, varid_omt, nt_omt_src)
        call get_var_time_len(ncid_omqv, varid_omqv, nt_omqv_src)
        call get_var_time_len(ncid_omqw, varid_omqw, nt_omqw_src)
        call get_var_time_len(ncid_omqr, varid_omqr, nt_omqr_src)
        call get_var_time_len(ncid_omqi, varid_omqi, nt_omqi_src)
        call get_var_time_len(ncid_omqs, varid_omqs, nt_omqs_src)
        call get_var_time_len(ncid_omqg, varid_omqg, nt_omqg_src)

        nt_read_src = nt_src
        nt_read_src = min(nt_read_src, nt_omega_src)
        nt_read_src = min(nt_read_src, nt_qv_src)
        nt_read_src = min(nt_read_src, nt_qw_src)
        nt_read_src = min(nt_read_src, nt_qr_src)
        nt_read_src = min(nt_read_src, nt_qi_src)
        nt_read_src = min(nt_read_src, nt_qs_src)
        nt_read_src = min(nt_read_src, nt_qg_src)
        nt_read_src = min(nt_read_src, nt_omt_src)
        nt_read_src = min(nt_read_src, nt_omqv_src)
        nt_read_src = min(nt_read_src, nt_omqw_src)
        nt_read_src = min(nt_read_src, nt_omqr_src)
        nt_read_src = min(nt_read_src, nt_omqi_src)
        nt_read_src = min(nt_read_src, nt_omqs_src)
        nt_read_src = min(nt_read_src, nt_omqg_src)

        if (use_period_filter .and. nt_read_src /= nt_src) then
            write(skip_detail,'(A, I0, A, I0)') 'using common cadence nt=', nt_read_src, ' from source nt=', nt_src
            call log_source_incident(incident_log_unit, source_date_current, 'cadence', 'source_set', skip_detail)
        end if

        if (.not. verify_time_axis_matches(ncid_omega, 'ptend_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omega_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'ptend_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qv, 'sphum_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qv_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'sphum_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qw, 'liq_wat_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qw_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'liq_wat_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qr, 'rainwat_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qr_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'rainwat_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qi, 'ice_wat_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qi_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'ice_wat_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qs, 'snowwat_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qs_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'snowwat_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_qg, 'graupel_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_qg_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'graupel_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omt, 'omT_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omt_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omT_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqv, 'omqv_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqv_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omqv_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqw, 'omql_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqw_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omql_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqr, 'omqr_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqr_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omqr_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqi, 'omqi_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqi_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omqi_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqs, 'omqs_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqs_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omqs_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if
        if (.not. verify_time_axis_matches(ncid_omqg, 'omqg_coarse', time_vals(1:nt_src), time_units, time_calendar, nt_src, nt_omqg_src, source_date_current)) then
            call log_skip_incident(incident_log_unit, source_date_current, 'time_axis', 'omqg_coarse', 'time axis mismatch')
            call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
            sources_skipped = sources_skipped + 1
            if (.not. continue_on_source_error) stop 1
            cycle source_loop
        end if

        nt_read_src = min(nt_read_src, nt_src)

    ! PR bins follow the same piecewise-log spacing used in threshold workflows.
    call init_piecewise_log_bins(pr_edges)
    do iedge = 1, nwork_edges
        work_edges(iedge) = -2500.0d0 + dble(iedge - 1)
    end do

    nchunks = ny / chunk_size
    if (mod(ny, chunk_size) /= 0) nchunks = nchunks + 1

    !$omp parallel default(shared) private(t,yc,ibuf,istate,ystart,y,x,p,work_acc,lift_acc,geometric_thickness,ilat,hist_y_start,hist_y_end,y_local_start,y_local_end,month_norm,hist_area_chunk,hist_work_chunk,hist_lift_chunk,hist2d_work_chunk,hist2d_lift_chunk,matches_period)
    !$omp master

    do t = 1, nt_read_src
        print *, 'timestep', t, 'of', nt_read_src
        if (use_period_filter) then
            call time_matches_target_period(time_vals(t), time_units, time_calendar, aggregation_mode, target_period_key, matches_period)
            if (.not. matches_period) cycle
            iday = 1
            month_step_count(1) = month_step_count(1) + 1
            if (.not. first_time_set) then
                month_time_vals(1) = time_vals(t)
                first_time_set = .true.
            end if
        else
            iday = month_index_of_t(t)
        end if
        buf_state = 0
        buf_ystart = 0
        fetch_ready = 0

        do yc = 1, nchunks
            ibuf = mod(yc - 1, nbuf) + 1
            ystart = (yc - 1) * chunk_size + 1
            if (ystart > ny) cycle

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

            print *, 'queuing y-chunk', yc, 'into pipeline (ibuf=', ibuf, ')'

            ! Read task: serialized NetCDF access fills one reusable buffer slot.
            !$omp task depend(out: fetch_ready(ibuf)) firstprivate(ibuf, ystart, t)
                !$omp critical(netcdf_io)
                    call read_chunk(ncid_dz, varid_dz, dz_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_temp, varid_temp, temp_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omega, varid_omega, omega_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qv, varid_qv, qv_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qw, varid_qw, qw_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qr, varid_qr, qr_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qi, varid_qi, qi_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qs, varid_qs, qs_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_qg, varid_qg, qg_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omt, varid_omt, omt_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqv, varid_omqv, omqv_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqw, varid_omqw, omqw_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqr, varid_omqr, omqr_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqi, varid_omqi, omqi_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqs, varid_omqs, omqs_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_omqg, varid_omqg, omqg_buffer(:,:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), nplev, t, ystart)
                    call read_chunk(ncid_pr, varid_pr, pr_buffer(:,:,:,ibuf), nx, min(chunk_size, ny - ystart + 1), 1, t, ystart)
                !$omp end critical(netcdf_io)
            !$omp end task

            ! Compute task: consumes a filled slot, accumulates local histograms,
            ! then reduces into shared monthly outputs.
            !$omp task depend(in: fetch_ready(ibuf)) firstprivate(ibuf, ystart, t, iday) &
            !$omp&    private(y,x,p,work_acc,lift_acc,geometric_thickness,ilat,iseg,hist_y_start,hist_y_end,y_local_start,y_local_end,band_has_data,month_norm,hist_area_chunk,hist_work_chunk,hist_lift_chunk,hist2d_work_chunk,hist2d_lift_chunk)
                allocate(hist_area_chunk(npr_edges - 1))
                allocate(hist_work_chunk(npr_edges - 1))
                allocate(hist_lift_chunk(npr_edges - 1))
                allocate(hist2d_work_chunk(npr_edges - 1, nwork_edges - 1))
                allocate(hist2d_lift_chunk(npr_edges - 1, nwork_edges - 1))

                do ilat = 1, nlat_bands
                    hist_area_chunk = 0.0d0
                    hist_work_chunk = 0.0d0
                    hist_lift_chunk = 0.0d0
                    hist2d_work_chunk = 0.0d0
                    hist2d_lift_chunk = 0.0d0

                    band_has_data = .false.
                    do iseg = 1, lat_band_segment_count(ilat)
                        if (lat_band_start(iseg, ilat) > lat_band_end(iseg, ilat)) cycle
                        hist_y_start = max(ystart, lat_band_start(iseg, ilat))
                        hist_y_end = min(ystart + min(chunk_size, ny - ystart + 1) - 1, lat_band_end(iseg, ilat))
                        if (hist_y_start > hist_y_end) cycle

                        band_has_data = .true.
                        y_local_start = hist_y_start - ystart + 1
                        y_local_end = hist_y_end - ystart + 1

                        do y = y_local_start, y_local_end
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
                                work_buffer(x,y,1,ibuf) = work_acc
                                lift_buffer(x,y,1,ibuf) = lift_acc
                            end do
                        end do

                        call accumulate_band_histograms(pr_buffer(:,y_local_start:y_local_end,1,ibuf), work_buffer(:,y_local_start:y_local_end,1,ibuf), lift_buffer(:,y_local_start:y_local_end,1,ibuf), &
                                                        cell_area(:,hist_y_start:hist_y_end), pr_edges, work_edges, hist_area_chunk, hist_work_chunk, hist_lift_chunk, hist2d_work_chunk, hist2d_lift_chunk)
                    end do

                    if (band_has_data) then
                        !$omp critical(histogram_accumulation)
                            hist_area_out(iday,:,ilat) = hist_area_out(iday,:,ilat) + hist_area_chunk
                            hist_work_out(iday,:,ilat) = hist_work_out(iday,:,ilat) + hist_work_chunk
                            hist_lift_out(iday,:,ilat) = hist_lift_out(iday,:,ilat) + hist_lift_chunk
                            hist2d_work_out(iday,:,:,ilat) = hist2d_work_out(iday,:,:,ilat) + transpose(hist2d_work_chunk)
                            hist2d_lift_out(iday,:,:,ilat) = hist2d_lift_out(iday,:,:,ilat) + transpose(hist2d_lift_chunk)
                        !$omp end critical(histogram_accumulation)
                    end if
                end do

                deallocate(hist_area_chunk)
                deallocate(hist_work_chunk)
                deallocate(hist_lift_chunk)
                deallocate(hist2d_work_chunk)
                deallocate(hist2d_lift_chunk)

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

        call close_source_files(ncid_temp, ncid_omega, ncid_qv, ncid_qw, ncid_qr, ncid_qi, ncid_qs, ncid_qg, ncid_omt, ncid_omqv, ncid_omqw, ncid_omqr, ncid_omqi, ncid_omqs, ncid_omqg, ncid_pr)
        sources_processed = sources_processed + 1
    end do source_loop

    if (incident_log_unit /= -1) then
        write(incident_log_unit,'(A,I0)') '#sources_total'//char(9), n_sources
        write(incident_log_unit,'(A,I0)') '#sources_processed'//char(9), sources_processed
        write(incident_log_unit,'(A,I0)') '#sources_skipped'//char(9), sources_skipped
        close(incident_log_unit)
    end if

    if (sources_processed <= 0) then
        write(error_unit,*) 'No valid sources processed. sources_total/sources_skipped=', n_sources, sources_skipped
        stop 1
    end if

    if (use_period_filter .and. month_step_count(1) <= 0) then
        write(error_unit,*) 'No timesteps matched target_period_key: ', trim(target_period_key)
        stop 1
    end if

    ! Convert grouped sums to grouped means over native timesteps.
    do iday = 1, nmonths
        if (month_step_count(iday) > 0) then
            month_norm = 1.0d0 / dble(month_step_count(iday))
            hist_area_out(iday,:,:) = hist_area_out(iday,:,:) * month_norm
            hist_work_out(iday,:,:) = hist_work_out(iday,:,:) * month_norm
            hist_lift_out(iday,:,:) = hist_lift_out(iday,:,:) * month_norm
            hist2d_work_out(iday,:,:,:) = hist2d_work_out(iday,:,:,:) * month_norm
            hist2d_lift_out(iday,:,:,:) = hist2d_lift_out(iday,:,:,:) * month_norm
        end if
    end do

    call check(nf90_create(trim(adjustl(path_hist_out)), nf90_clobber, ncid_hist_out))
    call check(nf90_def_dim(ncid_hist_out, 'time', nmonths, dimid_hist_time))
    call check(nf90_def_dim(ncid_hist_out, 'nbin_pr', npr_edges - 1, dimid_nbin_pr))
    call check(nf90_def_dim(ncid_hist_out, 'nbin_work', nwork_edges - 1, dimid_nbin_work))
    call check(nf90_def_dim(ncid_hist_out, 'nedges_pr', npr_edges, dimid_nedges_pr))
    call check(nf90_def_dim(ncid_hist_out, 'nedges_work', nwork_edges, dimid_nedges_work))
    call check(nf90_def_dim(ncid_hist_out, 'nlat', nlat_bands, dimid_hist_nlat))
    call check(nf90_def_dim(ncid_hist_out, 'nsegment_max', max_segments_per_band, dimid_hist_nsegment))

    call check(nf90_def_var(ncid_hist_out, 'time', nf90_double, (/dimid_hist_time/), varid_out_time))
    call check(nf90_put_att(ncid_hist_out, varid_out_time, 'long_name', 'time'))
    call check(nf90_put_att(ncid_hist_out, varid_out_time, 'standard_name', 'time'))
    call check(nf90_put_att(ncid_hist_out, varid_out_time, 'units', trim(time_units)))
    call check(nf90_put_att(ncid_hist_out, varid_out_time, 'calendar', trim(time_calendar)))
    call check(nf90_put_att(ncid_hist_out, varid_out_time, 'axis', 'T'))

    call check(nf90_def_var(ncid_hist_out, 'lat_band', nf90_double, (/dimid_hist_nlat/), varid_lat_band))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band, 'long_name', 'latitude band index'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band, 'units', '1'))
    call check(nf90_def_var(ncid_hist_out, 'lat_band_south', nf90_double, (/dimid_hist_nlat/), varid_lat_band_south))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_south, 'long_name', 'latitude band south boundary'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_south, 'units', 'degrees_north'))
    call check(nf90_def_var(ncid_hist_out, 'lat_band_north', nf90_double, (/dimid_hist_nlat/), varid_lat_band_north))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_north, 'long_name', 'latitude band north boundary'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_north, 'units', 'degrees_north'))
    call check(nf90_def_var(ncid_hist_out, 'lat_band_segment_count', nf90_int, (/dimid_hist_nlat/), varid_lat_band_segment_count))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_count, 'long_name', 'number of active latitude segments in each band'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_count, 'units', '1'))
    call check(nf90_def_var(ncid_hist_out, 'lat_band_segment_south', nf90_double, (/dimid_hist_nlat, dimid_hist_nsegment/), varid_lat_band_segment_south))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_south, 'long_name', 'south boundary of each latitude-band segment'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_south, 'units', 'degrees_north'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_south, '_FillValue', fill_value))
    call check(nf90_def_var(ncid_hist_out, 'lat_band_segment_north', nf90_double, (/dimid_hist_nlat, dimid_hist_nsegment/), varid_lat_band_segment_north))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_north, 'long_name', 'north boundary of each latitude-band segment'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_north, 'units', 'degrees_north'))
    call check(nf90_put_att(ncid_hist_out, varid_lat_band_segment_north, '_FillValue', fill_value))

    call check(nf90_def_var(ncid_hist_out, 'pr_edges', nf90_double, (/dimid_nedges_pr/), varid_out_lon))
    call check(nf90_def_var(ncid_hist_out, 'work_edges', nf90_double, (/dimid_nedges_work/), varid_out_lat))
    call check(nf90_put_att(ncid_hist_out, varid_out_lon, 'long_name', 'precipitation bin edges'))
    call check(nf90_put_att(ncid_hist_out, varid_out_lon, 'bin_count', npr_edges - 1))
    call check(nf90_put_att(ncid_hist_out, varid_out_lat, 'long_name', 'work/lift bin edges'))
    call check(nf90_put_att(ncid_hist_out, varid_out_lat, 'bin_count', nwork_edges - 1))

    ! NetCDF Fortran interface reverses on-disk dimension order relative to the
    ! dimid list. Pass dimids in reverse so ncdump shows time-first ordering.
    call check(nf90_def_var(ncid_hist_out, 'hist_area', nf90_double, (/dimid_hist_nlat, dimid_nbin_pr, dimid_hist_time/), varid_hist_area))
    call check(nf90_def_var(ncid_hist_out, 'hist_work', nf90_double, (/dimid_hist_nlat, dimid_nbin_pr, dimid_hist_time/), varid_hist_work))
    call check(nf90_def_var(ncid_hist_out, 'hist_lift', nf90_double, (/dimid_hist_nlat, dimid_nbin_pr, dimid_hist_time/), varid_hist_lift))
    call check(nf90_def_var(ncid_hist_out, 'hist2d_work', nf90_double, (/dimid_hist_nlat, dimid_nbin_pr, dimid_nbin_work, dimid_hist_time/), varid_hist2d_work))
    call check(nf90_def_var(ncid_hist_out, 'hist2d_lift', nf90_double, (/dimid_hist_nlat, dimid_nbin_pr, dimid_nbin_work, dimid_hist_time/), varid_hist2d_lift))
    call check(nf90_put_att(ncid_hist_out, varid_hist_area, 'long_name', 'area histogram by precipitation bin and latitude band'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_area, 'time_stat', trim(to_lower(aggregation_mode))//'_mean'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_area, 'time_aggregation', 'mean over native timesteps within each aggregation period'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_work, 'long_name', 'area-weighted work sum by precipitation bin and latitude band'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_work, 'time_stat', trim(to_lower(aggregation_mode))//'_mean'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_work, 'time_aggregation', 'mean over native timesteps within each aggregation period'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_lift, 'long_name', 'area-weighted lift sum by precipitation bin and latitude band'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_lift, 'time_stat', trim(to_lower(aggregation_mode))//'_mean'))
    call check(nf90_put_att(ncid_hist_out, varid_hist_lift, 'time_aggregation', 'mean over native timesteps within each aggregation period'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_work, 'long_name', 'joint area histogram of work and precipitation by latitude band'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_work, 'time_stat', trim(to_lower(aggregation_mode))//'_mean'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_work, 'time_aggregation', 'mean over native timesteps within each aggregation period'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_lift, 'long_name', 'joint area histogram of lift and precipitation by latitude band'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_lift, 'time_stat', trim(to_lower(aggregation_mode))//'_mean'))
    call check(nf90_put_att(ncid_hist_out, varid_hist2d_lift, 'time_aggregation', 'mean over native timesteps within each aggregation period'))

    call check(nf90_put_att(ncid_hist_out, nf90_global, 'Conventions', 'CF-1.8'))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'aggregation', trim(to_lower(aggregation_mode))//' mean over native timesteps'))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'calendar_assumption', 'standard/gregorian/julian'))
    if (use_custom_lat_band_segments) then
        lat_mode = 'custom_segments'
    else if (use_custom_lat_band_bounds) then
        lat_mode = 'custom_boundaries'
    else
        lat_mode = 'equal_width_default'
    end if
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'latitude_band_mode', trim(lat_mode)))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'nlat_bands', nlat_bands))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'sources_total', n_sources))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'sources_processed', sources_processed))
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'sources_skipped', sources_skipped))
    if (continue_on_source_error) then
        continue_mode_attr = 'true '
    else
        continue_mode_attr = 'false'
    end if
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'continue_on_source_error', trim(continue_mode_attr)))
    if (sources_skipped > 0) then
        call check(nf90_put_att(ncid_hist_out, nf90_global, 'skipped_sources_log', trim(incident_log_path)))
    end if
    call check(nf90_put_att(ncid_hist_out, nf90_global, 'threshold_file', ''))
    if (use_period_filter) then
        call check(nf90_put_att(ncid_hist_out, nf90_global, 'target_period_key', trim(target_period_key)))
    end if

    call check(nf90_enddef(ncid_hist_out))

    do ilat = 1, nlat_bands
        do iseg = lat_band_segment_count(ilat) + 1, max_segments_per_band
            lat_band_segment_south(iseg, ilat) = fill_value
            lat_band_segment_north(iseg, ilat) = fill_value
        end do
    end do

    write(*,'(A)') 'write var: time'
    call check(nf90_put_var(ncid_hist_out, varid_out_time, month_time_vals))
    write(*,'(A)') 'write var: lat_band'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band, [(dble(iedge), iedge=1,nlat_bands)]))
    write(*,'(A)') 'write var: lat_band_south'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band_south, lat_band_south(1:nlat_bands)))
    write(*,'(A)') 'write var: lat_band_north'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band_north, lat_band_north(1:nlat_bands)))
    write(*,'(A)') 'write var: lat_band_segment_count'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band_segment_count, lat_band_segment_count(1:nlat_bands)))
    write(*,'(A)') 'write var: lat_band_segment_south'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band_segment_south, transpose(lat_band_segment_south(:,1:nlat_bands))))
    write(*,'(A)') 'write var: lat_band_segment_north'
    call check(nf90_put_var(ncid_hist_out, varid_lat_band_segment_north, transpose(lat_band_segment_north(:,1:nlat_bands))))
    write(*,'(A)') 'write var: pr_edges'
    call check(nf90_put_var(ncid_hist_out, varid_out_lon, pr_edges))
    write(*,'(A)') 'write var: work_edges'
    call check(nf90_put_var(ncid_hist_out, varid_out_lat, work_edges))

    ! Write histogram slabs explicitly in variable-definition order
    ! (nlat, nbin_pr, [nbin_work], time) to avoid implicit rank-order mismatches.
    do iday = 1, nmonths
        do ilat = 1, nlat_bands
            write(*,'(A,2(I0,1X))') 'write var: hist_area iday ilat=', iday, ilat
            call check(nf90_put_var(ncid_hist_out, varid_hist_area, hist_area_out(iday,:,ilat), &
                                    start=(/ilat, 1, iday/), count=(/1, npr_edges - 1, 1/)))
            write(*,'(A,2(I0,1X))') 'write var: hist_work iday ilat=', iday, ilat
            call check(nf90_put_var(ncid_hist_out, varid_hist_work, hist_work_out(iday,:,ilat), &
                                    start=(/ilat, 1, iday/), count=(/1, npr_edges - 1, 1/)))
            write(*,'(A,2(I0,1X))') 'write var: hist_lift iday ilat=', iday, ilat
            call check(nf90_put_var(ncid_hist_out, varid_hist_lift, hist_lift_out(iday,:,ilat), &
                                    start=(/ilat, 1, iday/), count=(/1, npr_edges - 1, 1/)))
            write(*,'(A,2(I0,1X))') 'write var: hist2d_work iday ilat=', iday, ilat
            call check(nf90_put_var(ncid_hist_out, varid_hist2d_work, transpose(hist2d_work_out(iday,:,:,ilat)), &
                                    start=(/ilat, 1, 1, iday/), count=(/1, npr_edges - 1, nwork_edges - 1, 1/)))
            write(*,'(A,2(I0,1X))') 'write var: hist2d_lift iday ilat=', iday, ilat
            call check(nf90_put_var(ncid_hist_out, varid_hist2d_lift, transpose(hist2d_lift_out(iday,:,:,ilat)), &
                                    start=(/ilat, 1, 1, iday/), count=(/1, npr_edges - 1, nwork_edges - 1, 1/)))
        end do
    end do

    call check(nf90_close(ncid_hist_out))

    call check(nf90_close(ncid_dz))
    if (allocated(source_dates)) deallocate(source_dates)
    if (allocated(source_roots)) deallocate(source_roots)

contains

    subroutine read_chunk(ncid, varid, buffer, nx_in, ny_in, nz_in, t_in, ystart_in)
        integer, intent(in) :: ncid, varid, nx_in, ny_in, nz_in, t_in, ystart_in
        double precision, intent(out) :: buffer(nx_in, ny_in, nz_in, 1)
        integer :: ncstatus_local, inq_status
        character(len=nf90_max_name) :: var_name

        var_name = '<unknown>'
        ncstatus_local = nf90_get_var(ncid, varid, buffer, start=(/1, ystart_in, 1, t_in/), count=(/nx_in, ny_in, nz_in, 1/))
        if (ncstatus_local /= nf90_noerr) then
            inq_status = nf90_inquire_variable(ncid, varid, name=var_name)
            if (inq_status /= nf90_noerr) var_name = '<inq_failed>'
            write(error_unit,*) 'NetCDF read_chunk failure: var=', trim(var_name), ' status=', trim(nf90_strerror(ncstatus_local))
            write(error_unit,*) '  start=(1,', ystart_in, ',1,', t_in, ') count=(', nx_in, ',', ny_in, ',', nz_in, ',1)'
            stop 1
        end if
    end subroutine read_chunk

    subroutine close_source_files(ncid_temp_in, ncid_omega_in, ncid_qv_in, ncid_qw_in, ncid_qr_in, ncid_qi_in, ncid_qs_in, ncid_qg_in, ncid_omt_in, ncid_omqv_in, ncid_omqw_in, ncid_omqr_in, ncid_omqi_in, ncid_omqs_in, ncid_omqg_in, ncid_pr_in)
        integer, intent(inout) :: ncid_temp_in, ncid_omega_in, ncid_qv_in, ncid_qw_in, ncid_qr_in, ncid_qi_in, ncid_qs_in, ncid_qg_in
        integer, intent(inout) :: ncid_omt_in, ncid_omqv_in, ncid_omqw_in, ncid_omqr_in, ncid_omqi_in, ncid_omqs_in, ncid_omqg_in, ncid_pr_in

        if (ncid_temp_in > 0) call check(nf90_close(ncid_temp_in))
        if (ncid_omega_in > 0) call check(nf90_close(ncid_omega_in))
        if (ncid_qv_in > 0) call check(nf90_close(ncid_qv_in))
        if (ncid_qw_in > 0) call check(nf90_close(ncid_qw_in))
        if (ncid_qr_in > 0) call check(nf90_close(ncid_qr_in))
        if (ncid_qi_in > 0) call check(nf90_close(ncid_qi_in))
        if (ncid_qs_in > 0) call check(nf90_close(ncid_qs_in))
        if (ncid_qg_in > 0) call check(nf90_close(ncid_qg_in))
        if (ncid_omt_in > 0) call check(nf90_close(ncid_omt_in))
        if (ncid_omqv_in > 0) call check(nf90_close(ncid_omqv_in))
        if (ncid_omqw_in > 0) call check(nf90_close(ncid_omqw_in))
        if (ncid_omqr_in > 0) call check(nf90_close(ncid_omqr_in))
        if (ncid_omqi_in > 0) call check(nf90_close(ncid_omqi_in))
        if (ncid_omqs_in > 0) call check(nf90_close(ncid_omqs_in))
        if (ncid_omqg_in > 0) call check(nf90_close(ncid_omqg_in))
        if (ncid_pr_in > 0) call check(nf90_close(ncid_pr_in))

        ncid_temp_in = -1
        ncid_omega_in = -1
        ncid_qv_in = -1
        ncid_qw_in = -1
        ncid_qr_in = -1
        ncid_qi_in = -1
        ncid_qs_in = -1
        ncid_qg_in = -1
        ncid_omt_in = -1
        ncid_omqv_in = -1
        ncid_omqw_in = -1
        ncid_omqr_in = -1
        ncid_omqi_in = -1
        ncid_omqs_in = -1
        ncid_omqg_in = -1
        ncid_pr_in = -1
    end subroutine close_source_files

    subroutine log_skip_incident(log_unit, source_date, stage, var_name, detail)
        integer, intent(in) :: log_unit
        character(len=*), intent(in) :: source_date, stage, var_name, detail

        write(error_unit,'(A)') 'Skipping source '//trim(source_date)//' stage='//trim(stage)//' var='//trim(var_name)//' detail='//trim(detail)
        if (log_unit /= -1) then
            write(log_unit,'(A)') trim(source_date)//char(9)//trim(stage)//char(9)//trim(var_name)//char(9)//trim(detail)
        end if
    end subroutine log_skip_incident

    subroutine log_source_incident(log_unit, source_date, stage, var_name, detail)
        integer, intent(in) :: log_unit
        character(len=*), intent(in) :: source_date, stage, var_name, detail

        write(error_unit,'(A)') 'Source incident '//trim(source_date)//' stage='//trim(stage)//' var='//trim(var_name)//' detail='//trim(detail)
        if (log_unit /= -1) then
            write(log_unit,'(A)') trim(source_date)//char(9)//trim(stage)//char(9)//trim(var_name)//char(9)//trim(detail)
        end if
    end subroutine log_source_incident

    logical function verify_time_axis_matches(ncid, var_label, reference_time, reference_units, reference_calendar, expected_nt, candidate_nt, source_date) result(matches)
        integer, intent(in) :: ncid, expected_nt, candidate_nt
        character(len=*), intent(in) :: var_label, reference_units, reference_calendar, source_date
        double precision, intent(in) :: reference_time(:)
        integer :: varid_time_local, ncstatus_local, i
        integer :: ndims_local
        integer :: dimids_local(nf90_max_var_dims)
        double precision, allocatable :: candidate_time(:)
        character(len=255) :: candidate_units, candidate_calendar

        matches = .false.

        if (candidate_nt /= expected_nt) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': nt differs. expected/got=', expected_nt, candidate_nt
            return
        end if

        ncstatus_local = nf90_inq_varid(ncid, 'time', varid_time_local)
        if (ncstatus_local /= nf90_noerr) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': missing time variable.'
            return
        end if

        call check(nf90_inquire_variable(ncid, varid_time_local, ndims=ndims_local, dimids=dimids_local))
        if (ndims_local < 1) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': time variable has no dimensions.'
            return
        end if

        if (dimids_local(ndims_local) <= 0) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': invalid time dimension.'
            return
        end if

        allocate(candidate_time(candidate_nt))
        call check(nf90_get_var(ncid, varid_time_local, candidate_time))

        candidate_units = ''
        candidate_calendar = ''
        ncstatus_local = nf90_get_att(ncid, varid_time_local, 'units', candidate_units)
        if (ncstatus_local /= nf90_noerr) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': missing time units.'
            deallocate(candidate_time)
            return
        end if
        ncstatus_local = nf90_get_att(ncid, varid_time_local, 'calendar', candidate_calendar)
        if (ncstatus_local /= nf90_noerr) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': missing time calendar.'
            deallocate(candidate_time)
            return
        end if

        if (trim(to_lower(candidate_units)) /= trim(to_lower(reference_units)) .or. trim(to_lower(candidate_calendar)) /= trim(to_lower(reference_calendar))) then
            write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': time units/calendar differ.'
            deallocate(candidate_time)
            return
        end if

        if (use_period_filter) then
            do i = 1, min(expected_nt, candidate_nt)
                if (.not. time_values_match(reference_time(i), candidate_time(i))) then
                    write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': time coordinate differs at index ', i
                    deallocate(candidate_time)
                    return
                end if
            end do
        else
            do i = 1, expected_nt
                if (.not. time_values_match(reference_time(i), candidate_time(i))) then
                    write(error_unit,*) 'Structural time-axis mismatch in source date ', trim(source_date), ' variable ', trim(var_label), ': time coordinate differs at index ', i
                    deallocate(candidate_time)
                    return
                end if
            end do
        end if

        deallocate(candidate_time)
        matches = .true.
    end function verify_time_axis_matches

    logical function check_var_shape(ncid, varid, var_label, nx_expected, ny_expected, nz_expected, nt_expected, source_date, nx_got, ny_got, nz_got, nt_got)
        integer, intent(in) :: ncid, varid, nx_expected, ny_expected, nz_expected, nt_expected
        character(len=*), intent(in) :: var_label, source_date
        integer, intent(out) :: nx_got, ny_got, nz_got, nt_got

        call get_var_shape(ncid, varid, nx_var, ny_var, nz_var, nt_var)
        nx_got = nx_var
        ny_got = ny_var
        nz_got = nz_var
        nt_got = nt_var

        if (nx_var /= nx_expected .or. ny_var /= ny_expected .or. nz_var /= nz_expected .or. (.not. use_period_filter .and. nt_var /= nt_expected)) then
            write(error_unit,*) 'Structural shape mismatch in source date ', trim(source_date), ' variable ', trim(var_label)
            write(error_unit,'(A,4(I0,1X))') '  expected nx ny nz nt: ', nx_expected, ny_expected, nz_expected, nt_expected
            write(error_unit,'(A,4(I0,1X))') '  got      nx ny nz nt: ', nx_var, ny_var, nz_var, nt_var
            check_var_shape = .false.
            return
        end if
        check_var_shape = .true.
    end function check_var_shape

    subroutine get_var_shape(ncid, varid, nx_out, ny_out, nz_out, nt_out)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: nx_out, ny_out, nz_out, nt_out
        integer :: ndims_local
        integer :: dimids_local(nf90_max_var_dims)

        call check(nf90_inquire_variable(ncid, varid, ndims=ndims_local, dimids=dimids_local))

        select case (ndims_local)
        case (4)
            call check(nf90_inquire_dimension(ncid, dimids_local(1), len=nx_out))
            call check(nf90_inquire_dimension(ncid, dimids_local(2), len=ny_out))
            call check(nf90_inquire_dimension(ncid, dimids_local(3), len=nz_out))
            call check(nf90_inquire_dimension(ncid, dimids_local(4), len=nt_out))
        case (3)
            call check(nf90_inquire_dimension(ncid, dimids_local(1), len=nx_out))
            call check(nf90_inquire_dimension(ncid, dimids_local(2), len=ny_out))
            nz_out = 1
            call check(nf90_inquire_dimension(ncid, dimids_local(3), len=nt_out))
        case default
            write(error_unit,*) 'Unsupported variable rank in get_var_shape. rank=', ndims_local
            stop 1
        end select
    end subroutine get_var_shape

    logical function time_values_match(a, b)
        double precision, intent(in) :: a, b
        double precision :: scale

        scale = max(1.0d0, abs(a), abs(b))
        time_values_match = abs(a - b) <= 1.0d-10 * scale
    end function time_values_match

    subroutine get_var_time_len(ncid, varid, ntime_out)
        integer, intent(in) :: ncid, varid
        integer, intent(out) :: ntime_out
        integer :: ndims_local
        integer :: dimids_local(nf90_max_var_dims)

        call check(nf90_inquire_variable(ncid, varid, ndims=ndims_local, dimids=dimids_local))
        call check(nf90_inquire_dimension(ncid, dimids_local(ndims_local), len=ntime_out))
    end subroutine get_var_time_len

    subroutine resolve_lat_band_boundaries(nlat_in, count_in, custom_in, bounds_in, south_out, north_out)
        integer, intent(in) :: nlat_in, count_in
        logical, intent(in) :: custom_in
        double precision, intent(in) :: bounds_in(:)
        double precision, intent(out) :: south_out(:), north_out(:)
        integer :: i

        if (custom_in) then
            do i = 1, nlat_in
                south_out(i) = bounds_in(i)
                north_out(i) = bounds_in(i + 1)
            end do
        else
            do i = 1, nlat_in
                south_out(i) = -90.0d0 + dble(i - 1) * (180.0d0 / dble(nlat_in))
                north_out(i) = -90.0d0 + dble(i) * (180.0d0 / dble(nlat_in))
            end do
        end if
    end subroutine resolve_lat_band_boundaries

    subroutine resolve_lat_band_segments(nlat_in, count_in, custom_bounds_in, bounds_in, custom_segments_in, seg_count_inout, seg_south_inout, seg_north_inout, south_out, north_out)
        integer, intent(in) :: nlat_in, count_in
        logical, intent(in) :: custom_bounds_in, custom_segments_in
        double precision, intent(in) :: bounds_in(:)
        integer, intent(inout) :: seg_count_inout(:)
        double precision, intent(inout) :: seg_south_inout(:,:), seg_north_inout(:,:)
        double precision, intent(out) :: south_out(:), north_out(:)
        integer :: i, s

        if (custom_segments_in) then
            do i = 1, nlat_in
                south_out(i) = 1.0d99
                north_out(i) = -1.0d99
                do s = 1, seg_count_inout(i)
                    south_out(i) = min(south_out(i), seg_south_inout(s, i))
                    north_out(i) = max(north_out(i), seg_north_inout(s, i))
                end do
            end do
        else
            call resolve_lat_band_boundaries(nlat_in, count_in, custom_bounds_in, bounds_in, south_out, north_out)
            seg_count_inout = 0
            do i = 1, nlat_in
                seg_count_inout(i) = 1
                seg_south_inout(1, i) = south_out(i)
                seg_north_inout(1, i) = north_out(i)
            end do
        end if
    end subroutine resolve_lat_band_segments

    subroutine init_piecewise_log_bins(edges)
        double precision, intent(out) :: edges(:)
        integer :: i

        if (size(edges) /= npr_edges) then
            error stop 'init_piecewise_log_bins: unexpected edges size'
        end if

        do i = 1, 200
            edges(i) = 1.0d-6 * (1.0d-4 / 1.0d-6) ** (dble(i - 1) / 200.0d0)
        end do

        do i = 201, 800
            edges(i) = 1.0d-4 * (1.0d-2 / 1.0d-4) ** (dble(i - 201) / 600.0d0)
        end do

        do i = 801, 1201
            edges(i) = 1.0d-2 * (2.0d-1 / 1.0d-2) ** (dble(i - 801) / 400.0d0)
        end do
    end subroutine init_piecewise_log_bins

    subroutine resolve_lat_band_indices(lat_in, lat_south_in, lat_north_in, lat_start, lat_end)
        double precision, intent(in) :: lat_in(:)
        double precision, intent(in) :: lat_south_in, lat_north_in
        integer, intent(out) :: lat_start, lat_end
        integer :: i
        logical :: found_start, found_end

        found_start = .false.
        found_end = .false.
        lat_start = 1
        lat_end = size(lat_in)

        do i = 1, size(lat_in)
            if (.not. found_start .and. lat_in(i) >= lat_south_in) then
                lat_start = i
                found_start = .true.
            end if
            if (lat_in(i) < lat_north_in .or. (lat_in(i) >= 90.0d0 - 1.0d-8 .and. lat_north_in >= 90.0d0)) then
                lat_end = i
                found_end = .true.
            end if
        end do

        if (.not. found_start .or. .not. found_end) then
            lat_start = 1
            lat_end = 0
        end if
    end subroutine resolve_lat_band_indices

    ! Group contiguous timesteps by calendar month; month_times stores each
    ! group's representative time value (the first timestep in the month).
    subroutine build_month_groups(time_in, units, calendar, month_of_t, month_first_idx, month_times, month_step_count, nmonths)
        double precision, intent(in) :: time_in(:)
        character(len=*), intent(in) :: units, calendar
        integer, allocatable, intent(out) :: month_of_t(:), month_first_idx(:), month_step_count(:)
        double precision, allocatable, intent(out) :: month_times(:)
        integer, intent(out) :: nmonths
        integer :: i, year, month, day
        integer, allocatable :: month_key(:)
        double precision :: scale_days
        integer :: ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second

        call parse_time_units(units, ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second, scale_days)
        if (.not. is_supported_calendar(calendar)) then
            error stop 'build_month_groups: unsupported calendar'
        end if

        allocate(month_key(size(time_in)))
        do i = 1, size(time_in)
            call time_value_to_year_month(time_in(i), scale_days, ref_year, ref_month, ref_day, calendar, year, month)
            month_key(i) = year * 12 + month
        end do

        nmonths = 1
        do i = 2, size(month_key)
            if (month_key(i) /= month_key(i-1)) nmonths = nmonths + 1
        end do

        allocate(month_of_t(size(time_in)))
        allocate(month_first_idx(nmonths))
        allocate(month_times(nmonths))
        allocate(month_step_count(nmonths))
        month_step_count = 0

        nmonths = 1
        month_of_t(1) = 1
        month_first_idx(1) = 1
        month_times(1) = time_in(1)
        month_step_count(1) = 1
        do i = 2, size(month_key)
            if (month_key(i) /= month_key(i-1)) then
                nmonths = nmonths + 1
                month_first_idx(nmonths) = i
                month_times(nmonths) = time_in(i)
            end if
            month_of_t(i) = nmonths
            month_step_count(nmonths) = month_step_count(nmonths) + 1
        end do

        deallocate(month_key)
    end subroutine build_month_groups

    ! Group contiguous timesteps by calendar day; day_times stores each
    ! group's representative time value (the first timestep in the day).
    subroutine build_day_groups(time_in, units, calendar, day_of_t, day_first_idx, day_times, day_step_count, ndays)
        double precision, intent(in) :: time_in(:)
        character(len=*), intent(in) :: units, calendar
        integer, allocatable, intent(out) :: day_of_t(:), day_first_idx(:), day_step_count(:)
        double precision, allocatable, intent(out) :: day_times(:)
        integer, intent(out) :: ndays
        integer :: i, year, month, day
        integer, allocatable :: day_key(:)
        double precision :: scale_days
        integer :: ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second

        call parse_time_units(units, ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second, scale_days)
        if (.not. is_supported_calendar(calendar)) then
            error stop 'build_day_groups: unsupported calendar'
        end if

        allocate(day_key(size(time_in)))
        do i = 1, size(time_in)
            call time_value_to_year_month_day(time_in(i), scale_days, ref_year, ref_month, ref_day, calendar, year, month, day)
            day_key(i) = year * 10000 + month * 100 + day
        end do

        ndays = 1
        do i = 2, size(day_key)
            if (day_key(i) /= day_key(i-1)) ndays = ndays + 1
        end do

        allocate(day_of_t(size(time_in)))
        allocate(day_first_idx(ndays))
        allocate(day_times(ndays))
        allocate(day_step_count(ndays))
        day_step_count = 0

        ndays = 1
        day_of_t(1) = 1
        day_first_idx(1) = 1
        day_times(1) = time_in(1)
        day_step_count(1) = 1
        do i = 2, size(day_key)
            if (day_key(i) /= day_key(i-1)) then
                ndays = ndays + 1
                day_first_idx(ndays) = i
                day_times(ndays) = time_in(i)
            end if
            day_of_t(i) = ndays
            day_step_count(ndays) = day_step_count(ndays) + 1
        end do

        deallocate(day_key)
    end subroutine build_day_groups

    subroutine time_value_to_year_month(time_value, scale_days, ref_year, ref_month, ref_day, calendar, out_year, out_month)
        double precision, intent(in) :: time_value, scale_days
        integer, intent(in) :: ref_year, ref_month, ref_day
        character(len=*), intent(in) :: calendar
        integer, intent(out) :: out_year, out_month
        integer :: days_elapsed, year, month, day

        days_elapsed = int(floor(time_value * scale_days + 1.0d-9))
        year = ref_year
        month = ref_month
        day = ref_day
        call advance_date_by_days(year, month, day, days_elapsed, calendar)
        out_year = year
        out_month = month
    end subroutine time_value_to_year_month

    subroutine time_value_to_year_month_day(time_value, scale_days, ref_year, ref_month, ref_day, calendar, out_year, out_month, out_day)
        double precision, intent(in) :: time_value, scale_days
        integer, intent(in) :: ref_year, ref_month, ref_day
        character(len=*), intent(in) :: calendar
        integer, intent(out) :: out_year, out_month, out_day
        integer :: days_elapsed, year, month, day

        days_elapsed = int(floor(time_value * scale_days + 1.0d-9))
        year = ref_year
        month = ref_month
        day = ref_day
        call advance_date_by_days(year, month, day, days_elapsed, calendar)
        out_year = year
        out_month = month
        out_day = day
    end subroutine time_value_to_year_month_day

    subroutine advance_date_by_days(year, month, day, ndays, calendar)
        integer, intent(inout) :: year, month, day
        integer, intent(in) :: ndays
        character(len=*), intent(in) :: calendar
        integer :: remaining, dim

        remaining = ndays
        do while (remaining > 0)
            dim = days_in_month(year, month, calendar)
            if (day + remaining <= dim) then
                day = day + remaining
                remaining = 0
            else
                remaining = remaining - (dim - day + 1)
                day = 1
                call increment_month(year, month)
            end if
        end do
    end subroutine advance_date_by_days

    subroutine increment_month(year, month)
        integer, intent(inout) :: year, month
        month = month + 1
        if (month > 12) then
            month = 1
            year = year + 1
        end if
    end subroutine increment_month

    integer function days_in_month(year, month, calendar)
        integer, intent(in) :: year, month
        character(len=*), intent(in) :: calendar
        integer, dimension(12) :: mdays

        mdays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        if (is_supported_calendar(calendar)) then
            if (month == 2 .and. is_leap_year(year, calendar)) then
                days_in_month = 29
            else
                days_in_month = mdays(month)
            end if
        else
            error stop 'days_in_month: unsupported calendar'
        end if
    end function days_in_month

    logical function is_leap_year(year, calendar)
        integer, intent(in) :: year
        character(len=*), intent(in) :: calendar

        if (trim(to_lower(calendar)) == 'julian') then
            is_leap_year = (mod(year, 4) == 0)
        else
            is_leap_year = (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0))
        end if
    end function is_leap_year

    logical function is_supported_calendar(calendar)
        character(len=*), intent(in) :: calendar
        character(len=255) :: c

        c = trim(to_lower(calendar))
        is_supported_calendar = (c == 'standard' .or. c == 'gregorian' .or. c == 'julian')
    end function is_supported_calendar

    subroutine parse_time_units(units, ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second, scale_to_days)
        character(len=*), intent(in) :: units
        integer, intent(out) :: ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second
        double precision, intent(out) :: scale_to_days
        character(len=255) :: lower, rest
        integer :: idx

        lower = to_lower(trim(units))
        idx = index(lower, 'since')
        if (idx > 0) then
            rest = adjustl(units(idx + 5:))
        else
            rest = '1900-01-01 00:00:00'
        end if

        ref_year = 1900
        ref_month = 1
        ref_day = 1
        ref_hour = 0
        ref_minute = 0
        ref_second = 0

        if (len_trim(rest) >= 10) then
            read(rest(1:4), *) ref_year
            read(rest(6:7), *) ref_month
            read(rest(9:10), *) ref_day
            if (len_trim(rest) >= 19) then
                read(rest(12:13), *) ref_hour
                read(rest(15:16), *) ref_minute
                read(rest(18:19), *) ref_second
            end if
        end if

        if (index(lower, 'hour') == 1) then
            scale_to_days = 1.0d0 / 24.0d0
        else if (index(lower, 'minute') == 1) then
            scale_to_days = 1.0d0 / 1440.0d0
        else if (index(lower, 'second') == 1) then
            scale_to_days = 1.0d0 / 86400.0d0
        else if (index(lower, 'day') == 1) then
            scale_to_days = 1.0d0
        else
            scale_to_days = 1.0d0
        end if
    end subroutine parse_time_units

    ! Area-weighted accumulation for a single latitude-band chunk.
    subroutine accumulate_band_histograms(pr_slice, work_slice, lift_slice, area_slice, pr_edges_in, work_edges_in, hist_area_out, hist_work_out, hist_lift_out, hist2d_work_out, hist2d_lift_out)
        double precision, intent(in) :: pr_slice(:,:)
        double precision, intent(in) :: work_slice(:,:)
        double precision, intent(in) :: lift_slice(:,:)
        double precision, intent(in) :: area_slice(:,:)
        double precision, intent(in) :: pr_edges_in(:)
        double precision, intent(in) :: work_edges_in(:)
        double precision, intent(inout) :: hist_area_out(:)
        double precision, intent(inout) :: hist_work_out(:)
        double precision, intent(inout) :: hist_lift_out(:)
        double precision, intent(inout) :: hist2d_work_out(:,:)
        double precision, intent(inout) :: hist2d_lift_out(:,:)
        integer :: ix, iy, ipr, iwork
        double precision :: area_val, pr_val, work_val, lift_val

        do iy = 1, size(pr_slice, 2)
            do ix = 1, size(pr_slice, 1)
                area_val = area_slice(ix, iy)
                if (area_val <= 0.0d0) cycle
                pr_val = pr_slice(ix, iy)
                work_val = work_slice(ix, iy)
                lift_val = lift_slice(ix, iy)
                ipr = find_bin_index(pr_val, pr_edges_in)
                if (ipr > 0) then
                    hist_area_out(ipr) = hist_area_out(ipr) + area_val
                    hist_work_out(ipr) = hist_work_out(ipr) + work_val * area_val
                    hist_lift_out(ipr) = hist_lift_out(ipr) + lift_val * area_val
                end if
                iwork = find_bin_index(work_val, work_edges_in)
                if (ipr > 0 .and. iwork > 0) then
                    hist2d_work_out(iwork, ipr) = hist2d_work_out(iwork, ipr) + area_val
                end if
                iwork = find_bin_index(lift_val, work_edges_in)
                if (ipr > 0 .and. iwork > 0) then
                    hist2d_lift_out(iwork, ipr) = hist2d_lift_out(iwork, ipr) + area_val
                end if
            end do
        end do
    end subroutine accumulate_band_histograms

    ! Return 1-based bin index for half-open bins [edge_i, edge_{i+1}), with
    ! the top edge included in the final bin.
    integer function find_bin_index(value, edges)
        double precision, intent(in) :: value
        double precision, intent(in) :: edges(:)
        integer :: lo, hi, mid

        if (value < edges(1) .or. value > edges(size(edges))) then
            find_bin_index = 0
            return
        end if
        if (value == edges(size(edges))) then
            find_bin_index = size(edges) - 1
            return
        end if

        lo = 1
        hi = size(edges) - 1
        do while (lo <= hi)
            mid = (lo + hi) / 2
            if (value < edges(mid)) then
                hi = mid - 1
            else if (value >= edges(mid + 1)) then
                lo = mid + 1
            else
                find_bin_index = mid
                return
            end if
        end do
        find_bin_index = 0
    end function find_bin_index

    subroutine split_date_root_line(line_in, date_out, root_out)
        character(len=*), intent(in) :: line_in
        character(len=*), intent(out) :: date_out, root_out
        integer :: sep

        date_out = ''
        root_out = ''
        sep = index(line_in, '|')
        if (sep <= 1 .or. sep >= len_trim(line_in)) then
            error stop 'split_date_root_line: expected format date|source_root'
        end if
        date_out = adjustl(line_in(1:sep-1))
        root_out = adjustl(line_in(sep+1:len_trim(line_in)))
    end subroutine split_date_root_line

    subroutine time_matches_target_period(time_value, units, calendar, mode_in, target_key, matches)
        double precision, intent(in) :: time_value
        character(len=*), intent(in) :: units, calendar, mode_in, target_key
        logical, intent(out) :: matches
        integer :: year, month, day
        integer :: ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second
        double precision :: scale_days
        character(len=8) :: key_daily
        character(len=6) :: key_monthly

        call parse_time_units(units, ref_year, ref_month, ref_day, ref_hour, ref_minute, ref_second, scale_days)
        call time_value_to_year_month_day(time_value, scale_days, ref_year, ref_month, ref_day, calendar, year, month, day)

        write(key_daily, '(I4.4,I2.2,I2.2)') year, month, day
        write(key_monthly, '(I4.4,I2.2)') year, month

        select case (trim(to_lower(mode_in)))
        case ('daily')
            matches = (trim(target_key) == trim(key_daily))
        case ('monthly')
            matches = (trim(target_key) == trim(key_monthly))
        case default
            matches = .false.
        end select
    end subroutine time_matches_target_period

    function to_lower(text) result(out_text)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out_text
        integer :: i
        do i = 1, len(text)
            if (text(i:i) >= 'A' .and. text(i:i) <= 'Z') then
                out_text(i:i) = achar(iachar(text(i:i)) + 32)
            else
                out_text(i:i) = text(i:i)
            end if
        end do
    end function to_lower

end program compute_work_async_histograms_by_lat_band
