program compute_prate_thresholds_by_lat_band

    use netcdf
    use nc
    use hist
    use interp
    use iso_fortran_env, only: error_unit

    implicit none

    integer, parameter :: nbins = 1200
    integer, parameter :: nedges = nbins + 1
    integer, parameter :: npercentiles = 11
    integer, parameter :: nlat_bands = 18
    ! Latitude bands: 18 × 10-degree bands from -90 to +90
    ! Band 1: [-90, -80], Band 2: [-80, -70], ..., Band 18: [80, 90]

    character(len=1024) :: filenml, msg
    character(len=1024) :: history_root_part1, history_root_part2
    character(len=1024) :: date_list_file_part1, date_list_file_part2
    character(len=1024) :: output_file
    character(len=1024) :: input_path

    character(len=64), allocatable :: dates(:)
    real(8), allocatable :: bin_edges(:)
    real(8), allocatable :: hist_area_out(:,:)
    real(8), allocatable :: cell_area(:,:)
    real(8), allocatable :: cdf_edges(:)
    real(8), allocatable :: percentile_values(:)
    real(8), allocatable :: thresholds(:,:)

    logical :: have_grid
    logical :: have_clip
    integer :: ndates
    integer :: idate
    integer :: tmp_unit, ios
    integer :: lat_clip_start, lat_clip_end, lon_clip_start, lon_clip_end
    integer :: ilat

    namelist /config/ history_root_part1, history_root_part2, date_list_file_part1, date_list_file_part2, output_file

    history_root_part1 = ''
    history_root_part2 = ''
    date_list_file_part1 = ''
    date_list_file_part2 = ''
    output_file = 'thresholds_control_by_lat_band.txt'

    call get_command_argument(1, filenml)
    if (len_trim(filenml) == 0) then
        write(error_unit,*) 'Usage: compute_prate_thresholds_by_lat_band <config.nml>'
        stop 1
    end if

    open(newunit=tmp_unit, file=trim(adjustl(filenml)), status='old', iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to open configuration namelist, iomsg='//trim(msg)
        stop 1
    end if

    read(unit=tmp_unit, nml=config, iostat=ios, iomsg=msg)
    if (ios /= 0) then
        write(error_unit,*) 'Failed to read namelist, iomsg='//trim(msg)
        stop 1
    end if
    close(tmp_unit)

    if (len_trim(history_root_part1) == 0) then
        write(error_unit,*) 'Namelist item history_root_part1 is required.'
        stop 1
    end if

    if (len_trim(history_root_part2) == 0) then
        write(error_unit,*) 'Namelist item history_root_part2 is required.'
        stop 1
    end if

    if (len_trim(date_list_file_part1) == 0) then
        write(error_unit,*) 'Namelist item date_list_file_part1 is required.'
        stop 1
    end if

    if (len_trim(date_list_file_part2) == 0) then
        write(error_unit,*) 'Namelist item date_list_file_part2 is required.'
        stop 1
    end if

    allocate(bin_edges(nedges))
    call init_piecewise_log_bins(bin_edges)

    ! Allocate histogram for all latitude bands: (nbins, nlat_bands)
    allocate(hist_area_out(nbins, nlat_bands))
    hist_area_out = 0.0d0

    have_grid = .false.
    have_clip = .false.

    ! Pass 1: parse all precipitation files from part1 and part2 and accumulate area-weighted histogram per lat band.
    call process_history_list_pair(trim(adjustl(history_root_part1)), trim(adjustl(date_list_file_part1)), &
                                   bin_edges, hist_area_out, cell_area, have_grid, have_clip, &
                                   lat_clip_start, lat_clip_end, lon_clip_start, lon_clip_end)
    call process_history_list_pair(trim(adjustl(history_root_part2)), trim(adjustl(date_list_file_part2)), &
                                   bin_edges, hist_area_out, cell_area, have_grid, have_clip, &
                                   lat_clip_start, lat_clip_end, lon_clip_start, lon_clip_end)

    allocate(cdf_edges(nedges))
    allocate(percentile_values(npercentiles))
    percentile_values = (/ 0.5d0, 0.75d0, 0.9d0, 0.95d0, 0.97d0, 0.99d0, 0.995d0, 0.999d0, 0.9995d0, 0.9999d0, 0.99999d0 /)

    allocate(thresholds(nlat_bands, npercentiles))

    ! Pass 2: compute thresholds for each latitude band independently
    do ilat = 1, nlat_bands
        ! Guard against completely dry latitude bands (e.g., polar regions)
        if (sum(hist_area_out(:, ilat)) <= 0.0d0) then
            thresholds(ilat, :) = 0.0d0
            cycle
        end if

        call build_normalized_cdf(hist_area_out(:, ilat), cdf_edges)
        call compute_percentile_thresholds_loglog(bin_edges, cdf_edges, percentile_values, thresholds(ilat, :))
    end do

    call write_thresholds_ascii_2d(trim(adjustl(output_file)), percentile_values, thresholds)

contains

    subroutine process_history_list_pair(root, list_file, edges, hist_accum, area_ref, have_grid_ref, have_clip_ref, &
                                         lat_clip_start_ref, lat_clip_end_ref, lon_clip_start_ref, lon_clip_end_ref)
        character(len=*), intent(in) :: root, list_file
        real(8), intent(in) :: edges(:)
        real(8), intent(inout) :: hist_accum(:,:)
        real(8), allocatable, intent(inout) :: area_ref(:,:)
        logical, intent(inout) :: have_grid_ref
        logical, intent(inout) :: have_clip_ref
        integer, intent(inout) :: lat_clip_start_ref, lat_clip_end_ref, lon_clip_start_ref, lon_clip_end_ref

        integer :: i

        call read_date_list(list_file, dates, ndates)
        if (ndates <= 0) then
            write(error_unit,*) 'No dates found in list file: ', trim(adjustl(list_file))
            stop 1
        end if

        do i = 1, ndates
            call build_prate_path(trim(adjustl(root)), trim(adjustl(dates(i))), input_path)
            call accumulate_hist_from_file_by_lat_band(trim(adjustl(input_path)), edges, hist_accum, area_ref, have_grid_ref, have_clip_ref, &
                                           lat_clip_start_ref, lat_clip_end_ref, lon_clip_start_ref, lon_clip_end_ref)
        end do

        if (allocated(dates)) deallocate(dates)
    end subroutine process_history_list_pair

    subroutine read_date_list(path_list, dates_out, n_out)
        character(len=*), intent(in) :: path_list
        character(len=64), allocatable, intent(out) :: dates_out(:)
        integer, intent(out) :: n_out

        integer :: unit_list, ios_list, count, i
        character(len=256) :: line

        count = 0
        open(newunit=unit_list, file=trim(adjustl(path_list)), status='old', iostat=ios_list, iomsg=msg)
        if (ios_list /= 0) then
            write(error_unit,*) 'Failed to open date list file, iomsg='//trim(msg)
            stop 1
        end if

        do
            read(unit_list, '(A)', iostat=ios_list) line
            if (ios_list /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            count = count + 1
        end do

        rewind(unit_list)
        if (count > 0) then
            allocate(dates_out(count))
            i = 0
            do
                read(unit_list, '(A)', iostat=ios_list) line
                if (ios_list /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#' .or. line(1:1) == '!') cycle
                i = i + 1
                dates_out(i) = trim(line)
            end do
            n_out = i
        else
            allocate(dates_out(0))
            n_out = 0
        end if

        close(unit_list)
    end subroutine read_date_list

    subroutine build_prate_path(root, date_str, out_path)
        character(len=*), intent(in) :: root, date_str
        character(len=*), intent(out) :: out_path

        integer :: root_len

        root_len = len_trim(root)
        if (root_len <= 0) then
            out_path = ''
            return
        end if

        if (root(root_len:root_len) == '/') then
            out_path = trim(root)//trim(date_str)//'/PRATEsfc_coarse_C3072_1440x720.fre.nc'
        else
            out_path = trim(root)//'/'//trim(date_str)//'/PRATEsfc_coarse_C3072_1440x720.fre.nc'
        end if
    end subroutine build_prate_path

    subroutine init_piecewise_log_bins(edges)
        real(8), intent(out) :: edges(:)
        integer :: i

        if (size(edges) /= nedges) then
            error stop 'init_piecewise_log_bins: unexpected edges size'
        end if

        ! Tier 1: 1.0E-6 to 1.0E-4 (200 bins)
        do i = 1, 200
            edges(i) = 1.0d-6 * (1.0d-4 / 1.0d-6) ** (dble(i - 1) / 200.0d0)
        end do

        ! Tier 2: 1.0E-4 to 1.0E-2 (600 bins)
        do i = 201, 800
            edges(i) = 1.0d-4 * (1.0d-2 / 1.0d-4) ** (dble(i - 201) / 600.0d0)
        end do

        ! Tier 3: 1.0E-2 to 2.0E-1 (400 bins)
        do i = 801, 1201
            edges(i) = 1.0d-2 * (2.0d-1 / 1.0d-2) ** (dble(i - 801) / 400.0d0)
        end do
    end subroutine init_piecewise_log_bins

    subroutine accumulate_hist_from_file_by_lat_band(path_in, edges, hist_accum, area_ref, have_grid_ref, have_clip_ref, &
                                         lat_clip_start_ref, lat_clip_end_ref, lon_clip_start_ref, lon_clip_end_ref)
        character(len=*), intent(in) :: path_in
        real(8), intent(in) :: edges(:)
        real(8), intent(inout) :: hist_accum(:,:)
        real(8), allocatable, intent(inout) :: area_ref(:,:)
        logical, intent(inout) :: have_grid_ref
        logical, intent(inout) :: have_clip_ref
        integer, intent(inout) :: lat_clip_start_ref, lat_clip_end_ref, lon_clip_start_ref, lon_clip_end_ref

        integer :: ncid_in, varid_pr
        integer :: varid_lon, varid_lat
        integer :: ndims, natts, xtype
        integer :: dimids(nf90_max_var_dims)
        integer :: nx, ny, n3, n4
        integer :: dlen
        integer :: k3, k4, ilat, jlat

        real(8), allocatable :: lon(:), lat(:)
        real(8), allocatable :: pr3d(:,:,:)
        real(8), allocatable :: hist_area_chunk(:)
        real(8) :: lat_south_band, lat_north_band
        integer :: lat_start_band, lat_end_band, lat_start_clip, lat_end_clip

        call check(nf90_open(trim(path_in), nf90_nowrite, ncid_in))
        call check(nf90_inq_varid(ncid_in, 'PRATEsfc_coarse', varid_pr))
        call check(nf90_inquire_variable(ncid_in, varid_pr, xtype=xtype, ndims=ndims, dimids=dimids, natts=natts))

        if (ndims < 2 .or. ndims > 4) then
            write(error_unit,*) 'Unsupported PRATEsfc_coarse rank in file: ', trim(path_in)
            stop 1
        end if

        call check(nf90_inquire_dimension(ncid_in, dimids(1), len=nx))
        call check(nf90_inquire_dimension(ncid_in, dimids(2), len=ny))

        if (ndims >= 3) then
            call check(nf90_inquire_dimension(ncid_in, dimids(3), len=n3))
        else
            n3 = 1
        end if

        if (ndims == 4) then
            call check(nf90_inquire_dimension(ncid_in, dimids(4), len=n4))
        else
            n4 = 1
        end if

        if (.not. have_grid_ref) then
            call check(nf90_inq_varid(ncid_in, 'grid_xt_coarse', varid_lon))
            call check(nf90_inq_varid(ncid_in, 'grid_yt_coarse', varid_lat))
            allocate(lon(nx), lat(ny))
            call check(nf90_get_var(ncid_in, varid_lon, lon))
            call check(nf90_get_var(ncid_in, varid_lat, lat))
            call compute_cell_areas(lon, lat, area_ref)
            ! Global longitude clip to 0:360 (full domain)
            lon_clip_start_ref = 1
            lon_clip_end_ref = nx
            deallocate(lon, lat)
            have_grid_ref = .true.
        else
            if (.not. allocated(area_ref)) then
                write(error_unit,*) 'Internal error: area_ref should be allocated when have_grid_ref=true'
                stop 1
            end if
            if (size(area_ref, 1) /= nx .or. size(area_ref, 2) /= ny) then
                write(error_unit,*) 'Grid mismatch for file: ', trim(path_in)
                stop 1
            end if
        end if

        allocate(pr3d(nx, ny, 1))

        ! Read latitude grid if not already done
        if (.not. have_clip_ref) then
            allocate(lat(ny))
            call check(nf90_inq_varid(ncid_in, 'grid_yt_coarse', varid_lat))
            call check(nf90_get_var(ncid_in, varid_lat, lat))
            have_clip_ref = .true.
        end if

        ! Process each latitude band separately
        do ilat = 1, nlat_bands
            ! Band ilat covers latitude from lat_south_band to lat_north_band (in 10-degree increments)
            lat_south_band = -90.0d0 + dble(ilat - 1) * 10.0d0
            lat_north_band = lat_south_band + 10.0d0

            ! Find indices in the native latitude grid that fall within this band
            call resolve_lat_band_indices(lat, lat_south_band, lat_north_band, lat_start_band, lat_end_band)

            if (lat_start_band > lat_end_band) then
                cycle  ! No data in this band
            end if

            ! Accumulate histogram for this latitude band only
            do k4 = 1, n4
                do k3 = 1, n3
                    select case (ndims)
                    case (2)
                        call check(nf90_get_var(ncid_in, varid_pr, pr3d(:,:,1), start=(/1,1/), count=(/nx,ny/)))
                    case (3)
                        call check(nf90_get_var(ncid_in, varid_pr, pr3d(:,:,1), start=(/1,1,k3/), count=(/nx,ny,1/)))
                    case (4)
                        call check(nf90_get_var(ncid_in, varid_pr, pr3d(:,:,1), start=(/1,1,k3,k4/), count=(/nx,ny,1,1/)))
                    end select

                    ! Accumulate histogram only inside this latitude band.
                    call bin_sumarea(pr3d(lon_clip_start_ref:lon_clip_end_ref, lat_start_band:lat_end_band, :), &
                                     area_ref(lon_clip_start_ref:lon_clip_end_ref, lat_start_band:lat_end_band), &
                                     edges, hist_area_chunk)
                    hist_accum(:, ilat) = hist_accum(:, ilat) + hist_area_chunk
                    if (allocated(hist_area_chunk)) deallocate(hist_area_chunk)
                end do
            end do
        end do

        if (allocated(pr3d)) deallocate(pr3d)
        if (allocated(lat)) deallocate(lat)

        call check(nf90_close(ncid_in))

        dlen = len_trim(path_in)
        if (dlen > 0) then
            write(*,*) 'Processed file: ', trim(path_in)
        end if
    end subroutine accumulate_hist_from_file_by_lat_band

    subroutine resolve_lat_band_indices(lat, lat_south, lat_north, lat_start, lat_end)
        real(8), intent(in) :: lat(:)
        real(8), intent(in) :: lat_south, lat_north
        integer, intent(out) :: lat_start, lat_end

        integer :: i
        logical :: found_start, found_end

        found_start = .false.
        found_end = .false.
        lat_start = 1
        lat_end = size(lat)

        do i = 1, size(lat)
            if (.not. found_start .and. lat(i) >= lat_south) then
                lat_start = i
                found_start = .true.
            end if

            ! Prevent double-counting boundaries. Include the north pole (lat >= 90)
            ! only for the northernmost band where lat_north >= 90. Use a small
            ! tolerance to guard against floating-point representation of 90 deg.
            if (lat(i) < lat_north .or. (lat(i) >= 90.0d0 - 1.0d-8 .and. lat_north >= 90.0d0)) then
                lat_end = i
                found_end = .true.
            end if
        end do

        if (.not. found_start .or. .not. found_end) then
            lat_start = 1
            lat_end = 0  ! Signal no valid indices
        end if
    end subroutine resolve_lat_band_indices

    subroutine build_normalized_cdf(hist_counts, cdf_out)
        real(8), intent(in) :: hist_counts(:)
        real(8), intent(out) :: cdf_out(:)

        integer :: i
        real(8) :: total_area, running

        if (size(cdf_out) /= size(hist_counts) + 1) then
            error stop 'build_normalized_cdf: output size mismatch'
        end if

        total_area = sum(hist_counts)
        if (total_area <= 0.0d0) then
            error stop 'build_normalized_cdf: histogram has non-positive total area'
        end if

        cdf_out(1) = 0.0d0
        running = 0.0d0
        do i = 1, size(hist_counts)
            running = running + hist_counts(i)
            cdf_out(i + 1) = running / total_area
        end do

        cdf_out(size(cdf_out)) = 1.0d0
    end subroutine build_normalized_cdf

    subroutine compute_percentile_thresholds_loglog(edges, cdf_edges, percentiles, thresholds_out)
        real(8), intent(in) :: edges(:)
        real(8), intent(in) :: cdf_edges(:)
        real(8), intent(in) :: percentiles(:)
        real(8), intent(out) :: thresholds_out(:)

        integer :: i, nkeep
        real(8), parameter :: eps = 1.0d-14
        real(8), allocatable :: x_data(:), y_data(:)
        real(8), allocatable :: x_log(:), y_log(:)
        real(8), allocatable :: x_query(:), y_query(:)
        real(8) :: p_clamped, max_requested
        integer :: n_interp

        if (size(edges) /= size(cdf_edges)) then
            error stop 'compute_percentile_thresholds_loglog: edges/cdf size mismatch'
        end if

        if (size(thresholds_out) /= size(percentiles)) then
            error stop 'compute_percentile_thresholds_loglog: output size mismatch'
        end if

        allocate(x_data(size(cdf_edges)))
        allocate(y_data(size(edges)))

        nkeep = 0
        do i = 1, size(cdf_edges)
            if (cdf_edges(i) > 0.0d0) then
                if (nkeep == 0) then
                    nkeep = 1
                    x_data(nkeep) = cdf_edges(i)
                    y_data(nkeep) = edges(i)
                else
                    if (cdf_edges(i) > x_data(nkeep) + eps) then
                        nkeep = nkeep + 1
                        x_data(nkeep) = cdf_edges(i)
                        y_data(nkeep) = edges(i)
                    end if
                end if
            end if
        end do

        if (nkeep < 2) then
            error stop 'compute_percentile_thresholds_loglog: insufficient monotonic CDF support points'
        end if

        max_requested = maxval(percentiles)
        if (x_data(nkeep) < max_requested - eps) then
            write(error_unit, '(A,1X,ES16.8,1X,A,1X,ES16.8)') &
                'compute_percentile_thresholds_loglog: max valid CDF is', x_data(nkeep), &
                'but max requested percentile is', max_requested
            error stop 'compute_percentile_thresholds_loglog: requested percentiles exceed CDF coverage'
        end if

        n_interp = nkeep - 1

        allocate(x_log(n_interp), y_log(n_interp))
        do i = 1, n_interp
            x_log(i) = log10(1.0d0 - x_data(i))
        end do
        y_log = log10(y_data(1:n_interp))

        allocate(x_query(1), y_query(1))
        do i = 1, size(percentiles)
            p_clamped = min(max(percentiles(i), x_data(1)), x_data(n_interp))
            x_query(1) = log10(1.0d0 - p_clamped)
            call interp1(y_query, x_log, y_log, x_query, y_log(n_interp))
            thresholds_out(i) = 10.0d0 ** y_query(1)
        end do

        deallocate(x_data, y_data, x_log, y_log, x_query, y_query)
    end subroutine compute_percentile_thresholds_loglog

    subroutine write_thresholds_ascii_2d(path_out, percentiles, thresholds_in)
        character(len=*), intent(in) :: path_out
        real(8), intent(in) :: percentiles(:)
        real(8), intent(in) :: thresholds_in(:,:)

        integer :: unit_out, i, j, ios_out

        if (size(percentiles) /= size(thresholds_in, 2)) then
            error stop 'write_thresholds_ascii_2d: percentile/threshold size mismatch'
        end if

        if (size(thresholds_in, 1) /= nlat_bands) then
            error stop 'write_thresholds_ascii_2d: lat band size mismatch'
        end if

        open(newunit=unit_out, file=trim(path_out), status='replace', action='write', iostat=ios_out, iomsg=msg)
        if (ios_out /= 0) then
            write(error_unit,*) 'Failed to open output file, iomsg='//trim(msg)
            stop 1
        end if

        ! Header with percentile values
        write(unit_out, '(A)', advance='no') '# lat_band'
        do j = 1, size(percentiles)
            write(unit_out, '(A,F10.5)', advance='no') ' p', percentiles(j)
        end do
        write(unit_out, '(A)') ''

        ! Data rows: one per latitude band
        do i = 1, nlat_bands
            write(unit_out, '(I3)', advance='no') i
            do j = 1, size(percentiles)
                write(unit_out, '(1X,ES16.8)', advance='no') thresholds_in(i, j)
            end do
            write(unit_out, '(A)') ''
        end do

        close(unit_out)
    end subroutine write_thresholds_ascii_2d

end program compute_prate_thresholds_by_lat_band
