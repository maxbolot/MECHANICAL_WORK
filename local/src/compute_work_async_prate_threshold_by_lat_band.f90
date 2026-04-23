! async fork of compute_work_async.f90 with latitude-band-specific precipitation-threshold masking
!
! This variant reads thresholds in 18 latitude bands (10-degree wide) x 11 percentiles
! and applies per-band thresholds based on the latitude of each grid point.
!
! Key structural behavior is unchanged from the async implementation:
! NetCDF reads are launched in !$omp task blocks (serialized inside
! !$omp critical(netcdf_io)) so I/O and compute can overlap.
!
! Additional behavior in this fork:
! - Thresholds are read from an ASCII file with 18 x 11 entries (lat_band, percentile).
! - For each grid point, the latitude is used to determine the band index.
! - work/lift are computed once per grid point and then masked per percentile using
!   latitude-band-specific thresholds and PRATEsfc_coarse >= threshold(band, percentile).
! - Outputs include the leading percentile dimension plus global attributes with all thresholds.
!
! Buffer state convention:
!   0 = free
!   2 = compute finished for current chunk
!   3 = in use (being read OR computed)

program compute_work_async_prate_threshold_by_lat_band

    use omp_lib
    use netcdf
    use nc
    use iso_fortran_env, only: error_unit

    implicit none

    character(len=255) :: filenml, msg

    integer, parameter :: chunk_size = 144  ! chunk size along y axis
    integer, parameter :: nbuf = 2          ! number of pipeline buffers
    integer, parameter :: nlat_bands = 18
    integer, parameter :: npercentiles = 11

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
    double precision, dimension(:,:,:,:,:), allocatable :: pr_buffer         ! precipitation rate

    double precision, dimension(:,:,:,:,:), allocatable :: work_out_buffer
    double precision, dimension(:,:,:,:,:), allocatable :: lift_out_buffer

    double precision, dimension(:,:,:,:), allocatable :: work_accum_full     ! (npercentiles, nx, ny, nt)
    double precision, dimension(:,:,:,:), allocatable :: lift_accum_full
    double precision, dimension(:), allocatable :: lat_axis
    real(8), allocatable :: prate_thresholds_2d(:,:)                        ! (nlat_bands, npercentiles)
    real(8), allocatable :: prate_percentiles(:)

    integer :: ncid_hist, ncid_work_out, ncid_lift_out
    integer :: dimid_time, dimid_lev, dimid_lat, dimid_lon, dimid_out_pct
    integer :: varid_pr, varid_work_out, varid_lift_out
    integer :: nx, ny, nlev, nt, ndays
    integer :: i, j, k, iday, ibuf, ystart, x, y, p, ip, ilat
    integer :: fetch_ready(nbuf), compute_ready(nbuf), state(nbuf)
    integer :: iday_in_file, iday_out

    character(len=1024) :: path_hist, path_work_out, path_lift_out, path_thresholds
    character(len=255) :: date_str
    character(len=512) :: attr_name
    real(8) :: lat_south_band, lat_north_band
    real(8) :: work_acc, lift_acc, geometric_thickness, prate_val
    real(8) :: dt_inv

    namelist /config/ path_hist, path_work_out, path_lift_out, path_thresholds

    path_hist = ''
    path_work_out = ''
    path_lift_out = ''
    path_thresholds = ''

    call get_command_argument(1, filenml)
    if (len_trim(filenml) == 0) then
        write(error_unit,*) 'Usage: compute_work_async_prate_threshold_by_lat_band <config.nml>'
        stop 1
    end if

    open(newunit=i, file=trim(adjustl(filenml)), status='old', action='read', iostat=j, iomsg=msg)
    if (j /= 0) then
        write(error_unit,*) 'Failed to open configuration namelist, iomsg='//trim(msg)
        stop 1
    end if

    read(unit=i, nml=config, iostat=j, iomsg=msg)
    if (j /= 0) then
        write(error_unit,*) 'Failed to read namelist, iomsg='//trim(msg)
        stop 1
    end if
    close(i)

    if (len_trim(path_hist) == 0 .or. len_trim(path_work_out) == 0 .or. &
        len_trim(path_lift_out) == 0 .or. len_trim(path_thresholds) == 0) then
        write(error_unit,*) 'Error: all config items required (path_hist, path_work_out, path_lift_out, path_thresholds)'
        stop 1
    end if

    ! Read 2D thresholds (nlat_bands x npercentiles)
    call read_thresholds_ascii_2d(trim(adjustl(path_thresholds)), prate_percentiles, prate_thresholds_2d)

    call check(nf90_open(trim(adjustl(path_hist)), nf90_nowrite, ncid_hist))
    call check(nf90_inquire(ncid_hist, nDimensions=i, nVariables=j, nAttributes=k))
    call check(nf90_inquire_dimension(ncid_hist, 1, len=nx))
    call check(nf90_inquire_dimension(ncid_hist, 2, len=ny))
    call check(nf90_inquire_dimension(ncid_hist, 3, len=nlev))
    call check(nf90_inquire_dimension(ncid_hist, 4, len=nt))

    allocate(lat_axis(ny))
    call check(nf90_inq_varid(ncid_hist, 'grid_yt_coarse', varid_pr))
    call check(nf90_get_var(ncid_hist, varid_pr, lat_axis))

    ! Allocate buffers
    allocate(dz_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(temp_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omega_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qv_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qw_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qr_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qi_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qs_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(qg_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omt_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqv_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqw_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqr_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqi_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqs_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(omqg_buffer(nx, chunk_size, nlev, 1, nbuf))
    allocate(pr_buffer(nx, chunk_size, 1, 1, nbuf))

    allocate(work_out_buffer(npercentiles, nx, chunk_size, 1, nbuf))
    allocate(lift_out_buffer(npercentiles, nx, chunk_size, 1, nbuf))

    ! Count days in the file (for time dimension)
    ndays = nt
    allocate(work_accum_full(npercentiles, nx, ny, nt))
    allocate(lift_accum_full(npercentiles, nx, ny, nt))

    work_accum_full = 0.0d0
    lift_accum_full = 0.0d0

    state = 0
    fetch_ready = 0
    compute_ready = 0

    dt_inv = 1.0d0 / 3.0d0  ! time step in hours: assume 3 hours

    !$omp parallel num_threads(2) &
    !$omp default(shared)
    !$omp do collapse(1) schedule(static, 1)
    do ystart = 1, ny, chunk_size
        do ibuf = 1, nbuf
            !$omp critical(netcdf_io)
            state(ibuf) = 3
            !$omp end critical(netcdf_io)

            call read_chunk_from_hist(ncid_hist, ystart, min(ystart + chunk_size - 1, ny), 1, nt, &
                                      dz_buffer(:,:,:,:,ibuf), temp_buffer(:,:,:,:,ibuf), &
                                      omega_buffer(:,:,:,:,ibuf), &
                                      qv_buffer(:,:,:,:,ibuf), qw_buffer(:,:,:,:,ibuf), &
                                      qr_buffer(:,:,:,:,ibuf), qi_buffer(:,:,:,:,ibuf), &
                                      qs_buffer(:,:,:,:,ibuf), qg_buffer(:,:,:,:,ibuf), &
                                      omt_buffer(:,:,:,:,ibuf), omqv_buffer(:,:,:,:,ibuf), &
                                      omqw_buffer(:,:,:,:,ibuf), omqr_buffer(:,:,:,:,ibuf), &
                                      omqi_buffer(:,:,:,:,ibuf), omqs_buffer(:,:,:,:,ibuf), &
                                      omqg_buffer(:,:,:,:,ibuf), pr_buffer(:,:,:,:,ibuf))

            !$omp critical(netcdf_io)
            state(ibuf) = 2
            !$omp end critical(netcdf_io)
        end do

        do iday = 1, nt
            do ibuf = 1, nbuf
                if (state(ibuf) /= 2) cycle
                !$omp critical(netcdf_io)
                state(ibuf) = 3
                !$omp end critical(netcdf_io)

                ! Compute work/lift with latitude-band-specific threshold masking
                do y = 1, min(chunk_size, ny - ystart + 1)
                    do x = 1, nx
                        ! Determine latitude band for this grid point
                        ilat = get_lat_band_index(lat_axis(ystart + y - 1))

                        work_acc = 0.0d0
                        lift_acc = 0.0d0

                        do p = nlev, 1, -1
                            geometric_thickness = -dz_buffer(x,y,p,1,ibuf)
                            work_acc = work_acc - &
                                (omt_buffer(x,y,p,1,ibuf) / temp_buffer(x,y,p,1,ibuf) + 1.61d0 * (omqv_buffer(x,y,p,1,ibuf) + omega_buffer(x,y,p,1,ibuf) * qv_buffer(x,y,p,1,ibuf))) * geometric_thickness
                            lift_acc = lift_acc - &
                                (omega_buffer(x,y,p,1,ibuf) * (qv_buffer(x,y,p,1,ibuf) + qw_buffer(x,y,p,1,ibuf) + qr_buffer(x,y,p,1,ibuf) + qi_buffer(x,y,p,1,ibuf) + qs_buffer(x,y,p,1,ibuf) + qg_buffer(x,y,p,1,ibuf)) &
                                + omqv_buffer(x,y,p,1,ibuf) + omqw_buffer(x,y,p,1,ibuf) + omqr_buffer(x,y,p,1,ibuf) + omqi_buffer(x,y,p,1,ibuf) + omqs_buffer(x,y,p,1,ibuf) + omqg_buffer(x,y,p,1,ibuf)) * geometric_thickness
                        end do

                        prate_val = pr_buffer(x,y,1,1,ibuf)

                        ! Threshold masking per percentile rank using latitude-band-specific thresholds
                        do ip = 1, npercentiles
                            if (prate_val >= prate_thresholds_2d(ilat, ip)) then
                                work_out_buffer(ip,x,y,1,ibuf) = work_acc
                                lift_out_buffer(ip,x,y,1,ibuf) = lift_acc
                            else
                                work_out_buffer(ip,x,y,1,ibuf) = 0.0d0
                                lift_out_buffer(ip,x,y,1,ibuf) = 0.0d0
                            end if
                        end do
                    end do
                end do

                ! Accumulate into daily bins
                work_accum_full(:,:,ystart:ystart+min(chunk_size, ny-ystart+1)-1,iday) = &
                    work_accum_full(:,:,ystart:ystart+min(chunk_size, ny-ystart+1)-1,iday) + &
                    work_out_buffer(:,:,1:min(chunk_size, ny-ystart+1),1,ibuf)
                lift_accum_full(:,:,ystart:ystart+min(chunk_size, ny-ystart+1)-1,iday) = &
                    lift_accum_full(:,:,ystart:ystart+min(chunk_size, ny-ystart+1)-1,iday) + &
                    lift_out_buffer(:,:,1:min(chunk_size, ny-ystart+1),1,ibuf)

                !$omp critical(netcdf_io)
                state(ibuf) = 0
                !$omp end critical(netcdf_io)
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    ! Write output files
    call check(nf90_create(trim(adjustl(path_work_out)), nf90_clobber, ncid_work_out))
    call check(nf90_create(trim(adjustl(path_lift_out)), nf90_clobber, ncid_lift_out))

    call check(nf90_def_dim(ncid_work_out, 'percentile', npercentiles, dimid_out_pct))
    call check(nf90_def_dim(ncid_work_out, 'lon', nx, dimid_lon))
    call check(nf90_def_dim(ncid_work_out, 'lat', ny, dimid_lat))
    call check(nf90_def_dim(ncid_work_out, 'time', ndays, dimid_time))

    call check(nf90_def_dim(ncid_lift_out, 'percentile', npercentiles, dimid_out_pct))
    call check(nf90_def_dim(ncid_lift_out, 'lon', nx, dimid_lon))
    call check(nf90_def_dim(ncid_lift_out, 'lat', ny, dimid_lat))
    call check(nf90_def_dim(ncid_lift_out, 'time', ndays, dimid_time))

    ! Define variables
    call check(nf90_def_var(ncid_work_out, 'work', nf90_double, (/dimid_out_pct, dimid_lon, dimid_lat, dimid_time/), varid_work_out))
    call check(nf90_def_var(ncid_lift_out, 'lift', nf90_double, (/dimid_out_pct, dimid_lon, dimid_lat, dimid_time/), varid_lift_out))

    ! Store latitude-band-specific thresholds as global attributes for traceability
    do ilat = 1, nlat_bands
        do ip = 1, npercentiles
            write(attr_name, '(A,I2.2,A,I2.2,A)') 'threshold_band', ilat, '_p', nint(prate_percentiles(ip)*1000), '_prate'
            call check(nf90_put_att(ncid_work_out, nf90_global, trim(adjustl(attr_name)), prate_thresholds_2d(ilat, ip)))
            call check(nf90_put_att(ncid_lift_out, nf90_global, trim(adjustl(attr_name)), prate_thresholds_2d(ilat, ip)))
        end do
    end do

    call check(nf90_enddef(ncid_work_out))
    call check(nf90_enddef(ncid_lift_out))

    ! Write data
    call check(nf90_put_var(ncid_work_out, varid_work_out, work_accum_full(:,:,:,1:ndays), &
                            start=(/1,1,1,1/), count=(/npercentiles,nx,ny,ndays/)))
    call check(nf90_put_var(ncid_lift_out, varid_lift_out, lift_accum_full(:,:,:,1:ndays), &
                            start=(/1,1,1,1/), count=(/npercentiles,nx,ny,ndays/)))

    call check(nf90_close(ncid_work_out))
    call check(nf90_close(ncid_lift_out))
    call check(nf90_close(ncid_hist))

    if (allocated(dz_buffer)) deallocate(dz_buffer)
    if (allocated(pr_buffer)) deallocate(pr_buffer)
    if (allocated(work_out_buffer)) deallocate(work_out_buffer)
    if (allocated(lift_out_buffer)) deallocate(lift_out_buffer)
    if (allocated(work_accum_full)) deallocate(work_accum_full)
    if (allocated(lift_accum_full)) deallocate(lift_accum_full)
    if (allocated(prate_thresholds_2d)) deallocate(prate_thresholds_2d)
    if (allocated(prate_percentiles)) deallocate(prate_percentiles)
    if (allocated(lat_axis)) deallocate(lat_axis)

contains

    integer function get_lat_band_index(lat_value)
        real(8), intent(in) :: lat_value
        integer :: ilat

        ! Map latitude to 10-degree band index (1-18)
        ! Band 1: [-90, -80], Band 2: [-80, -70], ..., Band 18: [80, 90]
        ilat = int((lat_value + 90.0d0) / 10.0d0) + 1
        ilat = max(1, min(nlat_bands, ilat))
        get_lat_band_index = ilat
    end function get_lat_band_index

    subroutine read_chunk_from_hist(ncid_in, ystart, yend, tstart, tend, &
                                    dz, temp, omega, qv, qw, qr, qi, qs, qg, &
                                    omt, omqv, omqw, omqr, omqi, omqs, omqg, pr)
        integer, intent(in) :: ncid_in, ystart, yend, tstart, tend
        real(8), intent(out) :: dz(:,:,:,:)
        real(8), intent(out) :: temp(:,:,:,:)
        real(8), intent(out) :: omega(:,:,:,:)
        real(8), intent(out) :: qv(:,:,:,:)
        real(8), intent(out) :: qw(:,:,:,:)
        real(8), intent(out) :: qr(:,:,:,:)
        real(8), intent(out) :: qi(:,:,:,:)
        real(8), intent(out) :: qs(:,:,:,:)
        real(8), intent(out) :: qg(:,:,:,:)
        real(8), intent(out) :: omt(:,:,:,:)
        real(8), intent(out) :: omqv(:,:,:,:)
        real(8), intent(out) :: omqw(:,:,:,:)
        real(8), intent(out) :: omqr(:,:,:,:)
        real(8), intent(out) :: omqi(:,:,:,:)
        real(8), intent(out) :: omqs(:,:,:,:)
        real(8), intent(out) :: omqg(:,:,:,:)
        real(8), intent(out) :: pr(:,:,:,:)

        integer :: varid
        integer :: ylen

        !$omp critical(netcdf_io)
        ylen = yend - ystart + 1

        call check(nf90_inq_varid(ncid_in, 'dz', varid))
        call check(nf90_get_var(ncid_in, varid, dz, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'temp', varid))
        call check(nf90_get_var(ncid_in, varid, temp, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omega', varid))
        call check(nf90_get_var(ncid_in, varid, omega, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qv', varid))
        call check(nf90_get_var(ncid_in, varid, qv, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qw', varid))
        call check(nf90_get_var(ncid_in, varid, qw, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qr', varid))
        call check(nf90_get_var(ncid_in, varid, qr, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qi', varid))
        call check(nf90_get_var(ncid_in, varid, qi, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qs', varid))
        call check(nf90_get_var(ncid_in, varid, qs, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'qg', varid))
        call check(nf90_get_var(ncid_in, varid, qg, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omT', varid))
        call check(nf90_get_var(ncid_in, varid, omt, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqv', varid))
        call check(nf90_get_var(ncid_in, varid, omqv, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqw', varid))
        call check(nf90_get_var(ncid_in, varid, omqw, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqr', varid))
        call check(nf90_get_var(ncid_in, varid, omqr, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqi', varid))
        call check(nf90_get_var(ncid_in, varid, omqi, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqs', varid))
        call check(nf90_get_var(ncid_in, varid, omqs, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'omqg', varid))
        call check(nf90_get_var(ncid_in, varid, omqg, start=(/1,ystart,1,tstart/), count=(/nx,ylen,nlev,1/)))

        call check(nf90_inq_varid(ncid_in, 'PRATEsfc_coarse', varid))
        call check(nf90_get_var(ncid_in, varid, pr, start=(/1,ystart,1,tstart/), count=(/nx,ylen,1,1/)))

        !$omp end critical(netcdf_io)
    end subroutine read_chunk_from_hist

    subroutine read_thresholds_ascii_2d(path_in, percentiles, thresholds)
        character(len=*), intent(in) :: path_in
        real(8), allocatable, intent(out) :: percentiles(:)
        real(8), allocatable, intent(out) :: thresholds(:,:)

        integer :: unit_in, io, ilat, ip, npercentiles_read
        character(len=256) :: line
        real(8), allocatable :: row_data(:)

        ! First pass: count data rows (should be nlat_bands)
        open(newunit=unit_in, file=trim(path_in), status='old', action='read', iostat=io, iomsg=msg)
        if (io /= 0) then
            write(error_unit,*) 'Failed to open threshold file, iomsg='//trim(msg)
            stop 1
        end if

        ilat = 0
        do
            read(unit_in, '(A)', iostat=io) line
            if (io /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            ilat = ilat + 1
        end do

        if (ilat /= nlat_bands) then
            write(error_unit,*) 'Error: expected', nlat_bands, 'lat bands, found', ilat
            stop 1
        end if

        rewind(unit_in)

        ! Count percentiles from first data row
        npercentiles_read = 0
        do
            read(unit_in, '(A)', iostat=io) line
            if (io /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            ! Parse first data row to count fields
            allocate(row_data(50))
            read(line, *, iostat=io) row_data
            if (io == 0) then
                npercentiles_read = count(row_data /= 0.0d0) - 1  ! Subtract band index column
            end if
            deallocate(row_data)
            exit
        end do

        if (npercentiles_read <= 0) then
            error stop 'Could not determine percentile count from threshold file'
        end if

        if (.not. allocated(percentiles)) then
            allocate(percentiles(npercentiles_read))
            percentiles = [(dble(i)/npercentiles_read*0.999d0, i=1, npercentiles_read)]  ! Placeholder
        end if

        allocate(thresholds(nlat_bands, npercentiles_read))
        allocate(row_data(npercentiles_read + 1))

        rewind(unit_in)
        ilat = 0
        do
            read(unit_in, '(A)', iostat=io) line
            if (io /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            ilat = ilat + 1
            read(line, *, iostat=io) row_data
            if (io /= 0) cycle
            thresholds(ilat, :) = row_data(2:npercentiles_read+1)
        end do

        close(unit_in)
    end subroutine read_thresholds_ascii_2d

end program compute_work_async_prate_threshold_by_lat_band
