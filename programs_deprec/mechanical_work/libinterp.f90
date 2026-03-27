module libinterp

    implicit none
    
    contains
    
    subroutine interp1(yinterp, xv, yv, x, ydefault)

        ! linear interpolation subroutine
    
        double precision, dimension(:), intent(out) :: yinterp
        double precision, dimension(:), intent(in) :: xv, yv, x
        double precision, intent(in) :: ydefault
        integer :: nrowsinterp, nrowsdata, i, j
    
        nrowsinterp = size(x)
        nrowsdata = size(xv)
    
        do i=1,nrowsinterp
            if ((x(i) < xv(1)) .or. (x(i) > xv(nrowsdata))) then
                yinterp(i) = ydefault
            else
                do j=2,nrowsdata
                    if (x(i) <= xv(j)) then
                        yinterp(i) = (x(i)-xv(j-1)) / (xv(j)-xv(j-1)) * (yv(j)-yv(j-1)) + yv(j-1)
                        exit
                    end if
                enddo
            end if
        enddo
    
    
    end subroutine interp1
    
    
    subroutine check(istatus)

        ! catch and throw error according to NF90 library
    
        use netcdf
    
        integer, intent(in) :: istatus
    
        if (istatus /= nf90_noerr) then
            write(*,*) trim(adjustl(nf90_strerror(istatus)))
        end if
    
    
    end subroutine check
    
end module libinterp
    