subroutine ou_generator(tau,h,ou) ! in :ou_generator:ou_proc.f90
    
    implicit none
    real(8), intent(in) :: tau, h
    real(8), intent(inout), dimension(:,:) :: ou
    integer :: tt, kk, NN, len
    real(8) :: u
    
    NN = size(ou,2)
    len = size(ou,1)

    do kk = 1, NN
            call r8_normal(u)
            ou(1,kk) = sqrt(1./(2.*tau)) * u            
        do tt = 2,len
            call r8_normal(u)
            ou(tt,kk) = ou(tt-1,kk)*exp(-h/tau)+sqrt((1.-exp(-2.*h/tau))/(2.*tau)) * u
        end do
    end do

end subroutine


subroutine r8_normal(uno) ! in :r8_normal:ou_proc.f90

    implicit none
    real(8), intent(inout) :: uno!, due
    real(8) :: u1, u2
    real(8), parameter :: pi = 3.141592653
    call random_number(u1)
    call random_number(u2)
    uno = sqrt(-2. * log(u1)) * cos(2. * pi * u2)
    !due = sqrt(-2. * log(u1)) * sin(2. * pi * u2)
    
end subroutine



subroutine compute_struct_mixed(struct,db) ! in :compute_struct_mixed:ou_proc.f90
    
    implicit none
    real(8), intent(inout), dimension(:,:) :: struct
    real(8), intent(in), dimension(:,:,:) :: db
    real(8), dimension(:), allocatable :: diffx, diffy, diffz
    integer(8) :: tt, npart, times, tau_index, n_comp !, kk
    integer(8), dimension(34) :: taus
    real(8) :: x2y2, x2z2, y2z2, x2, y2, z2
    
    n_comp = size(db,1)
    if (n_comp .NE. 3) then
        write(*,*) "PROBLEM, DETECTED N_COMPONENTS DIFFERENT FROM 3"
        write(*,*) n_comp
    end if
    npart = size(db,3)
    times = size(db,2)
    !write(*,*) "From fortran: npart=", npart, " times=", times
    !write(*,*) "size_struct", size(struct,1), size(struct,2)

    allocate(diffx(npart))
    allocate(diffy(npart))
    allocate(diffz(npart))

    taus = (/ 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 18, 22, 26, 31, 38,   &
            46, 55, 66, 79, 95, 114, 137, 164, 197, 237, 284, 341, 410, &
            492, 590, 708, 850, 1020 /)

    do tau_index = 1, 34
        !write(*,*) "From fortran: ", taus(tau_index)
        write(*,*) "From fortran: tau_index=", tau_index
        diffx(:) = 0.
        diffy(:) = 0.
        diffz(:) = 0.
        x2y2 = 0.
        y2z2 = 0.
        x2z2 = 0.
        x2 = 0.
        y2 = 0.
        z2 = 0.
        do tt = 1, times - taus(tau_index)
            diffx(:) = db(1,tt+taus(tau_index),:) - db(1,tt,:)
            diffy(:) = db(2,tt+taus(tau_index),:) - db(2,tt,:)
            diffz(:) = db(3,tt+taus(tau_index),:) - db(3,tt,:)
            diffx = diffx**2
            diffy = diffy**2
            diffz = diffz**2
            x2 = x2 + sum(diffx)
            y2 = y2 + sum(diffy)
            z2 = z2 + sum(diffz)
            x2y2 = x2y2 + sum(diffx*diffy)
            y2z2 = y2z2 + sum(diffy*diffz)
            x2z2 = x2z2 + sum(diffx*diffz)
        end do
        !write(*,*) "From fortran: p2:", p2
        !write(*,*) "From fortran: p4:", p4
        !write(*,*) "From fortran: p6:", p6
        struct(tau_index,2) = x2 / ((times-taus(tau_index))*npart)
        struct(tau_index,3) = y2 / ((times-taus(tau_index))*npart)
        struct(tau_index,4) = z2 / ((times-taus(tau_index))*npart)

        struct(tau_index,5) = x2y2 / ((times-taus(tau_index))*npart)
        struct(tau_index,6) = y2z2 / ((times-taus(tau_index))*npart)
        struct(tau_index,7) = x2z2 / ((times-taus(tau_index))*npart)
        !write(*,*) "From fortran: time - taus(index) )", &
        !    times-taus(tau_index)
        !write(*,*) "From fortran: s2:", struct(tau_index,2) 
        !write(*,*) "From fortran: s4:", struct(tau_index,3) 
        !write(*,*) "From fortran: s6:", struct(tau_index,4) 
    end do
    !write(*,*) "taus still not written"
    struct(:,1) = taus
    !write(*,*) "taus written"
    !write(*,*) struct
    deallocate(diffx)
    deallocate(diffy)
    deallocate(diffz)
    write(*,*) "Done."

end subroutine



subroutine compute_struct(struct,db) ! in :compute_struct:ou_proc.f90
    
    implicit none
    real(8), intent(inout), dimension(:,:) :: struct
    real(8), intent(in), dimension(:,:) :: db
    real(8), dimension(:), allocatable :: diff
    integer(8) :: tt, npart, times, tau_index !, kk
    integer(8), dimension(34) :: taus
    real(8) :: p2, p4, p6, p8, p10, p12
    
    npart = size(db,2)
    times = size(db,1)
    !write(*,*) "From fortran: npart=", npart, " times=", times
    !write(*,*) "size_struct", size(struct,1), size(struct,2)

    allocate(diff(npart))

    taus = (/ 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 18, 22, 26, 31, 38,   &
            46, 55, 66, 79, 95, 114, 137, 164, 197, 237, 284, 341, 410, &
            492, 590, 708, 850, 1020 /)

    do tau_index = 1, 34
        !write(*,*) "From fortran: ", taus(tau_index)
        write(*,*) "From fortran: tau_index=", tau_index
        diff(:) = 0.
        p2 = 0.
        p4 = 0.
        p6 = 0.
        p8 = 0.
        p10 = 0.
        p12 = 0.
        do tt = 1, times - taus(tau_index)
           diff(:) = db(tt+taus(tau_index),:) - db(tt,:)
           p2 = p2 + sum(diff**2)
           p4 = p4 + sum(diff**4)
           p6 = p6 + sum(diff**6)
           p8 = p8 + sum(diff**8)
           p10 = p10 + sum(diff**10)
           p12 = p12 + sum(diff**12)
        end do
        !write(*,*) "From fortran: p2:", p2
        !write(*,*) "From fortran: p4:", p4
        !write(*,*) "From fortran: p6:", p6
        struct(tau_index,2) = p2 / ((times-taus(tau_index))*npart)
        struct(tau_index,3) = p4 / ((times-taus(tau_index))*npart)
        struct(tau_index,4) = p6 / ((times-taus(tau_index))*npart)
        struct(tau_index,5) = p8 / ((times-taus(tau_index))*npart)
        struct(tau_index,6) = p10 / ((times-taus(tau_index))*npart)
        struct(tau_index,7) = p12 / ((times-taus(tau_index))*npart)
        !write(*,*) "From fortran: time - taus(index) )", &
        !    times-taus(tau_index)
        !write(*,*) "From fortran: s2:", struct(tau_index,2) 
        !write(*,*) "From fortran: s4:", struct(tau_index,3) 
        !write(*,*) "From fortran: s6:", struct(tau_index,4) 
    end do
    !write(*,*) "taus still not written"
    struct(:,1) = taus
    !write(*,*) "taus written"
    !write(*,*) struct
    deallocate(diff)
    write(*,*) "Done."

end subroutine

subroutine cor_func(db,acf) ! in :corr_func:ou_proc.f90
    
    implicit none
    real(8), intent(inout), dimension(:,:) :: db
    real(8), intent(inout), dimension(:) :: acf
    real(8) :: mean
    integer :: npart, times, tt, tau, kk

    npart = size(db,2)
    times = size(db,1)
    if (size(acf) .ne. times) then
        call exit()
    end if

    ! CALCOLO FUNZIONE DI AUTOCORRELAZIONE 
    do kk = 1, npart
        acf(:) = 0.
        mean = 0.
        do tau = 0, times-1
            mean = mean + db(1+tau,kk)
            do tt = 1, times - tau
                acf(1+tau) = acf(1+tau) + db(tt+tau,kk)*db(tt,kk) 
            end do
            acf(1 + tau) = acf(1 + tau) / (times - tau)
        end do
        mean = mean / times
        db(:,kk) = ( acf(:) - mean**2. ) / ( acf(1) - mean**2. )
    end do

end subroutine

!subroutine cor_func(db,acf) ! in :corr_func:ou_proc.f90
!    
!    implicit none
!    real(8), intent(in), dimension(:,:) :: db
!    real(8), intent(inout), dimension(:) :: acf
!    real(8) :: std, mean
!    integer :: npart, times, tt, tau, kk
!
!    std = 0.
!    mean = 0.
!    acf(:) = 0.
!
!    npart = size(db,2)
!    times = size(db,1)
!    if (size(acf) .ne. times) then
!        call exit()
!    end if
!
!    ! CALCOLO FUNZIONE DI AUTOCORRELAZIONE 
!    do kk = 1, npart
!        do tau = 0, times-1
!            do tt = 1, times - tau
!                acf(1+tau) = acf(1+tau) + db(tt+tau,kk)*db(tt,kk) 
!            end do
!        end do
!    end do
!    do tau = 0, times-1
!        acf(1 + tau) = acf(1 + tau) / ((times - tau)*npart)
!    end do
!
!    ! CALCOLO MEDIA E VARIANZA
!    do kk = 1, npart
!        do tt = 1, times
!            mean = mean + db(tt,kk)
!        end do
!    end do
!
!    mean = mean / (npart * times)
!    std = acf(1) - mean**2.
!
!    acf(:) = ( acf(:) - mean**2. ) / std
!
!end subroutine



subroutine correlate(db,acf) ! in :correlate:ou_proc.f90
    
    implicit none
    real(8), intent(in), dimension(:) :: db
    real(8), intent(inout), dimension(:) :: acf
    real(8) :: mean
    integer :: times, tt, tau

    mean = 0.

    times = size(db)
    if (size(acf) .ne. times) then
        call exit()
    end if

    ! CALCOLO FUNZIONE DI AUTOCORRELAZIONE 
    do tau = 0, times-1
        do tt = 1, times - tau
            acf(1+tau) = acf(1+tau) + db(tt+tau)*db(tt)
        end do
        acf(1 + tau) = acf(1 + tau) / (times - tau)
    end do

    ! CALCOLO MEDIA E VARIANZA
    do tt = 1, times
        mean = mean + db(tt)
    end do

    mean = mean / times

    acf = acf - mean**2.

    acf = acf / acf(1)

end subroutine


subroutine histogram(xx, hist, bins) ! in :histogram:ou_proc.f90
    
    ! xx    : array with samples
    ! hist  : array to store histogram
    ! bins  : array to store bin edges

    implicit none
    real(8), intent(in), dimension(:) :: xx
    real(8), intent(inout), dimension(:) :: bins
    integer, intent(inout), dimension(:) :: hist
    integer :: nbins, nsamples, ii, kk
    real(8) :: xmax, xmin, binw

    ! INITIALIZATION AND PREPARATION
    ! nbins is the number of intervals

    nbins = size(bins)
    nsamples = size(xx)

    xmax = maxval(xx)
    xmin = minval(xx)
    binw = ( xmax - xmin ) / nbins
    xmax = xmax + binw ! aggiungo un bin in più agli estremi così non cado mai
    ! fuori dagli estremi nel calcolare kk
    xmin = xmin - binw
    binw = ( xmax - xmin ) / nbins

    do ii=1,nbins
        bins(ii) = xmin + (ii-0.5) * binw
    end do

    ! COMPUTATION OF HISTOGRAM
    
    do ii = 1, nsamples
        kk = int( ( xx(ii) - xmin ) / binw ) + 1
        hist(kk) = hist(kk) + 1
    end do

end subroutine

