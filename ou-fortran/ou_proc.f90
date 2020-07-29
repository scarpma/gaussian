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


subroutine compute_struct(struct,db) ! in :compute_struct:ou_proc.f90
    
    implicit none
    real(8), intent(inout), dimension(:,:) :: struct
    real(8), intent(in), dimension(:,:) :: db
    real(8), dimension(:), allocatable :: diff
    integer :: tt, npart, times, tau_index !, kk
    integer, dimension(23) :: taus
    real(8) :: p2, p4, p6
    
    npart = size(db,2)
    times = size(db,1)

    allocate(diff(npart))

    taus = (/ 2,3,4,5,6,8,11,14,19,25,34,45,59,79,104,138,184,244,323,429,568,754,1000 /)

    do tau_index = 1, 23
        write(*, fmt="(i0,2X)", advance="no") taus(tau_index)
        diff(:) = 0.
        p2 = 0.
        p4 = 0.
        p6 = 0.
        do tt = 1, times - taus(tau_index)
           diff(:) = db(tt+taus(tau_index),:) - db(tt,:)
           p2 = p2 + sum(diff**2)
           p4 = p4 + sum(diff**4)
           p6 = p6 + sum(diff**6)
        end do
        struct(tau_index,2) = p2 / ((times-taus(tau_index))*npart)
        struct(tau_index,3) = p4 / ((times-taus(tau_index))*npart)
        struct(tau_index,4) = p6 / ((times-taus(tau_index))*npart)
    end do
    struct(:,1) = taus
    deallocate(diff)

end subroutine
