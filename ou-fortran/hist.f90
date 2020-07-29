subroutine hist(db,bins) ! in :hist:hist.f90
    
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
