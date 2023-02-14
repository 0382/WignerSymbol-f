program test
    use iso_fortran_env, only: stderr => error_unit
    use testdrive
    use WignerSymbol
    implicit none

    type :: Moshinsky_case
        integer :: N, L, nr, lr, n1, l1, n2, l2, lambda
    end type
    type(Moshinsky_case), dimension(13) :: Moshinsky_test_set
    real(kind=8), dimension(13) :: Moshinsky_test_set_result

    integer :: stat, i
    type(testsuite_type), dimension(:), allocatable :: testsuites
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0
    call wigner_init(100, "nmax", 3)
    call init_mosh_test_set()

    testsuites = [ &
                 new_testsuite("test Moshinsky", collect_Moshinsky) &
                 ]
    do i = 1, size(testsuites)
        write (stderr, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, stderr, stat)
    end do
    if (stat > 0) then
        write (stderr, '(I0,1x,a)') stat, "test(s) failed!"
        error stop
    end if
contains
    pure function new_mosh(N, L, nr, lr, n1, l1, n2, l2, lambda) result(ans)
        integer, intent(in) :: N, L, nr, lr, n1, l1, n2, l2, lambda
        type(Moshinsky_case) :: ans
        ans%N = N
        ans%L = L
        ans%nr = nr
        ans%lr = lr
        ans%n1 = n1
        ans%l1 = l1
        ans%n2 = n2
        ans%l2 = l2
        ans%lambda = lambda
    end function new_mosh

    subroutine init_mosh_test_set()
        Moshinsky_test_set = &
            [new_mosh(0, 0, 0, 0, 0, 0, 0, 0, 0), &
             new_mosh(0, 1, 0, 0, 0, 0, 0, 1, 1), &
             new_mosh(0, 0, 0, 1, 0, 0, 0, 1, 1), &
             new_mosh(0, 2, 0, 0, 0, 0, 0, 2, 2), &
             new_mosh(0, 1, 0, 1, 0, 0, 0, 2, 2), &
             new_mosh(0, 0, 0, 2, 0, 0, 0, 2, 2), &
             new_mosh(1, 0, 0, 0, 0, 1, 0, 1, 0), &
             new_mosh(0, 1, 0, 1, 0, 1, 0, 1, 0), &
             new_mosh(0, 0, 1, 0, 0, 1, 0, 1, 0), &
             new_mosh(0, 1, 0, 1, 0, 1, 0, 1, 1), &
             new_mosh(0, 2, 0, 0, 0, 1, 0, 1, 2), &
             new_mosh(0, 1, 0, 1, 0, 1, 0, 1, 2), &
             new_mosh(0, 0, 0, 2, 0, 1, 0, 1, 2)]
    end subroutine init_mosh_test_set

    subroutine calc_exact_mosh_result(tan_beta)
        real(kind=8), intent(in) :: tan_beta
        real(kind=8) :: sec_beta, sin_beta, cos_beta
        sec_beta = sqrt(1.0d0 + tan_beta*tan_beta)
        cos_beta = 1.0d0/sec_beta
        sin_beta = tan_beta/sec_beta
        Moshinsky_test_set_result(1) = 1.0d0
        Moshinsky_test_set_result(2) = cos_beta
        Moshinsky_test_set_result(3) = -sin_beta
        Moshinsky_test_set_result(4) = cos_beta*cos_beta
        Moshinsky_test_set_result(5) = -sqrt(2.0d0)*sin_beta*cos_beta
        Moshinsky_test_set_result(6) = sin_beta*sin_beta
        Moshinsky_test_set_result(7) = sqrt(2.0d0)*sin_beta*cos_beta
        Moshinsky_test_set_result(8) = cos_beta*cos_beta - sin_beta*sin_beta
        Moshinsky_test_set_result(9) = -sqrt(2.0d0)*sin_beta*cos_beta
        Moshinsky_test_set_result(10) = -1.0d0
        Moshinsky_test_set_result(11) = sqrt(2.0d0)*sin_beta*cos_beta
        Moshinsky_test_set_result(12) = cos_beta*cos_beta - sin_beta*sin_beta
        Moshinsky_test_set_result(13) = -sqrt(2.0d0)*sin_beta*cos_beta
    end subroutine calc_exact_mosh_result

    subroutine collect_Moshinsky(testsuite)
        type(unittest_type), dimension(:), allocatable, intent(out) :: testsuite
        testsuite = &
            [new_unittest("moshinsky: tan(beta) = 1/3", test_Moshinsky_1), &
             new_unittest("moshinsky: tan(beta) = 1/2", test_Moshinsky_2), &
             new_unittest("moshinsky: tan(beta) = 1", test_Moshinsky_3), &
             new_unittest("moshinsky: tan(beta) = 2", test_Moshinsky_4), &
             new_unittest("moshinsky: tan(beta) = 3", test_Moshinsky_5) &
             ]
    end subroutine collect_Moshinsky

    subroutine test_Moshinsky_1(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=8) :: my, exact, diff
        type(Moshinsky_case) :: m
        real(kind=8), parameter :: tan_beta = 1.0d0/3.0d0
        integer :: idx
        call calc_exact_mosh_result(tan_beta)
        diff = 0.0d0
        do idx = 1, 13
            m = Moshinsky_test_set(idx)
            my = Moshinsky(m%N, m%L, m%nr, m%lr, m%n1, m%l1, m%n2, m%l2, m%lambda, tan_beta)
            exact = Moshinsky_test_set_result(idx)
            diff = diff + abs(my - exact)
        end do
        if (diff < 1d-12) diff = 0.0d0
        call check(error, diff, 0.0d0)
    end subroutine test_Moshinsky_1

    subroutine test_Moshinsky_2(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=8) :: my, exact, diff
        type(Moshinsky_case) :: m
        real(kind=8), parameter :: tan_beta = 1.0d0/2.0d0
        integer :: idx
        call calc_exact_mosh_result(tan_beta)
        diff = 0.0d0
        do idx = 1, 13
            m = Moshinsky_test_set(idx)
            my = Moshinsky(m%N, m%L, m%nr, m%lr, m%n1, m%l1, m%n2, m%l2, m%lambda, tan_beta)
            exact = Moshinsky_test_set_result(idx)
            diff = diff + abs(my - exact)
        end do
        if (diff < 1d-12) diff = 0.0d0
        call check(error, diff, 0.0d0)
    end subroutine test_Moshinsky_2

    subroutine test_Moshinsky_3(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=8) :: my, exact, diff
        type(Moshinsky_case) :: m
        real(kind=8), parameter :: tan_beta = 1.0d0
        integer :: idx
        call calc_exact_mosh_result(tan_beta)
        diff = 0.0d0
        do idx = 1, 13
            m = Moshinsky_test_set(idx)
            my = Moshinsky(m%N, m%L, m%nr, m%lr, m%n1, m%l1, m%n2, m%l2, m%lambda)
            exact = Moshinsky_test_set_result(idx)
            diff = diff + abs(my - exact)
        end do
        if (diff < 1d-12) diff = 0.0d0
        call check(error, diff, 0.0d0)
    end subroutine test_Moshinsky_3

    subroutine test_Moshinsky_4(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=8) :: my, exact, diff
        type(Moshinsky_case) :: m
        real(kind=8), parameter :: tan_beta = 2.0d0
        integer :: idx
        call calc_exact_mosh_result(tan_beta)
        diff = 0.0d0
        do idx = 1, 13
            m = Moshinsky_test_set(idx)
            my = Moshinsky(m%N, m%L, m%nr, m%lr, m%n1, m%l1, m%n2, m%l2, m%lambda, tan_beta)
            exact = Moshinsky_test_set_result(idx)
            diff = diff + abs(my - exact)
        end do
        if (diff < 1d-12) diff = 0.0d0
        call check(error, diff, 0.0d0)
    end subroutine test_Moshinsky_4

    subroutine test_Moshinsky_5(error)
        type(error_type), allocatable, intent(out) :: error
        real(kind=8) :: my, exact, diff
        type(Moshinsky_case) :: m
        real(kind=8), parameter :: tan_beta = 3.0d0
        integer :: idx
        call calc_exact_mosh_result(tan_beta)
        diff = 0.0d0
        do idx = 1, 13
            m = Moshinsky_test_set(idx)
            my = Moshinsky(m%N, m%L, m%nr, m%lr, m%n1, m%l1, m%n2, m%l2, m%lambda, tan_beta)
            exact = Moshinsky_test_set_result(idx)
            diff = diff + abs(my - exact)
        end do
        if (diff < 1d-12) diff = 0.0d0
        call check(error, diff, 0.0d0)
    end subroutine test_Moshinsky_5
end program test
