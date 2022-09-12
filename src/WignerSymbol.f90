module WignerSymbol
    implicit none
    private
    integer(kind=8), dimension(:), allocatable :: i64_binominal_data
    integer :: m_nmax = 0
    real(kind=8), dimension(:), allocatable :: m_binomial_data

    public :: wigner_init
    public :: binomial
    public :: CG
    public :: wigner3j
    public :: wigner6j
    public :: Racah
    public :: wigner9j

contains
    ! judge if a number is a odd number
    integer pure function isodd(n) result(ans)
        integer, intent(in) :: n
        ans = iand(n, 1)
    end function isodd
    ! judge if a number is a even number
    integer pure function iseven(n) result(ans)
        integer, intent(in) :: n
        ans = 1 - isodd(n)
    end function iseven
    ! judge if two number are same odd or same even
    logical pure function is_same_parity(m, n) result(ans)
        integer, intent(in) :: m, n
        ans = mod(ieor(m, n), 2) == 0
    end function is_same_parity
    ! return (-1)^n
    integer pure function iphase(n) result(ans)
        integer, intent(in) :: n
        ans = iseven(n) - isodd(n)
    end function iphase
    ! check if m-quantum number if one of the components of a the j-quantum number
    logical pure function check_jm(dj, dm) result(ans)
        integer, intent(in) :: dj, dm
        ans = is_same_parity(dj, dm) .and. abs(dm) <= dj
    end function check_jm
    ! judge if three angular momentum can couple
    logical pure function check_couple(dj1, dj2, dj3) result(ans)
        integer, intent(in) :: dj1, dj2, dj3
        ans = dj1 >= 0 .and. dj2 >= 0 &
              .and. is_same_parity(dj1 + dj2, dj3) &
              .and. dj3 <= (dj1 + dj2) .and. dj3 >= abs(dj1 - dj2)
    end function check_couple

    ! stroe binomial to nmax = n, the total size
    integer pure function binomial_data_size(n) result(ans)
        integer, intent(in) :: n
        integer :: x
        x = n/2 + 1
        ans = x*(x + isodd(n))
    end function binomial_data_size
    ! get the index of C(n, k) in m_binomial_data
    integer pure function binomial_index(n, k) result(ans)
        integer, intent(in) :: n, k
        integer :: x
        x = n/2 + 1
        ans = x*(x - iseven(n)) + k
    end function binomial_index

    ! checked binomial
    real(kind=8) function binomial(n, k) result(ans)
        integer, intent(in) :: n, k
        if (n < 0 .or. n > m_nmax .or. k < 0 .or. k > n) then
            ans = 0.d0
        else
            ans = m_binomial_data(binomial_index(n, min(k, n - k)))
        end if
    end function binomial

    ! do not check bound, it is always safe in Wigner Symbols calculation
    real(kind=8) function unsafe_binomial(n, k) result(ans)
        integer, intent(in) :: n, k
        ans = m_binomial_data(binomial_index(n, min(k, n - k)))
    end function unsafe_binomial

    real(kind=8) function CG(dj1, dj2, dj3, dm1, dm2, dm3) result(ans)
        integer, intent(in) :: dj1, dj2, dj3, dm1, dm2, dm3
        integer :: J, jm1, jm2, jm3, j1mm1, j2mm2, j3mm3, j2pm2
        integer :: low, high, z
        real(kind=8) :: A, B
        if (.not. (check_jm(dj1, dm1) .and. check_jm(dj2, dm2) .and. check_jm(dj3, dm3))) then
            ans = 0.d0
            return
        end if
        if (.not. check_couple(dj1, dj2, dj3)) then
            ans = 0.d0
            return
        end if
        if (dm1 + dm2 /= dm3) then
            ans = 0.d0
            return
        end if
        J = (dj1 + dj2 + dj3)/2
        jm1 = J - dj1
        jm2 = J - dj2
        jm3 = J - dj3
        j1mm1 = (dj1 - dm1)/2
        j2mm2 = (dj2 - dm2)/2
        j3mm3 = (dj3 - dm3)/2
        j2pm2 = (dj2 + dm2)/2
        A = sqrt(unsafe_binomial(dj1, jm2)*unsafe_binomial(dj2, jm3)/ &
                 (unsafe_binomial(J + 1, jm3)*unsafe_binomial(dj1, j1mm1)* &
                  unsafe_binomial(dj2, j2mm2)*unsafe_binomial(dj3, j3mm3)))
        B = 0.d0
        low = max(0, j1mm1 - jm2, j2pm2 - jm1)
        high = min(jm3, j1mm1, j2pm2)
        do z = low, high
            B = -B + unsafe_binomial(jm3, z)* &
                unsafe_binomial(jm2, j1mm1 - z)*unsafe_binomial(jm1, j2pm2 - z)
        end do
        ans = iphase(high)*A*B
    end function CG

    real(kind=8) function wigner3j(dj1, dj2, dj3, dm1, dm2, dm3) result(ans)
        integer, intent(in) :: dj1, dj2, dj3, dm1, dm2, dm3
        ans = iphase((dj1 + (dj3 + dm3)/2))*sqrt(1.d0/(dj3 + 1))* &
              CG(dj1, dj2, dj3, -dm1, -dm2, dm3)
    end function wigner3j

    real(kind=8) function wigner6j(dj1, dj2, dj3, dj4, dj5, dj6) result(ans)
        integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6
        integer :: j123, j156, j426, j453
        integer :: jpm123, jpm132, jpm231, jpm156, jpm426, jpm453
        integer :: low, high, x
        real(kind=8) :: A, B
        if (.not. (check_couple(dj1, dj2, dj3) .and. check_couple(dj1, dj5, dj6) .and. &
                   check_couple(dj4, dj2, dj6) .and. check_couple(dj4, dj5, dj3))) then
            ans = 0.d0
            return
        end if
        j123 = (dj1 + dj2 + dj3)/2
        j156 = (dj1 + dj5 + dj6)/2
        j426 = (dj4 + dj2 + dj6)/2
        j453 = (dj4 + dj5 + dj3)/2
        jpm123 = (dj1 + dj2 - dj3)/2
        jpm132 = (dj1 + dj3 - dj2)/2
        jpm231 = (dj2 + dj3 - dj1)/2
        jpm156 = (dj1 + dj5 - dj6)/2
        jpm426 = (dj4 + dj2 - dj6)/2
        jpm453 = (dj4 + dj5 - dj3)/2
        A = sqrt(unsafe_binomial(j123 + 1, dj1 + 1)*unsafe_binomial(dj1, jpm123)/ &
                 (unsafe_binomial(j156 + 1, dj1 + 1)*unsafe_binomial(dj1, jpm156)* &
                  unsafe_binomial(j453 + 1, dj4 + 1)*unsafe_binomial(dj4, jpm453)* &
                  unsafe_binomial(j426 + 1, dj4 + 1)*unsafe_binomial(dj4, jpm426)))
        B = 0.0d0
        low = max(j123, j156, j426, j453)
        high = min(jpm123 + j453, jpm132 + j426, jpm231 + j156)
        do x = low, high
            B = -B + unsafe_binomial(x + 1, j123 + 1)*unsafe_binomial(jpm123, x - j453)* &
                unsafe_binomial(jpm132, x - j426)*unsafe_binomial(jpm231, x - j156)
        end do
        ans = iphase(high)*A*B/(dj4 + 1)
    end function wigner6j

    real(kind=8) function Racah(dj1, dj2, dj3, dj4, dj5, dj6) result(ans)
        integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6
        ans = iphase((dj1 + dj2 + dj3 + dj4)/2)*wigner6j(dj1, dj2, dj5, dj4, dj3, dj6)
    end function Racah

    real(kind=8) function wigner9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) result(ans)
        integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9
        integer :: j123, j456, j789, j147, j258, j369
        integer :: pm123, pm132, pm231, pm456, pm465, pm564, pm789, pm798, pm897
        real(kind=8) :: P0_nu, P0_de, P0
        real(kind=8) :: PABC, At, Bt, Ct, Pt_de
        integer :: j19t, j26t, j48t
        integer :: dt, dtl, dth, x, xl, xh, y, yl, yh, z, zl, zh
        if (.not. (check_couple(dj1, dj2, dj3) .and. check_couple(dj4, dj5, dj6) .and. check_couple(dj7, dj8, dj9) &
                   .and. check_couple(dj1, dj4, dj7) .and. check_couple(dj2, dj5, dj8) .and. check_couple(dj3, dj6, dj9))) then
            ans = 0.d0
            return
        end if
        j123 = (dj1 + dj2 + dj3)/2
        j456 = (dj4 + dj5 + dj6)/2
        j789 = (dj7 + dj8 + dj9)/2
        j147 = (dj1 + dj4 + dj7)/2
        j258 = (dj2 + dj5 + dj8)/2
        j369 = (dj3 + dj6 + dj9)/2
        pm123 = (dj1 + dj2 - dj3)/2
        pm132 = (dj1 + dj3 - dj2)/2
        pm231 = (dj2 + dj3 - dj1)/2
        pm456 = (dj4 + dj5 - dj6)/2
        pm465 = (dj4 + dj6 - dj5)/2
        pm564 = (dj5 + dj6 - dj4)/2
        pm789 = (dj7 + dj8 - dj9)/2
        pm798 = (dj7 + dj9 - dj8)/2
        pm897 = (dj8 + dj9 - dj7)/2
        P0_nu = unsafe_binomial(j123 + 1, dj1 + 1)*unsafe_binomial(dj1, pm123)* &
                unsafe_binomial(j456 + 1, dj5 + 1)*unsafe_binomial(dj5, pm456)* &
                unsafe_binomial(j789 + 1, dj9 + 1)*unsafe_binomial(dj9, pm798)
        P0_de = unsafe_binomial(j147 + 1, dj1 + 1)*unsafe_binomial(dj1, (dj1 + dj4 - dj7)/2)* &
                unsafe_binomial(j258 + 1, dj5 + 1)*unsafe_binomial(dj5, (dj2 + dj5 - dj8)/2)* &
                unsafe_binomial(j369 + 1, dj9 + 1)*unsafe_binomial(dj9, (dj3 + dj9 - dj6)/2)
        P0 = sqrt(P0_nu/P0_de)
        dtl = max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
        dth = min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
        PABC = 0.d0
        do dt = dtl, dth, 2
            j19t = (dj1 + dj9 + dt)/2
            j26t = (dj2 + dj6 + dt)/2
            j48t = (dj4 + dj8 + dt)/2
            Pt_de = unsafe_binomial(j19t + 1, dt + 1)*unsafe_binomial(dt, (dj1 + dt - dj9)/2)* &
                    unsafe_binomial(j26t + 1, dt + 1)*unsafe_binomial(dt, (dj2 + dt - dj6)/2)* &
                    unsafe_binomial(j48t + 1, dt + 1)*unsafe_binomial(dt, (dj4 + dt - dj8)/2)
            Pt_de = Pt_de*(dt + 1)*(dt + 1)
            xl = max(j123, j369, j26t, j19t)
            xh = min(pm123 + j369, pm132 + j26t, pm231 + j19t)
            At = 0.d0
            do x = xl, xh
                At = -At + unsafe_binomial(x + 1, j123 + 1)*unsafe_binomial(pm123, x - j369)* &
                     unsafe_binomial(pm132, x - j26t)*unsafe_binomial(pm231, x - j19t)
            end do
            yl = max(j456, j26t, j258, j48t)
            yh = min(pm456 + j26t, pm465 + j258, pm564 + j48t)
            Bt = 0.d0
            do y = yl, yh
                Bt = -Bt + unsafe_binomial(y + 1, j456 + 1)*unsafe_binomial(pm456, y - j26t)* &
                     unsafe_binomial(pm465, y - j258)*unsafe_binomial(pm564, y - j48t)
            end do
            zl = max(j789, j19t, j48t, j147)
            zh = min(pm789 + j19t, pm798 + j48t, pm897 + j147)
            Ct = 0.d0
            do z = zl, zh
                Ct = -Ct + unsafe_binomial(z + 1, j789 + 1)*unsafe_binomial(pm789, z - j19t)* &
                     unsafe_binomial(pm798, z - j48t)*unsafe_binomial(pm897, z - j147)
            end do
            PABC = PABC + iphase(xh + yh + zh)*At*Bt*Ct/Pt_de
        end do
        ans = iphase(dth)*P0*PABC
    end function wigner9j

    subroutine wigner_init(num, type, rank)
        integer, intent(in) :: num
        character(len=*), intent(in) :: type
        integer, intent(in) :: rank
        integer :: n, k
        ! int64 can store to binomial(66, 33)
        ! use int64 to improve the precision of small binomials
        if (.not. allocated(i64_binominal_data)) then
            allocate (i64_binominal_data(0:1155))
            i64_binominal_data(:) = 1_8
            do n = 1, 66
                do k = 1, n/2
                    i64_binominal_data(binomial_index(n, k)) = &
                        i64_binominal_data(binomial_index(n - 1, min(k, n - 1 - k))) + &
                        i64_binominal_data(binomial_index(n - 1, k - 1))
                end do
            end do
        end if
        select case (type)
        case ("Jmax")
            select case (rank)
            case (3)
                call fill_binomial_data(3*num + 1)
            case (6)
                call fill_binomial_data(4*num + 1)
            case (9)
                call fill_binomial_data(5*num + 1)
            case default
                error stop "rank for Jmax must be 3, 6, 9"
            end select
        case ("2bjmax")
            select case (rank)
            case (3)
                call fill_binomial_data(2*num + 1)
            case (6)
                call fill_binomial_data(3*num + 1)
            case (9)
                call fill_binomial_data(4*num + 1)
            case default
                error stop "rank for 2bjmax must be 3, 6, 9"
            end select
        case ("nmax")
            call fill_binomial_data(num)
        case default
            error stop "type must Jmax, 2bjmax, or nmax"
        end select
    end subroutine wigner_init

    subroutine fill_binomial_data(nmax)
        integer, intent(in) :: nmax
        integer :: i, n, k, reserve_size
        real(kind=8), dimension(:), allocatable :: old_data
        if (nmax <= m_nmax) return
        if (nmax > 92679) error stop "nmax is too large" ! binomial_data_size > huge(integer)
        if (nmax <= 66) then ! just the sotred i64 binomials
            reserve_size = binomial_data_size(nmax)
            if (allocated(m_binomial_data)) deallocate (m_binomial_data)
            allocate (m_binomial_data(0:reserve_size - 1))
            do i = 0, reserve_size - 1
                m_binomial_data(i) = dble(i64_binominal_data(i))
            end do
        else
            reserve_size = binomial_data_size(nmax)
            ! if it's not the first time calling this function,
            ! and the m_binomial_data contains some data, use the old data
            if (allocated(m_binomial_data)) then
                allocate (old_data(0:ubound(m_binomial_data, dim=1)))
                old_data(:) = m_binomial_data(:)
                deallocate (m_binomial_data)
                allocate (m_binomial_data(0:reserve_size - 1))
                m_binomial_data(0:ubound(old_data, dim=1)) = old_data(:)
                do n = m_nmax + 1, nmax
                    do k = 0, n/2
                        m_binomial_data(binomial_index(n, k)) = binomial(n - 1, k) + binomial(n - 1, k - 1)
                    end do
                    m_nmax = m_nmax + 1
                end do
                deallocate (old_data)
            else
                ! if it's the first time calling this function
                ! use the i64_binomial_data for the first 66 lines
                ! and calculate the rest
                allocate (m_binomial_data(0:reserve_size - 1))
                do i = 0, 1155
                    m_binomial_data(i) = dble(i64_binominal_data(i))
                end do
                m_nmax = 66
                do n = 67, nmax
                    do k = 0, n/2
                        m_binomial_data(binomial_index(n, k)) = binomial(n - 1, k) + binomial(n - 1, k - 1)
                    end do
                    m_nmax = m_nmax + 1
                end do
            end if
        end if
        m_nmax = nmax
    end subroutine fill_binomial_data
end module WignerSymbol
