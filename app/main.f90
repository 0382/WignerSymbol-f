program test_wigner
    use WignerSymbol
    implicit none
    call wigner_init(20, "Jmax", 6)
    print*, wigner6j(6,6,6,6,6,6)
    print*, wigner9j(2,4,6,8,10,12,6,12,18)
end program test_wigner