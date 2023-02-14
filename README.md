# WignerSymbol-f

[中文](README_zh.md)

Calculate CG coefficient, Racah coefficient, and Wigner 3j, 6j, 9j coefficient. Calculation formula please see [CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl).

## Usage

The simple way to use this package is copy the file `src/WignerSymbol.f90` into your project.

```fortran
program test_wigner
    use WignerSymbol
    implicit none
    call wigner_init(20, "Jmax", 9)
    print*, wigner6j(6,6,6,6,6,6)
    print*, wigner9j(2,4,6,8,10,12,6,12,18)
end program test_wigner
```

Or you can use fpm
```toml
[dependencies]
WignerSymbol-f = { git="https://github.com/0382/WignerSymbol-f.git" }
```

## API
```fortran
! binomial
real(kind=8) pure function binomial(n, k)
    integer, intent(in) :: n, k
end function binomial
! CG coefficient
real(kind=8) pure function CG(dj1, dj2, dj3, dm1, dm2, dm3)
    integer, intent(in) :: dj1, dj2, dj3, dm1, dm2, dm3
end function CG
! Wigner 3j symbol
real(kind=8) pure function wigner3j(dj1, dj2, dj3, dm1, dm2, dm3)
    integer, intent(in) :: dj1, dj2, dj3, dm1, dm2, dm3
end function wigner3j
! Wigner 6j symbol
real(kind=8) pure function wigner6j(dj1, dj2, dj3, dj4, dj5, dj6)
    integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6
end function wigner6j
! Racah coefficient
real(kind=8) pure function Racah(dj1, dj2, dj3, dj4, dj5, dj6)
    integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6
end function Racah
! Wigner 9j symbol
real(kind=8) pure function wigner9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    integer, intent(in) :: dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9
end function wigner9j
! wigner d-function
real(kind=8) pure function dfunc(dj, dm1, dm2, beta) result(ans)
    integer, intent(in) :: dj, dm1, dm2
    real(kind=8), intent(in) :: beta
end function dfunc
! Buck et al. Nuc. Phys. A 600 (1996) 387-402
real(kind=8) pure function Moshinsky(N, L, nr, lr, n1, l1, n2, l2, lambda, tan_beta) result(ans)
    integer, intent(in) :: N, L, nr, lr, n1, l1, n2, l2, lambda
    real(kind=8), intent(in), optional :: tan_beta
end function Moshinsky
```
Apart from `binomial`, all the functions use double of the real angular momentum quantum number to avoid half integers. So if you want to calculate `<10|1/2,1/2;1/2,-1/2>`, you should call like this,
```fortran
use WignerSymbol
real(kind=8) :: x
x = CG(1,1,2,1,-1,0);
```

## The `wigner_init` subroutine

We calculate the Wigner Symbols with `binomial`s, and we will store some binomials first, then when we need one binomial, we just load it. In this package, the `binomial` function is only valid in the stored range. If you call a `binomial` function out of the range, it just gives you `0`.

The `wigner_init` subroutine is used to init the `binomial`s' table.
```fortran
subroutine wigner_init(num, type, rank)
    integer, intent(in) :: num
    character(len=*), intent(in) :: type
    integer, intent(in) :: rank
end subroutine wigner_init
```
and its parameters means

|                                       |    Calculate range    |   CG & 3j   | 6j & Racah  |     9j      |
| :-----------------------------------: | :-------------------: | :---------: | :---------: | :---------: |
|          meaning of `type`            | `type`\\\\`rank`      |      3      |      6      |      9      |
|        max angular momentum           |    `"Jmax"`           | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
| max two-body coupled angular momentum |   `"2bjmax"`          | `2*jmax+1`  | `3*jmax+1`  | `4*jmax+1`  |
|            max binomial               |    `"nmax"`           |   `nmax`    |   `namx`    |   `nmax`    |

The value in the table means the minimum binomial range to guarantee the Wigner Symbol calculation. You do not need to rememmber those values. You just need to find the maximum angular momentum in you canculation, then call the `wigner_init` subroutine before you calculate the Wigner Symbols.

### `"2bjmax"`

For example
```fortran
call wigner_init(21, "2bjmax", 6);
```

This means the maximum single particle angular momentum is `21/2`, and thus the maximum two-body coupled angular momentum is `21`, and the `rank = 6` means you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbol. In this example, the `n` of the stored maximum binomial coefficient is
```fortran
nmax = 3*jmax + 1 = 3 * 21 + 1 = 64
```

The `"2bjmax"` mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is just a single particle angular momentum, thus in this example no larger than `21/2`. With this assumption, `"2bjmax"` mode will use less memory than `"Jmax"` mode.

### `"Jmax"`

`"Jmax"` means the global maximum angular momentum, for every parameters. For another example,
```fortran
call wigner_init(Jmax, "Jmax", 6);
```

This means in this system, `Jmax = 21`, and we calculate to 6j symbols. Here
```fortran
nmax = 4*Jmax + 1 = 85
```

For quantum many body calculation with only two-body coupling, if single particle `jmax = 21/2`, then `Jmax = 21`. The `"2bjmax"` mode will cost less storage.

In the `"Jmax"` mode, it is always safe with out any assumption. Even having three-body coupling, you just need to use the maximum three body coupled angular momentum as `Jmax`, although it will cost more memory.

Actually, the memory used for store the binomials is not very large. A simple estimate is
```fortran
2 * nmax * nmax ! Byte
```
Even for `nmax = 1000` (`Jmax = 200` in 9j calculation, which is absolutly enough for most calculation), the memory cost is just 2MB, which is not very large.

### `"nmax"`

The `"nmax"` mode directly set `nmax`, and the `rank` parameter is ignored. This maybe useful when you only want to calculate `binomial`s using this package.

### Thread safety

You can call `wigner_init` several times in you program, if later call requires larger `nmax`, it will extent the stored binomials' table.

However, the `wigner_init` is **not** thread safe. So you shuld not call `wigner_init` dymanically in a multi-threading program. The recommended way to use this package is find the maximum angular momentum quantum number in you system, and call `wigner_init` at the beginning of the code, and then don't call it any more.
