# WignerSymbol-f

[English](README.md)

计算CG系数，Racah系数，Wigner 3j, 6j, 9j 系数。公式参考自[CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl)。

## 使用方法

最简单的使用方法是复制`src/WignerSymbol.f90`到你的项目里。

```fortran
program test_wigner
    use WignerSymbol
    implicit none
    call wigner_init(20, "Jmax", 9)
    print*, wigner6j(6,6,6,6,6,6)
    print*, wigner9j(2,4,6,8,10,12,6,12,18)
end program test_wigner
```

你也可以用fpm来安装
```toml
[dependencies]
WignerSymbol-f = { git="https://github.com/0382/WignerSymbol-f.git" }
```

## 提供的函数
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
其中，除了`binomial`函数之外，其余函数均使用真实角动量量子数的两倍作为参数，这是为了处理半整数角动量的情况。所以要计算`<10|1/2,1/2;1/2,-1/2>`这个CG系数，你需要调用的是
```fortran
use WignerSymbol
real(kind=8) :: x
x = CG(1,1,2,1,-1,0);
```

## `wigner_init`子程序

所有的系数都是使用二项式系数`binomial`来计算的，这个库会先存储下部分二项式系数，之后用到的时候直接取。所以这个库调用的`binomial`函数只有参数的范围在存储范围之内才获得正确的结果，否则直接返回零。

这个`wigner_init`子程序就是用来初始化二项式系数表的
```fortran
subroutine wigner_init(num, type, rank)
    integer, intent(in) :: num
    character(len=*), intent(in) :: type
    integer, intent(in) :: rank
end subroutine wigner_init
```
其含义如下

|                       |    计算范围    |   CG & 3j   | 6j & Racah  |     9j      |
| :-------------------: | :------------: | :---------: | :---------: | :---------: |
|       `type`的意义    | `type`\\`rank` |      3      |      6      |      9      |
|      最大角动量模式     |    `"Jmax"`    | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
|    最大两体角动量模式    |   `"2bjmax"`    | `2*jmax+1` | `3*jmax+1` | `4*jmax+1` |
|    最大二项式系数模式    |    `"nmax"`    |   `nmax`    |   `namx`    |   `nmax`    |

表格中的数据是最极端条件下保证不溢出最少要存多少二项式系数。不过使用的时候你不需要记住这些数值，只需要根据你的程序中出现的最大角动量和计算范围调用一下`wigner_init`子程序就好了。

### `"2bjmax"`

例如
```fortran
call wigner_init(21, "2bjmax", 6);
```

表示的意义是，我们体系中最大的单粒子轨道的角动量为`21/2`，于是最大可能的两体耦合角动量为`21`，同时代码中仅计算CG系数和6j系数，而不会计算9j系数。在这个例子中，保存的`binomial(n, k)`系数中，最大的`n`为
```fortran
nmax = 3*jmax + 1 = 3 * 21 + 1 = 64
```

`"2bjmax"`意味着你的程序仅需要计算两体耦合而不需要计算三体耦合或更高的耦合。这个模式假定在所有的Wigner系数计算中，至少有一个角动量是单粒子角动量，也即在这个例子中不超过`21/2`。在这个假定下，`"2bjmax"`模式能够少存一些二项式系数，占用内存更少。

### `"Jmax"`

`"Jmax"`最大角动量模式则不区分角动量的来源，表示全局的最大角动量。比如如果上述例子改成
```fortran
call wigner_init(Jmax, "Jmax", 6);
```

表示我们体系中最大可能的角动量为`Jmax = 21`，只计算到6j系数。那么对应的
```fortran
nmax = 4*Jmax + 1 = 85
```

对于量子多体计算来说，如果单粒子轨道最大角动量为`21/2`，那么两体耦合的最大角动量为`21`。采用`"2bjmax"`模式能够更加节省内存。

`"Jmax"`模式没有任何假定，它总是安全的。即使考虑到三体耦合，只要使用三体的总`Jmax`来`reserve`，就能够保证不会溢出，尽管这可能会浪费一些空间。

不过实际上存下二项式系数占用的内存并不多，简单的估算其占用内存约为
```fortran
2 * nmax * nmax ! Byte
```
即使`nmax = 1000`（对应计算到9j系数，`Jmax = 200`，对大多数计算绰绰有余），需要的内存也只有2MB，相对于科学计算的其他内存占用来说实在不算什么。

### `"nmax"`

`"nmax"`模式就是直接设置`nmax`，此时`rank`参数不起作用。可能只有你不想算Wigner系数而只想用这个库来计算二项式系数时会有用。

### 线程安全

你可以多次调用`wigner_init`子程序，后面的调用如果需要更大的`nmax`，那么它会拓展存储的二项式系数表。

不过请注意：`wigner_init`子程序**不是**线程安全的，如果你的程序是并行的，不要动态地调用`wigner_init`子程序。本库的推荐使用方法是，先计算出体系最大角动量，然后在程序开始时调用一次`wigner_init`子程序，之后就不应该继续调用它了。
