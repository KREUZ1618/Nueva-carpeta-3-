
program ADM1NEQ
    implicit none
    integer::EQ
    integer::VAR
    !número de ecuaciones
    parameter (EQ=2)
    !número de variables dependientes del tiempo
    parameter (VAR=2)
    !la variable extra es el tiempo VAR+1
    !vector de valores xn
    real*8,dimension(VAR+1)::x,x0
    real*8,dimension(EQ)::matrix
    !matriz de coeficientes
    real*8,dimension(VAR+1,EQ)::coef
    !h tamaño de paso,t tiempo inicial,tf tiempo final
    real*8::h,t,tf,N,h0
    real::y2,x2
    !N número de iteraciones
    integer::i,j

    !lectura de archivos
    OPEN(1,FILE="valoresinciales.dat",STATUS="OLD")
    OPEN(2,FILE="coeficientes.dat",STATUS="OLD")

    !vector de valores iniciales
    READ(1,*) x

    !matriz transpuesta,Fortran lee la matriz transpuesta de coeficientes
    READ(2,*)coef
    close(1)
    close(2)

    !prueba para ver como se realizan las operaciones
    matrix=matmul(x,coef)
    write(*,*) x,coef,matrix


    write(*,*) "escriba el tamaño de paso h"
    read(*,*)h

    write(*,*) "escriba el intervalo: tinicial,tfinal"
    read(*,*)t,tf


    open(3,FILE="tabla4.dat",STATUS="old")

    !número de iteraciones
    N=(tf-t)/h
    write(*,*)N
    !RK4 para N ecuaciones

    do i=1,N
    !subrutina RK4



    x0=x
    call RK4(h,x,coef,t,EQ)
    call ADM1(h,x,x0,coef,t,EQ)
    !incremento en el tiempo
    t=t+h

    !se guardan los valores almacenados en el vector x
    write(3,*)x

    end do
    close(3)



end program ADM1NEQ
subroutine EULER(h,x,x0,coef,t)

    integer::EQ
    integer::VAR
    parameter (EQ=2)
    !número de variables dependientes del tiempo
    parameter (VAR=2)
    !x vector de valores de xn,
    !xx vector de valores evaluados en la multiplicación matricial
    real*8,dimension(VAR+1)::x,xx,x0
    real*8,dimension(EQ)::k1,k2,k3,k4
    real*8,dimension(VAR+1,EQ)::coef
    real*8::h,t
    integer::i


    !X
    do i=1,EQ
    k1=matmul(x,coef)
    x(i)=x(i)+h*k1(i)
    end do


end subroutine
subroutine ADM1(h,x,x0,coef,t)

    integer::EQ
    integer::VAR
    parameter (EQ=2)
    !número de variables dependientes del tiempo
    parameter (VAR=2)
    !x vector de valores de xn,
    !xx vector de valores evaluados en la multiplicación matricial
    real*8,dimension(VAR+1)::x,xx,x0
    real*8,dimension(EQ)::k1,k2,k3,k4
    real*8,dimension(VAR+1,EQ)::coef
    real*8::h,t
    integer::i


    !X
    do i=1,EQ
    k1=matmul(x,coef)
    x(i)=x0(i)+h*k1(i)
    end do


end subroutine
subroutine RK4(h,x,coef,t)

    integer::EQ
    integer::VAR
    parameter (EQ=2)
    !número de variables dependientes del tiempo
    parameter (VAR=2)
    !x vector de valores de xn,
    !xx vector de valores evaluados en la multiplicación matricial(ec)
    real*8,dimension(VAR+1)::x,xx
    real*8,dimension(EQ)::k1,k2,k3,k4
    real*8,dimension(VAR+1,EQ)::coef
    real*8::h,t
    integer::i


    !K1
    k1=matmul(x,coef)
    write(*,*)"k1", k1

    !K2
    do i=1,EQ+1
    if(i.lt.EQ+1) then
    xx(i)=x(i)+k1(i)*h/2
    else
    xx(i)=t+h*0.5
    x(i)=xx(i)
    end if
    end do

    k2=matmul(xx,coef)
    write(*,*)"k2", k2

   !K3
    do i=1,EQ+1
    if(i.lt.EQ+1) then
    xx(i)=x(i)+k2(i)*h/2
    else
    xx(i)=t+h*0.5
    x(i)=xx(i)
    end if
    end do

    k3=matmul(xx,coef)
    write(*,*)"k3", k3

   !K4
    do i=1,EQ+1
    if(i.lt.EQ+1) then
    xx(i)=x(i)+k3(i)*h
    else
    xx(i)=t+h
    x(i)=xx(i)
    end if
    end do

    k4=matmul(xx,coef)
    write(*,*)"k4", k4


    !X
    do i=1,EQ
    x(i)=x(i)+(k1(i)+2*k2(i)+2*k3(i)+k4(i))*h/6
    end do


end subroutine
!solucion analítica

real function x2(t)
implicit none
real::t
x2=14*exp(-2*t)-8*exp(-800*t)
end

real function y2(t)
implicit none
real::t
y2=2*exp(-2*t)-8*exp(-800*t)
end
