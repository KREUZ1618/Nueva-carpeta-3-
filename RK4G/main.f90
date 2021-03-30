program RK4NEQ
    implicit none
    integer::EQ
    integer::VAR
    !número de ecuaciones
    parameter (EQ=5)
    !número de variables dependientes del tiempo 
    parameter (VAR=5)
    !los coeficientes extras son del tiempo y una constante
    !vector de valores xn
    real*8,dimension(VAR+2)::x
    real*8,dimension(EQ)::matrix
    !matriz de coeficientes
    real*8,dimension(VAR+2,EQ)::coef
    !h tamaño de paso,t tiempo inicial,tf tiempo final
    real*8::h,t,tf
    !N número de iteraciones
    integer::i,N

    !lectura de archivos
    OPEN(1,FILE="valoresinciales.dat",STATUS="OLD")
    OPEN(2,FILE="coeficientes.dat",STATUS="OLD")

    !vector de valores iniciales
    READ(1,*)x

    !matriz transpuesta,Fortran lee la matriz transpuesta de coeficientes
    READ(2,*)coef
    close(1)
    close(2)

    !prueba para ver como se realizan las operaciones
    matrix=matmul(x,coef)
    write(*,*) matrix


    write(*,*) "escriba el tamaño de paso h"
    read(*,*)h

    write(*,*) "escriba el intervalo: tinicial,tfinal"
    read(*,*)t,tf


    open(3,FILE="grafica1.dat",STATUS="old")

    !número de iteraciones
    N=(tf-t)/h
    !RK4 para N ecuaciones
    do i=1,N
    !subrutina RK4
    call RK4(h,x,coef,t,EQ)

    !incremento en el tiempo
    t=t+h

    !se guardan los valores almacenados en el vector x
    write(3,*)x
    end do
    close(3)



end program RK4NEQ

subroutine RK4(h,x,coef,t)

    integer::EQ
    integer::VAR
    parameter (EQ=5)
    !número de variables dependientes del tiempo
    parameter (VAR=5)
    !x vector de valores de xn,
    !xx vector de valores evaluados en la multiplicación matricial
    real*8,dimension(VAR+2)::x,xx
    real*8,dimension(EQ)::k1,k2,k3,k4
    real*8,dimension(VAR+2,EQ)::coef
    real*8::h,t
    integer::i

    !PARTE DEL VECTOR EVALUADO QUE SE MULTIPLICA POR LA PARTE CONSTANTE DE LA MATRIZ DE COEFCIENTES 
    xx(VAR+2)=1

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

