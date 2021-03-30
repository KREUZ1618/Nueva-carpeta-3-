program RK4
implicit none

integer::i,n
real::h,x,y,f,y1,x1,R,C,V,xp,dy,k1,k2,k3,k4,d,error1,error2,t,tf,y2,x2
real::c1,c2,c3,c4

write(*,*) "ti,tf,yo,yo'"
read(*,*) t,tf,y,x

write(*,*)'tamaño de h'
read(*,*)h


!se calcula el número de pasos
n=int((tf-t)/h)

!archivo tabla
open(1, file = 'tabla2.dat', status='old')
write(1,'(1x,T12,A,T15,A,T33,A,T50,A,T70,A,T90,A)') 't','x','y','x(t)','y(t)','y2-y1','x2-x1'

!archivo gráfica solucion numérica
open(2, file = 'grafica.dat', status='old')

!archivo gráfica solucion analítica
open(4, file = 'grafica2.dat', status='old')

!archivo gráfica error
open(3, file = 'graficaError.dat', status='old')


error1=abs(y2(t)-y)
error2=abs(x2(t)-x)
call d_grafica(2,t,x,x2(t))
call d_grafica(4,t,y,y2(t))
call d_tabla(1,t,x,y,x2(t),y2(t),error2,error1)
do i=1,n
 

 k1=y1(t,x,y)
 c1=x1(t,x,y)
 write(*,*) "k1",k1,c1


 k2=y1(t+0.5*h,x+0.5*c1*h,y+0.5*k1*h)
 c2=x1(t+0.5*h,x+0.5*c1*h,y+0.5*k1*h)
 write(*,*) "k2",k2,c2

 k3=y1(t+0.5*h,x+0.5*c2*h,y+0.5*k2*h)
 c3=x1(t+0.5*h,x+0.5*c2*h,y+0.5*k2*h)
 write(*,*) "k3",k3,c3

 k4=y1(t+h,x+c3*h,y+k3*h)
 c4=x1(t+h,x+c3*h,y+k3*h)
 write(*,*) "k4",k4,c4

 y=(y+(h*(k1+2*k2+2*k3+k4))/6)
 x=(x+(h*(c1+2*c2+2*c3+c4))/6)
 write(*,*) "x",y,x
 
 


t=t+h

error1=abs(y2(t)-y)
error2=abs(x2(t)-x)
call d_grafica(2,t,x,x2(t))
call d_grafica(4,t,y,y2(t))
call d_tabla(1,t,x,y,x2(t),y2(t),error2,error1)
write(*,*)"analítica",i,t,x2(t),y2(t),error2,error1
write(*,*)"númerica",i,t,x,y,error2,error1


end do

 close(4)
 close(3)
 close(2)
 close(1)
end program RK4




real function x1(t,x,y)
implicit none
real::x,y,t
x1=-3*x+4*y
end

real function y1(t,x,y)
implicit none
real::x,y,t
y1=-2*x+3*y

end

!solucion analítica

real function x2(t)
implicit none 
real::t
x2=3*exp(t)-2*exp(-t)
end

real function y2(t)
implicit none
real::t
y2=3*exp(t)-exp(-t)
end


subroutine d_grafica(archivo,x,y,y1)

 integer archivo
 real x,y,y1
 
 write(archivo,*) x,y,y1
end 


subroutine d_tabla(archivo,t,x,y,x2,y2,error1,error2)


 integer archivo
 real t,x,y,x2,y2,error1,error2
 write(archivo,*) t,x,y,x2,y2,error1,error2
 
 end 

