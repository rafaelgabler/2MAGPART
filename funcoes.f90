module funcoes

use variaveis

contains

!********************************** Subrotinas Utilizadas no codigo ************************!
subroutine resvel(a,b,c,d,e,f)
real a                      ! componente de velocidade em questao
real b			    ! passo de tempo
real c                      ! numero de stokes
real d                      ! somatorio das demais forcas na direcao considerada
real e                      ! vale 0 ou 1, se for zero nao considera empuxo liquido, se for 1 considera 
real f			    ! termo que multiplica U na edo
real k1,k2,k3,k4            ! variaveis internas utilizadas para o runge-kutta de quarta ordem

k1=b*(-a*f-e+d)/c
k2=b*(f*(-a-0.5*k1)-e+d)/c
k3=b*(f*(-a-0.5*k2)-e+d)/c
k4=b*(f*(-a-k3)-e+d)/c

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine resvel

!*******************************************************************************************!

subroutine resvel_sem_inercia(a,b,c)
real a                      ! componente de velocidade em questao
real b			    ! somatorio das demais forcas na direcao considerada
real c                      ! vale 0 ou 1, se for zero nao considera empuxo liquido, se for 1 considera

a=-c+b 
end subroutine resvel_sem_inercia

!*******************************************************************************************!
subroutine resomega(a,b,c,d,e)
real a                      ! componente de velocidade em questao
real b			    ! passo de tempo
real c                      ! numero de stokes rotacional
real d                      ! somatorio dos demais torques na direcao considerada
real e                      ! termo que multiplica omega na edo
real k1,k2,k3,k4            ! variaveis internas utilizadas para o runge-kutta de quarta ordem

k1=b*(-a*e+d)/c
k2=b*(e*(-a-0.5*k1)+d)/c
k3=b*(e*(-a-0.5*k2)+d)/c
k4=b*(e*(-a-k3)+d)/c

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine resomega

!*******************************************************************************************!

subroutine resomega_sem_inercia(a,b)
real a                      ! componente de velocidade angular em questao
real b			    ! somatorio dos demais torques na direcao considerada

a=b 
end subroutine resomega_sem_inercia

!*******************************************************************************************!

subroutine evoldip(a,b,c,d,e,f)
real a			    ! Momento de dipolo atual da particula na direcao considerada
real b			    ! Momento de dipolo atual da particula na direcao j
real c			    ! Momento de dipolo atual da particula na direcao k
real d                      ! Velocidade angular da particula na direcao j
real e                      ! Velocidade angular da particula na direcao k
real f                      ! Passo de tempo numerico

a=a+(d*c -e*b)*f


end subroutine evoldip  

!*******************************************************************************************!

subroutine respos(a,b,c)
real a                      ! posicao
real b			    ! passo de tempo
real c                      ! componente de velocidade em questao
real k1,k2,k3,k4            ! variaveis internas utilizadas para o runge-kutta de quarta ordem

k1=b*(c)
k2=b*((0.5*k1)+c)
k3=b*((0.5*k2)+c)
k4=b*((k3)+c)

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine respos

!*******************************************************************************************!

subroutine randomica(a,b,c,n,d)
real a,b                    ! a,b = range do numero randomico  
integer n, m
real c(n)                   ! c = sequencia randomica gerada
integer d,i,e
integer f(8)
integer, allocatable :: seed(:)

 call random_seed(size = m)
allocate (seed(m))
 
 CALL DATE_AND_TIME(values=f)
 CALL SYSTEM_CLOCK(count=e)

do i = 1,m
seed(i) =  37*d + f(8)*i*d + e*(3*d+i)
end do

 call random_seed(put = seed)

 call random_number(c)

 c = a+(b-a)*c

end subroutine randomica

!*******************************************************************************************!
subroutine media_direcao(X,N,rea,media,direcao)
integer N,rea
integer i,j,k
integer direcao
real X(rea,N,3)
real media


media=sum(X(:,:,direcao))/(N*rea)


end subroutine media_direcao
!*******************************************************************************************!


end module funcoes
