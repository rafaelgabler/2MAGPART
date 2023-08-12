subroutine principal

use variaveis
use funcoes

write(*,*) 'Processando...'

comentarios=FALSE

! Determinando o inverso do parametro de polidispersidade, importante para determinar 
! as funcoes resistencia da esfera 2, X22a e Y22a

alfa=1.0/beta


! Determinando a posicao inicial das particulas

do j=1,rea
X(j,1,1)=X_11
X(j,1,2)=X_12
X(j,1,3)=X_13
X(j,2,1)=X_21
X(j,2,2)=X_22
X(j,2,3)=X_23
end do

if(pos_inicial) then
do j=1,rea
X(j,1,1)=X_11
X(j,1,2)=X_12
X(j,1,3)=X_13
X(j,2,1)=X_21
X(j,2,2)=X_22
X(j,2,3)=X_23
end do

else

if(quad.eq.1) then
do j=1,(rea**0.5)
do i=1,(rea**0.5)
 X(((j-1)*(rea**0.5)+i),2,1)=  X_11 + i*(5.0/(rea**0.5))
 X(((j-1)*(rea**0.5)+i),2,2)=  X_12 + j*(5.0/(rea**0.5))
end do
end do
end if

if(quad.eq.2) then
do j=1,(rea**0.5)
do i=1,(rea**0.5)
 X(((j-1)*(rea**0.5)+i),2,1)=  X_11 +  i*(5.0/(rea**0.5))
 X(((j-1)*(rea**0.5)+i),2,2)=  X_12 -  j*(5.0/(rea**0.5))
end do
end do
end if

if(quad.eq.3) then
do j=1,(rea**0.5)
do i=1,(rea**0.5)
 X(((j-1)*(rea**0.5)+i),2,1)=  X_11 -  i*(5.0/(rea**0.5))
 X(((j-1)*(rea**0.5)+i),2,2)=  X_12 -  j*(5.0/(rea**0.5))
end do
end do
end if

if(quad.eq.4) then
do j=1,(rea**0.5)
do i=1,(rea**0.5)
 X(((j-1)*(rea**0.5)+i),2,1)=  X_11 -  i*(5.0/(rea**0.5))
 X(((j-1)*(rea**0.5)+i),2,2)=  X_12 +  j*(5.0/(rea**0.5))
end do
end do
end if


do j=1,rea
posicao_inicial(j,1)=X(j,2,1) 
posicao_inicial(j,2)=X(j,2,2)
end do

open(14*rea,file='diagrama_da_condicao_inicial.plt')
write(14*rea,*)'Variables="X","Y"'
do j=1,rea
write(14*rea,*)  X(j,2,1),X(j,2,2)
end do 


end if

! Determinando a distancia na qual se desliga a forca de interacao magnetica
! para evitar overlaps de particulas

!if(St.ge.0.005) then
!if(St.le.0.7) then
!dist = 1.6853*((Arm)**0.0975)
!end if
!end if

!if(St.ge.0.701) then
!if(St.le.5.0) then
!dist = 1.8401*((Arm)**0.1399)
!end if
!end if

!if(St.ge.5.01) then
!if(St.le.10.01) then
!dist = 2.3409*((Arm)**0.1628)
!end if
!end if

dist=0.0

! Determinando o vetor momento de dipolo das duas particulas ao longo de todas as realizacoes

if(dipolo)then

do j=1,rea
Di(j,1,1)=Di_11
Di(j,1,2)=Di_12
Di(j,1,3)=Di_13
Di(j,2,1)=Di_21
Di(j,2,2)=Di_22
Di(j,2,3)=Di_23
end do

! Normalizando o vetor momento de dipolo inicial das particulas

do j=1,rea
do i=1,2
modip=((Di(j,i,1)**2.0)+(Di(j,i,2)**2.0)+(Di(j,i,3)**2.0))**0.5
Di(j,i,1)=Di(j,i,1)/modip
Di(j,i,2)=Di(j,i,2)/modip
Di(j,i,3)=Di(j,i,3)/modip
modip=0.0
end do
end do

else

! Distribuindo os momentos de dipolo iniciais das particulas

 call randomica(-1.0,1.0,nr,(3*2*rea),2)

do j=1,rea
do i=1,N
Di(j,i,1)= nr((i*2+(i-2)+(N*3*(j-1))))
Di(j,i,2)= nr((i*2+(i-1)+(N*3*(j-1))))
Di(j,i,3)= nr((i*2+(i)+(N*3*(j-1))))
end do
end do

! Normalizando estes vetores

do j=1,rea
do i=1,N
modip=((Di(j,i,1)**2.0)+(Di(j,i,2)**2.0)+(Di(j,i,3)**2.0))**0.5
Di(j,i,1)=Di(j,i,1)/modip
Di(j,i,2)=Di(j,i,2)/modip
Di(j,i,3)=Di(j,i,3)/modip
modip=0.0
end do
end do

end if

if(comentarios) then
write(*,*) 'Momento de dipolo inicial:'
write(*,*) Di

pause

write(*,*) 'Modulo dos momentos de dipolo iniciais:'
do j=1,rea
do i=1,N
modip=((Di(j,i,1)**2.0)+(Di(j,i,2)**2.0)+(Di(j,i,3)**2.0))**0.5
write(*,*) modip
end do
end do

pause
end if

! Zerando o vetor tempo

T(1)=0.0


! Criando os arquivos correspondentes a cada realizacao

if(arquivo)then

do i=1,rea
write(rea_char, '(I3)') i
open (i,file='posicao'//rea_char//'.plt')
write(i,*) 'Variables="X","Y","Z"'
end do

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt')
write(rea+i,*) 'Variables="U","V","W"'
end do


do i=1,rea
write(rea_char, '(I3)') i
open (rea+rea+rea+i,file='dipolo'//rea_char//'.plt')
write(rea+rea+rea+i,*) 'Variables="Di1","Di2","Di3"'
end do

open(5*rea,file='mag_tempo.plt')
write(5*rea,*)'Variables="M","t"'

do i=1,rea
write(rea_char, '(I3)') i
open (rea+rea+i,file='trajetoria'//rea_char//'.plt')
write(rea+rea+i,*) 'Variables="X1","X2","X3"'
end do

end if


open (10*rea,file='distancia.plt')
write(10*rea,*) 'Variables="D","t"'

open(100*rea,file='analise_rotacional.plt')
write(100*rea,*) 'Variables="t","Z","theta","W"'

open(204*rea,file='variacao_dipolo.plt')
write(204*rea,*) 'Variables="X","Y","Z","Dix","Diy","Diz"'

open(666*rea,file='variacao_dipolo.gnu')

734 FORMAT(12F12.5)

open(5*rea,file='mag_tempo.plt')

! Escrevendo nos arquivos criados as informacoes referentes a condicao inicial

509 FORMAT(F30.4,F30.4,F30.4)

k=1

write(204*rea,'(A12,I6,A1)') 'zone t="',k,'"' 
write(204*rea,*) X(1,1,1),X(1,1,2),X(1,1,3),Di(1,1,1),Di(1,1,2),Di(1,1,3)
write(204*rea,*) X(1,2,1),X(1,2,2),X(1,2,3),Di(1,2,1),Di(1,2,2),Di(1,2,3)

write(666*rea,734) X(1,1,1),X(1,1,2),X(1,1,3),Di(1,1,1),Di(1,1,2),Di(1,1,3),X(1,2,1),X(1,2,2),X(1,2,3),Di(1,2,1),Di(1,2,2),Di(1,2,3)

if(arquivo) then

do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"' 
do i=1,N
write(j,509)X(j,i,1),X(j,i,2),   &
X(j,i,3)
end do
end do

do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"' 
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do

end if


write(*,*) 'Iteracao:', k

modomega=0.0
theta=acos((Di(1,1,1)*Di(1,2,1))+(Di(1,1,2)*Di(1,2,2))+(Di(1,1,3)*Di(1,2,3)))

!************** A PARTIR DESTE MOMENTO A SIMULACAO COMECA **********************************! 

!Iniciando o processo iterativo

do k=2,npast

write(*,*) 'Iteracao:', k

! Zerando as forcas e torques

FORCAS=0.0
TORQUES=0.0
T(k)=T(k-1) + dt

!************* Determinando as forcas que agem sobre cada particula ************************!


!************************** FORCAS BROWNIANAS **********************************************!

if(browniano) then

 call randomica(-1.0,1.0,nr,(3*2*rea),2+k)
 !$OMP PARALLEL DO
do j=1,rea

i=1

nr1= nr((i*2+(i-2)+(2*3*(j-1))))
nr2= nr((i*2+(i-1)+(2*3*(j-1))))
nr3= nr((i*2+(i)+(2*3*(j-1))))

modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5

nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand

FORCAS(1,j,i,1)= ((6.0/(Pe*dt))**0.5)*nr1
FORCAS(1,j,i,2)= ((6.0/(Pe*dt))**0.5)*nr2
FORCAS(1,j,i,3)= ((6.0/(Pe*dt))**0.5)*nr3

i=2

nr1= nr((i*2+(i-2)+(2*3*(j-1))))
nr2= nr((i*2+(i-1)+(2*3*(j-1))))
nr3= nr((i*2+(i)+(2*3*(j-1))))

modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5

nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand
 
FORCAS(1,j,i,1)= alfa**(1.0/2.0)*((6.0/(Pe*dt))**0.5)*nr1
FORCAS(1,j,i,2)= alfa**(1.0/2.0)*((6.0/(Pe*dt))**0.5)*nr2
FORCAS(1,j,i,3)= alfa**(1.0/2.0)*((6.0/(Pe*dt))**0.5)*nr3
end do
 !$OMP END PARALLEL DO

else

 !$OMP PARALLEL DO
do j=1,rea
FORCAS(1,j,1,1)= 0.0
FORCAS(1,j,1,2)= 0.0
FORCAS(1,j,1,3)= 0.0

FORCAS(1,j,2,1)= 0.0
FORCAS(1,j,2,2)= 0.0
FORCAS(1,j,2,3)= 0.0
end do
 !$OMP END PARALLEL DO
end if

if(comentarios) then
write(*,*) 'Forcas brownianas:'
do j=1,rea
do i=1,2
write(*,*) FORCAS(1,j,i,1),FORCAS(1,j,i,2),FORCAS(1,j,i,3)
end do
end do
pause
end if

!*** Determinando as forcas magneticas devido as interacoes entre momentos de dipolo ***!
if(magnetico)then
 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
do q=1,2
if(i.ne.q) then
! Calcula-se a distancia entre a particula em questao e as outras particulas
r=( ((X(j,i,1)-X(j,q,1))**2.0) + ((X(j,i,2)-X(j,q,2))**2.0) + ((X(j,i,3)-X(j,q,3))**2.0))**0.5
if(r.le.dist) then
aux1(j,q)=0.0
aux2(j,q)=0.0
aux2(j,q)=0.0
else
! Calculando o vetor R_{ij} que liga uma particula i a uma particula j

rij(1)=X(j,i,1)-X(j,q,1)
rij(2)=X(j,i,2)-X(j,q,2)
rij(3)=X(j,i,3)-X(j,q,3)

! Normalizando o vetor R_{ij}
modrij=((rij(1)**2.0)+(rij(2)**2.0)+(rij(3)**2.0))**0.5
rij(1)=rij(1)/modrij
rij(2)=rij(2)/modrij
rij(3)=rij(3)/modrij

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(1)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,1)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,1)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))   &
+(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(1)

aux1(j,q)=(Stm/(r**4.0))*(termo1+termo2+termo3-5*termo4)  

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(2)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,2)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,2)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))   &
+(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(2)

aux2(j,q)=(Stm/(r**4.0))*(termo1+termo2+termo3-5*termo4)

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(3)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,3)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,3)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+   &
(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(3)

aux3(j,q)=(Stm/(r**4.0))*(termo1+termo2+termo3-5*termo4)

end if
end if
end do

FORCAS(2,j,i,1)=sum(aux1(j,:))
FORCAS(2,j,i,2)=sum(aux2(j,:))
FORCAS(2,j,i,3)=sum(aux3(j,:))

aux1=0.0
aux2=0.0
aux3=0.0
end do
end do
 !$OMP END PARALLEL DO
end if

if(comentarios) then
write(*,*) 'Forcas por interacao magnetica:'
do j=1,rea
do i=1,2
write(*,*) FORCAS(2,j,i,1),FORCAS(2,j,i,2),FORCAS(2,j,i,3)
end do
end do
pause
end if

!*******************FORCA POR INTERACAO HIDRODINAMICA **************************************!


 !$OMP PARALLEL DO
do j=1,rea

r=( ((X(j,1,1)-X(j,2,1))**2.0) + ((X(j,1,2)-X(j,2,2))**2.0) + ((X(j,1,3)-X(j,2,3))**2.0))**0.5
qui=abs(r-1-beta)
s=((2.0*r)*beta/(1.0+beta))


! Montando as funcoes de beta para calculo das funcoes resistencia X11a,X12a,X21a e determinacao 
! das forcas por interacoes hidrodinamicas 

f0=1.0
f1=3.0*beta
f2=9.0*beta
f3=(-4.0*beta)+(27.0*(beta**2.0))-(4.0*(beta)**3.0)
f4=(-24.0*beta)+(81.0*(beta)**2.0)+(36.0*(beta)**3.0)
f5=(72.0*(beta)**2.0)+(243.0*(beta)**3.0)+(72.0*(beta)**4.0)
f6=(16.0*beta)+(108.0*(beta)**2.0)+(281.0*(beta)**3.0)   &
+(648.0*(beta)**4.0)+(144.0*(beta)**5.0)
f7=(288.0*(beta)**2.0)+(1620.0*(beta)**3.0)+(1515.0*(beta)**4.0)   &
+(1620.0*(beta)**5.0)+(288.0*(beta)**6.0)
f8=(576.0*(beta)**2.0)+(4848.0*(beta)**3.0)+(5409.0*(beta)**4.0)   &
+(4524.0*(beta)**5.0)+(3888.0*(beta)**6.0)+(576.0*(beta)**7.0)
f9=(1152.0*(beta)**2.0)+(9072.0*(beta)**3.0)+(14752.0*(beta)**4.0)   &
+(26163.0*(beta)**5.0)+(14752.0*(beta)**6.0)+(9072.0*(beta)**7.0)   &
+(1152.0*(beta)**8.0)
f10=(2304.0*(beta)**2.0)+(20736.0*(beta)**3.0)+(42804.0*(beta)**4.0)   &
+(115849.0*(beta)**5.0)+(76176.0*(beta)**6.0)+(39264.0*(beta)**7.0)   &
+(20736.0*(beta)**8.0)+(2304.0*(beta)**9.0)
f11=(4608.0*(beta)**2.0)+(46656.0*(beta)**3.0)+(108912.0*(beta)**4.0)   &
+(269100.0*(beta)**5.0)+(319899.0*(beta)**6.0)+(269100.0*(beta)**7.0)   &
+(108912.0*(beta)**8.0)+(46656.0*(beta)**9.0)+(4608.0*(beta)**10.0)

! Montando as funcoes g1,g2,g3

g1=(2.0*(beta**2.0))/((1.0+beta)**3.0) 
g2=(0.2*beta*(1.0+7.0*beta+beta**2.0))/((1.0+beta)**3.0)
g3=(1.0/42.0)*(1.0+(18.0*beta)-(29.0*beta**2.0)+(18.0*beta**3.0)+(beta**4.0))/((1.0+beta)**3.0)

! Montando as funcoes m1

m1(1)=-1.0
m1(2)=-2.0
m1(3)=1.0
m1(4)=2.0
m1(5)=3.0
m1(6)=4.0
m1(7)=5.0
m1(8)=6.0
m1(9)=7.0
m1(10)=8.0
m1(11)=9.0

! Montando as funcoes resistencia X11a,X12a,X21a

X11a=(g1/(1.0-(4.0/(s**2.0)))) - g2*log(1.0-(4.0/(s**2.0)))   &
-(g3*(1.0-(4.0/(s**2.0))))*log(1.0-(4.0/(s**2.0))) + f0 -g1   &
+(((1.0/(2.0**2.0))/((1.0+beta)**2.0))*f2 - g1 - g2 + (4.0/(m1(2)*2.0))*g3)*((2.0/s)**2.0)   &
+(((1.0/(2.0**4.0))/((1.0+beta)**4.0))*f4 - g1 - 0.5*g2 + (4.0/(m1(4)*4.0))*g3)*((2.0/s)**4.0)   &
+(((1.0/(2.0**6.0))/((1.0+beta)**6.0))*f6 - g1 - 0.333*g2 + (4.0/(m1(6)*6.0))*g3)*((2.0/s)**6.0)   &
+(((1.0/(2.0**8.0))/((1.0+beta)**8.0))*f8 - g1 - 0.25*g2 + (4.0/(m1(8)*8.0))*g3)*((2.0/s)**8.0)   &
+(((1.0/(2.0**10.0))/((1.0+beta)**10.0))*f10 - g1 - 0.2*g2 + (4.0/(m1(10)*10.0))*g3)*((2.0/s)**10.0)  


X12a=(2.0/s)*(g1/(1.0-(4.0/(s**2.0)))) + g2*log((s+2.0)/(s-2.0))   & 
+ (g3*(1.0-(4.0/(s**2.0))))*log((s+2.0)/(s-2.0)) + 4.0*(g3/s)   &
+(((1.0/(2.0**1.0))/((1.0+beta)**1.0))*f1 - g1 - (2.0/1.0)*g2 + (4.0/(m1(1)*1.0))*g3)*((2.0/s)**1.0)   &
+(((1.0/(2.0**3.0))/((1.0+beta)**3.0))*f3 - g1 - (2.0/3.0)*g2 + (4.0/(m1(3)*3.0))*g3)*((2.0/s)**3.0)   &
+(((1.0/(2.0**5.0))/((1.0+beta)**5.0))*f5 - g1 - (2.0/5.0)*g2 + (4.0/(m1(5)*5.0))*g3)*((2.0/s)**5.0)   &
+(((1.0/(2.0**7.0))/((1.0+beta)**7.0))*f7 - g1 - (2.0/7.0)*g2 + (4.0/(m1(7)*7.0))*g3)*((2.0/s)**7.0)   &
+(((1.0/(2.0**9.0))/((1.0+beta)**9.0))*f9 - g1 - (2.0/9.0)*g2 + (4.0/(m1(9)*9.0))*g3)*((2.0/s)**9.0)   &
+(((1.0/(2.0**11.0))/((1.0+beta)**11.0))*f11 - g1 - (2.0/11.0)*g2 + (4.0/(m1(11)*11.0))*g3)*((2.0/s)**11.0)   


X12a=(-2.0/(1.0+beta))*X12a

X21a=X12a


! Montando as funcoes de beta para calculo das funcoes resistencia X22a 

f0=1.0
f2=9.0*alfa
f4=(-24.0*alfa)+(81.0*(alfa)**2.0)+(36.0*(alfa)**3.0)
f6=(16.0*alfa)+(108.0*(alfa)**2.0)+(281.0*(alfa)**3.0)   &
+(648.0*(alfa)**4.0)+(144.0*(alfa)**5.0)
f8=(576.0*(alfa)**2.0)+(4848.0*(alfa)**3.0)+(5409.0*(alfa)**4.0)   &
+(4524.0*(alfa)**5.0)+(3888.0*(alfa)**6.0)+(576.0*(alfa)**7.0)
f10=(2304.0*(alfa)**2.0)+(20736.0*(alfa)**3.0)+(42804.0*(alfa)**4.0)   &
+(115849.0*(alfa)**5.0)+(76176.0*(alfa)**6.0)+(39264.0*(alfa)**7.0)   &
+(20736.0*(alfa)**8.0)+(2304.0*(alfa)**9.0)

g1=(2.0*(alfa**2.0))/((1.0+alfa)**3.0) 
g2=(0.2*alfa*(1.0+7.0*alfa+alfa**2.0))/((1.0+alfa)**3.0)
g3=(1.0/42.0)*(1.0+(18.0*alfa)-(29.0*alfa**2.0)+(18.0*alfa**3.0)+(alfa**4.0))/((1.0+alfa)**3.0)


! Montando a funcao resistencia X22a

X22a=(g1/(1.0-(4.0/(s**2.0)))) - g2*log(1.0-(4.0/(s**2.0)))   &
-(g3*(1.0-(4.0/(s**2.0))))*log(1.0-(4.0/(s**2.0))) + f0 -g1   &
+(((1.0/(2.0**2.0))/((1.0+alfa)**2.0))*f2 - g1 - g2 + (4.0/(m1(2)*2.0))*g3)*((2.0/s)**2.0)   &
+(((1.0/(2.0**4.0))/((1.0+alfa)**4.0))*f4 - g1 - 0.5*g2 + (4.0/(m1(4)*4.0))*g3)*((2.0/s)**4.0)   &
+(((1.0/(2.0**6.0))/((1.0+alfa)**6.0))*f6 - g1 - 0.333*g2 + (4.0/(m1(6)*6.0))*g3)*((2.0/s)**6.0)   &
+(((1.0/(2.0**8.0))/((1.0+alfa)**8.0))*f8 - g1 - 0.25*g2 + (4.0/(m1(8)*8.0))*g3)*((2.0/s)**8.0)   &
+(((1.0/(2.0**10.0))/((1.0+alfa)**10.0))*f10 - g1 - 0.2*g2 + (4.0/(m1(10)*10.0))*g3)*((2.0/s)**10.0)  

! Montando as funcoes de beta para calculo das funcoes resistencia Y11a,Y12a,Y21a 

f0=1.0
f1=(3.0/2.0)*beta
f2=(9.0/4.0)*beta
f3=(2.0*beta)+(27.0/8.0)*(beta**2.0)+(2.0*beta**3.0)
f4=(96.0/16.0)*beta+(81.0/16.0)*(beta**2.0)+(288.0/16.0)*(beta**3.0)
f5=(1008.0/32.0)*(beta**2.0)+(243.0/32.0)*(beta**3.0)+(1008.0/32.0)*(beta**4.0)
f6=(4.0*beta)+(3456.0/64.0)*(beta**2.0)+(1241.0/64.0)*(beta**3.0)   &
+(5184.0/64.0)*(beta**4.0)+(4608.0/64.0)*(beta**5.0)
f7=(18432.0/128.0)*(beta**2.0)+(16848.0/128.0)*(beta**3.0)   &
+(19083.0/128.0)*(beta**4.0)+(16848.0/128.0)*(beta**5.0)   &
+(18432.0/128.0)*(beta**6.0)
f8=(71424.0/256.0)*(beta**2.0)+(136352.0/256.0)*(beta**3.0)   &
+(126369.0/256.0)*(beta**4.0)-(3744.0/256.0)*(beta**5.0)   &
+(165888.0/256.0)*(beta**6.0)+(73728.0/256.0)*(beta**7.0)
f9=(294912.0/512.0)*(beta**2.0)+(580608.0/512.0)*(beta**3.0)   &
+(967088.0/512.0)*(beta**4.0)+(766179.0/512.0)*(beta**5.0)   &
+(967088.0/512.0)*(beta**6.0)+(580608.0/512.0)*(beta**7.0)   &
+(294912.0/512.0)*(beta**8.0)
f10=(1179648.0/1024.0)*(beta**2.0)+(2011392.0/1024.0)*(beta**3.0)   &
+(6303168.0/1024.0)*(beta**4.0)+(10548393.0/1024.0)*(beta**5.0)   &
+(8654976.0/1024.0)*(beta**6.0)-(179712.0/1024.0)*(beta**7.0)   &
+(3981312.0/1024.0)*(beta**8.0)+(1179648.0/1024.0)*(beta**9.0)
f11=(4718592.0/2048.0)*(beta**2.0)+(14598144.0/2048.0)*(beta**3.0)   &
+(22600704.0/2048.0)*(beta**4.0)+(43912080.0/2048.0)*(beta**5.0)   &
+(95203835.0/2048.0)*(beta**6.0)+(43912080.0/2048.0)*(beta**7.0)   &
+(22600704.0/2048.0)*(beta**8.0)+(14598144.0/2048.0)*(beta**9.0)   &
+(4718592.0/2048.0)*(beta**10.0)

! Montando as funcoes g2,g3

g2=(0.2667*beta*(2.0+beta+2.0*beta**2.0))/((1.0+beta)**3.0)
g3=(2.0/375.0)*(16.0-(45.0*beta)+(58.0*beta**2.0)-(45.0*beta**3.0)+(16*beta**4.0))/((1.0+beta)**3.0)


! Montando as funcoes resistencia Y11a,Y12a,Y21a

Y11a=-g2*log(1.0-(4.0/(s**2.0)))-g3*(1.0-(4.0/(s**2.0)))*log(1.0-(4.0/(s**2.0)))+f0   &
+(((1.0/(2.0**2.0))/((1.0+beta)**2.0))*f2 - g2 + (4.0/(m1(2)*2.0))*g3)*((2.0/s)**2.0)   &
+(((1.0/(2.0**4.0))/((1.0+beta)**4.0))*f4 - 0.5*g2 + (4.0/(m1(4)*4.0))*g3)*((2.0/s)**4.0)   &
+(((1.0/(2.0**6.0))/((1.0+beta)**6.0))*f6 - 0.333*g2 + (4.0/(m1(6)*6.0))*g3)*((2.0/s)**6.0)   &
+(((1.0/(2.0**8.0))/((1.0+beta)**8.0))*f8 - 0.25*g2 + (4.0/(m1(8)*8.0))*g3)*((2.0/s)**8.0)   &
+(((1.0/(2.0**10.0))/((1.0+beta)**10.0))*f10 - 0.2*g2 + (4.0/(m1(10)*10.0))*g3)*((2.0/s)**10.0)  

Y12a=g2*log((s+2.0)/(s-2.0))+g3*(1.0-(4.0/(s**2.0)))*log((s+2.0)/(s-2.0))+(4.0*g3/s)   &
+(((1.0/(2.0**1.0))/((1.0+beta)**1.0))*f1 - (2.0/1.0)*g2 + (4.0/(m1(1)*1.0))*g3)*((2.0/s)**1.0)   &
+(((1.0/(2.0**3.0))/((1.0+beta)**3.0))*f3 - (2.0/3.0)*g2 + (4.0/(m1(3)*3.0))*g3)*((2.0/s)**3.0)   &
+(((1.0/(2.0**5.0))/((1.0+beta)**5.0))*f5 - (2.0/5.0)*g2 + (4.0/(m1(5)*5.0))*g3)*((2.0/s)**5.0)   &
+(((1.0/(2.0**7.0))/((1.0+beta)**7.0))*f7 - (2.0/7.0)*g2 + (4.0/(m1(7)*7.0))*g3)*((2.0/s)**7.0)   &
+(((1.0/(2.0**9.0))/((1.0+beta)**9.0))*f9 - (2.0/9.0)*g2 + (4.0/(m1(9)*9.0))*g3)*((2.0/s)**9.0)   &
+(((1.0/(2.0**11.0))/((1.0+beta)**11.0))*f11 - (2.0/11.0)*g2 + (4.0/(m1(11)*11.0))*g3)*((2.0/s)**11.0)   
Y12a=(-2.0/(1.0+beta))*Y12a

Y21a=Y12a

! Montando as funcoes de beta para calculo das funcoes resistencia Y22a 


f0=1.0
f2=(9.0/4.0)*alfa
f4=(96.0/16.0)*alfa+(81.0/16.0)*(alfa**2.0)+(288.0/16.0)*(alfa**3.0)
f6=(4.0*alfa)+(3456.0/64.0)*(alfa**2.0)+(1241.0/64.0)*(alfa**3.0)   &
+(5184.0/64.0)*(alfa**4.0)+(4608.0/64.0)*(alfa**5.0)
f8=(71424.0/256.0)*(alfa**2.0)+(136352.0/256.0)*(alfa**3.0)   &
+(126369.0/256.0)*(alfa**4.0)-(3744.0/256.0)*(alfa**5.0)   &
+(165888.0/256.0)*(alfa**6.0)+(73728.0/256.0)*(alfa**7.0)
f10=(1179648.0/1024.0)*(alfa**2.0)+(2011392.0/1024.0)*(alfa**3.0)   &
+(6303168.0/1024.0)*(alfa**4.0)+(10548393.0/1024.0)*(alfa**5.0)   &
+(8654976.0/1024.0)*(alfa**6.0)-(179712.0/1024.0)*(alfa**7.0)   &
+(3981312.0/1024.0)*(alfa**8.0)+(1179648.0/1024.0)*(alfa**9.0)

g2=(0.2667*alfa*(2.0+alfa+2.0*alfa**2.0))/((1.0+alfa)**3.0)
g3=(2.0/375.0)*(16.0-(45.0*alfa)+(58.0*alfa**2.0)-(45.0*alfa**3.0)+(16*alfa**4.0))/((1.0+alfa)**3.0)

! Montando as funcoes resistencia

Y22a=-g2*log(1.0-(4.0/(s**2.0)))-g3*(1.0-(4.0/(s**2.0)))*log(1.0-(4.0/(s**2.0)))+f0   &
+(((1.0/(2.0**2.0))/((1.0+alfa)**2.0))*f2 - g2 + (4.0/(m1(2)*2.0))*g3)*((2.0/s)**2.0)   &
+(((1.0/(2.0**4.0))/((1.0+alfa)**4.0))*f4 - 0.5*g2 + (4.0/(m1(4)*4.0))*g3)*((2.0/s)**4.0)   &
+(((1.0/(2.0**6.0))/((1.0+alfa)**6.0))*f6 - 0.333*g2 + (4.0/(m1(6)*6.0))*g3)*((2.0/s)**6.0)   &
+(((1.0/(2.0**8.0))/((1.0+alfa)**8.0))*f8 - 0.25*g2 + (4.0/(m1(8)*8.0))*g3)*((2.0/s)**8.0)   &
+(((1.0/(2.0**10.0))/((1.0+alfa)**10.0))*f10 - 0.2*g2 + (4.0/(m1(10)*10.0))*g3)*((2.0/s)**10.0)  

! Calculando as forcas hidrodinamicas 

! Primeiramente vamos determinar o vetor normalizado d=r/|r| para a esfera 1:

d(j,1)=(X(j,1,1)-X(j,2,1))/r
d(j,2)=(X(j,1,2)-X(j,2,2))/r
d(j,3)=(X(j,1,3)-X(j,2,3))/r

modd(j)=((d(j,1)**2.0)+(d(j,2)**2.0)+(d(j,3)**2.0))**0.5

d(j,1)=d(j,1)/modd(j)
d(j,2)=d(j,2)/modd(j)
d(j,3)=d(j,3)/modd(j)


FORCAS(3,j,1,1)=X11a*(d(j,1)*d(j,1)*U(j,1,1)+d(j,1)*d(j,2)*U(j,1,2)+d(j,1)*d(j,3)*U(j,1,3))   &
+Y11a*((1-d(j,1)*d(j,1))*U(j,1,1)-d(j,1)*d(j,2)*U(j,1,2)-d(j,1)*d(j,3)*U(j,1,3))   &
+X12a*(d(j,1)*d(j,1)*U(j,2,1)+d(j,1)*d(j,2)*U(j,2,2)+d(j,1)*d(j,3)*U(j,2,3))   &
+Y12a*((1-d(j,1)*d(j,1))*U(j,2,1)-d(j,1)*d(j,2)*U(j,2,2)-d(j,1)*d(j,3)*U(j,2,3))   

FORCAS(3,j,1,2)=X11a*(d(j,2)*d(j,1)*U(j,1,1)+d(j,2)*d(j,2)*U(j,1,2)+d(j,2)*d(j,3)*U(j,1,3))   &
+Y11a*((1-d(j,2)*d(j,2))*U(j,1,2)-d(j,2)*d(j,1)*U(j,1,1)-d(j,2)*d(j,3)*U(j,1,3))   &
+X12a*(d(j,2)*d(j,1)*U(j,2,1)+d(j,2)*d(j,2)*U(j,2,2)+d(j,2)*d(j,3)*U(j,2,3))   &
+Y12a*((1-d(j,2)*d(j,2))*U(j,2,2)-d(j,2)*d(j,1)*U(j,2,1)-d(j,2)*d(j,3)*U(j,2,3))   

FORCAS(3,j,1,3)=X11a*(d(j,3)*d(j,1)*U(j,1,1)+d(j,3)*d(j,2)*U(j,1,2)+d(j,3)*d(j,3)*U(j,1,3))   &
+Y11a*((1-d(j,3)*d(j,3))*U(j,1,3)-d(j,3)*d(j,1)*U(j,1,1)-d(j,3)*d(j,2)*U(j,1,2))   &
+X12a*(d(j,3)*d(j,1)*U(j,2,1)+d(j,3)*d(j,2)*U(j,2,2)+d(j,3)*d(j,3)*U(j,2,3))   &
+Y12a*((1-d(j,3)*d(j,3))*U(j,2,3)-d(j,3)*d(j,1)*U(j,2,1)-d(j,3)*d(j,2)*U(j,2,2))   


! Agora vamos determinar o vetor normalizado d=r/|r| para a esfera 2:

d(j,1)=(X(j,2,1)-X(j,1,1))/r
d(j,2)=(X(j,2,2)-X(j,1,2))/r
d(j,3)=(X(j,2,3)-X(j,1,3))/r

modd(j)=((d(j,1)**2.0)+(d(j,2)**2.0)+(d(j,3)**2.0))**0.5

d(j,1)=d(j,1)/modd(j)
d(j,2)=d(j,2)/modd(j)
d(j,3)=d(j,3)/modd(j)

FORCAS(3,j,2,1)=X21a*(d(j,1)*d(j,1)*U(j,1,1)+d(j,1)*d(j,2)*U(j,1,2)+d(j,1)*d(j,3)*U(j,1,3))   &
+Y21a*((1-d(j,1)*d(j,1))*U(j,1,1)-d(j,1)*d(j,2)*U(j,1,2)-d(j,1)*d(j,3)*U(j,1,3))   &
+X22a*(d(j,1)*d(j,1)*U(j,2,1)+d(j,1)*d(j,2)*U(j,2,2)+d(j,1)*d(j,3)*U(j,2,3))   &
+Y22a*((1-d(j,1)*d(j,1))*U(j,2,1)-d(j,1)*d(j,2)*U(j,2,2)-d(j,1)*d(j,3)*U(j,2,3))   

FORCAS(3,j,2,2)=X21a*(d(j,2)*d(j,1)*U(j,1,1)+d(j,2)*d(j,2)*U(j,1,2)+d(j,2)*d(j,3)*U(j,1,3))   &
+Y21a*((1-d(j,2)*d(j,2))*U(j,1,2)-d(j,2)*d(j,1)*U(j,1,1)-d(j,2)*d(j,3)*U(j,1,3))   &
+X22a*(d(j,2)*d(j,1)*U(j,2,1)+d(j,2)*d(j,2)*U(j,2,2)+d(j,2)*d(j,3)*U(j,2,3))   &
+Y22a*((1-d(j,2)*d(j,2))*U(j,2,2)-d(j,2)*d(j,1)*U(j,2,1)-d(j,2)*d(j,3)*U(j,2,3))   

FORCAS(3,j,2,3)=X21a*(d(j,3)*d(j,1)*U(j,1,1)+d(j,3)*d(j,2)*U(j,1,2)+d(j,3)*d(j,3)*U(j,1,3))   &
+Y21a*((1-d(j,3)*d(j,3))*U(j,1,3)-d(j,3)*d(j,1)*U(j,1,1)-d(j,3)*d(j,2)*U(j,1,2))   &
+X22a*(d(j,3)*d(j,1)*U(j,2,1)+d(j,3)*d(j,2)*U(j,2,2)+d(j,3)*d(j,3)*U(j,2,3))   &
+Y22a*((1-d(j,3)*d(j,3))*U(j,2,3)-d(j,3)*d(j,1)*U(j,2,1)-d(j,3)*d(j,2)*U(j,2,2))   

FORCAS(3,j,1,1)=-FORCAS(3,j,1,1)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))
FORCAS(3,j,1,2)=-FORCAS(3,j,1,2)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))
FORCAS(3,j,1,3)=-FORCAS(3,j,1,3)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))

FORCAS(3,j,2,1)=-FORCAS(3,j,2,1)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))
FORCAS(3,j,2,2)=-FORCAS(3,j,2,2)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))
FORCAS(3,j,2,3)=-FORCAS(3,j,2,3)*(beta+1.0)/(2.0*beta)!*(1/(6.0*acos(-1.0)))



end do
 !$OMP END PARALLEL DO

if(comentarios) then
write(*,*) 'Forcas por interacao hidrodinamica:'
do j=1,rea
do i=1,2
write(*,*) FORCAS(3,j,i,1),FORCAS(3,j,i,2),FORCAS(3,j,i,3)
end do
end do
pause
end if



!********Determinando as forcas magneticas devido a um campo externo aplicado***************!
!if(magnetico)then
!if(externo) then
!do j=1,rea
!do i=1,2
!if(X(j,i,3).le.(h-2.0)) then
!FORCAS(4,j,i,3)= 2.0*(alfa_campo/Pe)*(Di(j,i,3))/((100.0-X(j,i,3))**3.0)
!end if
!end do
!end do
!end if
!end if

if(comentarios) then
write(*,*) 'Forcas por campo externo:'
do j=1,rea
do i=1,2
write(*,*) FORCAS(4,j,i,1),FORCAS(4,j,i,2),FORCAS(4,j,i,3)
end do
end do
pause
end if

!*******************************************************************************************!
! Calculando as forcas totais que agem sobre cada particula

 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
Ft(j,i,1)= FORCAS(1,j,i,1) + FORCAS(2,j,i,1) + FORCAS(3,j,i,1)
Ft(j,i,2)= FORCAS(1,j,i,2) + FORCAS(2,j,i,2) + FORCAS(3,j,i,2)
Ft(j,i,3)= FORCAS(1,j,i,3) + FORCAS(2,j,i,3) + FORCAS(3,j,i,3) !+ FORCAS(4,j,i,3)
end do
end do
 !$OMP END PARALLEL DO

!*****Iniciando o processo de solucao da equacao do momento da quantidade de movimento *****!


! Este processo serve na verdade para fornecer uma evolucao temporal dos vetores momento de 
! dipolo magnetico de cada particula a partir da equacao do momento angular


! Vendo a distancia entre as particulas em todas as realizacoes para contar a quantidade de trajetorias
! agregativas e difusivas

 !$OMP PARALLEL DO
do j=1,rea
r=(((X(j,1,1)-X(j,2,1))**2.0) + ((X(j,1,2)-X(j,2,2))**2.0) + ((X(j,1,3)-X(j,2,3))**2.0))**0.5
contagem(j)=r
if(r.le.(1.0+(1.0/beta)))then
agregativa=agregativa+1
end if
end do
 !$OMP END PARALLEL DO



if(torque) then

!*********** Determinando os torques brownianos sobre as particulas ************************!

if(browniano) then

 call randomica(-1.0,1.0,nr,(3*2*rea),npast+k)
 !$OMP PARALLEL DO
do j=1,rea
do i=1,2

nr1= nr((i*2+(i-2)+(2*3*(j-1))))
nr2= nr((i*2+(i-1)+(2*3*(j-1))))
nr3= nr((i*2+(i)+(2*3*(j-1))))
modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5
nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand
 
TORQUES(1,j,i,1) = ((6.0/(Per*dt))**0.5)*nr1
TORQUES(1,j,i,2) = ((6.0/(Per*dt))**0.5)*nr2
TORQUES(1,j,i,3) = ((6.0/(Per*dt))**0.5)*nr3
end do
end do
 !$OMP END PARALLEL DO
else
 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
TORQUES(1,j,i,1) = 0.0
TORQUES(1,j,i,2) = 0.0
TORQUES(1,j,i,3) = 0.0
end do
end do
 !$OMP END PARALLEL DO
end if

if(comentarios) then
write(*,*) 'Torques brownianos:'
do j=1,rea
do i=1,2
write(*,*) TORQUES(1,j,i,1),TORQUES(1,j,i,2),TORQUES(1,j,i,3)
end do
end do
pause
end if

!************************ Determinando os torques magneticos sobre as particulas ***********!
if(magnetico) then
 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
do q=1,2
if(i.ne.q) then
! Calcula-se a distancia entre a particula em questao e as outras particulas
r=( ((X(j,i,1)-X(j,q,1))**2.0) + ((X(j,i,2)-X(j,q,2))**2.0) + ((X(j,i,3)-X(j,q,3))**2.0))**0.5
if(r.le.dist) then
aux1(j,q)=0.0
aux2(j,q)=0.0
aux3(j,q)=0.0
else

! Calculando o vetor R_{ij} que liga uma particula i a uma particula j

rij(1)=X(j,i,1)-X(j,q,1)
rij(2)=X(j,i,2)-X(j,q,2)
rij(3)=X(j,i,3)-X(j,q,3)

! Normalizando o vetor R_{ij}
modrij=((rij(1)**2.0)+(rij(2)**2.0)+(rij(3)**2.0))**0.5
rij(1)=rij(1)/modrij
rij(2)=rij(2)/modrij
rij(3)=rij(3)/modrij

termo1=((Di(j,i,2)*Di(j,q,3))-(Di(j,i,3)*Di(j,q,2)))
termo2=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))
termo3=((Di(j,i,2)*rij(3))-(Di(j,i,3)*rij(2)))

aux1(j,q)=(Stmr/(r**3.0))*(((-1.0/3.0)*termo1)+(termo2*termo3))  

termo1=((Di(j,i,3)*Di(j,q,1))-(Di(j,i,1)*Di(j,q,3)))
termo2=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))
termo3=((Di(j,i,3)*rij(1))-(Di(j,i,1)*rij(3)))

aux2(j,q)=(Stmr/(r**3.0))*(((-1.0/3.0)*termo1)+(termo2*termo3))  

termo1=((Di(j,i,1)*Di(j,q,2))-(Di(j,i,2)*Di(j,q,1)))
termo2=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))
termo3=((Di(j,i,1)*rij(2))-(Di(j,i,2)*rij(1)))

aux3(j,q)=(Stmr/(r**3.0))*(((-1.0/3.0)*termo1)+(termo2*termo3))  

end if
end if

end do
TORQUES(2,j,i,1)=sum(aux1(j,:))
TORQUES(2,j,i,2)=sum(aux2(j,:))
TORQUES(2,j,i,3)=sum(aux3(j,:))

aux1=0.0
aux2=0.0
aux3=0.0
end do
end do
 !$OMP END PARALLEL DO
end if

if(comentarios) then
write(*,*) 'Torques por interacao magnetica:'
do j=1,rea
do i=1,2
write(*,*) TORQUES(2,j,i,1),TORQUES(2,j,i,2),TORQUES(2,j,i,3)
end do
end do
pause
end if

!***Determinando o torque magnetico devido a aplicacao de um campo externo****!
if(magnetico)then
if(externo) then
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
TORQUES(3,j,i,1)=(alfa_campo/Per)*Di(j,i,2)
TORQUES(3,j,i,2)=-(alfa_campo/Per)*Di(j,i,1)
TORQUES(3,j,i,3)=0.0
end do
end do
 !$OMP END PARALLEL DO
end if
end if

if(comentarios) then
write(*,*) 'Torques devido a um campo externo:'
do j=1,rea
do i=1,2
write(*,*) TORQUES(3,j,i,1),TORQUES(3,j,i,2),TORQUES(3,j,i,3)
end do
end do
pause
end if

!***********Determinando o somatorio dos torques que atuam sobre a particula****************!

 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) + TORQUES(3,j,i,1)
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) + TORQUES(3,j,i,1)
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) + TORQUES(3,j,i,1)
end do
end do
 !$OMP END PARALLEL DO

! Calculando a velocidade angular das particulas (Utilizando Runge-Kutta de 4 ordem)

 !$OMP PARALLEL DO
do j=1,rea
 call resomega(W(j,1,1),dt,Str,Tt(j,1,1),1.0)
 call resomega(W(j,1,2),dt,Str,Tt(j,1,2),1.0)
 call resomega(W(j,1,3),dt,Str,Tt(j,1,3),1.0)

 call resomega(W(j,2,1),dt,(Str/(beta**4.0)),Tt(j,2,1),beta**3.0)
 call resomega(W(j,2,2),dt,(Str/(beta**4.0)),Tt(j,2,2),beta**3.0)
 call resomega(W(j,2,3),dt,(Str/(beta**4.0)),Tt(j,2,3),beta**3.0)
end do
 !$OMP END PARALLEL DO

if(comentarios) then
write(*,*) 'Velocidade angular na iteracao atual:'
do j=1,rea
do i=1,2
write(*,*) W(j,i,1),W(j,i,2),W(j,i,3)
end do
end do
pause
end if

! Evoluindo o vetor momento de dipolo magneticos das particulas

 !$OMP PARALLEL DO
do j=1,rea
do i=1,2
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt)
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt)
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt)
end do
end do
 !$OMP END PARALLEL DO

end if

if(comentarios) then
write(*,*) 'Momentos de dipolo na iteracao atual:'
do j=1,rea
do i=1,2
write(*,*) Di(j,i,1),Di(j,i,2),Di(j,i,3)
end do
end do
pause
end if

! Resolvendo a velocidade de cada particula

 !$OMP PARALLEL DO
do j=1,rea
 call resvel(U(j,1,1),dt,St,Ft(j,1,1),0.0,0.0)
 call resvel(U(j,1,2),dt,St,Ft(j,1,2),0.0,0.0)
 call resvel(U(j,1,3),dt,St,Ft(j,1,3),1.0,0.0)
  
 call resvel(U(j,2,1),dt,(St/((beta**3.0)*lambda)),Ft(j,2,1),0.0,0.0)
 call resvel(U(j,2,2),dt,(St/((beta**3.0)*lambda)),Ft(j,2,2),0.0,0.0)
 call resvel(U(j,2,3),dt,(St/((beta**3.0)*lambda)),Ft(j,2,3),(1.0/(gama*(beta**3.0))),0.0)
end do
 !$OMP END PARALLEL DO

if(comentarios) then
write(*,*) 'Velocidade na iteracao atual:'
do j=1,rea
do i=1,2
write(*,*) U(j,i,1),U(j,i,2),U(j,i,3)
end do
end do
pause
end if

! Calculando a velocidade relativa entre as particulas

 !$OMP PARALLEL DO
do j=1,rea
do q=1,3
Urel(j,q)=-U(j,1,q)+U(j,2,q)
end do
end do
 !$OMP END PARALLEL DO

! Resolvendo a posicao de cada particula

 !$OMP PARALLEL DO
do j=1,rea
 call respos(X(j,1,1),dt,U(j,1,1))
 call respos(X(j,1,2),dt,U(j,1,2))
 call respos(X(j,1,3),dt,U(j,1,3))

 call respos(X(j,2,1),dt,U(j,2,1))
 call respos(X(j,2,2),dt,U(j,2,2))
 call respos(X(j,2,3),dt,U(j,2,3))
end do
 !$OMP END PARALLEL DO


! VERIFICANDO QUANTAS TRAJETORIAS TIVERAM COLISAO

! Determinando a distancia entre as particulas


do j=1,rea
i=1
q=2
! Calcula-se a distancia entre a particula em questao e as outras particulas
r=( ((X(j,i,1)-X(j,q,1))**2.0) + ((X(j,i,2)-X(j,q,2))**2.0) + ((X(j,i,3)-X(j,q,3))**2.0))**0.5
! Se houve colisao entao eu vou marcar em qual realizacao ela ocorreu, para nao plotar 
! essa realizacao no diagrama de reversibilidade
if(r.le.(1.0+alfa+0.2)) then
 colisao(j)=1
end if
end do



! Resolvendo a trajetoria relativa 

 !$OMP PARALLEL DO
do j=1,rea
do q=1,3
 Xrel(j,q)=X(j,2,q)-X(j,1,q)
end do
end do
 !$OMP END PARALLEL DO

if(mag_tempo) then
 call media_direcao(Di,N,rea,magtempo(k),3)
write(5*rea,*)magtempo(k),k*dt
end if
 

mod1=(((Di(1,1,1)**2.0)+(Di(1,1,2)**2.0)+(Di(1,1,3)**2.0))**0.5)+0.00001
mod2=(((Di(1,2,1)**2.0)+(Di(1,2,2)**2.0)+(Di(1,2,3)**2.0))**0.5)+0.00001


thetaux=((Di(1,1,1)*Di(1,2,1))+(Di(1,1,2)*Di(1,2,2))+(Di(1,1,3)*Di(1,2,3)))/(mod1*mod2)
theta=acos(thetaux)
modomega=((W(1,1,1)**2.0)+(W(1,1,2)**2.0)+(W(1,1,3)**2.0))**0.5

!write(*,*) thetaux,theta
!pause

write(100*rea,*) k*dt,(X(1,1,3)-X(1,2,3)),theta,modomega

! Escrevendo em um arquivo de saida a posicao e velocidade de cada particula em cada 
! instante de tempo alem da trajetoria relativa para cada realizacao


teste1=k/n3
teste2=k/n2

if(teste1.eq.teste2) then
write(204*rea,'(A12,I6,A1)') 'zone t="',k,'"' 
write(204*rea,*) X(1,1,1),X(1,1,2),X(1,1,3),Di(1,1,1),Di(1,1,2),Di(1,1,3)
write(204*rea,*) X(1,2,1),X(1,2,2),X(1,2,3),Di(1,2,1),Di(1,2,2),Di(1,2,3)

write(666*rea,734) X(1,1,1),X(1,1,2),X(1,1,3),Di(1,1,1),Di(1,1,2),Di(1,1,3),X(1,2,1),X(1,2,2),X(1,2,3),Di(1,2,1),Di(1,2,2),Di(1,2,3)
end if

if(teste1.eq.teste2) then
if(arquivo) then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"' 
do i=1,N
write(j,509)X(j,i,1),X(j,i,2),   &
X(j,i,3)
end do
end do

do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"' 
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do

do j=1,rea
write(rea+rea+j,509)Xrel(j,1),Xrel(j,2),   &
Xrel(j,3)
end do

end if
end if



do j=1,rea
distancia(j,k)=(((X(j,1,1)-X(j,2,1))**2.0) + ((X(j,1,2)-X(j,2,2))**2.0) + ((X(j,1,3)-X(j,2,3))**2.0))**0.5
end do

end do

! Escrevendo em um arquivo a distancia entre as esferas ao longo do calculo para todas as realizacoes

do j=1,rea
write(10*rea,*) 'zone t="',j,'"'

do k=1,npast
teste1=k/n3
teste2=k/n2
if(teste1.eq.teste2) then
write(10*rea,*) distancia(j,k), k*dt 
end if
end do
end do



if(arquivo) then
do j=1,rea
write(rea+rea+rea+j,509)Di(j,1,1),Di(j,1,2),   &
Di(j,1,3)
write(rea+rea+rea+j,509)Di(j,2,1),Di(j,2,2),   &
Di(j,2,3)
end do
end if




!********************** DETERMINANDO O DIAGRAMA DE REVERSIBILIDADE *************************!
if(diagrama) then
open(6*rea,file='diagrama_de_reversibilidade.plt')
write(6*rea,*)'Variables="X","Y"'


do j=1,rea
if(colisao(j).ne.1) then
write(6*rea,*)  X(j,2,1),X(j,2,2)
end if
end do


!do j=1,rea
!write(*,*) j, colisao(j)
!end do

reaok=sum(colisao)

!write(*,*) reaok

!pause

!************************** DETERMINANDO OS COEFICIENTES DE DIFUSAO ************************!

do j=1,rea
if(colisao(j).eq.0) then
posicao_final(j,1)=X(j,1,1) 
posicao_final(j,2)=X(j,1,2)

desloc_self(j,1)=(posicao_final(j,1)-posicao_inicial(j,1))**2.0
desloc_self(j,2)=(posicao_final(j,2)-posicao_inicial(j,2))**2.0

desloc_col(j,1)=(posicao_final(j,1)-posicao_inicial(j,1))*Xrel(j,1)
desloc_col(j,2)=(posicao_final(j,2)-posicao_inicial(j,2))*Xrel(j,2)
end if
end do



Dself1=((3.0/8.0)*acos(-1.0))*((0.5*(1.0+alfa))**3.0)*sum(desloc_self(:,1))
Dself2=((3.0/8.0)*acos(-1.0))*((0.5*(1.0+alfa))**3.0)*sum(desloc_self(:,2))

Dcol1=((3.0/8.0)*acos(-1.0))*((0.5*(1.0+alfa))**3.0)*sum(desloc_col(:,1))
Dcol2=((3.0/8.0)*acos(-1.0))*((0.5*(1.0+alfa))**3.0)*sum(desloc_col(:,2))

Dglobal1=Dself1+Dcol1
Dglobal2=Dself2+Dcol2

open(7*rea,file='difusao.plt')
write(7*rea,*)'Variables="Ds1","Ds2","Dc1","Dc2","D1","D2"'

write(7*rea,*)  Dself1,Dself2,Dcol1,Dcol2,Dglobal1,Dglobal2

do j=1,rea
posicao_final(j,1)=X(j,2,1)-posicao_inicial(j,1) 
posicao_final(j,2)=X(j,2,2)-posicao_inicial(j,2)
deltaquadrado(j,1)=posicao_final(j,1)**2.0
deltaquadrado(j,2)=posicao_final(j,2)**2.0
difusao_frc1=sum(deltaquadrado(:,1))/rea
difusao_frc2=sum(deltaquadrado(:,2))/rea
end do

write(*,*) 'A razao entre os deslocamentos ao quadrado medios e:',difusao_frc1/difusao_frc2

do j=1,8*rea
 close(j)
end do

else

do j=1,6*rea
 close(j)
end do

end if

difusiva=rea-agregativa

razao3=agregativa
razao4=difusiva
razao5=rea

razao1=(razao3/razao5)
razao2=(razao4/razao5)

!******* CHECAR SE razao1 e razao2 são realmente números reais ou estão definidos como inteiros *******!

write(*,*) 'A frequencia de trajetorias agregativas e:', razao1 , agregativa, rea
write(*,*) 'A frequencia de trajetorias difusivas e:', razao2 , difusiva, rea


! Desalocando as matrizes da memoria

if(browniano) then
deallocate(nr)
end if
deallocate(X)
deallocate(U) 
deallocate(W)
deallocate(Urel)
deallocate(Xrel)
deallocate(FORCAS)
deallocate(Ft)
deallocate(TORQUES)
deallocate(Tt)
deallocate(Di)
deallocate(aux1)
deallocate(aux2)
deallocate(aux3)
deallocate(T)
deallocate(magtempo)
deallocate(d)


end subroutine principal
