subroutine zeratudo

use variaveis
! Aloca as variaveis e zera tudo para o inicio da simulacao


! Determinando o passo de tempo e o numero de passos de tempo para poder alocar alguns vetores

dt=min(St/10.0,Pe/10.0,0.001)

npast = tempo/dt

! Quantidade de numeros randomicos necessarios em cada timestep

if(browniano) then
nnr=3*2*rea
end if

N=2

if(browniano) then
allocate(nr(nnr))
end if

allocate(X(rea,2,3))
allocate(colisao(rea))
allocate(posicao_inicial(rea,2))
allocate(posicao_final(rea,2))
allocate(desloc_self(rea,2))
allocate(desloc_col(rea,2))
allocate(U(rea,2,3)) 
allocate(W(rea,2,3))
allocate(Urel(rea,3))
allocate(Xrel(rea,3))
allocate(FORCAS(4,rea,2,3))
allocate(Ft(rea,2,3))
allocate(TORQUES(3,rea,2,3))
allocate(Tt(rea,2,3))
allocate(Di(rea,2,3))
allocate(aux1(rea,2))
allocate(aux2(rea,2))
allocate(aux3(rea,2))
allocate(deltaquadrado(rea,2))
allocate(T(npast))
allocate(magtempo(npast))
allocate(d(rea,3))
allocate(modd(rea))
allocate(distancia(rea,npast))
allocate(contagem(rea))



! Zerando tudo

razao1=0.0
razao2=0.0
razao3=0.0
razao4=0.0
razao5=0.0
parteativa=0.0
posicao_inicial=0.0
posicao_final=0.0
desloc_self=0.0
desloc_col=0.0
Dself1=0.0
Delf2=0.0
Dcol1=0.0
Dcol2=0.0
Dglobal1=0.0
Dglobal2=0.0
difusao_frc1=0.0
difusao_frc2=0.0
distancia=0.0
deltaquadrado=0.0
contagem=0.0
auxflut=0.0
difusaoaux=0.0
errocor=0.0
autocor=0.0
auxcor=0.0
errodif=0.0
difaux=0.0
dif=0.0
modd=0.0
funcaor=0.0
errofmedia=0.0
errovmedia=0.0
fmedia=0.0
vmedia=0.0
velmedia=0.0
flutmag=0.0
flut=0.0
ULINHA=0.0
d=0.0
aux_erro_vel=0.0
aux_erro_var=0.0
aux_erro_cor=0.0
variancia=0.0
var=0.0
errovar=0.0
rrea=0.0
U=0.0
W=0.0
X=0.0
V=0.0
Urel=0.0
Xrelaux=0.0
Xrelmedio=0.0
Xrel=0.0
Di=0.0
Ft=0.0
Tt=0.0
FORCAS=0.0
TORQUES=0.0
aux1=0.0
aux2=0.0
aux3=0.0
T=0.0
colisao=0
Urel=0.0
magtempo=0.0
UMEDIA=0.0
SIGMA=0.0
agregativa=0
difusiva=0
r=0.0
s=0.0
qui=0.0
f0=0.0
f1=0.0
f2=0.0
f3=0.0
f4=0.0
f5=0.0
f6=0.0
f7=0.0
f8=0.0
f9=0.0
f10=0.0
f11=0.0
g1=0.0
g2=0.0
g3=0.0
m1=0.0
X11a=0.0
X12a=0.0 
X21a=0.0 
X22a=0.0
Y11a=0.0  
Y12a=0.0 
Y21a=0.0 
Y22a=0.0
A11x=0.0 
A12x=0.0
A22x=0.0 
A21x=0.0
A11y=0.0 
A12y=0.0 
A22y=0.0 
A21y=0.0
Wx=0.0
modip=0.0
teste1=0.0
magnetizacao_media=0.0
nr1=0.0
nr2=0.0
nr3=0.0
dist=0.0
nr=0.0
termo1=0.0
termo2=0.0
termo3=0.0
termo4=0.0
modrij=0.0			
rij(3)=0.0			
modrand=0.0
theta=0.0 
modomega=0.0 
thetaux=0.0
mod1=0.0
mod2=0.0

end subroutine zeratudo
