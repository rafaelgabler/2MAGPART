subroutine saida
use variaveis
use funcoes

write(*,*) 'Iniciando o modulo de processamento de dados...'
write(*,*) ''
! Alocando algumas variaveis

npast=npast/n2
N=2

allocate(U(rea,N,3)) 
allocate(V(rea,N,3)) 
allocate(X(rea,N,3)) 
allocate(Xrelaux(npast,rea,3)) 
allocate(Xrelmedio(npast,3)) 
allocate(Di(rea,N,3))
allocate(flut(N,npast,rea,3))
allocate(velmedia(npast,rea,3))
allocate(vmedia(npast,3))
allocate(errovmedia(npast,3))
allocate(auxcor(rea,N,npast,3))
allocate(funcaor(npast,3))
allocate(dif(npast,3))
allocate(difaux(rea,npast,3))
allocate(errodif(npast,3))
allocate(autocor(rea,N,npast,3))
allocate(errocor(npast,3))
allocate(aux_erro_vel(npast,rea,3))
allocate(aux_erro_var(npast,rea,3))
allocate(aux_erro_cor(rea,npast,3))
allocate(variancia(npast,rea,3))
allocate(var(npast,3))
allocate(errovar(npast,3))
allocate(rrea(rea,npast,3))

nativa=0
particulas_ativas=0
509 FORMAT(F30.4,F30.4,F30.4)

!*************** TRACANDO A MEDIA DAS TRAJETORIAS RELATIVAS ***************************!

do i=1,rea

write(rea_char, '(I3)') i
open (rea+rea+i,file='trajetoria'//rea_char//'.plt',STATUS='OLD')
read (rea+rea+i,'(A)') linha1 
do k=1,npast-1
read (rea+rea+i,509) Xrelaux(k,i,1),Xrelaux(k,i,2),Xrelaux(k,i,3)
end do

end do



do k=1,npast-1
Xrelmedio(k,1)=sum(Xrelaux(k,:,1))/rea
Xrelmedio(k,2)=sum(Xrelaux(k,:,2))/rea
Xrelmedio(k,3)=sum(Xrelaux(k,:,3))/rea
end do

! Escrevendo em um arquivo de saida a trajetoria relativa media

open(6*rea,file='trajetoria_final.plt')
write(6*rea,*)'Variables="X1","X2","X3"'
do k=1,npast-1
write(6*rea,*) Xrelmedio(k,1),Xrelmedio(k,2),Xrelmedio(k,3)
end do


if(arquivo) then
!*************************************** VELOCIDADE MEDIA **********************************!
do i=1,rea
write(rea_char, '(I3)') i
open (i,file='posicao'//rea_char//'.plt',STATUS='OLD')
end do

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do

do j=1,rea
read (j,'(A)') linha1 
read (rea+j,'(A)') linha1 
do k=1,npast-1
read (j,'(A)') linha2 
read(rea+j,'(A)') linha2
do i=1,N
read(j,509)X(j,i,1),X(j,i,2),X(j,i,3)
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)
end do

! Determinando a velocidade media em cada time-step e em cada realizacao (para cada direcao)
! considerando somente as particulas da parte ativa da suspensao

 call media_direcao(U,N,rea,velmedia(k,j,1),1)
 call media_direcao(U,N,rea,velmedia(k,j,2),2)
 call media_direcao(U,N,rea,velmedia(k,j,3),3)

end do
end do

! Calculando agora a velocidade media em cada time-step, tirando uma media em cima das
! realizacoes do valor de velmedia(k,j,direcao)

do k=1,npast-1
vmedia(k,1)=sum(velmedia(k,:,1))/rea
vmedia(k,2)=sum(velmedia(k,:,2))/rea
vmedia(k,3)=sum(velmedia(k,:,3))/rea
end do


! Calculando o erro-bar da velocidade media

do k=1,npast-1
do j=1,rea
aux_erro_vel(k,j,1)=(velmedia(k,j,1)-vmedia(k,1))**2.0
aux_erro_vel(k,j,2)=(velmedia(k,j,2)-vmedia(k,2))**2.0
aux_erro_vel(k,j,3)=(velmedia(k,j,3)-vmedia(k,3))**2.0
end do
end do

do k=1,npast-1
errovmedia(k,1)=((1.0/(rea))*sum(aux_erro_vel(k,:,1)))**0.5
errovmedia(k,2)=((1.0/(rea))*sum(aux_erro_vel(k,:,2)))**0.5
errovmedia(k,3)=((1.0/(rea))*sum(aux_erro_vel(k,:,3)))**0.5
end do

! Fechando os arquivos
do j=1,2*rea
 close(j)
end do

! Escrevendo em um arquivo de saida a velocidade media e seu errobar

510 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)

i=2*rea+1
open (i,file='velocidade_media.plt')
write(i,*) 'Variables="U","V","W","UMEDIA","DU","DV","DW","T"'

do k=1,npast-1
write(i,510)vmedia(k,1),vmedia(k,2),vmedia(k,3),   &
(((vmedia(k,1)**2.0)+(vmedia(k,2)**2.0)+(vmedia(k,3)**2.0))**0.5),   &
errovmedia(k,1),errovmedia(k,2),errovmedia(k,3),k*dt*n2
end do

write(*,*) 'Analise estatistica em cima da velocidade media - OK'

!*******************************************************************************************!

!********************************************** VARIANCIA **********************************!
do i=1,rea
write(rea_char, '(I3)') i
open (i,file='posicao'//rea_char//'.plt',STATUS='OLD')
end do

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do


do j=1,rea
read (j,'(A)') linha1 
read (rea+j,'(A)') linha1 
do k=1,npast-1
read (j,'(A)') linha2 
read(rea+j,'(A)') linha2
do i=1,N
read(j,509)X(j,i,1),X(j,i,2),X(j,i,3)
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

! Determinando a variancia em cada time-step e em cada realizacao (para cada direcao)
! considerando somente as particulas da parte ativa da suspensao

V(j,i,1) = (U(j,i,1)-velmedia(k,j,1))**2.0
V(j,i,2) = (U(j,i,2)-velmedia(k,j,2))**2.0
V(j,i,3) = (U(j,i,3)-velmedia(k,j,3))**2.0
end do

 call media_direcao(V,N,rea,variancia(k,j,1),1)
 call media_direcao(V,N,rea,variancia(k,j,2),2)
 call media_direcao(V,N,rea,variancia(k,j,3),3)

V=0.0

end do
end do


! Calculando agora a variancia em cada time-step, tirando uma media em cima das
! realizacoes do valor de velmedia(k,j,direcao)

do k=1,npast-1
var(k,1)=sum(variancia(k,:,1))/rea
var(k,2)=sum(variancia(k,:,2))/rea
var(k,3)=sum(variancia(k,:,3))/rea
end do


! Calculando o erro-bar da variancia

do k=1,npast-1
do j=1,rea
aux_erro_var(k,j,1)=(variancia(k,j,1)-var(k,1))**2.0
aux_erro_var(k,j,2)=(variancia(k,j,2)-var(k,2))**2.0
aux_erro_var(k,j,3)=(variancia(k,j,3)-var(k,3))**2.0
end do
end do

do k=1,npast-1
errovar(k,1)=((1.0/(rea))*sum(aux_erro_var(k,:,1)))**0.5
errovar(k,2)=((1.0/(rea))*sum(aux_erro_var(k,:,2)))**0.5
errovar(k,3)=((1.0/(rea))*sum(aux_erro_var(k,:,3)))**0.5
end do

! Fechando os arquivos
do j=1,2*rea
 close(j)
end do


! Escrevendo em um arquivo de saida a variancia e seu errobar

i=2*rea+2
open (i,file='variancia.plt')
write(i,*) 'Variables="Var1","Var2","Var3","VarMedia","DVar1","DVar2","DVar3","T"'


do k=1,npast-1
write(i,510)var(k,1),var(k,2),var(k,3),   &
(((var(k,1)**2.0)+(var(k,2)**2.0)+(var(k,3)**2.0))**0.5),   &
errovar(k,1),errovar(k,2),errovar(k,3),k*dt*n2
end do


write(*,*) 'Analise estatistica em cima da variancia - OK'
!*******************************************************************************************!


!************************************* MAGNETIZACAO DE EQUILIBRIO **************************!

write(*,*)''
write(*,*) 'Determinando a magnetizacao de equilibrio'
write(*,*)''

do i=1,rea
write(rea_char, '(I3)') i
open (rea+rea+i,file='dipolo'//rea_char//'.plt',STATUS='OLD')
end do

do j=1,rea
read (rea+rea+j,'(A)') linha3 
do i=1,N
read(rea+rea+j,509)Di(j,i,1),Di(j,i,2),Di(j,i,3)
end do
end do


! Determinando a magnetizacao de equilibrio da suspensao
if(magnetizacao) then
 call media_direcao(Di,N,rea,magnetizacao_media,3)
write(*,*) 'A magnetizacao de equilibrio da suspensao e:',magnetizacao_media,'(Mo/Ms)'
end if

! Fechando os arquivos dipolo___.plt

do j=1,rea
 close((2*rea)+j)
end do

!*******************************************************************************************!

! Gerando um arquivo relatorio para a analise short-time

difusao=((dif(npast-3,1)**2.0)+(dif(npast-3,2)**2.0)+(dif(npast-3,3)**2.0))**0.5

UMEDIA(1)=sum(vmedia(:,1))/(npast-1.0)
UMEDIA(2)=sum(vmedia(:,2))/(npast-1.0)
UMEDIA(3)=sum(vmedia(:,3))/(npast-1.0)

SIGMA(1)=sum(var(:,1))/(npast-1.0)
SIGMA(2)=sum(var(:,2))/(npast-1.0)
SIGMA(3)=sum(var(:,3))/(npast-1.0)

umsobreene=1.0/N

open (rea+rea+rea+1,file='relatorio_short_time.plt')
write(rea+rea+rea+1,*)'variables="1/N","N","U","V","W","UMEDIA","SIGMA_1","SIGMA_2","SIGMA_3","SIGMA_MEDIA","M"'
write(rea+rea+rea+1,*) umsobreene,N,UMEDIA(1),UMEDIA(2),UMEDIA(3),   &
(((UMEDIA(1)**2.0)+(UMEDIA(2)**2.0)+(UMEDIA(3)**2.0))**0.5), SIGMA(1),SIGMA(2),SIGMA(3),   &
(((SIGMA(1)**2.0)+(SIGMA(2)**2.0)+(SIGMA(3)**2.0))**0.5), magnetizacao_media

write(*,*) 'Geracao de arquivo relatorio para a analise short-time - OK'

write(*,*) ''

stop

!***************************************** FUNCAO AUTOCORRELACAO ***************************!

! Abrindo novamente os arquivos

do i=1,rea
write(rea_char, '(I3)') i
open (i,file='posicao'//rea_char//'.plt',STATUS='OLD')
end do

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do

! Calculando a flutuacao em cada time-step e realizacao para cada particula

do j=1,rea
read (j,'(A)') linha1 
read (rea+j,'(A)') linha1 
do k=1,npast-1
read (j,'(A)') linha2 
read(rea+j,'(A)') linha2

do i=1,N
read(j,509)X(j,i,1),X(j,i,2),X(j,i,3)
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

flut(i,k,j,1)= U(j,i,1)-velmedia(k,j,1)
flut(i,k,j,2)= U(j,i,2)-velmedia(k,j,2)
flut(i,k,j,3)= U(j,i,3)-velmedia(k,j,3)
end do

end do
end do


! Montando a funcao auto-correlacao das flutuacoes de velocidade para cada particula e cada realizacao

write(*,*) 'Montando a funcao auto-correlacao das flutuacoes de velocidade para cada particula e cada realizacao'

do j=1,rea

 !$OMP PARALLEL DO
do q=1,N
do i=1,npast-1

 !$OMP PARALLEL DO
do k=1,npast-i
auxcor(j,q,k,1)=(flut(q,k,j,1)*flut(q,k+i-1,j,1))/(flut(q,k,j,1)**2.0)
auxcor(j,q,k,2)=(flut(q,k,j,2)*flut(q,k+i-1,j,2))/(flut(q,k,j,2)**2.0)
auxcor(j,q,k,3)=(flut(q,k,j,3)*flut(q,k+i-1,j,3))/(flut(q,k,j,3)**2.0)
end do
 !$OMP END PARALLEL DO

autocor(j,q,i,1)=sum(auxcor(j,q,:,1))/(npast-i)
autocor(j,q,i,2)=sum(auxcor(j,q,:,2))/(npast-i)
autocor(j,q,i,3)=sum(auxcor(j,q,:,3))/(npast-i)
auxcor=0.0

end do

end do
 !$OMP END PARALLEL DO
end do

! Tirando a media da funcao auto-correlacao em cima das particulas

write(*,*) 'Tirando a media da funcao auto-correlacao em cima das particulas'

do j=1,rea
do k=1,npast-1
rrea(j,k,1)= sum(autocor(j,:,k,1))/N
rrea(j,k,2)= sum(autocor(j,:,k,2))/N
rrea(j,k,3)= sum(autocor(j,:,k,3))/N
end do
end do


! Tirando uma media em cima das realizacoes

write(*,*) 'Tirando uma media em cima das realizacoes'

do k=1,npast-1
funcaor(k,1) = sum(rrea(:,k,1))/rea
funcaor(k,2) = sum(rrea(:,k,2))/rea
funcaor(k,3) = sum(rrea(:,k,3))/rea
end do


! Calculando o erro-bar da funcao autocorrelacao


do j=1,rea
do k=1,npast-2
aux_erro_cor(j,k,1)=(rrea(j,k,1)-funcaor(k,1))**2.0
aux_erro_cor(j,k,2)=(rrea(j,k,2)-funcaor(k,2))**2.0
aux_erro_cor(j,k,3)=(rrea(j,k,3)-funcaor(k,3))**2.0
end do
end do

! Calculando o erro-bar da funcao autocorrelacao
do k=1,npast-2
errocor(k,1)=((1.0/(rea))*sum(aux_erro_cor(:,k,1)))**0.5
errocor(k,2)=((1.0/(rea))*sum(aux_erro_cor(:,k,2)))**0.5
errocor(k,3)=((1.0/(rea))*sum(aux_erro_cor(:,k,3)))**0.5
end do


! Colocando um filtro numerico para consertar um bug no valor de R(0)

funcaor(1,1)=1.0
funcaor(1,2)=1.0
funcaor(1,3)=1.0

! Fechando os arquivos
do j=1,2*rea
 close(j)
end do

! Escrevendo em um arquivo a funcao autocorrelacao das flutuacoes e seu errobar

i=2*rea+3
open (i,file='autocorrelacao.plt')
write(i,*) 'Variables="R1","R2","R3","RMEDIA","DR1","DR2","DR3","T"'

do k=1,npast-1
write(i,510)funcaor(k,1),funcaor(k,2),funcaor(k,3),   &
(((funcaor(k,1)**2.0)+(funcaor(k,2)**2.0)+(funcaor(k,3)**2.0))**0.5),   &
errocor(k,1),errocor(k,2),errocor(k,3),k*dt*n2
end do

write(*,*) 'Analise estatistica em cima da funcao autocorrelacao - OK'

!*******************************************************************************************!

!**************************************** COEFICIENTE DE DIFUSAO ***************************!

do q=1,rea
difaux(q,1,1)=0.0
difaux(q,1,2)=0.0
difaux(q,1,3)=0.0
end do

do q=1,rea
do k=2,npast-3
difaux(q,k,1)=difaux(q,k-1,1)+((rrea(q,k,1)+rrea(q,k-1,1))*(n2*dt)/2.0)
difaux(q,k,2)=difaux(q,k-1,2)+((rrea(q,k,2)+rrea(q,k-1,2))*(n2*dt)/2.0)
difaux(q,k,3)=difaux(q,k-1,3)+((rrea(q,k,3)+rrea(q,k-1,3))*(n2*dt)/2.0)
end do
end do


! Determinando o coeficiente de difusao

dif(1,1)=0.0
dif(1,2)=0.0
dif(1,3)=0.0

do k=2,npast-3
dif(k,1)=sum(difaux(:,k,1))/rea
dif(k,2)=sum(difaux(:,k,2))/rea
dif(k,3)=sum(difaux(:,k,3))/rea
end do


! Iniciando o calculo do erro do coeficiente de difusao

do k=1,npast-2
do q=1,rea
difaux(q,k,1)=(difaux(q,k,1)-dif(k,1))**2.0
difaux(q,k,2)=(difaux(q,k,2)-dif(k,2))**2.0
difaux(q,k,3)=(difaux(q,k,3)-dif(k,3))**2.0
end do
end do


! Calculando o erro-bar do coeficiente de difusao
do k=1,npast-2
errodif(k,1)=((1.0/(rea))*sum(difaux(:,k,1)))**0.5
errodif(k,2)=((1.0/(rea))*sum(difaux(:,k,2)))**0.5
errodif(k,3)=((1.0/(rea))*sum(difaux(:,k,3)))**0.5
end do


! Escrevendo em um arquivo o coeficiente de difusao e seu erro-bar

513 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
i=2*rea+4
open (i,file='difusao.plt')
write(i,*) 'Variables="D1","D2","D3","DMEDIO","DD1","DD2","DD3","T"'

do k=1,npast-1
write(i,510)dif(k,1),dif(k,2),dif(k,3),   &
(((dif(k,1)**2.0)+(dif(k,2)**2.0)+(dif(k,3)**2.0))**0.5),   &
errodif(k,1),errodif(k,2),errodif(k,3),k*dt*n2
end do


write(*,*) 'Analise estatistica em cima do coeficiente de difusao - OK'

!*******************************************************************************************!


! Contando o numero de dimeros, trimeros e agregados com mais particulas.
! Essa subrotina ja traca a fdp do numero de particulas por agregado de todas
! as realizacoes. 

!do  j=1,rea
! call analise_de_agregados(X,N,j,rea)
!end do

end if

write(*,*) 'Geracao dos arquivos de saida - OK'

! Desalocando as variaveis para liberar espaco na memoria

write(*,*) ''
write(*,*) 'Desalocando as matrizes criadas no modulo de pos-processamento...'
write(*,*) ''

deallocate(U) 
write(*,*) 'Desalocando matriz 1 - OK'
deallocate(V) 
write(*,*) 'Desalocando matriz 2 - OK'
deallocate(X)
write(*,*) 'Desalocando matriz 3 - OK'
deallocate(Xrelaux)
write(*,*) 'Desalocando matriz 4 - OK'
deallocate(Xrelmedio)
write(*,*) 'Desalocando matriz 5 - OK'
deallocate(Di)
write(*,*) 'Desalocando matriz 6 - OK'
deallocate(flut)
write(*,*) 'Desalocando matriz 7 - OK'
deallocate(flut)
write(*,*) 'Desalocando matriz 8 - OK'
deallocate(velmedia)
write(*,*) 'Desalocando matriz 9 - OK'
deallocate(vmedia)
write(*,*) 'Desalocando matriz 10 - OK'
deallocate(errovmedia)
write(*,*) 'Desalocando matriz 11 - OK'
deallocate(auxcor)
write(*,*) 'Desalocando matriz 12 - OK'
deallocate(funcaor)
write(*,*) 'Desalocando matriz 13 - OK'
deallocate(dif)
write(*,*) 'Desalocando matriz 14 - OK'
deallocate(difaux)
write(*,*) 'Desalocando matriz 15 - OK'
deallocate(errodif)
write(*,*) 'Desalocando matriz 16- OK'
deallocate(autocor)
write(*,*) 'Desalocando matriz 17- OK'
deallocate(errocor)
write(*,*) 'Desalocando matriz 18- OK'
deallocate(aux_erro_vel)
write(*,*) 'Desalocando matriz 19- OK'
deallocate(aux_erro_var)
write(*,*) 'Desalocando matriz 20- OK'
deallocate(aux_erro_cor)
write(*,*) 'Desalocando matriz 21- OK'
deallocate(variancia)
write(*,*) 'Desalocando matriz 22- OK'
deallocate(var)
write(*,*) 'Desalocando matriz 23- OK'
deallocate(errovar)
write(*,*) 'Desalocando matriz 24- OK'
deallocate(rrea)
write(*,*) 'Desalocando matriz 25- OK'

end subroutine saida
