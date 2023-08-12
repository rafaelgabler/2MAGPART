module variaveis

! Definindo variaveis
integer i,j,k,q, npast
integer adimen, quad
real dt, razao1, razao2, razao3, razao4, razao5
real, allocatable :: parteativa(:)
real, allocatable :: posicao_inicial(:,:)
real, allocatable :: posicao_final(:,:)
real, allocatable :: desloc_self(:,:)
real, allocatable :: desloc_col(:,:)
real Dself1,Dself2,Dcol1,Dcol2
real Dglobal1, Dglobal2
real difusao_frc1, difusao_frc2
real, allocatable :: distancia(:,:)
real, allocatable :: deltaquadrado(:,:)
real, allocatable :: contagem(:)
real, allocatable :: auxflut(:,:,:)
real, allocatable :: difusaoaux(:,:)
real, allocatable :: errocor(:,:)
real, allocatable :: autocor(:,:,:,:)
real, allocatable :: auxcor(:,:,:,:)
real, allocatable :: errodif(:,:)
real, allocatable :: difaux(:,:,:)
real, allocatable :: dif(:,:)
real, allocatable :: modd(:)
real, allocatable :: funcaor(:,:)
real, allocatable :: errofmedia(:,:)
real, allocatable :: errovmedia(:,:)
real, allocatable :: fmedia(:,:)
real, allocatable :: vmedia(:,:)
real, allocatable :: velmedia(:,:,:)
real, allocatable :: flutmag(:,:)
real, allocatable :: flut(:,:,:,:)
real, allocatable :: ULINHA(:,:,:)
real, allocatable :: d(:,:)
real, allocatable :: aux_erro_vel(:,:,:)
real, allocatable :: aux_erro_var(:,:,:)
real, allocatable :: aux_erro_cor(:,:,:)
real, allocatable :: variancia(:,:,:)
real, allocatable :: var(:,:)
real, allocatable :: errovar(:,:)
real, allocatable :: rrea(:,:,:)
real, allocatable :: U(:,:,:)
real, allocatable :: W(:,:,:)
real, allocatable :: X(:,:,:)
real, allocatable :: V(:,:,:)
real, allocatable :: Urel(:,:)
real, allocatable :: Xrelaux(:,:,:)
real, allocatable :: Xrelmedio(:,:)
real, allocatable :: Xrel(:,:)
real, allocatable :: Di(:,:,:)
real, allocatable :: Ft(:,:,:)
real, allocatable :: Tt(:,:,:)
real, allocatable :: FORCAS(:,:,:,:)
real, allocatable :: TORQUES(:,:,:,:)
real, allocatable :: aux1(:,:)
real, allocatable :: aux2(:,:)
real, allocatable :: aux3(:,:)
real, allocatable :: T(:)
integer, allocatable :: colisao(:)
real, allocatable :: magtempo(:)
real UMEDIA(3)
real SIGMA(3)
real St, Stm, Pe, Pa, Pre, Ppsi, beta, lambda, gama, alfa
real alfa_campo
real Str,Per,Stmr
real r
real s, qui
real X_11,X_12,X_13,X_21,X_22,X_23
real Di_11,Di_12,Di_13,Di_21,Di_22,Di_23
real f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
real g1,g2,g3,m1(11)
real X11a, X12a, X21a, X22a
real Y11a, Y12a, Y21a, Y22a
real A11x, A12x, A22x, A21x
real A11y, A12y, A22y, A21y
real Wx
real tempo
real ti,tf,tpros
logical inercia
logical torque
logical browniano
logical magnetico
logical estatistica
logical dipolo
logical pos_inicial
logical externo
logical mag_tempo
logical comentarios
logical diagrama
logical arquivo
real modip
integer nnr
integer rea
integer n2
integer reaok
integer agregativa,difusiva
real n3
integer N
real teste1
real magnetizacao_media
integer teste2
real nr1,nr2,nr3
real dist
 character(3) rea_char
 character(25) linha1
 character(19) linha2
real, allocatable :: nr(:)  	
real termo1,termo2		
real termo3,termo4		
real modrij			
real rij(3)			
real modrand	
real theta, modomega, thetaux,mod1,mod2
 		

!********************* LEGENDA DAS VARIAVEIS UTILIZADAS NESTE CODIGO ***********************!

! VARIAVEIS INTEIRAS

! i,j,k	  - representam variaveis inteiras utilizadas em loops;

! MATRIZES UTILIZADAS

! U(2,3)  - e uma matriz que contem os valores de velocidade das particulas 1 e 2 nas direcoes 1, 2 e 3;
! Urel(3) - e um vetor que contem a velocidade relativa das particulas nas direcoes 1, 2 e 3;
! X(2,3)  - e uma matriz que contem os valores da posicao das particulas 1 e 2 nas direcoes 1, 2 e 3;
! Xrel(3) - e um vetor que contem a trajetoria relativa das particulas nas direcoes 1, 2 e 3;
! Ft(2,3) - matriz que contem a forca total sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;
! Fb(2,3) - matriz que contem a forca browniana sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;
! Fm(2,3) - matriz que contem a forca magnetica sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;
! Fe(2,3) - matriz que contem a forca de repulsao eletrostatica sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;
! Fa(2,3) - matriz que contem a forca atrativa de van der Waals sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;
! Fh(2,3) - matriz que contem a forca de interacao hidrodinamica sobre as particulas 1 e 2 nas direcoes 1, 2 e 3;


! PARAMETROS ADIMENSIONAIS

! St     - Numero de Stokes;
! Stm    - Numero de Stokes magnetico;
! Pe	 - Numero de Peclet;
! Pa	 - Parametro adimensional da forca atrativa de van der Waals;
! Pre	 - Parametro adimensional da forca repulsiva eletrostatica;
! Ppsi	 - Razao entre as funcoes psi2/psi1 utilizadas no calculo da forca de repulsao eletrostatica;
! beta   - Razao entre os raios das particulas a2/a1 (e um parametro de polidispersidade);
! lambda - Razao entre as massas especificas das particulas 2 e 1;
! gama	 - Razao entre delta rho 2 e delta rho 1;

! VARIAVEIS ALOCAVEIS

! T	- vetor que armazena o instante de tempo adimensional de cada iteracao

! OUTRAS VARIAVEIS

! nr - vetor que armazena numeros randomicos necessarios para expressar movimento browniano, 
! posicao inical e momento de dipolo inicial

! termo1, termo2, termo3, termo4 - termos utilizados para implementacao de forcas magneticas;

! dt		    - e o passo de tempo da simulacao;
! tempo		    - e o tempo total da simulacao;
! r		    - distancia entre os centros das particulas;
! qui		    - distancia entre as superficies das particulas;	
! f1,f2,f3, etc...  - funcoes utilizadas para o calculo das interacoes hidrodinamicas (usadas para calculo das funcoes X..a);
! X11a,X12a,etc...  - funcoes resistencia para calculo de interacoes hidrodinamicas de longo alcance;
! A11x,A12x,etc...  - funcoes resistencia para calculo de interacoes hidrodinamicas de curto alcance;
! Wx		    - funcao utilizada no calculo das funcoes resistencia para regimes de curto alcance;

! VARIAVEIS LOGICAS

! inercia	    - verifica que a particula possui ou nao inercia;
! torque	    - verifica se ha necessidade de se resolver a equacao do torque;

!*******************************************************************************************!



end module variaveis
