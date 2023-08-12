subroutine entrada

use variaveis

!******************************** Adquirindo informacoes da simulacao **********************!


! Programa de leitura de arquivos de entrada para processamento com o programa two_particles
! Autor: Rafael Gabler Gontijo


505   FORMAT(1X,A40,1X,E10.4E2)
507   FORMAT(1X,A40,1X,I6)
508   FORMAT(1X,A40,L10)
509   FORMAT(1X,A40,1X,I1)

open(6,file='entrada.dat')
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,508) texto,inercia
      READ (6,508) texto,torque
      READ (6,508) texto,browniano
      READ (6,508) texto,magnetico
      READ (6,508) texto,estatistica
      READ (6,508) texto,dipolo
      READ (6,508) texto,pos_inicial
      READ (6,508) texto,externo
      READ (6,508) texto,mag_tempo
      READ (6,508) texto,diagrama
      READ (6,508) texto,arquivo
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ(6,505) texto,St
      READ(6,505) texto,beta
      READ(6,505) texto,lambda
      READ(6,505) texto,gama
      READ(6,505) texto,Pe
      READ(6,505) texto,Str
      READ(6,505) texto,Per
      READ(6,505) texto,Stm
      READ(6,505) texto,Stmr
      READ(6,505) texto,alfa_campo
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ(6,505) texto,X_11
      READ(6,505) texto,X_12
      READ(6,505) texto,X_13
      READ(6,505) texto,X_21
      READ(6,505) texto,X_22
      READ(6,505) texto,X_23
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ(6,505) texto,Di_11
      READ(6,505) texto,Di_12
      READ(6,505) texto,Di_13
      READ(6,505) texto,Di_21
      READ(6,505) texto,Di_22
      READ(6,505) texto,Di_23
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ (6,'(A)') texto
      READ(6,505) texto,tempo
      READ(6,507) texto,n2
      READ(6,507) texto,rea
      READ(6,509) texto,quad
n3=n2
 close(6)

end subroutine entrada


