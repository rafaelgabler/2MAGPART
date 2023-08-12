program deux_particles
use variaveis
! Titulo do programa e apresentacao do mesmo

print *,'******************************************************************************'
print *,'*                     Programa Two Particles - v.2.1 - 2011                  *'
print *,'*                                                                            *'
print *,'* Autor: Msc. Rafael Gabler Gontijo                                          *'
print *,'*                                                                            *'
print *,'*                                                                            *'
print *,'******************************************************************************'

print *,''
print *,'Interacao entre duas particulas magneticas'
print *,''
print *,''

 call cpu_time(ti)

 call entrada
 call zeratudo
 call principal
if(estatistica)then
 call saida
end if
 call cpu_time(tf)

 tpros=tf-ti

write(*,*) ''
write(*,*) 'O tempo total de processamento foi de:',tpros,'segundos'
write(*,*) ''


end program deux_particles
