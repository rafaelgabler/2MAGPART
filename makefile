#--------------------------------------------------------#
# Diretivas de compilação para o programa Deux Particles #
#                                                     	 #
#             Escrito inicialmente por                   #
#               Rafael Gabler Gontijo                    #
#                        em                              #
#                    05/05/2010                          #
#                                                        #
# Objetivo: Compilar uma determinada versão do código    #
#                                                        #
# Versão atual: - Aceita compilação para GNU/Linux ou    #
#               - Compilador: GNU gfortran.              #
#                                                        #
#                                                        #
# Última modificação por: Rafael Gabler Gontijo          #
# Última modificação em: 12/08/2010                      #
#--------------------------------------------------------#
#
# Arquivos fonte padrão (linux) e outros sistemas
#
SRC-LNX = funcoes.f90 variaveis.f90 entrada.f90 zeratudo.f90  \
          principal.f90 saida.f90 two_particles.f90  \
          
 OBJ-LNX = $(SRC-LNX:.f90=.o)

#
# Definição do compilador específico e das flags

#ENGINE-LNX = gfortran
ENGINE-LNX = ifort
FLAGS-LNX = -m64 -O2 -openmp

#
# Objetivo principal: geração de executável em linux

#cavmag2d.ex : $(SRC-LNX) 
#	$(ENGINE-LNX) -o two_particles.ex $(SRC-LNX)

deux_particles.ex : $(SRC-LNX) 
	$(ENGINE-LNX) $(FLAGS-LNX) -o two_particles.ex $(SRC-LNX)


#
# Limpeza
clean :
	rm two_particles.ex
