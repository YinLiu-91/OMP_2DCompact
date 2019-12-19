
include Makefile.in

OBJECTS  = OMP_2DCompact.o AbstractSingleBlockMesh.cpp Utils.o Pade6.o Compact10Filter.o WideGaussianFilter.o Options.o CurvilinearCSolver_Core.o CurvilinearCSolver_BC.o SpongeBC.o CurvilinearInterpolator.o LADKawai.o
#POSTPROOBJ = Utils.o Derivatives.o PostProcess.o VisitWriter.o

all: OMP_2DCOMPACT.exe

OMP_2DCompact.o: OMP_2DCompact.cpp Macros.hpp Options.hpp TimeStepping.hpp Domain.hpp BC.hpp AbstractCSolver.hpp AbstractRK.hpp TVDRK3.hpp LSLDDRK4.hpp AbstractSingleBlockMesh.hpp AlgebraicSingleBlockMesh.hpp AbstractDerivatives.hpp Pade6.hpp CurvilinearCSolver.hpp CurvilinearInterpolator.hpp 
	$(CC) $(CFLAGS) -c $< 

CurvilinearCSolver_Core.o: CurvilinearCSolver_Core.cpp CurvilinearCSolver.hpp Macros.hpp Options.hpp Utils.hpp AbstractCSolver.hpp Pade6.hpp Compact10Filter.hpp
	$(CC) $(CFLAGS) -c $< 

CurvilinearCSolver_BC.o: CurvilinearCSolver_BC.cpp CurvilinearCSolver.hpp Macros.hpp Options.hpp Utils.hpp AbstractCSolver.hpp Pade6.hpp Compact10Filter.hpp
	$(CC) $(CFLAGS) -c $< 

SpongeBC.o: SpongeBC.cpp SpongeBC.hpp Macros.hpp Options.hpp Domain.hpp IdealGas.hpp BC.hpp Utils.hpp AbstractSingleBlockMesh.hpp 
	$(CC) $(CFLAGS) -c $< 

AbstractSingleBlockMesh.o: AbstractSingleBlockMesh.cpp AbstractSingleBlockMesh.hpp AbstractCSolver.hpp Domain.cpp Pade6.hpp Adt.hpp Options.hpp
	$(CC) $(CFLAGS) -c $< 

Pade6.o: Pade6.cpp Pade6.hpp AbstractDerivatives.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp Options.hpp 
	$(CC) $(CFLAGS) -c $< 

Compact10Filter.o: Compact10Filter.cpp Compact10Filter.hpp AbstractFilter.hpp AbstractDerivatives.hpp Options.hpp Macros.hpp Domain.hpp Utils.hpp BC.hpp
	$(CC) $(CFLAGS) -c $< 

WideGaussianFilter.o: WideGaussianFilter.cpp WideGaussianFilter.hpp AbstractFilter.hpp AbstractDerivatives.hpp Options.hpp Macros.hpp Domain.hpp Utils.hpp BC.hpp
	$(CC) $(CFLAGS) -c $< 

LADKawaii.o: LADKawai.cpp LADKawai.hpp AbstractLAD.hpp AbstractDerivatives.hpp AbstractFilter.hpp Options.hpp Macros.hpp Utils.hpp WideGaussianFilter.hpp Pade6.hpp 
	$(CC) $(CFLAGS) -c $< 

Utils.o: Utils.cpp Domain.hpp Utils.hpp 
	$(CC) $(CFLAGS) -c $< 

Options.o: Options.cpp Options.hpp
	$(CC) $(CFLAGS) -c $< 

CurvilinearInterpolator.o: CurvilinearInterpolator.cpp CurvilinearInterpolator.hpp Macros.hpp Utils.hpp AbstractCSolver.hpp AlgebraicSingleBlockMesh.hpp BC.hpp Domain.hpp Options.hpp
	$(CC) $(CFLAGS) -c $< 



OMP_2DCOMPACT.exe:  $(OBJECTS)
	$(CC) $(CFLAGS) -I$(INC) $(OBJECTS) -o $@ -L$(LIB) $(LIBF)

clean: 
	rm -rf   *.o 


