
include Makefile.in

OBJECTS  = OMP_2DCompact.o AbstractSingleBlockMesh.cpp Utils.o Pade6.o Compact10Filter.o Options.o  
#POSTPROOBJ = Utils.o Derivatives.o PostProcess.o VisitWriter.o

all: OMP_2DCOMPACT.exe

OMP_2DCompact.o: OMP_2DCompact.cpp Macros.hpp Options.hpp TimeStepping.hpp Domain.hpp BC.hpp AbstractCSolver.hpp AbstractRK.hpp TVDRK3.hpp AbstractSingleBlockMesh.hpp AlgebraicSingleBlockMesh.hpp AbstractDerivatives.hpp Pade6.hpp
	$(CC) $(CFLAGS) -c $< 

AbstractSingleBlockMesh.o: AbstractSingleBlockMesh.cpp AbstractSingleBlockMesh.hpp AbstractCSolver.hpp Domain.cpp Pade6.hpp Adt.hpp Options.hpp

Pade6.o: Pade6.cpp Pade6.hpp AbstractDerivatives.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp Options.hpp 
	$(CC) $(CFLAGS) -c $< 

Compact10Filter.o: Compact10Filter.cpp Compact10Filter.hpp AbstractFilter.hpp AbstractDerivatives.hpp Options.hpp Macros.hpp Domain.hpp Utils.hpp BC.hpp

Utils.o: Utils.cpp Domain.hpp Utils.hpp 
	$(CC) $(CFLAGS) -c $< 

Options.o: Options.cpp Options.hpp
	$(CC) $(CFLAGS) -c $< 

OMP_2DCOMPACT.exe:  $(OBJECTS)
	$(CC) $(CFLAGS) -I$(INC) $(OBJECTS) -o $@ -L$(LIB) $(LIBF)

clean: 
	rm -rf   *.o 


