# Declare compiler, required flags, and dependencies
CC = g++
CXXFLAGS = -std=c++17 -pedantic -O2 -g
LDLIBS = -lboost_program_options
HDRS = Euler.h FluxUpdate.h NumericalFlux.h ConvertVars.h EquationofStates.h ExactRPSolver.h TypeDefs.h GhostFluidMethod.h interpolation.hpp
TARGET = Solve
# SRCS = EulerSolver.cpp Euler.cpp FluxUpdate.cpp NumericalFlux.cpp ConvertVars.cpp EquationofStates.cpp
OBJS = EulerSolver.o Euler.o FluxUpdate.o NumericalFlux.o ConvertVars.o EquationofStates.o ExactRPSolver.o TypeDefs.o GhostFluidMethod.o

# Define make target
default: $(TARGET)
all: $(TARGET)

# Generate object files
%.cpp:

%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

# Link object file and produce executable
$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: testToro1 testToro2 testToro3 testFedkiwTestB testWangTestB testWangTestB2Dx testWangTestB2Dy testWangTestB2D35 testHeliumBubble clean cleanPlots cleanData

# Test cases for the code
# N = 100, 200, 400   0.01, 0.005, 0.0025
testToro1:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=1 --name="Data/Toro1N400.dat"

testToro2:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=2 --name="Data/Toro2N400.dat"

testToro3:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=3 --name="Data/Toro3N400.dat"

testFedkiwTestB:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=4 --name="Data/FedkiwTestBN400.dat"

testWangTestB:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=5 --name="Data/WangTestBN400Shock.dat"

# 2D one resolution
testWangTestB2Dx:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=400 --Ny=400 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=5 --name="Data/WangTestB2Dx.dat"

testWangTestB2Dy:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=400 --Ny=400 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=6 --name="Data/WangTestB2Dy.dat"

testWangTestB2D35:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=100 --Ny=100 --nVar=4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=7 --name="Data/WangTestB2D35.dat"

testHeliumBubble:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=325 --yBottomDomain=-44.5 --yTopDomain=44.5 --Nx=325 --Ny=89 --nVar=4 --Cfl=0.8 --NxGhost=2 --NyGhost=2 --testNumber=8 --name="Data/HeliumBubble.dat"

# Clean commands
clean:
	rm -f $(TARGET) *.o

cleanPlots:
	rm -f Plots/*

cleanData:
	rm -f Data/*