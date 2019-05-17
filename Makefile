CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = Solver.o Particle.o ParticleToy.o RodConstraint.o GravityForce.o SpringForce.o CircularWireConstraint.o imageio.o

project1: $(OBJS)
	$(CXX) -o bin/$@ $^ -Llib -lfreeglut -lglu32 -lopengl32 -lpng12
	
clean:
	del $(OBJS) bin\project1.exe
