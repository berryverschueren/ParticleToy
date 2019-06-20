CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = ParticleToy.o FluidSystem.o RigidBody.o

project1: $(OBJS)
	$(CXX) -o bin/$@ $^ -Llib -lfreeglut -lglu32 -lopengl32 -lpng12
	
clean:
	del $(OBJS) bin\project1.exe
