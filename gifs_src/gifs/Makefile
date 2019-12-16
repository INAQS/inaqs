
CPP=icpc
CPPFLAGS= -std=c++11 -O2 -Werror -Wall

INCLUDE= -I../../include

OBJ= qm.o properties.o gifs.o hallo.o

libgifs: $(OBJ)
	$(CPP) -shared -o ../../lib/libgifs.so $(OBJ)


%.o: %.cpp
	$(CPP) $(INCLUDE) $(CPPFLAGS) -fPIC -c $<
