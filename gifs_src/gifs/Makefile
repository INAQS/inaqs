CPP=g++
CPPFLAGS= -std=c++11 -O2 -Werror -Wall

INCLUDE= -I../../include

OBJ= qm.o properties.o gifs.o

libgifs: $(OBJ)
	mkdir -p ../../lib
	$(CPP) -shared -o ../../lib/libgifs.so $(OBJ)


%.o: %.cpp
	$(CPP) $(INCLUDE) $(CPPFLAGS) -fPIC -c $<
