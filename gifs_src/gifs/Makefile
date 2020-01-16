CPP = g++
CPPFLAGS = -std=c++11 -O0 -Werror -Wall -Wextra -g
#CPPFLAGS += -fsanitize=address

# The address sanitizer is an excellent way to find memory errors, but
# can be a bit of a pain to use. You will need to LD_PRELOAD the the
# asan shared library and will probably want to turn-off leak
# detection as below. For more details and documentation, see:
# https://github.com/google/sanitizers/wiki/AddressSanitizer

#LD_PRELOAD=/opt/gcc/7.3/lib64/libasan.so ASAN_OPTIONS=detect_leaks=0

INCLUDE = -I../../include

OBJ = gifs.o properties.o bomd.o gifs_implementation.o qm_interface.o qm_qchem.o

libgifs: $(OBJ)
	mkdir -p ../../lib
	$(CPP) -shared -o ../../lib/libgifs.so $(OBJ)


%.o: %.cpp
	$(CPP) $(INCLUDE) $(CPPFLAGS) -fPIC -c $<
