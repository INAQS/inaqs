CPP = g++
CPPFLAGS += -std=c++11 -O0 -Werror -Wall -Wextra -g -Wshadow=compatible-local
CPPFLAGS += -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5
#CPPFLAGS += -fanalyzer # only with gcc v. 10+
#CPPFLAGS += -fsanitize=address

# The address sanitizer is an excellent way to find memory errors, but
# can be a bit of a pain to use. You will need to LD_PRELOAD the the
# asan shared library and will probably want to turn-off leak
# detection as below. For more details and documentation, see:
# https://github.com/google/sanitizers/wiki/AddressSanitizer

#LD_PRELOAD=/opt/gcc/7.3/lib64/libasan.so ASAN_OPTIONS=detect_leaks=0

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:%.cpp=%.o)


TDIR := tests
TSRCS := $(wildcard $(TDIR)/*.cpp)
TOBJS := $(TSRCS:%.cpp=%.o)


INCLUDE := -I../../include -I../../armadillo/include -I../ -I./
LIBDIR := ../../lib
LIB := libgifs.so
LIBRARY := $(LIBDIR)/$(LIB)

# automatic dependency generation; see:
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
COMPILE.cpp = $(CPP) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -fPIC -c
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@
LINKLIBS = -lhdf5

$(LIBRARY): $(OBJS) | $(LIBDIR)
	$(CPP) -shared $(OBJS) $(LINKLIBS) -o $@

# no dependency generation for tests
$(TDIR)/%.o : $(TDIR)/%.cpp
	$(COMPILE.cpp) -c $(OUTPUT_OPTION) $<

%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR) $(TDEPDIR)
	$(COMPILE.cpp) $(DEPFLAGS) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

clean:
	rm *.o || true
	rm $(TDIR)/*.o || true
	rm -rf $(DEPDIR) || true
	rm $(LIBRARY) || true

tests: $(TDIR)/test

$(TDIR)/test: $(TOBJS) $(LIBRARY)
	$(CPP) $(TOBJS) -lgifs -L$(LIBDIR) -lopenblas $(LINKLIBS) -o $@
	LD_LIBRARY_PATH=$(LIBDIR):$(LD_LIBRARY_PATH) ./$@ --success

$(LIBDIR): ; mkdir -p $@
$(DEPDIR): ; mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
