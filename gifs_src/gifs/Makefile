CPP = g++
CPPFLAGS += -std=c++11 -O0 -Werror -Wall -Wextra -g -DARMA_DONT_USE_WRAPPER
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

INCLUDE := -I../../include -I../../armadillo/include
LIBDIR := ../../lib
LIB := libgifs.so
LIBRARY := $(LIBDIR)/$(LIB)

# automatic dependency generation; see:
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
COMPILE.cpp = $(CPP) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -fPIC -c
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

$(LIBRARY): $(OBJS) | $(LIBDIR)
	$(CPP) $(CPPFLAGS) -shared -o $@ $(OBJS)

%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

clean:
	rm *.o || true
	rm -rf $(DEPDIR) || true
	rm $(LIBRARY) || true


$(LIBDIR): ; mkdir -p $@
$(DEPDIR): ; mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
