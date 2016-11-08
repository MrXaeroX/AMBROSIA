#
# Makefile for x86 Linux
#

EXENAME=ambrosia

ARCH=i686
OPTIMIZE=2

#make sure this is the correct compiler for your system
CC=g++
CP=cp
STRIP=strip -x

EXE_SRCDIR=../src
EXE_DSTDIR=../bin/bin_x64
EXE_OBJDIR=./obj-$(ARCH)

SHARED_SRCDIR=../src/shared

BASE_CFLAGS= -DLINUX -DNDEBUG
WARNOPT=-Wall -Wextra -Wno-variadic-macros

#full optimization
SHCFLAGS=-m64 $(WARNOPT) -pedantic -O$(OPTIMIZE) -mtune=$(ARCH) -march=pentium \
		 -mmmx -msse -msse2 -msse3 -ffast-math -fno-strict-aliasing
CFLAGS=$(BASE_CFLAGS) $(SHCFLAGS)

INCLUDEDIRS= -I$(SHARED_SRCDIR)
LIBRARYDIRS=

LDFLAGS=-ldl -lstdc++ -static-libgcc -lpthread

DO_CC=$(CC) $(CFLAGS) $(INCLUDEDIRS) -o $@ -c $<

#############################################################################
# SETUP AND BUILD
#############################################################################

$(EXE_OBJDIR)/%.o: $(EXE_SRCDIR)/%.cpp
	$(DO_CC)

OBJ = \
	$(EXE_OBJDIR)/md_common.o \
	$(EXE_OBJDIR)/md_config.o \
	$(EXE_OBJDIR)/md_log.o \
	$(EXE_OBJDIR)/md_main.o \
	$(EXE_OBJDIR)/md_math.o \
	$(EXE_OBJDIR)/md_model.o \
	$(EXE_OBJDIR)/md_model_energy.o \
	$(EXE_OBJDIR)/md_model_nbp.o \
	$(EXE_OBJDIR)/md_model_qeq.o \
	$(EXE_OBJDIR)/md_model_top.o \
	$(EXE_OBJDIR)/md_pdb_format.o \
	$(EXE_OBJDIR)/md_scrti.o \
	$(EXE_OBJDIR)/md_threads.o
	
$(EXENAME) : neat $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LIBRARYDIRS) $(LDFLAGS)
	$(STRIP) $(EXENAME)
	$(CP) $(EXENAME) $(EXE_DSTDIR)/$(EXENAME)

neat:
	-mkdir -p $(EXE_OBJDIR)
	
clean:
	-rm -f $(OBJ)
	-rm -f $(EXENAME)
spotless: clean
	-rm -r $(EXE_OBJDIR)
