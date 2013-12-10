# Makefile for rocktools, a 3-D triangle mesh manipulation toolkit

# user-customizable (only used when "make install")
#BIN = /usr/local/bin
BIN = ~/bin
#BIN = .

CC=gcc
LINKER=gcc

CFLAGS=-Wall
ifdef DEBUG
  CFLAGS+=-g -p -ggdb -fbounds-check
else
  CFLAGS+=-O3 -funroll-loops -fomit-frame-pointer
  CFLAGS+=-mtune=generic
endif

ifeq ($(UNAME), Linux)
  LIBS=
endif
ifeq ($(UNAME), Darwin)
  LIBS=-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include -L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib
endif
LIBS+=-lm -lpng

HFILES = structs.h
CFILES = inout.c\
	utils.c

EXE = rockdetail\
	rockcreate\
	rocksmooth\
	rocktrim\
	rockerode\
	rockconvert\
	rockmarker\
	rockxray\
	rocksplit\
	rockpng\
	rockslice\
	rockinfo

all : $(EXE)
	@echo "Rocktools made"

install : $(EXE)
	cp rockcreate $(BIN)/rockcreate
	cp rockdetail $(BIN)/rockdetail
	cp rocksmooth $(BIN)/rocksmooth
	cp rocktrim $(BIN)/rocktrim
	cp rockerode $(BIN)/rockerode
	cp rockconvert $(BIN)/rockconvert
	cp rockmarker $(BIN)/rockmarker
	cp rockxray $(BIN)/rockxray
	cp rocksplit $(BIN)/rocksplit
	cp rockinfo $(BIN)/rockinfo
	cp rockpng $(BIN)/rockpng
	cp rockdice $(BIN)/rockdice
	cp rockslice $(BIN)/rockslice

# alternatively, make vort3d with the debug flags turned on
debug: CFLAGS = -pg -ggdb -Wall -ftrapv
debug: all

# ask that make not delete the object files
#.PRECIOUS : objects/%.o

#objects/%.o : %.c $(HFILES)
#	$(CC) $(CFLAGS) -c -o $@ $<

rockcreate: rockcreate.c createutil.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -DADJ_NODE -DCONN -o $@ rockcreate.c createutil.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockdetail: rockdetail.c detailutil.c smoothutil.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -DADJ_NODE -DCONN -DDETAIL -o $@ rockdetail.c detailutil.c smoothutil.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockerode : rockerode.c erodeutil.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -DADJ_NODE -DERODE -DCONN -o $@ rockerode.c erodeutil.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rocksmooth : rocksmooth.c smoothutil.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -DADJ_NODE -DCONN -o $@ rocksmooth.c smoothutil.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockconvert: rockconvert.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rockconvert.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rocktrim: rocktrim.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rocktrim.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockinfo: rockinfo.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rockinfo.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockpng: rockpng.c pngutil.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rockpng.c pngutil.c $(CFILES) $(LIBS)
#@echo "$(@F) made"

rocksplit: rocksplit.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rocksplit.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

rockslice: rockslice.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rockslice.c $(CFILES) $(LIBS)

# any others that don't need connectivity or adjacent nodes
rock% : rock%.c %util.c $(CFILES) $(HFILES) Makefile
	$(CC) $(CFLAGS) -o $@ rock$*.c $*util.c $(CFILES) $(LIBS)
#	@echo "$(@F) made"

clean : 
	rm -f $(EXE)
