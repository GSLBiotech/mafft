PREFIX = /usr/local
LIBDIR = $(PREFIX)/libexec/mafft
BINDIR = $(PREFIX)/bin
MANDIR = $(PREFIX)/share/man/man1
DESTDIR = 

ifeq ($(CC), cl)
	CFLAGS = /nologo /O2
	OBJ = obj
else
	OBJ = o
	#MNO_CYGWIN = -mno-cygwin

	ENABLE_MULTITHREAD = -Denablemultithread
	# Comment out the above line if your compiler 
	# does not support TLS (thread-local strage).

	#ENABLE_ATOMIC = -Denableatomic
	# Comment out the above line if your compiler 
	# does not support "atomic_int".

	CC = gcc
	#CC = icc
	CFLAGS = -O3
	#CFLAGS = -O3 -fPIC
	# add -fPIC when building .so files

	#CC = icc
	#CFLAGS = -fast
	# if you have icc, use this.

	#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3
	#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3  -fprofile-arcs -ftest-coverage 
	#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -DMALLOC_CHECK_=3 # for shark, valgrind
	#CFLAGS =  -O0  -fPIC -pedantic -Wall -std=c99 -g -DMALLOC_CHECK_=3 -lprofiler  # ?



	ifdef ENABLE_MULTITHREAD
	LIBS = -lm  -lpthread
	else
	LIBS = -lm
	endif

	ifdef ENABLE_ATOMIC
	STDF = -std=c11
	else
	STDF = -std=c99
	endif
endif

MYCFLAGS = $(MNO_CYGWIN) $(ENABLE_MULTITHREAD) $(ENABLE_ATOMIC) $(STDF) $(CFLAGS)

INSTALL = install

PROGS = dvtditr dndfast7 dndblast sextet5 mafft-distance pairlocalalign \
		multi2hat3s pairash addsingle maffttext2hex hex2maffttext \
        splittbfast disttbfast tbfast nodepair mafft-profile f2cl mccaskillwrap contrafoldwrap countlen \
		seq2regtable regtable2seq score getlag dndpre setcore replaceu restoreu setdirection makedirectionlist version
SOS = libdisttbfast.so
DLLS = libdisttbfast.dll
DYLIBS = libdisttbfast.dylib

PERLPROGS = mafftash_premafft.pl seekquencer_premafft.pl
SCRIPTS = mafft mafft-homologs.rb mafft-sparsecore.rb
OBJSETDIRECTION = mtxutl.$(OBJ) io.$(OBJ) setdirection.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJREPLACEU = mtxutl.$(OBJ) io.$(OBJ) replaceu.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJRESTOREU = mtxutl.$(OBJ) io.$(OBJ) restoreu.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJREGTABLE2SEQ = mtxutl.$(OBJ) io.$(OBJ) regtable2seq.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJSEQ2REGTABLE = mtxutl.$(OBJ) io.$(OBJ) seq2regtable.$(OBJ) defs.$(OBJ)
OBJCOUNTLEN = mtxutl.$(OBJ) io.$(OBJ) countlen.$(OBJ) defs.$(OBJ)
OBJF2CL = mtxutl.$(OBJ) io.$(OBJ) f2cl.$(OBJ) constants.$(OBJ) defs.$(OBJ)
OBJMCCASKILLWRAP = mtxutl.$(OBJ) io.$(OBJ) mccaskillwrap.$(OBJ) constants.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJCONTRAFOLDWRAP = mtxutl.$(OBJ) io.$(OBJ) contrafoldwrap.$(OBJ) constants.$(OBJ) defs.$(OBJ) mltaln9.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJMULTI2HAT3S = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) \
				 multi2hat3s.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJPAIRASH = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
			 Falign.$(OBJ) MSalignmm.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) MSalign11.$(OBJ) suboptalign11.$(OBJ) SAalignmm.$(OBJ) \
			 pairash.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJPAIRLOCALALIGN = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		            Falign.$(OBJ) MSalignmm.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) MSalign11.$(OBJ) suboptalign11.$(OBJ) SAalignmm.$(OBJ) \
					pairlocalalignmain.$(OBJ) pairlocalalign.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJDUMMY = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			disttbfast_dummy.$(OBJ) dummy.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJSPLITFROMALN = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ)  Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			splitfromaln.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJSPLITTBFAST = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ)  Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			splittbfast.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) 
OBJSPLITTBFASTP = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJDISTTBFAST = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			disttbfast.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) addfunctions.$(OBJ)
OBJMAKEDIRECTIONLIST = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ)  Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			makedirectionlist.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) addfunctions.$(OBJ)
OBJTBFAST = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) MSalignmm.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) \
			tbfast.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) addfunctions.$(OBJ) \
			pairlocalalign.$(OBJ) MSalign11.$(OBJ) 
OBJNODPAIR = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) MSalignmm.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) \
			nodepair.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) addfunctions.$(OBJ) \
			pairlocalalign.$(OBJ) MSalign11.$(OBJ)
OBJADDSINGLE = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) MSalignmm.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) \
			addsingle.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ) addfunctions.$(OBJ)
OBJSETCORE = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) \
			setcore.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJTDITR = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) nj.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ)  \
		   Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) fftFunctions.$(OBJ) fft.$(OBJ)  \
		   tditeration.$(OBJ) tditr.$(OBJ) defs.$(OBJ) SAalignmm.$(OBJ) treeOperation.$(OBJ)
OBJDVTDITR = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) nj.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
			Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) MSalignmm.$(OBJ) fftFunctions.$(OBJ) fft.$(OBJ)  \
			tditeration.$(OBJ) dvtditr.$(OBJ) defs.$(OBJ) SAalignmm.$(OBJ) treeOperation.$(OBJ) addfunctions.$(OBJ) 
OBJGETLAG = mtxutl.$(OBJ) io.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ)  \
		    Falign.$(OBJ) Falign_localhom.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) SAalignmm.$(OBJ) MSalignmm.$(OBJ) Dalignmm.$(OBJ) \
			getlag.$(OBJ) defs.$(OBJ) fft.$(OBJ) fftFunctions.$(OBJ)
OBJGAPFILL = mtxutl.$(OBJ) io.$(OBJ) constants.$(OBJ) gapfill.$(OBJ) defs.$(OBJ)
OBJDNDFAST5 = dndfast5.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ)
OBJDNDBLAST = dndblast.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJDNDFAST7 = dndfast7.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJDNDFAST6 = dndfast6.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ)
OBJDNDFAST4 = dndfast4.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ)
OBJDNDFAST6 = dndfast6.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) defs.$(OBJ)
OBJSEXTET5 = io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) sextet5.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJDISTANCE = io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) mafft-distance.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJTRIPLET6 = io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) triplet6.$(OBJ) defs.$(OBJ)
OBJTRIPLET5 = io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) triplet5.$(OBJ) defs.$(OBJ)
OBJOCTET4 = io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) octet4.$(OBJ) defs.$(OBJ)
OBJDNDPRE = dndpre.$(OBJ) io.$(OBJ) constants.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)
OBJGALN = io.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) tddis.$(OBJ) constants.$(OBJ) partSalignmm.$(OBJ) MSalignmm.$(OBJ) Lalignmm.$(OBJ) rna.$(OBJ) Salignmm.$(OBJ) Dalignmm.$(OBJ) \
          SAalignmm.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ) Falign.$(OBJ) Falign_localhom.$(OBJ) fftFunctions.$(OBJ) fft.$(OBJ) mafft-profile.$(OBJ) defs.$(OBJ)
OBJSCORE = io.$(OBJ) mtxutl.$(OBJ) mltaln9.$(OBJ) score.$(OBJ) constants.$(OBJ) defs.$(OBJ) Galign11.$(OBJ) Lalign11.$(OBJ) genalign11.$(OBJ)

HEADER = mltaln.h mtxutl.h mafft.h
FFTHEADER = fft.h

MANPAGES = mafft.1 mafft-homologs.1




all : $(PERLPROGS) $(PROGS) $(SCRIPTS)
	cp $(SCRIPTS) ../scripts
	chmod 755 ../scripts/*
	cp $(PERLPROGS) $(PROGS) ../binaries
	chmod 755 ../binaries/*
	cp $(MANPAGES) ../binaries
	@echo done.

sos : $(SOS)
dylibs : $(DYLIBS)
dlls : $(DLLS)

univscript: univscript.tmpl Makefile
	sed "s:_PROGS:$(PROGS):" univscript.tmpl  > univscript

mafft: mafft.tmpl mltaln.h
	sed "s:_LIBDIR:$(LIBDIR):" mafft.tmpl  > mafft

mafft-homologs.rb: mafft-homologs.tmpl
#	cp mafft-homologs.tmpl mafft-homologs.rb
	sed "s:_BINDIR:$(BINDIR):" mafft-homologs.tmpl  > mafft-homologs.rb

mafft-sparsecore.rb: mafft-sparsecore.tmpl
#	cp mafft-sparsecore.tmpl mafft-sparsecore.rb
	sed "s:_BINDIR:$(BINDIR):" mafft-sparsecore.tmpl  > mafft-sparsecore.rb

mltaln.h : functions.h
	touch mltaln.h

version : version.c mltaln.h
	$(CC) -o $@ version.c $(MYCFLAGS) $(LDFLAGS) $(LIBS) 

maffttext2hex : maffttext2hex.c
	$(CC) -o $@ maffttext2hex.c $(MYCFLAGS) $(LDFLAGS) $(LIBS)

hex2maffttext : hex2maffttext.c
	$(CC) -o $@ hex2maffttext.c $(MYCFLAGS) $(LDFLAGS) $(LIBS)

tbfast : $(OBJTBFAST)
	$(CC) -o $@ $(OBJTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

nodepair : $(OBJNODPAIR)
	$(CC) -o $@ $(OBJNODPAIR) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

addsingle : $(OBJADDSINGLE)
	$(CC) -o $@ $(OBJADDSINGLE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

disttbfast : $(OBJDISTTBFAST)
	$(CC) -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.so : $(OBJDISTTBFAST)
	$(CC) -shared -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.dylib : $(OBJDISTTBFAST)
	$(CC) -dynamiclib -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

libdisttbfast.dll : $(OBJDISTTBFAST)
	$(CC) -shared -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

makedirectionlist : $(OBJMAKEDIRECTIONLIST)
	$(CC) -o $@ $(OBJMAKEDIRECTIONLIST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

splittbfast : $(OBJSPLITTBFAST)
	$(CC) -o $@ $(OBJSPLITTBFAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

splitfromaln : $(OBJSPLITFROMALN)
	$(CC) -o $@ $(OBJSPLITFROMALN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dummy : $(OBJDUMMY)
	$(CC) -o $@ $(OBJDUMMY) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

setcore : $(OBJSETCORE)
	$(CC) -o $@ $(OBJSETCORE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

countlen : $(OBJCOUNTLEN)
	$(CC) -o $@ $(OBJCOUNTLEN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

seq2regtable : $(OBJSEQ2REGTABLE)
	$(CC) -o $@ $(OBJSEQ2REGTABLE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

regtable2seq : $(OBJREGTABLE2SEQ)
	$(CC) -o $@ $(OBJREGTABLE2SEQ) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

setdirection : $(OBJSETDIRECTION)
	$(CC) -o $@ $(OBJSETDIRECTION) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

replaceu : $(OBJREPLACEU)
	$(CC) -o $@ $(OBJREPLACEU) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

restoreu : $(OBJRESTOREU)
	$(CC) -o $@ $(OBJRESTOREU) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

f2cl : $(OBJF2CL)
	$(CC) -o $@ $(OBJF2CL) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

mccaskillwrap : $(OBJMCCASKILLWRAP)
	$(CC) -o $@ $(OBJMCCASKILLWRAP) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

contrafoldwrap : $(OBJCONTRAFOLDWRAP)
	$(CC) -o $@ $(OBJCONTRAFOLDWRAP) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

pairlocalalign : $(OBJPAIRLOCALALIGN)
	$(CC) -o $@ $(OBJPAIRLOCALALIGN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

pairash : $(OBJPAIRASH)
	$(CC) -o $@ $(OBJPAIRASH) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

multi2hat3s : $(OBJMULTI2HAT3S)
	$(CC) -o $@ $(OBJMULTI2HAT3S) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

getlag : $(OBJGETLAG)
	$(CC) -o $@ $(OBJGETLAG) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

tditr : $(OBJTDITR)
	$(CC) -o $@ $(OBJTDITR) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dvtditr : $(OBJDVTDITR)
	$(CC) -o $@ $(OBJDVTDITR) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

mafft-profile : $(OBJGALN)
	$(CC) -o $@ $(OBJGALN) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

gapfill : $(OBJGAPFILL)
	$(CC) -o $@ $(OBJGAPFILL) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndfast4 : $(OBJDNDFAST4)
	$(CC) -o $@ $(OBJDNDFAST4) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndfast5 : $(OBJDNDFAST5)
	$(CC) -o $@ $(OBJDNDFAST5) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndfast6 : $(OBJDNDFAST6)
	$(CC) -o $@ $(OBJDNDFAST6) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndfast7 : $(OBJDNDFAST7)
	$(CC) -o $@ $(OBJDNDFAST7) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndblast : $(OBJDNDBLAST)
	$(CC) -o $@ $(OBJDNDBLAST) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndfast3 : $(OBJDNDFAST3)
	$(CC) -o $@ $(OBJDNDFAST3) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

triplet : $(OBJTRIPLET)
	$(CC) -o $@ $(OBJTRIPLET) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

triplet3 : $(OBJTRIPLET3)
	$(CC) -o $@ $(OBJTRIPLET3) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

sextet3 : $(OBJSEXTET3)
	$(CC) -o $@ $(OBJSEXTET3) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

sextet4 : $(OBJSEXTET4)
	$(CC) -o $@ $(OBJSEXTET4) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

sextet5 : $(OBJSEXTET5)
	$(CC) -o $@ $(OBJSEXTET5) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

mafft-distance : $(OBJDISTANCE)
	$(CC) -o $@ $(OBJDISTANCE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

triplet5 : $(OBJTRIPLET5)
	$(CC) -o $@ $(OBJTRIPLET5) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

triplet6 : $(OBJTRIPLET6)
	$(CC) -o $@ $(OBJTRIPLET6) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

octet4 : $(OBJOCTET4)
	$(CC) -o $@ $(OBJOCTET4) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

dndpre : $(OBJDNDPRE)
	$(CC) -o $@ $(OBJDNDPRE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

score : $(OBJSCORE)
	$(CC) -o $@ $(OBJSCORE) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

genMtx : $(OBJGENMTX)
	$(CC) -o $@ $(OBJGENMTX) $(MYCFLAGS) $(LDFLAGS) $(LIBS)

mafftash_premafft.pl : mafftash_premafft.tmpl
	cp mafftash_premafft.tmpl mafftash_premafft.pl

seekquencer_premafft.pl : seekquencer_premafft.tmpl
	cp seekquencer_premafft.tmpl seekquencer_premafft.pl

gapfill.$(OBJ) : gapfill.c $(HEADER)
	$(CC) $(MYCFLAGS) -c gapfill.c

mltaln9.$(OBJ) : mltaln9.c $(HEADER)
	$(CC) $(MYCFLAGS) -c mltaln9.c

tddis.$(OBJ) : tddis.c $(HEADER)
	$(CC) $(MYCFLAGS) -c tddis.c

constants.$(OBJ) : constants.c miyata.h miyata5.h blosum.c DNA.h JTT.c $(HEADER)
	$(CC) $(MYCFLAGS) -c constants.c

defs.$(OBJ) : defs.c 
	$(CC) $(MYCFLAGS) -c defs.c

#A+++alignmm.$(OBJ) : SA+++alignmm.c $(HEADER)
#	$(CC) $(MYCFLAGS) -c SA+++alignmm.c -o A+++alignmm.$(OBJ)

Salignmm.$(OBJ) : Salignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Salignmm.c 

Dalignmm.$(OBJ) : Dalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Dalignmm.c 

MSalignmm.$(OBJ) : MSalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c MSalignmm.c 

partSalignmm.$(OBJ) : partSalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c partSalignmm.c  

Lalign11.$(OBJ) : Lalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Lalign11.c 

genalign11.$(OBJ) : genalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c genalign11.c 

suboptalign11.$(OBJ) : suboptalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c suboptalign11.c 

Galign11.$(OBJ) : Galign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Galign11.c 

MSalign11.$(OBJ) : MSalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c MSalign11.c 

SAalignmm.$(OBJ) : SAalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c SAalignmm.c -o SAalignmm.$(OBJ)

Lalignmm.$(OBJ) : Lalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Lalignmm.c 

rna.$(OBJ) : rna.c $(HEADER)
	$(CC) $(MYCFLAGS) -c rna.c 

disttbfast.$(OBJ) : disttbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c disttbfast.c

splitfromaln.$(OBJ) : splitfromaln.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splitfromaln.c

splittbfast.$(OBJ) : splittbfast.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splittbfast.c

splittbfast2.$(OBJ) : splittbfast2.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splittbfast2.c

makedirectionlist.$(OBJ) : makedirectionlist.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c makedirectionlist.c

disttbfast_dummy.$(OBJ) : disttbfast_dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c disttbfast_dummy.c

dummy.$(OBJ) : dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c dummy.c

tbfast.$(OBJ) : tbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c tbfast.c

nodepair.$(OBJ) : nodepair.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c nodepair.c

addsingle.$(OBJ) : addsingle.c $(HEADER) $(FFTHEADER)
addsingle.$(OBJ) : addsingle.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c addsingle.c

tbfast2.$(OBJ) : tbfast2.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c tbfast2.c

setcore.$(OBJ) : setcore.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c setcore.c

getlag.$(OBJ) : getlag.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c getlag.c

tditr.$(OBJ) : tditr.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c tditr.c

dvtditr.$(OBJ) : dvtditr.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c dvtditr.c

tditeration.$(OBJ) : tditeration.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c tditeration.c

mafft-profile.$(OBJ) : mafft-profile.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c mafft-profile.c

dndfast4.$(OBJ) : dndfast4.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast4.c

dndfast5.$(OBJ) : dndfast5.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast5.c

dndfast6.$(OBJ) : dndfast6.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast6.c

dndfast7.$(OBJ) : dndfast7.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast7.c

dndblast.$(OBJ) : dndblast.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndblast.c

dndfast3.$(OBJ) : dndfast3.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast3.c

dndpre.$(OBJ) : dndpre.c $(HEADER)
	$(CC) $(MYCFLAGS) -c dndpre.c

countlen.$(OBJ) : countlen.c $(HEADER)
	$(CC) $(MYCFLAGS) -c countlen.c

seq2regtable.$(OBJ) : seq2regtable.c $(HEADER)
	$(CC) $(MYCFLAGS) -c seq2regtable.c

regtable2seq.$(OBJ) : regtable2seq.c $(HEADER)
	$(CC) $(MYCFLAGS) -c regtable2seq.c

f2cl.$(OBJ) : f2cl.c $(HEADER)
	$(CC) $(MYCFLAGS) -c f2cl.c

setdirection.$(OBJ) : setdirection.c $(HEADER)
	$(CC) $(MYCFLAGS) -c setdirection.c

replaceu.$(OBJ) : replaceu.c $(HEADER)
	$(CC) $(MYCFLAGS) -c replaceu.c

restoreu.$(OBJ) : restoreu.c $(HEADER)
	$(CC) $(MYCFLAGS) -c restoreu.c

mccaskillwrap.$(OBJ) : mccaskillwrap.c $(HEADER)
	$(CC) $(MYCFLAGS) -c mccaskillwrap.c

contrafoldwrap.$(OBJ) : contrafoldwrap.c $(HEADER)
	$(CC) $(MYCFLAGS) -c contrafoldwrap.c

pairlocalalign.$(OBJ) : pairlocalalign.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pairlocalalign.c

pairlocalalignmain.$(OBJ) : pairlocalalignmain.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pairlocalalignmain.c

pairash.$(OBJ) : pairash.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pairash.c

multi2hat3s.$(OBJ) : multi2hat3s.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c multi2hat3s.c

io.$(OBJ) : io.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c io.c

nj.$(OBJ) : nj.c $(HEADER)
	$(CC) $(MYCFLAGS) -c nj.c

treeOperation.$(OBJ) : treeOperation.c $(HEADER)
	$(CC) $(MYCFLAGS) -c treeOperation.c

sextet5.$(OBJ) : sextet5.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c sextet5.c

mafft-distance.$(OBJ) : mafft-distance.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c mafft-distance.c

maffttext2hex.$(OBJ) : maffttext2hex.c
	$(CC) $(MYCFLAGS) -c maffttext2hex.c

hex2maffttext.$(OBJ) : hex2maffttext.c
	$(CC) $(MYCFLAGS) -c hex2maffttext.c

triplet6.$(OBJ) : triplet6.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c triplet6.c

fft.$(OBJ) : fft.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fft.c 

fftFunctions.$(OBJ) : fftFunctions.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fftFunctions.c

Falign.$(OBJ) : Falign.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c Falign.c

Falign_localhom.$(OBJ) : Falign_localhom.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c Falign_localhom.c

mtxutl.$(OBJ) : mtxutl.c 
	$(CC) $(MYCFLAGS) -c mtxutl.c

addfunctions.$(OBJ) : addfunctions.c $(HEADER)
	$(CC) $(MYCFLAGS) -c addfunctions.c

score.$(OBJ) : score.c $(HEADER)
	$(CC) $(MYCFLAGS) -c score.c

clean :
	rm -f *.$(OBJ) *.a *.exe *~ $(PERLPROGS) $(PROGS) $(SCRIPTS) $(SOS) $(DYLIBS) $(DLLS) *.gcda *.gcno
#	rm -f ../binaries/* ../scripts/*

install : all
	mkdir -p $(DESTDIR)$(LIBDIR)
	chmod 755 $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(BINDIR)
	chmod 755 $(DESTDIR)$(BINDIR)
	chmod 755 $(SCRIPTS)
	$(INSTALL) $(SCRIPTS)  $(DESTDIR)$(BINDIR)
	chmod 755 $(PROGS) ||:     # in MinGW, it's ok if this fails
	$(INSTALL) -s $(PROGS) $(DESTDIR)$(LIBDIR)
	$(INSTALL) $(PERLPROGS) $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m 644 $(MANPAGES) $(DESTDIR)$(LIBDIR)

	( cd $(DESTDIR)$(BINDIR); \
rm -f linsi ginsi einsi fftns fftnsi nwns nwnsi xinsi qinsi; \
rm -f mafft-linsi mafft-ginsi mafft-einsi mafft-fftns mafft-fftnsi mafft-nwns mafft-nwnsi mafft-xinsi mafft-qinsi mafft-randomcore.rb ; \
ln -s mafft linsi; ln -s mafft ginsi; ln -s mafft fftns; \
ln -s mafft fftnsi; ln -s mafft nwns; ln -s mafft nwnsi; \
ln -s mafft einsi; \
ln -s mafft mafft-linsi; ln -s mafft mafft-ginsi; ln -s mafft mafft-fftns; \
ln -s mafft mafft-fftnsi; ln -s mafft mafft-nwns; ln -s mafft mafft-nwnsi; \
ln -s mafft mafft-einsi; ln -s mafft mafft-xinsi; ln -s mafft mafft-qinsi;\
rm -f mafft-profile mafft-profile.exe; ln -s $(LIBDIR)/mafft-profile .; \
rm -f mafft-distance mafft-distance.exe; ln -s $(LIBDIR)/mafft-distance . )

	mkdir -p $(DESTDIR)$(MANDIR)
	chmod 755 $(DESTDIR)$(MANDIR)
	$(INSTALL) -m 644 $(MANPAGES) $(DESTDIR)$(MANDIR)
# remove incorrectly installed manpages by previous versions
#	rm -f /usr/local/man/man1/mafft.1 /usr/local/man/man1/mafft-homologs.1 
