#
#                  this is the makefile for titan
#                       release version 1.2
#                   its executable is titan.x
#
FFLAGS = -O -cg89 #-mips2 #-trapuv
CFLAG1 = $(FFLAGS)
CFLAG2 = $(FFLAGS) -c
CC = f77 
RM = /bin/rm
AR = ar crv
RL = ranlib              # ranlib not necessary on HP or SGI

OBJI =  start01.o   start02.o   start03.o   start04.o   start05.o  \
	start06.o   start07.o   start08.o   start09.o   start10.o  \
	writout.o   writsty.o     titan.o  gridinit.o 

OBJB =  advectc.o   advecti.o   archive.o    clzfil.o   contin.o   \
	diffuse.o    eddfac.o      eddf.o      eddg.o      eos.o   \
	gasmomh.o              gasmomrh.o                          \
	gasnrgh.o   gasnrgr.o  gasnrgrh.o                          \
	   grid.o     intrc.o     intrp.o      mass.o              \
 	 matgen.o    matslv.o    maxdel.o      opac.o   opnfil.o   \
	   peqh.o      peqr.o     peqrh.o       peq.o   pardoc.o   \
	radmomr.o  radmomrh.o   radnrgr.o  radnrgrh.o              \
	readeos.o  readopac.o     reset.o                          \
	  start.o      step.o    timstp.o    totnrg.o              \
	updateh.o   updater.o  updaterh.o    update.o  viscous.o

OBJL =     blas.o  craypack.o   linpack.o   slvpack.o

OBJP =  hdfinit.o    hstore.o

OBJT =     code.o

BOX = titan_box.a      # compiled titan source code
LIB = titan_lib.a      # compiled math routines

CODE = -o titan.x      # titan executable

codes:                  $(OBJI) $(OBJB) $(OBJL) $(OBJP)
	$(CC) $(CFLAG1) $(OBJI) $(OBJB) $(OBJL) $(OBJP) $(CODE)

boxes:                  $(OBJI) $(OBJB) $(OBJL) $(OBJP)
	$(CC) $(CFLAG2) $(OBJI) $(OBJB) $(OBJL) $(OBJP)
	$(AR) $(LIB)                    $(OBJL)
	$(AR) $(BOX)            $(OBJB)
	$(RL)                   $(BOX) $(LIB)
	$(CC) $(CFLAG1) $(OBJI) $(BOX) $(LIB) $(OBJP) $(CODE)

titan:	                $(OBJI) $(BOX) $(LIB) $(OBJP)
	$(CC) $(CFLAG1) $(OBJI) $(BOX) $(LIB) $(OBJP) $(CODE)

code1:	                $(OBJT)
	$(CC) $(CFLAG1) $(OBJT) $(CODE)

clean:	
	$(RM) *.o *.x *.a

trash:
	$(RM) *.trace core
