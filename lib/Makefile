# PORT phase 1: build libtitan.so plus standalone executables from the
# patched sources in ../fortran (read-only here; objects go to obj/).
# Usage:  make -j1 [STATIC=1]
SRC   = ../fortran
FC    = gfortran
# -finit-local-zero: only so writout's never-set phantom-row local odepth(ngrs-1)
# prints 0 as in the reference .out; solver state does not depend on it.
# -fno-automatic is NOT needed (verified: all dumps byte-identical without it).
BASE  = -O2 -std=legacy -fPIC -ffp-contract=off -finit-local-zero -I$(SRC) -Isrc
ifeq ($(STATIC),1)
BASE += -fno-automatic
endif
FFLAGS = $(BASE)

OBJI = start01 start02 start03 start04 start05 start06 start07 start08 \
       start09 start10 writout writsty gridinit
OBJB = advectc advecti archive clzfil contin diffuse eddfac eddf eddg eos \
       gasmomh gasmomrh gasnrgh gasnrgr gasnrgrh grid intrc intrp mass \
       matgen matslv maxdel opac opnfil peqh peqr peqrh peq pardoc \
       radmomr radmomrh radnrgr radnrgrh readeos readopac reset \
       start step timstp totnrg updateh updater updaterh update viscous
OBJL = blas craypack linpack slvpack
OBJP = hdfinit hstore
CORE = $(addprefix obj/,$(addsuffix .o,$(OBJI) $(OBJB) $(OBJL) $(OBJP)))

all: libtitan.so titan.x titan_trace.x

# src/ overrides (start03/08/09 with TITAN_NCELL hook) take precedence over ../fortran
obj/%.o: src/%.f | obj
	$(FC) $(FFLAGS) -c $< -o $@
obj/%.o: $(SRC)/%.f | obj
	$(FC) $(FFLAGS) -c $< -o $@
obj:
	mkdir -p obj

# library: eos.o/opac.o renamed to eos_f_/opac_f_ and dispatched via hooks.c
LIBCORE = $(filter-out obj/eos.o obj/opac.o,$(CORE)) obj/eos_h.o obj/opac_h.o obj/hooks.o
obj/eos_h.o: obj/eos.o
	objcopy --redefine-sym eos_=eos_f_ $< $@
obj/opac_h.o: obj/opac.o
	objcopy --redefine-sym opac_=opac_f_ $< $@
obj/hooks.o: src/hooks.c | obj
	gcc -O2 -fPIC -c $< -o $@
libtitan.so: $(LIBCORE) obj/titan_lib.o obj/port_msg.o
	$(FC) -shared $(FFLAGS) $^ -o $@
titan.x: $(CORE) obj/titan.o obj/port_msg.o
	$(FC) $(FFLAGS) $^ -o $@
titan_trace.x: $(CORE) obj/titan_trace.o obj/port_msg.o
	$(FC) $(FFLAGS) $^ -o $@
clean:
	rm -rf obj libtitan.so titan.x titan_trace.x
