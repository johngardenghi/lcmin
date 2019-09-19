# Prerequisites path
VPATH := $(SOLV):$(HSL):$(F90PPATH):$(INTER)/cutest:$(INTER)/cuter

# Executable name
EXEC := lcmin

# SIFDECODE generated fonts
SIF := ELFUN EXTER GROUP RANGE

# Checks for availability of the problem set by the user
SIFPROB := $(shell if [ -e $(MASTSIF)/$(PROBNAME).SIF ]; then echo true; fi)
F90PROB := $(shell if [ -e $(F90PPATH)/$(PROBNAME).f90 ]; then echo true; fi)

# Library flags
LFLAGS := -llcmin -lalgencan -lhsl_ma57 -lhsl_ma48 -lmetis -lblas -lm

# Objects for linking LCMIN
COMMONOBJ := algsolver.o constants.o constraints.o initial_point.o	\
             kktsys.o lcmin.o lcmin_solver.o lsys_ma57.o mc58ad.o	\
             seval.o testing.o user_subs.o util.o variables.o

# Main Rules
all: lib
cutest: lib sifdec-cutest cutest-link

# ----------------------------------------- LINKING RULES -----------------------------------------
ifeq ($(SIFPROB), true)

sifdec-cutest:
	$(SIFDECODE)/bin/sifdecoder -A $(MYARCH) $(PROBNAME).SIF
	$(FC) $(FFLAGS) -c $(addsuffix .f, $(SIF))

else

sifdec-cutest:
	$(error $(PROBNAME).SIF does not exist in MASTSIF=$(MASTSIF) directory)

endif

cutest-link: lcmin_cutest.o $(addsuffix .o, $(SIF))
	$(FC) $(FFLAGS) -L$(CUTEST)/objects/$(MYARCH)/double -L$(LIB) -o $(EXEC) $^ -lcutest $(LFLAGS)
# -------------------------------------------------------------------------------------------------

# ----------------------------------------- LIBRARY RULE ------------------------------------------
lib: $(COMMONOBJ)
	cp lcmin.mod $(INC)
	ar rcs $(LIB)/liblcmin.a $(COMMONOBJ)
# -------------------------------------------------------------------------------------------------

# --------------------------------------- COMPILATION RULES ---------------------------------------
algsolver.o: algsolver.f90 constants.o util.o variables.o
	$(FC) $(FFLAGS) -c -o $@ $<

constants.o: constants.f90
	$(FC) $(FFLAGS) -c -o $@ $<

constraints.o: constraints.f90 constants.o mc58ad.o util.o	\
               variables.o testing.o
	$(FC) $(FFLAGS) -I$(INC) -c -o $@ $<

initial_point.o: initial_point.f90 algsolver.o constants.o	\
                 lsys_ma57.o util.o variables.o
	$(FC) $(FFLAGS) -c -o $@ $<

kktsys.o: kktsys.f90 constants.o lsys_ma57.o util.o variables.o
	$(FC) $(FFLAGS) -c -o $@ $<

lcmin.o: lcmin.f90 constants.o constraints.o initial_point.o	\
         lcmin_solver.o seval.o testing.o user_subs.o util.o	\
         variables.o
	$(FC) $(FFLAGS) -c -o $@ $<

lcmin_cutest.o: lcmin_cutest.f90
	$(FC) $(FFLAGS) -I$(INC) -c -o $@ $<

lcmin_solver.o: lcmin_solver.f90 constants.o constraints.o		\
	        initial_point.o kktsys.o seval.o util.o variables.o
	$(FC) $(FFLAGS) -c -o $@ $<

lsys_ma57.o: lsys_ma57.f90 constants.o variables.o
	$(FC) $(FFLAGS) -I$(INC) -c -o $@ $<

mc58ad.o: mc58ad.f
	$(FC) $(FFLAGS) -c -o $@ $<

seval.o: seval.f90 constants.o user_subs.o util.o
	$(FC) $(FFLAGS) -c -o $@ $<

testing.o: testing.f90
	$(FC) $(FFLAGS) -c -o $@ $<

user_subs.o: user_subs.f90 constants.o
	$(FC) $(FFLAGS) -c -o $@ $<

util.o: util.f90 constants.o
	$(FC) $(FFLAGS) -c -o $@ $<

variables.o: variables.f90 constants.o
	$(FC) $(FFLAGS) -c -o $@ $<

# -------------------------------------------------------------------------------------------------

# --------------------------------------- INSTALLATION RULES --------------------------------------
install-cutest:
	mkdir -p $(BIN)/cutest
	mv -f $(EXEC) $(BIN)/cutest
	mv -f OUTSDIF.d $(BIN)/cutest
	mv -f AUTOMAT.d $(BIN)/cutest
# -------------------------------------------------------------------------------------------------

# ----------------------------------------- CLEANING RULES ----------------------------------------
clean:
	rm -f $(addsuffix .f, $(SIF))
	rm -f $(OBJ)/*.SIF
	rm -f $(OBJ)/*.o
	rm -f $(OBJ)/*.mod

distclean: clean
	rm -f $(BIN)/f90/$(EXEC)

	rm -f $(BIN)/cuter/$(EXEC)
	rm -f $(BIN)/cuter/OUTSDIF.d
	rm -f $(BIN)/cuter/AUTOMAT.d

	rm -f $(BIN)/cutest/$(EXEC)
	rm -f $(BIN)/cutest/OUTSDIF.d
	rm -f $(BIN)/cutest/AUTOMAT.d

	rm -f $(LIB)/liblcmin.a
	rm -f $(INC)/lcmin.mod
# -------------------------------------------------------------------------------------------------
