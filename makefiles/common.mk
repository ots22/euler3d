SRC_DIR = ../../src
SRC = $(wildcard $(SRC_DIR)/*.f90 $(SRC_DIR)/*.F90)
OBJ = $(notdir $(SRC:.f90=.o))
OBJ := $(notdir $(OBJ:.F90=.o))

FCFLAGS += -L../../lib
FCFLAGS += -I../../lib

ALL_OBJ = $(OBJ)

$(SRC_DIR)/depend: $(SRC)
	makedepf90 -W -m"%m.mod" -b"." $^ > $@

-include $(SRC_DIR)/depend

euler3d: $(SRC_DIR)/depend $(ALL_OBJ)
	$(FC) $(FCFLAGS) $(ALL_OBJ) -o euler3d $(LDFLAGS)

# these rules gain more dependencies from 'depend', some of which are
# .mod files (hence the filter on the dependencies)
%.mod:
	$(FC) $(FCFLAGS) $(FCSYNTAXONLY) $(filter %.f90 %.F90, $^)

%.o:
	$(FC) $(FCFLAGS) $(filter %.f90 %.F90, $^) -c -o $@
