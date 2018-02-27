PHONY: tests clean all
.DEFAULT_GOAL = debug

include $(PFUNIT)/include/base.mk

SRC_DIR=../src
TEST_DIR=.

include $(SRC_DIR)/make_include

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=$(OBJ_EXT))

SRC_SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRC_OBJECTS := $(patsubst %.f90,%.o,$(SRC_SRCS))
SRC_DIR_OBJS = $(filter-out ../src/main.o, $(SRC_OBJECTS))

testSuites.inc: $(SRCS)

INC_FLAGS = -I$(SRC_DIR) -I$(PFUNIT)/mod -I$(PFUNIT)/include
FCFLAGS += $(INC_FLAGS)
LIBS = -L$(PFUNIT)/lib -lpfunit
FPPFLAGS=

debug: FCFLAGS += $(DEBUGFLAGS)
debug: SRC_TARGET = debug
debug: tests

fast: FCFLAGS += $(FASTFLAGS)
fast: LIBS += -qopenmp
fast: SRC_TARGET=fast
fast: tests

bhtree:
	make -C $(SRC_DIR) $(SRC_TARGET)

tests: bhtree $(OBJS) test.x

test.x: testSuites.inc $(OBJS)
	$(FC) -o $@ $(INC_FLAGS) $(PFUNIT)/include/driver.F90 $(OBJS) $(SRC_DIR_OBJS) $(LIBS) $(FPPFLAGS)


%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

%$(OBJ_EXT): %.F90
	$(FC) -c $(FCFLAGS) $(FPPFLAGS) $<

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT) tests.xml
