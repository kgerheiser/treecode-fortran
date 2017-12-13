PHONY: tests clean all
include $(PFUNIT)/include/base.mk
.DEFAULT_GOAL = tests

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=$(OBJ_EXT))

SRC_DIR=../src
TEST_DIR=.

tests: $(OBJS)

testSuites.inc: $(SRCS)

FC=ifort
DEBUGFLAGS := -check all -traceback -g -fstack-protector -warn all -implicitnone
FFLAGS = $(DEBUGFLAGS)
FFLAGS += -I$(SRC_DIR) -I$(PFUNIT)/mod
LIBS = -L$(PFUNIT)/lib -lpfunit
FPPFLAGS=

bhtree:
	make -C $(SRC_DIR) test

tests: test.x

test.x: testSuites.inc $(OBJS)
	$(FC) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(SRC_DIR)/*$(OBJ_EXT) $(LIBS) $(FPPFLAGS)


%.F90: %.pf bhtree
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $<

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT) tests.xml
