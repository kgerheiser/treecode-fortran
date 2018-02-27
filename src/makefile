include make_include

# Source files to compile
SRC := $(wildcard *.f90)
OBJECTS := $(patsubst %.f90,%.o,$(SRC))

all: bhtree $(DEP_FILE)

test: FCFLAGS += $(DEBUGFLAGS)

bhtree: $(OBJECTS)
	$(FC) $(FCFLAGS) $^ -o $@

debug: FCFLAGS += $(DEBUGFLAGS)
debug: bhtree

fast: FCFLAGS += $(FASTFLAGS)
fast: bhtree

%.o: %.f90
	$(FC) -c $< $(FCFLAGS)

.PHONY: clean
clean:
	rm -f *.o *.mod 

# Make dependencies
.PHONY: depend
depend: $(DEP_FILE)

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(DEP_FILE): $(SRC)
	@echo "Making dependencies!"
	$(MAKEDEPEND) -w -o ./$(DEP_FILE) -f $(SRC)

include $(DEP_FILE)

