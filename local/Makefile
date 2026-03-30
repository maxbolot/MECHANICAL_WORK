FC := ifx

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

TARGET := $(BIN_DIR)/compute_work
SOURCES := $(SRC_DIR)/compute_work.f90
OBJECTS := $(OBJ_DIR)/compute_work.o

EXT_INC_DIR := $(HOME)/local/include
EXT_LIB_DIR := $(HOME)/local/lib

NETCDF_FFLAGS := $(shell command -v nf-config >/dev/null 2>&1 && nf-config --fflags)
NETCDF_FLIBS := $(shell command -v nf-config >/dev/null 2>&1 && nf-config --flibs)

FFLAGS := -O3 -qopenmp $(NETCDF_FFLAGS) -I$(EXT_INC_DIR)
LDFLAGS := -qopenmp -L$(EXT_LIB_DIR)
LDLIBS := -linterp -lnc -lhist $(NETCDF_FLIBS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(OBJ_DIR)/compute_work.o: $(SRC_DIR)/compute_work.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -f $(OBJECTS) $(TARGET)
