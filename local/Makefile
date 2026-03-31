FC := ifx
SHELL := /bin/bash

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

TARGET       := $(BIN_DIR)/compute_work
TARGET_ASYNC := $(BIN_DIR)/compute_work_async
SOURCES      := $(SRC_DIR)/compute_work.f90
SOURCES_ASYNC:= $(SRC_DIR)/compute_work_async.f90
OBJECTS      := $(OBJ_DIR)/compute_work.o
OBJECTS_ASYNC:= $(OBJ_DIR)/compute_work_async.o

EXT_INC_DIR := $(HOME)/local/include
EXT_LIB_DIR := $(HOME)/local/lib

MODULE_INIT := source /usr/share/Modules/init/bash >/dev/null 2>&1 || source /etc/profile.d/modules.sh >/dev/null 2>&1
MODULE_LOAD := module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

FFLAGS := -O3 -qopenmp -I$(EXT_INC_DIR)
LDFLAGS := -qopenmp -L$(EXT_LIB_DIR)
LDLIBS := -linterp -lnc -lhist

.PHONY: all clean

all: $(TARGET) $(TARGET_ASYNC)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(MODULE_INIT) && \
	$(MODULE_LOAD) && \
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS) $$(nf-config --flibs)

$(TARGET_ASYNC): $(OBJECTS_ASYNC) | $(BIN_DIR)
	$(MODULE_INIT) && \
	$(MODULE_LOAD) && \
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS_ASYNC) $(LDLIBS) $$(nf-config --flibs)

$(OBJ_DIR)/compute_work.o: $(SRC_DIR)/compute_work.f90 | $(OBJ_DIR)
	$(MODULE_INIT) && \
	$(MODULE_LOAD) && \
	$(FC) $(FFLAGS) $$(nf-config --fflags) -c $< -o $@

$(OBJ_DIR)/compute_work_async.o: $(SRC_DIR)/compute_work_async.f90 | $(OBJ_DIR)
	$(MODULE_INIT) && \
	$(MODULE_LOAD) && \
	$(FC) $(FFLAGS) $$(nf-config --fflags) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -f $(OBJECTS) $(OBJECTS_ASYNC) $(TARGET) $(TARGET_ASYNC)
