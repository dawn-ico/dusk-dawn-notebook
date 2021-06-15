DAWN_DIR=/opt/dawn
ECKIT_DIR=/opt/anaconda3/envs/ecmwf-lab/lib/metview-bundle/
ATLAS_DIR=/opt/anaconda3/envs/ecmwf-lab/lib/metview-bundle/
ATLASUTILS_DIR=/opt/AtlasUtils/

all: driver.cpp laplacian_fvm_cxx-naive.cpp
	@g++ \
		-o driver.cpp.o \
		-c driver.cpp \
		-DDAWN_ENABLE_BINDGEN \
		-DDAWN_PRECISION=1 \
		-I$(DAWN_DIR)/lib/python3.9/site-packages/dawn4py/_external_src \
		-isystem $(ATLAS_DIR)/include \
		-isystem $(ATLASUTILS_DIR)/include/ \
		-isystem $(ECKIT_DIR)/include \
		-isystem $(ECKIT_DIR)/include/eckit \
		-isystem $(ECKIT_DIR)/include/eckit/geometry \
		-isystem $(ECKIT_DIR)/include/eckit/linalg \
		-isystem $(ECKIT_DIR)/include/eckit/mpi \
		-isystem $(ECKIT_DIR)/include/eckit/option \
		-fopenmp \
		-std=gnu++17
	@g++\
		driver.cpp.o \
		-o runner \
		-Wl,-rpath,$(ATLAS_DIR)/lib:$(ECKIT_DIR)/lib: \
		$(ATLASUTILS_DIR)/lib/libatlasUtilsLib.a \
		$(ATLAS_DIR)/lib/libatlas.so \
		$(ECKIT_DIR)/lib/libeckit_geometry.so \
		$(ECKIT_DIR)/lib/libeckit_linalg.so \
		$(ECKIT_DIR)/lib/libeckit_mpi.so \
		$(ECKIT_DIR)/lib/libeckit_option.so \
		$(ECKIT_DIR)/lib/libeckit.so \
		-lm \
		-lpthread \
		-ldl \
		-lrt \
		-lpthread \
		-ldl
	@mkdir -p out

clean:
	@rm \
		driver.cpp.o \
		runner
