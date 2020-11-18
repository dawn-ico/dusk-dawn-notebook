all: driver.cpp diff_ops_cxx-naive.cpp
	@g++ \
		-o driver.cpp.o \
		-c driver.cpp \
		-DDAWN_ENABLE_BINDGEN \
		-DDAWN_PRECISION=1 \
		-I/usr/src/dawn/dawn/src/dawn4py/_external_src \
		-I/usr/src/dawn/build/_deps/gridtools-src/include \
		-isystem /usr/local/atlas/include \
		-isystem /usr/local/AtlasUtils/include/ \
		-isystem /usr/local/eckit/include \
		-isystem /usr/local/eckit/include/eckit \
		-isystem /usr/local/eckit/include/eckit/geometry \
		-isystem /usr/local/eckit/include/eckit/linalg \
		-isystem /usr/local/eckit/include/eckit/mpi \
		-isystem /usr/local/eckit/include/eckit/option \
		-fopenmp \
		-std=gnu++17
	@g++\
		driver.cpp.o \
		-o runner \
		-Wl,-rpath,/usr/local/atlas/lib:/usr/local/eckit/lib: \
		/usr/local/AtlasUtils/lib/libatlasUtilsLib.a \
		/usr/local/atlas/lib/libatlas.so \
		/usr/local/eckit/lib/libeckit_geometry.so \
		/usr/local/eckit/lib/libeckit_linalg.so \
		/usr/local/eckit/lib/libeckit_mpi.so \
		/usr/local/eckit/lib/libeckit_option.so \
		/usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so \
		/usr/lib/x86_64-linux-gnu/libpthread.so \
		/usr/local/eckit/lib/libeckit.so \
		/usr/lib/x86_64-linux-gnu/librt.so \
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
