cmake_minimum_required(VERSION 3.0.0)

# set the project name
project(Shock_Solver LANGUAGES CXX CUDA)

# set(CMAKE_CXX_COMPILER "/opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/icpc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-zmm-usage=high -qopenmp -qmkl -Wall -pedantic -g -p -O3")
set(CMAKE_CXX_STANDARD 17)

# CUDA stuff
set(CUDA_PATH "/usr/local/cuda")
set(CUDA_INC_PATH "${CUDA_PATH}/include")
set(CUDA_BIN_PATH "$(CUDA_PATH)/bin")
set(CUDA_LIB_PATH "$(CUDA_PATH)/lib64")
find_package(CUDAToolkit)

set(NVCC "$(CUDA_BIN_PATH)/nvcc")
set(NVCC_FLAGS "$(-m64)")
set(NVCC_FLAGS "${NVCC_FLAGS} -g -dc -Wno-deprectaed-gpu-targets --std=c++11 --expt-relaxed-constexpr")
set(NVCC_GENCODES "$(-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_35,code=sm_35 \
		-gencode arch=compute_50,code=sm_50 \
		-gencode arch=compute_52,code=sm_52 \
		-gencode arch=compute_60,code=sm_60 \
		-gencode arch=compute_61,code=sm_61)")
set(CUDA_LINK_FLAGS "$(-dlink -Wno-deprecated-gpu-targets)")
set(INCLUDE "$(-I${CUDA_INC_PATH})")
set(CUDA_LIBS "$(-L${CUDA_LIB_PATH} -lcudart -lcufft -lcublas -lcudnn -lcurand)")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INCLUDE} ${CUDA_LIBS}")


file(GLOB SRC_FILES Code/*.cpp)
file(GLOB COMMON_FILES 
    "Code/*.cpp"
    "Code/*.cu")



#add the executable
add_executable(Shock_solver_CPU ${SRC_FILES})
add_executable(Shock_solver_GPU ${COMMON_FILES})

