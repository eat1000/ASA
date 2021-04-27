OUTPUT = asa
CXX = g++
EIGEN_PATH = /usr/local/include/Eigen
MKL_PATH = /opt/intel/mkl
CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include
LIB += -lz -Wl,--start-group  $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.so $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.so $(MKL_PATH)/lib/intel64/libmkl_core.so -Wl,--end-group  
SRC = asa.cpp geno_calc.cpp

main : $(OUTPUT)
$(OUTPUT): 
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(SRC) $(LIB)

