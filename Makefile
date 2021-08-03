OUTPUT1 = asa
OUTPUT2 = psnps
CXX = g++
EIGEN_PATH = /usr/local/include/Eigen
MKL_PATH = /opt/intel/mkl
CXXFLAGS1 = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include
CXXFLAGS2 = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG
LIB1 += -lz -Wl,--start-group  $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.so $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.so $(MKL_PATH)/lib/intel64/libmkl_core.so -Wl,--end-group  
SRC1 = asa.cpp geno_calc.cpp
SRC2 = psnps.cpp

main : $(OUTPUT1) $(OUTPUT2)
$(OUTPUT1): 
	$(CXX) $(CXXFLAGS1) -o $(OUTPUT1) $(SRC1) $(LIB1)
$(OUTPUT2):        
	$(CXX) $(CXXFLAGS2) -o $(OUTPUT2) $(SRC2)
