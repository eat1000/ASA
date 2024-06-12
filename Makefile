OUTPUT1 = asa
OUTPUT2 = psnps
OUTPUT3 = lai
OUTPUT4 = kartag
CXX = g++
PG = -pg
LDHTS = -lhts -L/usr/local/lib
EIGEN_PATH = /usr/local/include/Eigen
MKL_PATH = /opt/intel/mkl
HTS_PATH = /usr/local/include/htslib
CXXFLAGS1 = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include
CXXFLAGS2 = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG
CXXFLAGS3 = -w -O3 -m64 -lz -lbz2 -llzma -lpthread -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include -I $(HTS_PATH)
LIB1 += -lz -Wl,--start-group  $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.so $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.so $(MKL_PATH)/lib/intel64/libmkl_core.so -Wl,--end-group
SRC1 = asa.cpp geno_calc.cpp
SRC2 = psnps.cpp
SRC3 = lai.cpp
SRC4 = kartag.cpp

main : $(OUTPUT1) $(OUTPUT2) $(OUTPUT3) $(OUTPUT4)
$(OUTPUT1): 
	$(CXX) $(CXXFLAGS1) -o $(OUTPUT1) $(SRC1) $(LIB1)
$(OUTPUT2):        
	$(CXX) $(CXXFLAGS2) -o $(OUTPUT2) $(SRC2)
$(OUTPUT3):
	$(CXX) $(CXXFLAGS3) -o $(OUTPUT3) $(SRC3) $(LDHTS) -lz
$(OUTPUT4):
	$(CXX) $(CXXFLAGS3) -o $(OUTPUT4) $(SRC4) $(LDHTS)

