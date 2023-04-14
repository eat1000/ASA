#ifndef _ASA_H_
#define _ASA_H_
#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <stdio.h>
#include <ctype.h>
#include <ctime>
#include <numeric> 
#include <math.h>
#include <omp.h>
#include <Eigen/Dense>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include <map>
#include <thread>
#include <iomanip>
#include <zlib.h>
#include <unordered_set>
#include <getopt.h>

using namespace boost;
using namespace Eigen;
using namespace std;

//#define SINGLE_PRECISION
#ifdef SINGLE_PRECISION
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef float eigenVar;
#else
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef double eigenVar;
#endif

class ASA
{
private:
//Genotype data input
ifstream bedfile;
int indi_num, snp_num, remain_snp_num;
const char* delimiter = " \t";
vector<string> famid, subid, snp, allele1, allele2;
vector<int> chrom, position, snp_to_delete;
//SNP information file
ifstream infofile;
vector<string> pop, info_pop, info_pop_MiA, pop_MiA, info_pop_MaA, pop_MaA;
vector<eigenVar> pop_maf, info_pop_maf;
int info_snp_num, not_in_info_num;
//GRM
eigenMatrix Z, Z_V, Z_D;
vector<vector<unsigned int>> Z_N;
//Calculating aiv, psv
eigenMatrix aiv, psv;
vector<string> diff_pop;
vector<int> info_in_bim;
int pop_num;

public:
ofstream writeLOG;
void readFAM(string);
void readBIM(string);
void make_X_Z(string, int, bool, bool, bool);
void grm_id_output(string);
void grm_output(string);
void grm_eigen();
void grm_eigvec_pc_eigval_output(string, string, string, int);
void make_X_psv(string, int);
void make_X_aiv(string, int);
void make_X_mle(string, int);
void psv_output(string);
void aiv_output(string);
void mle_output(string);
void snp_info_check(string);
void snp_info_input(string, bool);
void test(string);
};
#endif
