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
        int indi_num;
        int snp_num;
        const char* delimiter = " \t";
        vector<string> famid;
        vector<string> subid;
        vector<int> chrom;
        vector<string> snp;
        vector<int> position;
        vector<string> allele1;
        vector<string> allele2;
        vector<int> allele1_num;
        vector<int> allele2_num;
        vector<int> miss_num;
        vector<double> miss_fre;
        vector<double> allele1_fre;
        vector<double> allele2_fre;
        vector<double> maf;
        vector<string> minor_allele;
        vector<int> snp_to_delete;
        vector<int> min_alle_eqt_alle2;
        vector<pair<int, pair<int, int> > > chr_idx;
        int remain_snp_num;
        //SNP information file
        ifstream infofile;
        vector<string> pop, info_pop;
        vector<eigenVar> pop_maf, info_pop_maf;
        int info_snp_num, not_in_info_num;
	//GRM
        eigenMatrix Z;
        eigenMatrix Z_V;
        eigenMatrix Z_D;
	vector<vector<unsigned int>> Z_N;
        //Calculating aiv, psv
        eigenMatrix aiv, psv;
        vector<string> diff_pop;
	vector<int> info_in_bim;
        int pop_num;
		
    public:
        void allele_info_stat(string, double, double, double, double, int);
        void allele_info_output(string);
        ofstream writeLOG;
        template <typename T> void sort_indexes(const vector<T> & , vector<int> &);
        template <typename T1, typename T2, typename T3> static bool cmp_3arr_1(const pair<T1, pair<T2, T3> > &, const pair<T1, pair<T2, T3> > &);
        template <typename T1, typename T2, typename T3> static bool cmp_3arr_2(const pair<T1, pair<T2, T3> > &, const pair<T1, pair<T2, T3> > &);
        template <typename T1, typename T2> static bool cmp_2arr_1(const pair<T1, T2> &, const pair<T1, T2> &);
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
        void snp_info_input(string);
//	double cubic_fun(double, double, double, double, double);
//	double cubic_fun_der(double, double, double, double, double);
//	double cubic_eq_solver(double, double, double, double, double);
};
#endif
