#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cstring>
#include<algorithm>
#include <numeric>
#include<map>
#include<stdio.h>
#include<getopt.h>
#include<cmath>
#include<set>
#include <bitset>
#include<unordered_set>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <omp.h>
#include <thread>


using namespace std;
using namespace boost;

class Dataset {
private:
const char* delimiter = " \t";
vector<string> pop_name, nonpop_name, famid, subid, chrom, snp, allele1, allele2, MA, ref_id, ref_pop, ref_pop_name;
vector<int> position, ref_pop_cnt, id_pop, pop_pos, nonpop_pos, pop_snp_idx;
vector<bool> snp_pop;
vector<vector<int>> allele1_cnt, miss_cnt;
vector<vector<double>> miss_frq, maf;
int pop_num, nonpop_num, batch_size, indi_num, snp_num, remain_snp_num, ref_indi_num, ref_pop_num, pop_snp_num, out_snp_num, max_thread_num;
double maf_min, maf_max, miss_max;
bool error_flag = false;	
string bed_file, fam_file, bim_file, ref_file, out_file;	
ifstream bedfile, famfile, bimfile, reffile, outfile;
set<string> subid_set, ref_pop_set, ref_id_set;
unordered_multiset<string> ref_pop_multiset;
map<string, int> id_pop_map;

public:
int flag[25]{};
string arg[25];
ofstream logFile;
int thread_num;

void help()
{
	cout << "Population-Specific SNP Screener (PSNPS)" << endl;
	cout << "Options:" << endl;
	cout << "--bfile\t\t" << "Input genotype file in plink binary format." << endl;
	cout << "--ref-pop\t" << "Input file that describes reference populations. The file is expected to have two columns without headers:" << endl;
	cout << "\t\t" << "the first is individual ID and the second is the population that the individual belongs to." << endl;
	cout << "--pop\t\t" << "Screen SNPs that are specific to the specified population." << endl;
	cout << "\t\t" << "If multiple populations are specified, SNPs polymorphic in all specified populations and monmophic in unspecified populations will be found." << endl;
	cout << "--snp-num\t" << "Number of population-specific SNPs to be saved [default: all SNPs specific to the specified population(s)]." << endl;
	cout << "--freq\t\t" << "Calculate and output allele frequencies in the reference populations." << endl;
	cout << "--out\t\t" << "Output file for saving population-specific SNPs or allele frequencies [default: psnps]." << endl;
	cout << "--maf-min\t" << "Exclude SNPs with MAFs smaller than the specified value in the population(s) specified by --pop [default: 0]." << endl;
	cout << "\t\t" << "Minor alleles are determined by the total samples of the reference populations." << endl;
	cout << "--maf-max\t" << "Exclude SNPs with MAFs larger than the specified value in the population(s) specified by --pop [default: 0.5]." << endl;
	cout << "--miss-max\t" << "Exclude SNPs with missing rates larger than the specified value in the reference populations [default: 0]." << endl;
	cout << "--batch-size\t" << "Number of SNPs to be processed in a batch [default: 10000]." << endl;
	cout << "--thread-num\t" << "Number of threads on which the program will be running [default: thread number in your machine - 1]." << endl;
}

Dataset(int argc, char *argv[]) 
{
	int opt, longindex;
	const char *optstring = "";
	struct option long_options[] =
	{
		{ "bfile", required_argument,  NULL, 0},
		{ "ref-pop", required_argument,  NULL, 1},
		{ "out", required_argument, NULL, 2},
		{ "pop", required_argument, NULL, 3},
		{ "maf-min", required_argument, NULL, 4},
		{ "maf-max", required_argument, NULL, 5},
		{ "miss-max", required_argument, NULL, 6},
		{ "batch-size", required_argument, NULL, 7},
		{ "freq", no_argument, NULL, 8},
		{ "snp-num", required_argument, NULL, 9},
		{ "thread-num", required_argument, NULL, 10},
		{ "help", no_argument, NULL, 11},
		{0, 0, 0, 0}
	};

	while ((opt = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1)
	{
		switch (opt)
		{
		case 0: flag[0] = 1; arg[0] = optarg; break;
		case 1: flag[1] = 1; arg[1] = optarg; break;
		case 2: flag[2] = 1; arg[2] = optarg; break;
		case 3: 
		{
			flag[3] = 1;
			pop_name.push_back(optarg); 
			for (int i = optind; i < argc; i++) 
			{
			if (argv[i][0] == '-') break;
			else pop_name.push_back(argv[i]);
			}
			pop_num = pop_name.size();
			break;
		}
		case 4: flag[4] = 1; arg[4] = optarg; break;
		case 5: flag[5] = 1; arg[5] = optarg; break;
		case 6: flag[6] = 1; arg[6] = optarg; break;
		case 7: flag[7] = 1; arg[7] = optarg; break;
		case 8: flag[8] = 1; break;
		case 9: flag[9] = 1; arg[9] = optarg; break;
		case 10: flag[10] = 1; arg[10] = optarg; break;
		case 11:  help(); exit(0); break;
		default: break;
		}
	}

	if (flag[2] == 1) out_file = arg[2];
	else out_file = "psnps";
	logFile.open(out_file + ".log", ios::out);
	cout << "Population-Specific SNP Screener (PSNPS)" << endl;
	cout << "Options specified:" << endl;
	logFile << "Population-Specific SNP Screener (PSNPS)" << endl;
	logFile << "Options specified:" << endl;
	for (int i = 0; i < 3; i++) 
	{
		if (flag[i] == 1) 
		{
			cout << "--" << long_options[i].name << " " << arg[i] << endl;
			logFile << "--" << long_options[i].name << " " << arg[i] << endl;
		}
	}
	if (flag[3] == 1) 
	{
		cout << "--" << long_options[3].name;
		logFile << "--" << long_options[3].name;
		for (int i = 0; i < pop_num; i++) 
		{
			cout << " " << pop_name[i];
			logFile << " " << pop_name[i];
		}
		cout << endl;
		logFile << endl;
	}
	for (int i = 4; i < 8; i++)
	{
		if (flag[i] == 1) 
		{
			cout << "--" << long_options[i].name << " " << arg[i] << endl;
			logFile << "--" << long_options[i].name << " " << arg[i] << endl;
		}
	}
	if (flag[8] == 1) 
	{
		cout << "--" << long_options[8].name << endl;
		logFile << "--" << long_options[8].name << endl;
	}
	for (int i = 9; i < 11; i++)
	{
		if (flag[i] == 1) 
		{
			cout << "--" << long_options[i].name << " " << arg[i] << endl;
			logFile << "--" << long_options[i].name << " " << arg[i] << endl;
		}
	}

	if (flag[0] == 1)
	{
		bed_file = arg[0] + ".bed";
		bedfile.open(bed_file.c_str(), ios::in | ios::binary);
		if(!bedfile)
		{
			cout << "ERROR: " << bed_file << " was not found!" << endl;
			logFile << "ERROR: " << bed_file << " was not found!" << endl;
			error_flag = true;
		}
		else bedfile.close();
		fam_file = arg[0] + ".fam";
		famfile.open(fam_file.c_str(), ios::in);
		if(!famfile)
		{
			cout << "ERROR: " << fam_file << " was not found!" << endl;
			logFile << "ERROR: " << fam_file << " was not found!" << endl;
			error_flag = true;
		}
		else famfile.close();
		bim_file = arg[0] + ".bim";
		bimfile.open(bim_file.c_str(), ios::in);
		if(!bimfile)
		{
			cout << "ERROR: " << bim_file << " was not found!" << endl;
			logFile << "ERROR: " << bim_file << " was not found!" << endl;
			error_flag = true;
		}
		else bimfile.close();
	}
	else
	{
		cout << "ERROR: " << "use --bfile to specify genotype data." << endl;
		logFile << "ERROR: " << "use --bfile to specify genotype data." << endl;
		error_flag = true;
	}
	if (flag[1] == 1)
	{
		ref_file = arg[1];
		reffile.open(ref_file.c_str(), ios::in);
		if(!reffile)
		{
			cout << "ERROR: " << ref_file << " was not found!" << endl;
			logFile << "ERROR: " << ref_file << " was not found!" << endl;
			error_flag = true;
		}
		else reffile.close();
	}
	else
	{
		cout << "ERROR: " << "use --ref-pop to specify the file that describes reference populations." << endl;
		logFile << "ERROR: " << "use --ref-pop to specify the file that describes reference populations." << endl;
		error_flag = true;
	}
	if (flag[3] == 1)
	{
		for (int i = 0; i < pop_num; i++) 
		{
			for (int j = i + 1; j < pop_num; j++) 
			{
				if (pop_name[i] == pop_name[j])
				{
					cout << "ERROR: duplicated population " << pop_name[i] << " was specified by --pop." << endl;
					logFile << "ERROR: duplicated population " << pop_name[i] << " was specified by --pop." << endl;
					error_flag = true;
				}
			}
		}
	}
	if (flag[4] == 1)
	{
		maf_min =  atof(arg[4].c_str());
		if (maf_min < 0 || maf_min > 1)
		{
			cout << "ERROR: --maf-min should be between 0 and 1." << endl;
			logFile << "ERROR: --maf-min should be between 0 and 1." << endl;
			error_flag = true;
		}
	}
	else maf_min = 0.0;
	if (flag[5] == 1)
	{
		maf_max =  atof(arg[5].c_str());
		if (maf_max < 0 || maf_max > 1)
		{
			cout << "ERROR: --maf-max should be between 0 and 1." << endl;
			logFile << "ERROR: --maf-max should be between 0 and 1." << endl;
			error_flag = true;
		}
	}
	else maf_max = 0.5;
	if (flag[4] == 1 && flag[5] == 1 && maf_min >= maf_max)
	{
		cout << "ERROR: --maf-min should be smaller than --maf-max." << endl;
		logFile << "ERROR: --maf-min should be smaller than --maf-max." << endl;
		error_flag = true;
	}
	if (flag[6] == 1)
	{
		miss_max =  atof(arg[6].c_str());
		if (miss_max < 0 || miss_max > 1)
		{
			cout << "ERROR: --miss-max should be between 0 and 1." << endl;
			logFile << "ERROR: --miss-max should be between 0 and 1." << endl;
			error_flag = true;
		}
	}
	else miss_max = 0;
	if (flag[7] == 1)
	{
		batch_size =  atoi(arg[7].c_str());
		if (batch_size <= 0)
		{
			cout << "ERROR: --batch-size should be larger than 0." << endl;
			logFile << "ERROR: --batch-size should be larger than 0." << endl;
			error_flag = true;
		}
	}
	else batch_size = 10000;
	if (flag[9] == 1)
	{
		out_snp_num =  atoi(arg[9].c_str());
		if (out_snp_num <= 0)
		{
			cout << "ERROR: --snp-num should be larger than 0." << endl;
			logFile << "ERROR: --snp-num should be larger than 0." << endl;
			error_flag = true;
		}
	}
	max_thread_num = std::thread::hardware_concurrency();
	if (flag[10] == 1)
	{
		thread_num =  atoi(arg[10].c_str());
		if (thread_num < 1 || thread_num > max_thread_num)
		{
			cout << "ERROR: --thread-num should be between 1 and " << max_thread_num << "." << endl;
			logFile << "ERROR: --thread-num should be between 1 and " << max_thread_num << "." << endl;
			error_flag = true;
		}
	}
	else if (max_thread_num > 1) thread_num = max_thread_num - 1;
	else thread_num = 1;
	if (flag[3] == 0 && flag[8] == 0)
	{
		cout << "ERROR: use --pop to screen population-specific SNPs or --freq to calculate allele frequencies in the reference populations." << endl;
		logFile << "ERROR: use --pop to screen population-specific SNPs or --freq to calculate allele frequencies in the reference populations." << endl;
		error_flag = true;
	}

	if (error_flag)
	{
		logFile.close();
		exit(0);
	}
}

~Dataset() {logFile.close();}

void readFAM()
{
	string thisline;
	vector<string> parsedline;
	cout << "Loading FAM file [" + fam_file + "]... " << flush;
	logFile << "Loading FAM file [" + fam_file + "]... " << flush;
	famfile.open(fam_file.c_str(), ios::in);
	while (getline(famfile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		famid.push_back(parsedline[0]);
		subid.push_back(parsedline[1]);
		subid_set.insert(parsedline[1]);
	}
	famfile.close();
	indi_num = subid.size();
	if (indi_num != subid_set.size())
	{
		cout << "ERROR: duplicated individual ID found in [" + fam_file + "]." << endl;
		logFile << "ERROR: duplicated individual ID found in [" + fam_file + "]." << endl;
		exit(0);
	}
	cout << "done." << endl;
	logFile << "done." << endl;
	cout << indi_num << " individuals loaded from FAM file [" + fam_file + "]." << endl;
	logFile << indi_num << " individuals loaded from FAM file [" + fam_file + "]." << endl;
}

void readBIM()
{
	string thisline;
	vector<string> parsedline;
	cout << "Loading BIM file [" + bim_file + "]... " << flush;
	logFile << "Loading BIM file [" + bim_file + "]... " << flush;
	snp_num = 0;
	bimfile.open(bim_file.c_str(), ios::in);
	while (getline(bimfile, thisline, '\n')) snp_num ++;
	chrom.resize(snp_num);
	snp.resize(snp_num);
	position.resize(snp_num);
	allele1.resize(snp_num);
	allele2.resize(snp_num);
	int i = 0;
	bimfile.clear();
	bimfile.seekg(0, ios::beg);
	while (getline(bimfile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		chrom[i] = parsedline[0];
		snp[i] = parsedline[1];
		position[i] = atoi((parsedline[3]).c_str());
		allele1[i] = parsedline[4];
		allele2[i] = parsedline[5];
		i ++;
	}
	bimfile.close();
	snp_num = snp.size();
	cout << "done." << endl;
	logFile << "done." << endl;
	cout << snp_num << " SNPs loaded from BIM file [" + bim_file + "]." << endl;
	logFile << snp_num << " SNPs loaded from BIM file [" + bim_file + "]." << endl;
}

void readREF()
{
	string thisline;
	vector<string> parsedline;
	set<string>::iterator tmpiter;
	cout << "Loading reference-population file [" + ref_file + "]... " << flush;
	logFile << "Loading reference-population file [" + ref_file + "]... " << flush;
	reffile.open(ref_file.c_str(), ios::in);
	while (getline(reffile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		ref_id.push_back(parsedline[0]);
		ref_id_set.insert(parsedline[0]);
		ref_pop.push_back(parsedline[1]);
		ref_pop_set.insert(parsedline[1]);
		ref_pop_multiset.insert(parsedline[1]);
	}
	reffile.close();
	ref_indi_num = ref_id.size();
	if (ref_indi_num != ref_id_set.size())
	{
		cout << "ERROR: duplicated individual ID found in [" + ref_file + "]." << endl;
		logFile << "ERROR: duplicated individual ID found in [" + ref_file + "]." << endl;
		exit(0);
	}
	for (tmpiter = ref_pop_set.begin(); tmpiter != ref_pop_set.end(); tmpiter++) 
	{
		ref_pop_name.push_back(*tmpiter);
		ref_pop_cnt.push_back(ref_pop_multiset.count(*tmpiter));
	}
	ref_pop_num = ref_pop_name.size();
	cout << "done." << endl;
	logFile << "done." << endl;
	cout <<  ref_pop_num << " reference populations and " << ref_indi_num << " individuals loaded from the reference population file [" + ref_file + "]." << endl;
	logFile << ref_pop_num << " reference populations and " << ref_indi_num << " individuals loaded from the reference population file [" + ref_file + "]." << endl;
	cout << "Population (Sample Size):" << endl;
	logFile << "Population (Sample Size):" << endl;
	for (int i = 0; i < ref_pop_num; i++) 
	{
		cout << ref_pop_name[i] << " (" << ref_pop_cnt[i] << ")" << endl;
		logFile << ref_pop_name[i] << " (" << ref_pop_cnt[i] << ")" << endl;
	}
	for (int i = 0; i < ref_indi_num; i++) 
	{
		if (subid_set.find(ref_id[i]) == subid_set.end())
		{
			cout << "ERROR: individual " << ref_id[i] << " in the reference-population file [" + ref_file + "] was not found in FAM file [" << fam_file << "]." << endl;
			logFile << "ERROR: individual " << ref_id[i] << " in the reference-population file [" + ref_file + "] was not found in FAM file [" << fam_file << "]." << endl;
			exit(0);
		}
		int pop_tmp = -1;
		for (int j = 0; j < ref_pop_num; j++)
		{
			if (ref_pop[i] == ref_pop_name[j]) pop_tmp = j;
		}
		id_pop_map.insert(make_pair(ref_id[i], pop_tmp));
	}
	id_pop.resize(indi_num);
	map<string, int>::iterator it;
	for (int i = 0; i < indi_num; i++) 
	{
		id_pop[i] = -1;
		it =  id_pop_map.find(subid[i]);
		if (it != id_pop_map.end()) id_pop[i] = it->second;
	}
	pop_pos.resize(pop_num);
	for (int i = 0; i < ref_pop_num; i++) 
	{
		bool pop_found = false;
		for (int j = 0; j < pop_num; j++)
		{
			if (pop_name[j] == ref_pop_name[i]) 
			{
				pop_found = true;
				pop_pos[j] = i;
			}
		}
		if (!pop_found)
		{
			nonpop_name.push_back(ref_pop_name[i]);
			nonpop_pos.push_back(i);
		}
	}
	nonpop_num = nonpop_name.size();
	for (int i = 0; i < pop_num; i++)
	{
		if (pop_pos[i] == -1)
		{
			cout << "ERROR: population " << pop_name[i] << " specified by --pop was not found in the [" << ref_file << "]." << endl;
			logFile << "ERROR: population " << pop_name[i] << " specified by --pop was not found in the [" << ref_file << "]." << endl;
			exit(0);
		}
	}
	subid_set.clear();
	ref_id_set.clear();
	ref_pop_set.clear();
	ref_pop_multiset.clear();
}

void allele_stat()
{
	allele1_cnt.resize(snp_num);
	maf.resize(snp_num);
	MA.resize(snp_num);
	miss_cnt.resize(snp_num);
	miss_frq.resize(snp_num);
	char tmp1char[1];
	bedfile.open(bed_file.c_str(), ios::in | ios::binary);
	cout << "Calculating allele frequencies... " << flush;
	logFile << "Calculating allele frequencies... " << flush;
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	int line_byte_size = ceil(double(indi_num) / 4);
	int batch_snp_size = batch_size;
	int batch_byte_size = batch_snp_size * line_byte_size;
	int batch_num = ceil(double(snp_num) / batch_snp_size);
	char *batch_char = new char[batch_byte_size];
	string batch_str;
	for (int batch_idx = 1; batch_idx <= batch_num; batch_idx ++ )
	{
		int snp_start = (batch_idx - 1) * batch_snp_size;
		int snp_end = batch_idx * batch_snp_size;
		if (batch_idx == batch_num)
		{
			snp_end = snp_num;
			batch_byte_size = (snp_end - snp_start) * line_byte_size;
		}
		bedfile.read(batch_char, batch_byte_size);
		batch_str.assign(batch_char, batch_byte_size);
		dynamic_bitset<unsigned char> batch_bit(batch_str.begin(), batch_str.end());

		#pragma omp parallel for
		for (int i = snp_start; i < snp_end; ++ i)
		{
			allele1_cnt[i].resize(ref_pop_num);
			maf[i].resize(ref_pop_num);
			miss_cnt[i].resize(ref_pop_num);
			miss_frq[i].resize(ref_pop_num);
			for (int j = 0; j < ref_pop_num; j++)
			{
				allele1_cnt[i][j] = 0;
				miss_cnt[i][j] = 0;
			}
			int pos = (line_byte_size) * (i - snp_start) * 8;
			for (int j = 0; j < indi_num; ++ j)
			{
				bool tmp1bool = batch_bit[pos ++ ];
				bool tmp2bool = batch_bit[pos ++ ];
				if (!tmp1bool && !tmp2bool)
				{
					for (int k = 0; k < ref_pop_num; k++)
					{
						if (id_pop[j] == k) allele1_cnt[i][k] += 2;
					}
					continue;
				}
				if (!tmp1bool && tmp2bool)
				{
					for (int k = 0; k < ref_pop_num; k++)
					{
						if (id_pop[j] == k) allele1_cnt[i][k] ++;
					}
					continue;
				}
				if (tmp1bool && !tmp2bool)
				{
					for (int k = 0; k < ref_pop_num; k++)
					{
						if (id_pop[j] == k) miss_cnt[i][k] ++;
					}
				}
			}
			double allele1_total = accumulate(allele1_cnt[i].begin(), allele1_cnt[i].end(), 0);
			int miss_total = accumulate(miss_cnt[i].begin(), miss_cnt[i].end(), 0);
			double A1freq = 0.0;
			if (miss_total != ref_indi_num) A1freq = allele1_total / double(ref_indi_num - miss_total) / 2.0;
			for (int k = 0; k < ref_pop_num; k++)
			{
				if(miss_cnt[i][k] == ref_pop_cnt[k])
            				{
					maf[i][k] = 0.0;
					miss_frq[i][k] = 1.0;
				}
				else
				{
					maf[i][k] = double(allele1_cnt[i][k]) / double(ref_pop_cnt[k] - miss_cnt[i][k]) / 2.0;
					miss_frq[i][k] = double(miss_cnt[i][k]) / ref_pop_cnt[k];
				}
				if(A1freq > 0.5)
				{
					MA[i] = allele2[i];
					if(miss_cnt[i][k] < ref_pop_cnt[k]) maf[i][k] = 1.0 - maf[i][k];
				}
				else
				{
					MA[i] = allele1[i];
				}
			}
		}
	}
	delete [] batch_char;
	bedfile.close();

	int low_MAF_num = 0;
	int high_MAF_num = 0;
	int miss_frq_num = 0;
	remain_snp_num = 0;
	snp_pop.resize(snp_num);
	for(int i = 0; i < snp_num; ++ i)
	{
		snp_pop[i] = true;
		for (int j = 0; j < pop_num; j++)
		{
			if(maf[i][pop_pos[j]] < maf_min)
			{
				low_MAF_num ++;
				snp_pop[i] = false;
			}
			if(maf[i][pop_pos[j]] > maf_max)
			{
				high_MAF_num ++;
				snp_pop[i] = false;
			}
		}
		bool miss_frq_flag = false;
		for (int j = 0; j < ref_pop_num; j++)
		{
			if (miss_frq[i][j] > miss_max) 
			{
				miss_frq_flag = true;
				snp_pop[i] = false;
			}
		}
		if (miss_frq_flag) miss_frq_num ++;
		if(snp_pop[i])  remain_snp_num++;
	}
	cout << "done." << endl;
	logFile << "done." << endl;
	cout << low_MAF_num << " SNPs removed due to MAFs < " << maf_min << ".\n";
	logFile << low_MAF_num << " SNPs removed due to MAFs < " << maf_min << ".\n";
	cout << high_MAF_num << " SNPs removed due to MAFs > " << maf_max << ".\n";
	logFile << high_MAF_num << " SNPs removed due to MAFs > " << maf_max << ".\n";	
	cout << miss_frq_num << " SNPs removed due to missing rate > " << miss_max << ".\n";
	logFile << miss_frq_num << " SNPs removed due to missing rate > " << miss_max << ".\n";
	if(remain_snp_num == 0)
	{
		cout << "ERROR: 0 SNPs pass filters." << endl;
		logFile << "ERROR: 0 SNPs pass filters." << endl;
		exit(0);
	}
	cout << remain_snp_num << " SNPs passed filters." << endl;
	logFile << remain_snp_num << " SNPs passed filters." << endl;
}

void pop_snp()
{
	pop_snp_num = 0;
	for(int i = 0; i < snp_num; ++ i)
	{
		if (snp_pop[i])
		{
			for (int j = 0; j < pop_num; j++)
			{
				if (maf[i][pop_pos[j]] == 0) snp_pop[i] =false;
			}
			for (int j = 0; j < nonpop_num; j++)
			{
				if (maf[i][nonpop_pos[j]] > 0) snp_pop[i] = false;
			}
			if (snp_pop[i])
			{
				pop_snp_idx.push_back(i);
				pop_snp_num++;
			}
		}
	}
	cout << pop_snp_num << " SNPs are specific to " << pop_name[0];
	for (int i = 1; i < pop_num; i++) cout << ", " << pop_name[i];
	cout << "." << endl;
	logFile << pop_snp_num << " SNPs are specific to " << pop_name[0];
	for (int i = 1; i < pop_num; i++) logFile << ", " << pop_name[i];
	logFile << "." << endl;
}

void allele_freq_output()
{
	ofstream outfile;
	outfile.open((out_file + ".frq").c_str(), ios::out);
	cout << "Saving allele frequencies to [" + out_file + ".frq]... " << flush;
	logFile << "Saving allele frequencies to [" + out_file + ".frq]... ";
	outfile << "chrom\t" << "pos\t" << "SNP\t" << "A1\t" << "A2\t" << "MA";
	for (int i = 0; i < ref_pop_num; i++) outfile << "\t" << ref_pop_name[i] << "_MAF\t"  << ref_pop_name[i] << "_miss_rate";
	outfile << endl;
	for(int i = 0; i<snp_num; ++ i)
	{
		outfile << chrom[i] << "\t";
		outfile << position[i] << "\t";
		outfile << snp[i] << "\t";
		outfile << allele1[i] << "\t";
		outfile << allele2[i] << "\t";
		outfile << MA[i];
		for (int j = 0; j < ref_pop_num; j++) outfile << "\t" << maf[i][j] << "\t" << miss_frq[i][j];
		outfile << endl;
	}
	cout << "done." << endl;
	logFile << "done." << endl;
	outfile.close();
}

void pop_snp_output()
{
	ofstream snpfile;
	string snp_file;
	snp_file = out_file + "." + pop_name[0];
	for (int i = 1; i < pop_num; i++) snp_file = snp_file + "." + pop_name[i];
	snpfile.open(snp_file.c_str(), ios::out);
	if (flag[9] == 1 && pop_snp_num > out_snp_num)
	{
		random_shuffle(pop_snp_idx.begin(), pop_snp_idx.end());
		cout << "Saving randomly selected " << out_snp_num << " population-specific SNPs to [" + snp_file + "]... " << flush;
		logFile << "Saving randomly selected " << out_snp_num << " population-specific SNPs to [" + snp_file + "]... ";
	}
	else
	{
		out_snp_num = pop_snp_num;
		cout << "Saving the " << out_snp_num << " population-specific SNPs to [" + snp_file + "]... " << flush;
		logFile << "Saving the " << out_snp_num << " population-specific SNPs to [" + snp_file + "]... ";
	}

	for(int i = 0; i < out_snp_num; i++)
	{
		snpfile << snp[pop_snp_idx[i]];
		for (int j = 0; j < pop_num; j++) snpfile << "\t" << pop_name[j] << "\t" << maf[pop_snp_idx[i]][pop_pos[j]];
		snpfile << endl;
	}
	cout << "done." << endl;
	logFile << "done." << endl;
	snpfile.close();
}
};

int main(int argc, char ** argv)
{
	time_t current_time;
	char * current_time_str;
	Dataset dat(argc, argv);
	stringstream ss;
	ss << dat.thread_num;
	setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
	omp_set_num_threads(dat.thread_num);
	cout << "The program will be running on " << dat.thread_num << " threads." << endl;
	dat.logFile << "The program will be running on " << dat.thread_num << " threads." << endl;
	current_time = time(0);
	current_time_str = ctime(&current_time);
	cout << "Start time: " << current_time_str << flush;
	dat.logFile << "Start time: " << current_time_str;
	dat.readFAM();
	dat.readBIM();
	dat.readREF();
	if (dat.flag[8] == 1 || dat.flag[3] == 1) dat.allele_stat();
	if (dat.flag[8] == 1) dat.allele_freq_output();
	if (dat.flag[3] == 1)	
	{
		dat.pop_snp();
		dat.pop_snp_output();
	}
	current_time = time(0);
	current_time_str = ctime(&current_time);
	cout << "End time: " << current_time_str << endl;
	dat.logFile << "End time: " << current_time_str << endl;
	return 0;
}
