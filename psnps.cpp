﻿#include<iostream>
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
vector<string> pop_name, nonpop_name, famid, subid, chrom, snp, allele1, allele2, ref_id, ref_pop, ref_pop_name;
vector<int> position, ref_pop_cnt, id_pop, pop_pos, nonpop_pos, pop_snp_idx;
vector<bool> snp_pop;
vector<vector<int>> allele1_cnt, miss_cnt, pop_snp_selected;
vector<vector<double>> miss_frq, maf;
int pop_num, nonpop_num, batch_size, indi_num, snp_num, remain_snp_num, ref_indi_num, ref_pop_num, pop_snp_num, out_snp_num, max_thread_num, dist_min;
double maf_min, maf_max, miss_max;
bool error_flag = false;	
string bed_file, fam_file, bim_file, ref_file, out_file, panel_file;	
ifstream bedfile, famfile, bimfile, reffile, outfile;
set<string> subid_set, ref_pop_set, ref_id_set;
unordered_multiset<string> ref_pop_multiset;
map<string, int> id_pop_map;

public:
int flag[25]{};
string arg[25];
ofstream logFile;
int thread_num;

void cover()
{
	cout << "+================================================+" << endl;
	cout << "|                                                |" << endl;
	cout << "|    Population-Specific SNP Screener (PSNPS)    |" << endl;
	cout << "|    ASA version 1.2.0                           |" << endl;
	cout << "|                                                |" << endl;
	cout << "+================================================+" << endl;
}

void help()
{
	cover();
	cout << "Options:" << endl;
	cout << "--bfile\t\t" << "Input genotype file in plink binary format." << endl;
	cout << "--ref-pop\t" << "Input file that describes reference populations. The file is expected to have two columns without headers:" << endl;
	cout << "\t\t" << "the first is individual ID and the second is the population that the individual belongs to." << endl;
	cout << "--out\t\t" << "Output file for saving population-specific SNPs or allele frequencies [default: psnps]." << endl;
	cout << "--snp-panel\t" << "Screen SNPs that are specific to one reference population at a time and compose a SNP panel." << endl;
	cout << "--pop\t\t" << "Screen SNPs that are specific to the specified population(s)." << endl;
	cout << "\t\t" << "If multiple populations are specified, SNPs polymorphic in all specified populations and monmophic in unspecified populations will be found." << endl;
	cout << "--snp-num\t" << "Number of population-specific SNPs to be saved in each screening [default: all SNPs specific to the population(s)]." << endl;
	cout << "--random-seed\t" << "Set a random seed for selecting the population-specific SNPs to be saved." << endl;
	cout << "--freq\t\t" << "Calculate and output allele frequencies in the reference populations." << endl;
	cout << "--maf-min\t" << "Exclude SNPs with MAFs smaller than the specified value in the population(s) specified by --pop [default: 0]." << endl;
	cout << "\t\t" << "Minor alleles are determined by the total samples of the reference populations." << endl;
	cout << "--maf-max\t" << "Exclude SNPs with MAFs larger than the specified value in the population(s) specified by --pop [default: 0.5]." << endl;
	cout << "--miss-max\t" << "Exclude SNPs with missing rates larger than the specified value in the reference populations [default: 0]." << endl;
	cout << "--dist-min\t" << "Exclude SNPs with distances from previous ones less than or equal to the specified value [default: 0]." << endl;
	cout << "--tv-only\t" << "Keep SNPs that are transversions only." << endl;
	cout << "--batch-size\t" << "Number of SNPs to be processed in a batch [default: 10000]." << endl;
	cout << "--thread-num\t" << "Number of threads on which the program will be running [default: thread number in your machine - 1]." << endl;
	cout << "--sort\t\t" << "Sort the specified SNP panel file by chromosome and position." << endl;
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
		{ "dist-min", required_argument, NULL, 12},
		{ "random-seed", no_argument, NULL, 13},
		{ "tv-only", no_argument, NULL, 14},
		{ "snp-panel", no_argument, NULL, 15},
		{ "sort", required_argument, NULL, 16},
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
		case 12: flag[12] = 1; arg[12] = optarg; break;
		case 13: flag[13] = 1; break;
		case 14: flag[14] = 1; break;
		case 15: flag[15] = 1; break;
		case 16: flag[16] = 1; arg[16] = optarg; panel_file = optarg; break;
		default: break;
		}
	}

	cover();
	if (flag[2] == 1) out_file = arg[2];
	else out_file = "psnps";
	logFile.open(out_file + ".log", ios::out);
	cout << "Options specified:" << endl;
	logFile << "Options specified:" << endl;
	for (int i = 0; i < 3; i++) 
	{
		if (flag[i] == 1) 
		{
			cout << "--" << long_options[i].name << " " << arg[i] << endl;
			logFile << "--" << long_options[i].name << " " << arg[i] << endl;
		}
	}
	if (flag[15] == 1) 
	{
		cout << "--" << long_options[15].name << " " << endl;
		logFile << "--" << long_options[15].name << " " << endl;
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
	if (flag[16] == 1) 
	{
		cout << "--" << long_options[16].name << " " << arg[16] << endl;
		logFile << "--" << long_options[16].name << " " << arg[16] << endl;
	}
	for (int i = 4; i < 7; i++)
	{
		if (flag[i] == 1) 
		{
			cout << "--" << long_options[i].name << " " << arg[i] << endl;
			logFile << "--" << long_options[i].name << " " << arg[i] << endl;
		}
	}
	if (flag[7] == 1) 
	{
		cout << "--" << long_options[7].name << " " << arg[7] << endl;
		logFile << "--" << long_options[7].name << " " << arg[7] << endl;
	}
	if (flag[12] == 1) 
	{
		cout << "--" << long_options[12].name << " " << arg[12] << endl;
		logFile << "--" << long_options[12].name << " " << arg[12] << endl;
		dist_min =  atoi(arg[12].c_str());
		if (dist_min < 0)
		{
			cout << "ERROR: --dist-min should be larger or equal to 0." << endl;
			logFile << "ERROR: --dist-min should be larger or equal to 0." << endl;
			error_flag = true;
		}
	}
	else dist_min = 0;
	if (flag[14] == 1) 
	{
		cout << "--" << long_options[14].name << endl;
		logFile << "--" << long_options[14].name << endl;
	}
	if (flag[8] == 1) 
	{
		cout << "--" << long_options[8].name << endl;
		logFile << "--" << long_options[8].name << endl;
	}
	if (flag[9] == 1) 
	{
		cout << "--" << long_options[9].name << " " << arg[9] << endl;
		logFile << "--" << long_options[9].name << " " << arg[9] << endl;
	}
	if (flag[13] == 1) 
	{
		cout << "--" << long_options[13].name << endl;
		logFile << "--" << long_options[13].name << endl;
	}
	for (int i = 10; i < 11; i++)
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
	else if (flag[16] == 0)
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
	else if (flag[16] == 0)
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
	if (flag[3] == 0 && flag[8] == 0 && flag[15] == 0 && flag[16] == 0)
	{
		cout << "ERROR: use --pop to screen population-specific SNPs, --snp-panel to  screen population-specific SNPs and compose a SNP panel, or --freq to calculate allele frequencies in the reference populations." << endl;
		logFile << "ERROR: use --pop to screen population-specific SNPs, --snp-panel to  screen population-specific SNPs and compose a SNP panel, or --freq to calculate allele frequencies in the reference populations." << endl;
		error_flag = true;
	}
	if (flag[3] == 1 && flag[15] == 1)
	{
		cout << "ERROR: use --pop to screen population-specific SNPs or --snp-panel to  screen population-specific SNPs and compose a SNP panel." << endl;
		logFile << "ERROR: use --pop to screen population-specific SNPs or --snp-panel to  screen population-specific SNPs and compose a SNP panel." << endl;
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
	int chrom_this, pos_this;
	int i = 0, chrom_last = -1, pos_last = -1, sex_snp_num = 0, pos_del_num = 0, miss_allele_num = 0, ts_snp_num = 0;
	bool tv;
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
	snp_pop.resize(snp_num, true);
	bimfile.clear();
	bimfile.seekg(0, ios::beg);
	while (getline(bimfile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		chrom_this = atoi(parsedline[0].c_str());
		chrom[i] = parsedline[0];
		snp[i] = parsedline[1];
		pos_this = atoi(parsedline[3].c_str());
		position[i] = pos_this;
		allele1[i] = parsedline[4];
		allele2[i] = parsedline[5];
		if (allele1[i] == "*" || allele1[i] == "." || allele2[i] == "*" || allele2[i] == ".")
		{
			snp_pop[i] = false;
			miss_allele_num++;
		}
		if (flag[14] == 1)
		{
			tv = false;
			if (allele1[i] == "A" || allele1[i] == "G")
			{
				if (allele2[i] == "C" || allele2[i] == "T") tv = true;
			}
			else if (allele1[i] == "C" || allele1[i] == "T")
			{
				if (allele2[i] == "A" || allele2[i] == "G") tv = true;
			}
			if (!tv)
			{
				snp_pop[i] = false;
				ts_snp_num++;
			}
		}
		if (chrom_this < 1 || chrom_this > 22)
		{
			snp_pop[i] = false;
			sex_snp_num++;
		}
		if (chrom_this == chrom_last && pos_this == pos_last)
		{
			snp_pop[i] = false;
			pos_del_num++;
			if (snp_pop[i - 1])
			{
				snp_pop[i - 1] = false;
				pos_del_num++;
			}
			i ++;
			continue;
		}
		else if (chrom_this == chrom_last && pos_this - pos_last <= dist_min)
		{
			snp_pop[i] = false;
			pos_del_num++;
			i ++;
			continue;
		}
		chrom_last = chrom_this;
		pos_last = pos_this;
		i++;
	}
	bimfile.close();
	snp_num = snp.size();
	cout << "done." << endl;
	logFile << "done." << endl;
	cout << snp_num << " SNPs loaded from BIM file [" + bim_file + "]." << endl;
	logFile << snp_num << " SNPs loaded from BIM file [" + bim_file + "]." << endl;
	if (pos_del_num > 0)
	{
		cout << pos_del_num << " SNPs removed due to their distances from previous ones less then or equal to " << dist_min << "." << endl;
		logFile << pos_del_num << " SNPs removed due to their distances from previous ones less then or equal to " << dist_min << "." << endl;
	}
	if (sex_snp_num > 0)
	{
		cout << sex_snp_num << " SNPs on sex chromosomes removed." << endl;
		logFile << sex_snp_num << "SNPs on sex chromosomes removed." << endl;
	}
	if (miss_allele_num > 0)
	{
		cout << miss_allele_num << " SNPs removed due to missing allele information." << endl;
		logFile << miss_allele_num << " SNPs removed due to missing allele information." << endl;
	}
	if (flag[14] == 1 && ts_snp_num > 0) 
	{
		cout << ts_snp_num << " SNPs removed due to not being transversions." << endl;
		logFile << ts_snp_num << " SNPs removed due to not being transversions." << endl;
	}
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

	subid_set.clear();
	ref_id_set.clear();
	ref_pop_set.clear();
	ref_pop_multiset.clear();
}

void allele_stat()
{
	allele1_cnt.resize(snp_num);
	maf.resize(snp_num);
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
			string str_tmp;
			for (int j = 0; j < ref_pop_num; j++)
			{
				allele1_cnt[i][j] = 0;
				miss_cnt[i][j] = 0;
			}
			if (!snp_pop[i]) continue;
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
			if (miss_total != ref_indi_num) A1freq = double(allele1_total) / double(ref_indi_num - miss_total) / 2.0;
			if (A1freq > 0.5)
			{
				str_tmp = allele2[i];
				allele2[i] = allele1[i];
				allele1[i] = str_tmp;
			}
			for (int k = 0; k < ref_pop_num; k++)
			{
				if (miss_cnt[i][k] == ref_pop_cnt[k])
            				{
					maf[i][k] = 0.0;
					miss_frq[i][k] = 1.0;
				}
				else
				{
					maf[i][k] = double(allele1_cnt[i][k]) / double(ref_pop_cnt[k] - miss_cnt[i][k]) / 2.0;
					if (A1freq > 0.5) maf[i][k] = 1.0 - maf[i][k];
					miss_frq[i][k] = double(miss_cnt[i][k]) / double(ref_pop_cnt[k]);
				}
			}
		}
	}
	delete [] batch_char;
	bedfile.close();
	cout << "done." << endl;
	logFile << "done." << endl;
}


void allele_freq_output()
{
	ofstream outfile;
	outfile.open((out_file + ".frq").c_str(), ios::out);
	cout << "Saving allele frequencies to [" + out_file + ".frq]... " << flush;
	logFile << "Saving allele frequencies to [" + out_file + ".frq]... ";
	outfile << "chrom\t" << "pos\t" << "SNP\t" << "MinorAllele\t" << "MajorAllele";
	for (int i = 0; i < ref_pop_num; i++) outfile << "\t" << ref_pop_name[i] << "_MAF\t"  << ref_pop_name[i] << "_miss_rate";
	outfile << endl;
	for(int i = 0; i<snp_num; ++ i)
	{
		outfile << chrom[i] << "\t";
		outfile << position[i] << "\t";
		outfile << snp[i] << "\t";
		outfile << allele1[i] << "\t";
		outfile << allele2[i];
		for (int j = 0; j < ref_pop_num; j++) outfile << "\t" << maf[i][j] << "\t" << miss_frq[i][j];
		outfile << endl;
	}
	cout << "done." << endl;
	logFile << "done." << endl;
	outfile.close();
}

void pop_snp()
{
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

	int low_MAF_num = 0, high_MAF_num = 0, miss_frq_num = 0;
	remain_snp_num = 0;
	for(int i = 0; i < snp_num; ++ i)
	{
		for (int j = 0; j < pop_num; j++)
		{
			if (maf[i][pop_pos[j]] < maf_min)
			{
				low_MAF_num ++;
				snp_pop[i] = false;
			}
			if (maf[i][pop_pos[j]] > maf_max)
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
		if (snp_pop[i])  remain_snp_num++;
	}
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

	pop_snp_num = 0;
	for (int i = 0; i < snp_num; ++ i)
	{
		if (snp_pop[i])
		{
			for (int j = 0; j < pop_num; j++)
			{
				if (maf[i][pop_pos[j]] == 0) snp_pop[i] =false;
			}
			for (int j = 0; j < nonpop_num; j++)
			{
				if (maf[i][nonpop_pos[j]] > 0) snp_pop[i] = false; //==0
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

void pop_snp_output()
{
	ofstream snpfile;
	string snp_file;
	snp_file = out_file + "." + pop_name[0];
	for (int i = 1; i < pop_num; i++) snp_file = snp_file + "." + pop_name[i];
	snpfile.open(snp_file.c_str(), ios::out);
	if (flag[9] == 1 && pop_snp_num > out_snp_num)
	{
		if (flag[13] ==1) srand(time(0));
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

	// for(int i = 0; i < out_snp_num; i++)
	// {
	// 	snpfile << snp[pop_snp_idx[i]];
	// 	for (int j = 0; j < pop_num; j++) snpfile << "\t" << pop_name[j] << "\t" << maf[pop_snp_idx[i]][pop_pos[j]];
	// 	snpfile  << "\t" << allele1[pop_snp_idx[i]]  << "\t" << allele2[pop_snp_idx[i]] << endl;
	// }

		for(int i = 0; i < out_snp_num; i++)
	{
		snpfile << snp[pop_snp_idx[i]];
		for (int j = 0; j < pop_num; j++) snpfile << "\t" << pop_name[j] << "\t" << maf[pop_snp_idx[i]][pop_pos[j]];
		snpfile  << "\t" << allele1[pop_snp_idx[i]]  << "\t" << allele2[pop_snp_idx[i]];
		snpfile  << "\t" << chrom[pop_snp_idx[i]]  << "\t" << position[pop_snp_idx[i]] << endl;
	}

	cout << "done." << endl;
	logFile << "done." << endl;
	snpfile.close();
}

void  snp_panel()
{
	pop_snp_selected.resize(ref_pop_num);
	int miss_frq_num = 0;
	for(int i = 0; i < snp_num; ++ i)
	{
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
	}

	vector<bool> this_snp_pop;
	for (int pop_idx = 0; pop_idx < ref_pop_num; pop_idx++)
	{
	this_snp_pop = snp_pop;
	cout << "Screening SNPs specific to " << ref_pop_name[pop_idx] << "..." << flush;
	logFile << "Screening SNPs specific to " << ref_pop_name[pop_idx] << "..." << flush;
	int low_MAF_num = 0, high_MAF_num = 0;
	remain_snp_num = 0;
	for(int i = 0; i < snp_num; ++ i)
	{
		if (maf[i][pop_idx] < maf_min)
		{
			low_MAF_num ++;
			this_snp_pop[i] = false;
		}
		if (maf[i][pop_idx] > maf_max)
		{
			high_MAF_num ++;
			this_snp_pop[i] = false;
		}
		if (this_snp_pop[i])  remain_snp_num++;
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

	pop_snp_num = 0;
	pop_snp_idx.clear();
	for (int i = 0; i < snp_num; ++ i)
	{
		if (snp_pop[i])
		{
			if (maf[i][pop_idx] == 0) this_snp_pop[i] =false;
			for (int j = 0; j < ref_pop_num; j++)
			{
				if (j == pop_idx) continue;
				if (maf[i][j] > 0) this_snp_pop[i] = false; //==0
			}
			if (this_snp_pop[i])
			{
				pop_snp_idx.push_back(i);
				pop_snp_num++;
			}
		}
	}
	cout << pop_snp_num << " SNPs are specific to " << ref_pop_name[pop_idx] << "." << endl;
	logFile << pop_snp_num << " SNPs are specific to " << ref_pop_name[pop_idx] << "." << endl;
	if (flag[9] == 1 && pop_snp_num > out_snp_num)
	{
		if (flag[13] ==1) srand(time(0));
		random_shuffle(pop_snp_idx.begin(), pop_snp_idx.end());
		for (int i = 0; i < out_snp_num; ++ i) pop_snp_selected[pop_idx].push_back(pop_snp_idx[i]);
		cout << out_snp_num << " SNPs specific to " << ref_pop_name[pop_idx] << " were randomly selected." << endl;
		logFile << out_snp_num << " SNPs specific to " << ref_pop_name[pop_idx] << " were randomly selected." << endl;
	}
	else
	{
		for (int i = 0; i < pop_snp_num; ++ i) pop_snp_selected[pop_idx].push_back(pop_snp_idx[i]);
		cout << pop_snp_num << " SNPs specific to " << ref_pop_name[pop_idx] << " were selected." << endl;
		logFile << pop_snp_num << " SNPs specific to " << ref_pop_name[pop_idx] << " were selected." << endl;
	}
	}
}

void snp_panel_output()
{
	ofstream snpfile;
	string snp_file;
	int panel_snp_num = 0;
	snp_file = out_file + ".panel";
	snpfile.open(snp_file.c_str(), ios::out);
	cout << "Saving the SNP panel to [" + snp_file + "]... " << flush;
	logFile << "Saving the SNP panel to [" + snp_file + "]... " << flush;
	for (int pop_idx = 0; pop_idx < ref_pop_num; pop_idx++)
	{
		for (int i = 0; i < pop_snp_selected[pop_idx].size(); i++)
		{
			snpfile << snp[pop_snp_selected[pop_idx][i]];
			snpfile << "\t" << ref_pop_name[pop_idx] << "\t" << maf[pop_snp_selected[pop_idx][i]][pop_idx];
			snpfile  << "\t" << allele1[pop_snp_selected[pop_idx][i]]  << "\t" << allele2[pop_snp_selected[pop_idx][i]];
			snpfile  << "\t" << chrom[pop_snp_selected[pop_idx][i]] << "\t" << position[pop_snp_selected[pop_idx][i]] << endl;
			panel_snp_num++;
		}
	}
	cout << "done." << endl << panel_snp_num << " SNPs were saved." << endl;
	logFile << "done." << endl << panel_snp_num << " SNPs were saved." << endl;
	snpfile.close();
}

struct SNPinfo
{
 	string info_snp, info_pop, info_pop_MiA, info_pop_MaA;
	double info_pop_maf;
	int chrom, position;
};
struct cmp 
{
	bool operator() (const pair<int, int>& lhs, const pair<int, int>& rhs) const 
	{
		if (lhs.first == rhs.first) 
		{
			return lhs.second < rhs.second;
		} 
		else 
		{	return lhs.first < rhs.first;
		}
	}
};

void panel_sort()
{
	multimap<pair<int,int>,SNPinfo,cmp> snp_map;
	set<string> info_pop_set;
	unordered_multiset<string> info_pop_multiset;
	vector <string>  parsedline, info_diff_pop;
	vector <int> info_pop_snp_num;
	string thisline;
	ifstream infofile;
	infofile.open(panel_file.c_str(),ios::in);
	if(!infofile)
	{
		cout << "ERROR: " << panel_file << " was not found!" << endl;
		logFile << "ERROR: " << panel_file << " was not found!" << endl;
		exit (0);
	}
	int info_col_num = 0;
	getline(infofile, thisline, '\n');
	split(parsedline, thisline, is_any_of(delimiter));
	info_col_num = parsedline.size();
	infofile.close(); 
	if (info_col_num < 7)
	{
		cout << "ERROR: " << panel_file << " has " << info_col_num << " columns, 7 columns are expeted." << endl;
		logFile << "ERROR: " << panel_file << " has " << info_col_num << " columns, 7 columns are expeted." << endl;
		exit(0);
	}

	int info_snp_num = 0;
	infofile.open(panel_file.c_str(),ios::in);
	while (getline(infofile,thisline,'\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		SNPinfo snp_info;
		snp_info.info_snp = parsedline[0];
		snp_info.info_pop = parsedline[1];
		snp_info.info_pop_maf = atof(parsedline[2].c_str());
		snp_info.info_pop_MiA = parsedline[3];
		snp_info.info_pop_MaA = parsedline[4];
		snp_info.chrom = atoi(parsedline[5].c_str());
		snp_info.position = atoi(parsedline[6].c_str());
		snp_map.insert(make_pair(make_pair(snp_info.chrom, snp_info.position), snp_info));
		info_pop_set.insert(parsedline[1]);
		info_pop_multiset.insert(parsedline[1]);
		info_snp_num++;
	}

	infofile.close();

	cout << info_snp_num << " SNPs loaded from the panel file [" << panel_file << "]." << endl;
	logFile << info_snp_num << " SNPs loaded from the panel file [" << panel_file << "]." << endl;
	int info_pop_num = info_pop_set.size();
	cout << info_pop_num << " populations found, population ID and number of the population-specific SNPs are " << endl;
 	logFile << info_pop_num << " populations found, population ID and number of the population-specific SNPs are " << endl;
	set<string>::iterator it;
	for (it = info_pop_set.begin(); it != info_pop_set.end(); ++ it)
	{
		info_diff_pop.push_back(*it);
		info_pop_snp_num.push_back(info_pop_multiset.count(*it)); 
	}
	for(int i = 0; i < info_pop_num; ++ i)
	{
		cout << info_diff_pop[i] << ": " << info_pop_snp_num[i] << "\t";
		logFile << info_diff_pop[i] << ": " << info_pop_snp_num[i] << "\t";
	}
 	cout << endl;
	logFile << endl;

	string snp_file;
	ofstream snpfile;
	snp_file = out_file + ".panel";
	snpfile.open(snp_file.c_str(), ios::out);
	cout << "Saving the sorted SNP panel to [" + snp_file + "]... " << flush;
	logFile << "Saving the sorted SNP panel to [" + snp_file + "]... " << flush;
	info_snp_num = 0;


	for (auto& item : snp_map) 
	{
		snpfile << item.second.info_snp;
		snpfile << "\t" << item.second.info_pop << "\t" << item.second.info_pop_maf;
		snpfile  << "\t" << item.second.info_pop_MiA  << "\t" << item.second.info_pop_MaA;
		snpfile  << "\t" << item.second.chrom << "\t" << item.second.position << endl;
		info_snp_num++;
	}
	cout << "done." << endl << info_snp_num << " SNPs were saved." << endl;
	logFile << "done." << endl << info_snp_num << " SNPs were saved." << endl;
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
	if (dat.flag[16] == 1) 
	{
		dat.panel_sort();
		current_time = time(0);
		current_time_str = ctime(&current_time);
		cout << "End time: " << current_time_str << endl;
		dat.logFile << "End time: " << current_time_str << endl;
		return 0;
	}
	dat.readFAM();
	dat.readBIM();
	dat.readREF();
	if (dat.flag[8] == 1 || dat.flag[3] == 1 || dat.flag[15] == 1) dat.allele_stat();
	if (dat.flag[8] == 1) dat.allele_freq_output();
	if (dat.flag[3] == 1)	
	{
		dat.pop_snp();
		dat.pop_snp_output();
	}
	if (dat.flag[15] == 1)	
	{
		dat.snp_panel();
		dat.snp_panel_output();
	}
	current_time = time(0);
	current_time_str = ctime(&current_time);
	cout << "End time: " << current_time_str << endl;
	dat.logFile << "End time: " << current_time_str << endl;
	return 0;
}
