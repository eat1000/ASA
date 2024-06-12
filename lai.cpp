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
#include <Eigen/Dense>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include <map>
#include <set>
#include <iomanip>
#include <zlib.h>
#include <unordered_set>
#include <getopt.h>
#include <unistd.h>
#include <fcntl.h>
#include <random>
#include <sam.h>
#include <vcf.h>
#include <tbx.h>
#include <hts.h>
#include <synced_bcf_reader.h>
#include <unordered_map>

using namespace boost;
using namespace Eigen;
using namespace std;

//#define SINGLE_PRECISION
#ifdef SINGLE_PRECISION
typedef ArrayXf eigenArray;
typedef float eigenVar;
#else
typedef ArrayXXd eigenArray;
typedef double eigenVar;
#endif

class LAI{
private:
int indi_num, snp_num, remain_snp_num;
const char* delimiter = " \t";
vector<string> famid, subid, snp, allele1, allele2;
vector<int> chrom, position, snp_to_delete, pop_chrom, pop_position;
ifstream infofile;
vector<string> info_snp, pop, info_pop, info_diff_pop, info_pop_MiA, pop_MiA, info_pop_MaA, pop_MaA; 
vector<eigenVar> pop_maf, info_pop_maf;
int info_snp_num, info_pop_num, not_in_info_num;
vector<int> info_in_bim, info_chrom, info_position; 

struct SNPinfo
{
    string info_snp, info_pop, info_pop_MiA, info_pop_MaA;
    eigenVar info_pop_maf;
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
        {
            return lhs.first < rhs.first;
        }
    }
};
multimap<pair<int, int>, SNPinfo, cmp> snp_map;
multimap<string, int> info_pop_snp;
int find_position(const std::vector<std::string>& vec, const std::string& target)
{
    auto iter = std::find(vec.begin(), vec.end(), target);
    if (iter != vec.end())
    {
        return iter - vec.begin();
    }
    else
    {
        return -1;
    }
};

struct Window_info
{
	int diff_pop_num, snp_start, snp_end, prev_start, prev_end, half_window_size;
	vector<int> snp_cnt;

	void window_move(vector<int>& info_chrom, vector<int>& info_position, vector<string>&info_pop, vector<string>&info_diff_pop, int snp_index, int info_snp_num)
	{
		prev_start = snp_start;
		prev_end = snp_end;

		for (int k = prev_end + 1; k < info_snp_num; k++)
		{
			if (k == info_snp_num - 1)
			{
				if (info_position[k] - info_position[snp_index] <= half_window_size)
				{
					for (int j = 0; j < diff_pop_num; j++)
					{
						if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
					}
					snp_end = k;
				}
				break;
			}
			if (info_chrom[k] != info_chrom[k + 1])
			{
				if (info_position[k] - info_position[snp_index] <= half_window_size)
				{
					for (int j = 0; j < diff_pop_num; j++)
					{
						if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
					}
					snp_end = k;
				}
				break;
			}
			if (info_position[k] - info_position[snp_index] <= half_window_size)
			{
				for (int j = 0; j < diff_pop_num; j++)
				{
					if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
				}
				snp_end = k;
			}
			else break;
		} 

		for (int k = prev_start; k < snp_index; k++)
		{
			if (info_position[snp_index] - info_position[k] > half_window_size)
			{
				for (int j = 0; j < diff_pop_num; j++)
				{
					if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] --;
				}
				snp_start = k + 1;
			}
			else break;
		} 
	}

	void window_ini(vector<int>& info_chrom, vector<int>& info_position, vector<string>&info_pop, vector<string>&info_diff_pop, int window_size, int info_snp_num)
	{
		diff_pop_num = info_diff_pop.size();
		snp_cnt.resize(diff_pop_num);
		snp_cnt.assign(diff_pop_num, 0);
		snp_start = 0;
		snp_end = 0;
		prev_start = 0;
		prev_end = 0;
		half_window_size = window_size / 2;

		for (int k = 0; k < info_snp_num; k++)
		{
			if (k == info_snp_num - 1)
			{
				for (int j = 0; j < diff_pop_num; j++)
				{
					if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
				}
				snp_end = k;
				break;
			}
			if (info_chrom[k] != info_chrom[k + 1])
			{
				for (int j = 0; j < diff_pop_num; j++)
				{
					if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
				}
				snp_end = k;
				break;
			}
			if (info_position[k] - info_position[snp_start] <= window_size)
			{
				for (int j = 0; j < diff_pop_num; j++)
				{
					if (info_pop[k] == info_diff_pop[j]) snp_cnt[j] ++;
				}
				snp_end = k;
			}
			else break;
		} 
	}
};

public:
string ofile = "lai", ifile, vfile;
bool match_alleles_flag = false, outfile_flag = false, ifile_flag = false, analysis_flag = false;
bool pop_allele_flag = false, window_size_flag = false, lai_flag = false, laiv2lai_flag = false, laiv_flag = false, vfile_flag = false;
string pop_allele_id, lgzfile;
double laiv_min = 0.000001;
int chr = 1, window_size = 2000000, gap_max = 1000000000, allele_min = 2;

ofstream writeLOG;

void cover()
{
	cout << "+======================================+" << endl;
	cout << "|                                      |" << endl;
	cout << "|    Local Ancestry Inference (LAI)    |" << endl;
	cout << "|    ASA version 1.2.0                 |" << endl;
	cout << "|                                      |" << endl;
	cout << "+======================================+" << endl;
}

void help()
{
	cover();
	cout << "--vfile\t\t" << "Prefix of haplotype file in bgzipped VCF format [.vcf.gz] with the associated index file in tbi or csi format [.vcf.gz.tbi] or [.vcf.gz.csi]." << endl;
	cout << "--ifile\t\t" << "Input file that defines a panel of population-specific SNPs. The file is expected to have seven columns without headers, which are in the order of SNP ID, population that the SNP is specific to, MAF in the population, minor and major alleles in reference populations, chromsome and position." << endl;
	cout << "--out\t\t" << "Prefix of output file [default: lai]." << endl;
	cout << "--chr\t\t" << "Specify the chromosome for the analysis [default: 1]." << endl;
	cout << "--pop-allele-id\t" << "Output populaiton-specific alleles carried by the speficied individual." << endl;
	cout << "--lai\t\t" << "Calculate and output local ancestral inference." << endl;
	cout << "--laiv\t\t" << "Output local ancestral information vectors [default: no output]." << endl;
	cout << "--laiv2lai\t" << "Read the specifed local ancestral infomation vector file and output the local ancestral inference." << endl;
	cout << "--match-alleles\t" << "Both minor and major alleles of the population-specific SNPs have to match the two alleles in the haplotype file [default: as least one allele of the population-specific SNPs has to match one of the two alleles in the haplotype file]." << endl;
	cout << "--window-size\t" << "Window size in the local ancestral information calculation [default: 2000000] (bp)." << endl;
	cout << "--laiv-min\t" << "Minimum value of laiv for calling ancestries [default: 0.000001]." << endl;
	cout << "--allele-min\t" << "Minimum number of alleles in the window for calling ancestries [default: 2]." << endl;
//	cout << "--gap-max\t" << "Maximum gap between two population-specific SNPs that is allowed in a window [default: 1000000] (bp)." << endl;
}

LAI(int argc, char *argv[]) 
{
	int opt, longindex;
	const char *optstring = "";
	struct option long_options[] =
	{
		{ "vfile", required_argument, NULL, 0},
		{ "ifile", required_argument,  NULL, 1},
		{ "out", required_argument, NULL, 2},
		{ "help", no_argument, NULL, 3},
		{ "match-alleles", no_argument, NULL, 4},
		{ "pop-allele-id", required_argument, NULL, 5},
		{ "window-size", required_argument, NULL, 6},
		{ "lai", no_argument, NULL, 7},
		{ "laiv-min", required_argument, NULL, 8},
		{ "chr", required_argument, NULL, 9},
		{ "laiv2lai", required_argument, NULL, 10},
//		{ "gap-max", required_argument, NULL, 11},
		{ "laiv", no_argument, NULL, 12},
		{ "allele-min", required_argument, NULL, 13},
		{0, 0, 0, 0}
    	};

	while ((opt = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1)
	{
		switch (opt)
		{
		case 0: vfile = optarg; vfile_flag = true; break;
		case 1: ifile = optarg; ifile_flag = true; break;
		case 2: ofile = optarg; outfile_flag = true; break;
		case 3: help(); exit(0); break;
		case 4: match_alleles_flag = true; break;
		case 5: pop_allele_id = optarg; pop_allele_flag = true; analysis_flag = true; break;
		case 6: window_size = atoi(optarg); window_size_flag = true; break; 
		case 7: lai_flag = true; analysis_flag = true; break;
		case 8: laiv_min = atof(optarg); break;
		case 9: chr = atoi(optarg); break;
		case 10: lgzfile = optarg; laiv2lai_flag = true; analysis_flag = true; break;
		case 11: gap_max = atoi(optarg); break;
		case 12: laiv_flag = true; break;
		case 13: allele_min = atoi(optarg); break;
		default: break;
		}
	}

	cover();

    writeLOG.open((ofile + ".log").c_str(), ios::out);
	cout << "Options specified:" << endl;
	writeLOG << "Options specified:" << endl;
	if (!analysis_flag)
	{
		cout << "ERROR: no analysis was specified, use --lai, --laiv2lai, or --pop-allele-id for an analysis." << endl;
		writeLOG << "ERROR: no analysis was specified, use --lai, --laiv2lai, or --pop-allele-id for an analysis." << endl;
		exit(0);
	}
	if (window_size_flag && window_size <= 0)
	{
		cout << "Error: --window-size should be larger than 0." << endl;
		writeLOG << "Error: --window-size should be larger than 0." << endl;
		exit(0);
	}
	if (lai_flag)
	{
		cout << "--lai " << endl;
		if (vfile_flag) cout << "--vfile " << vfile << endl;
		else
		{
			cout << "ERROR: used --vfile to specify the VCF file" << endl;
			writeLOG << "ERROR: used --vfile to specify the VCF file" << endl;
			exit(0);
		}
		if (ifile_flag) cout << "--ifile " << ifile << endl;
		else
		{
			cout << "ERROR: used --ifile to specify the file that defines a panel of population-specific SNPs." << endl;
			writeLOG << "ERROR: used --ifile to specify the file that defines a panel of population-specific SNPs" << endl;
			exit(0);
		}
		cout << "--out " << ofile << endl;
		cout << "--chr " << chr << endl;
		cout << "--window-size " << window_size << endl;
		if (match_alleles_flag) cout << "--match-alleles " << endl;
//        		cout << "--gap-max " << gap_max << endl;
		cout << "--laiv-min " << laiv_min << endl;
		cout << "--allele-min " << allele_min << endl;
		if (laiv_flag) cout << "--laiv " << endl;
		writeLOG << "--lai" << endl;
		writeLOG << "--vfile " << vfile << endl;
		writeLOG << "--ifile " << ifile << endl;
		writeLOG << "--out " << ofile << endl;
		writeLOG << "--chr " << chr << endl;
		writeLOG << "--window-size " << window_size << endl;
		if (match_alleles_flag) writeLOG << "--match-alleles " << endl;
//        		writeLOG << "--gap-max " << gap_max << endl;
		writeLOG << "--laiv-min " << laiv_min << endl;
		writeLOG << "--allele-min " << allele_min << endl;
		if (laiv_flag) writeLOG << "--laiv " << endl;
	}
	if (laiv2lai_flag)
	{			
		cout << "--laiv2lai " << lgzfile << endl;
		cout << "--out " << ofile << endl;
//		cout << "--gap-max " << gap_max << endl;
		cout << "--laiv-min " << laiv_min << endl;
		writeLOG << "--laiv2lai " << lgzfile << endl;
		writeLOG << "--out " << ofile << endl;
//		writeLOG << "--gap-max " << gap_max << endl;
		writeLOG << "--laiv-min " << laiv_min << endl;
	}
	if (pop_allele_flag)
	{
		cout << "--pop-allele-id " << pop_allele_id << endl;
		if (vfile_flag) cout << "--vfile " << vfile << endl;
		else
		{
			cout << "ERROR: used --vfile to specify the VCF file" << endl;
			writeLOG << "ERROR: used --vfile to specify the VCF file" << endl;
			exit(0);
		}
		if (ifile_flag) cout << "--ifile " << ifile << endl;
		else
		{
			cout << "ERROR: used --ifile to specify the file that defines a panel of population-specific SNPs." << endl;
			writeLOG << "ERROR: used --ifile to specify the file that defines a panel of population-specific SNPs" << endl;
			exit(0);
		}
		cout << "--out " << ofile << endl;
		if (match_alleles_flag) cout << "--match-alleles " << endl;
		writeLOG << "--pop-allele-id " << pop_allele_id.c_str() << endl;
		writeLOG << "--vfile " << vfile << endl;
		writeLOG << "--ifile " << ifile << endl;
		writeLOG << "--out " << ofile << endl;
		if (match_alleles_flag) writeLOG << "--match-alleles " << endl;
	}
	cout << endl;
	writeLOG << endl;
}
~LAI() {writeLOG.close();}

void vcf_read()
{
	string vcfName = (vfile + ".vcf.gz");
	htsFile *fp = hts_open(vcfName.c_str(), "r");
	if (!fp) 
	{
		cout << "ERROR: Failed to open [" << vcfName << "]." << endl;
		writeLOG << "ERROR: Failed to open [" << vcfName << "]." << endl;
		exit (0);
	}

	//input and check format
	
	int is_vcf = 0, is_bcf = 0;
	char *version;
	version = hts_format_description(hts_get_format(fp));
	is_vcf = hts_get_format(fp)->format == vcf;
	is_bcf = hts_get_format(fp)->format == bcf;
	if (is_vcf) 
	{
		cout << "[" << vcfName << "] was detected as being " << version << "." << endl;
		writeLOG << "[" << vcfName << "] was detected as being " << version << "." << endl;
	} 
	else if (is_bcf) 
	{
		cout << "[" << vcfName << "] was detected as being " << version << "." << endl;
		writeLOG << "[" << vcfName << "] was detected as being " << version << "." << endl;
	} 
	else 
	{
		cout << "ERROR: [" << vcfName << "] is of unsupported file format." << endl;
		writeLOG << "ERROR: [" << vcfName << "] is of unsupported file format." << endl;
		exit (0);
	}
	free(version); 
	cout << "Scanning haplotype file [" << vcfName << "]... " << flush;
	writeLOG << "Scanning haplotype file [" << vcfName << "]... " << flush;

	bcf_hdr_t *hdr = bcf_hdr_read(fp); 
	bcf1_t *rec = bcf_init1();
	indi_num = bcf_hdr_nsamples(hdr);
	for (int i = 0; i < indi_num; ++i) 
	{
        		const char *sampleId = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
       		if (sampleId != NULL) 
		{
            			subid.push_back(sampleId);
       		}
	}

	bool genotype_check = true;
	while (bcf_read(fp, hdr, rec) == 0) 
	{
		if (bcf_unpack(rec, BCF_UN_ALL) != 0)
		{
			cout << endl << "ERROR: Failed to unpack haplotype file [" << vcfName << "]." << flush;
			writeLOG << endl << "ERROR: Failed to unpack haplotype file [" << vcfName << "]." << flush;
			bcf_close(fp);
			exit (0);
		}
		if (genotype_check)
		{
			int n_fmt = 0;
			int* fmt = NULL;
			if (bcf_get_genotypes(hdr, rec, &fmt, &n_fmt) <= 0) 
			{
				cout << endl << "ERROR: Failed to find haplotype data in [" << vcfName << "]." << endl;
				writeLOG << endl << "ERROR: Failed to find haplotype data in [" << vcfName << "]." << endl;
				bcf_close(fp);
				exit(0);
			}
			genotype_check = false;
		}
		const char *v_chrom = bcf_hdr_id2name(hdr, rec->rid);
		int v_chrom_int = atoi(v_chrom);
		if (v_chrom_int != chr) continue;
		chrom.push_back(v_chrom_int);
		const char *v_snp = rec->d.id ? rec->d.id : "NoID";
		snp.push_back(string(v_snp));
		int v_pos = rec->pos + 1;
		position.push_back(v_pos);
		const char *v_ref = rec->d.allele[0];
		allele1.push_back(string(v_ref));
		const char *v_alt = rec->d.allele[1];
		allele2.push_back(string(v_alt));
		snp_to_delete.push_back(0);
	}
	snp_num = snp.size();
	bcf_hdr_destroy(hdr);
	bcf_destroy(rec);
	cout << "done." << endl;
	writeLOG << "done." << endl;
	cout << indi_num << " individuals and " << snp_num << " SNPs on chromosome " << chr << " were found in haplotype file [" << vcfName << "]." << endl;
	writeLOG << indi_num << " individuals and " << snp_num << " SNPs on chromosome " << chr << " were found in haplotype file [" << vcfName << "]." << endl;
}

void snp_info_input()
{
	int record_num;
	string thisline, this_snp;
	string vcfName = (vfile + ".vcf.gz");
	vector <string> parsedline;
	vector <int> info_pop_snp_num;
	//multimap<string, int> info_pop_snp;
	set<string> pop_set, info_pop_set;
	unordered_multiset<string> info_pop_multiset;

	ifstream infofile;
	infofile.open(ifile.c_str(),ios::in);
	if (!infofile)
	{
		cout << "ERROR: " << ifile << " was not found!" << endl;
		writeLOG << "ERROR: " << ifile << " was not found!" << endl;
		exit (0);
	}
	int info_col_num = 0;
	getline(infofile, thisline, '\n');
	split(parsedline, thisline, is_any_of(delimiter));
	info_col_num = parsedline.size();
	infofile.close(); 
	if (info_col_num < 7)
	{
		cout << "ERROR: " << ifile << " has " << info_col_num << " columns, 7 columns are expected." << endl;
		writeLOG << "ERROR: " << ifile << " has " << info_col_num << " columns, 7 columns are expected." << endl;
		exit(0);
	}

	infofile.open(ifile.c_str(),ios::in);
	while (getline(infofile,thisline,'\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		if (stoi(parsedline[5]) != chr) continue;
		SNPinfo snp_info;
		snp_info.info_snp = parsedline[0];
		snp_info.info_pop = parsedline[1];
		snp_info.info_pop_maf = atof(parsedline[2].c_str());
		snp_info.info_pop_MiA = parsedline[3];
		snp_info.info_pop_MaA = parsedline[4];
		snp_info.chrom = atoi(parsedline[5].c_str());
		snp_info.position = atoi(parsedline[6].c_str());
		snp_map.insert(make_pair(make_pair(snp_info.chrom, snp_info.position), snp_info));
	}
	infofile.close();

	info_snp_num = 0;
	for (auto& item : snp_map) 
	{
		this_snp = item.second.info_snp;
		info_pop_snp.insert(multimap<string, int>::value_type(this_snp, info_snp_num));
		info_snp.push_back(item.second.info_snp);
		info_pop_MiA.push_back(item.second.info_pop_MiA);
		info_pop_MaA.push_back(item.second.info_pop_MaA);
		info_pop.push_back(item.second.info_pop);
		info_pop_set.insert(item.second.info_pop);
		info_pop_multiset.insert(item.second.info_pop);
		info_pop_maf.push_back(item.second.info_pop_maf);
		info_chrom.push_back(item.second.chrom);
		info_position.push_back(item.second.position);
		info_snp_num++;
	}

	cout << info_snp_num << " SNPs on chromosome " << chr << " were loaded from the SNP panel file [" << ifile << "]." << endl;
	writeLOG << info_snp_num << " SNPs on chromosome " << chr << " were loaded from the SNP panel file [" << ifile << "]." << endl;
	info_pop_num = info_pop_set.size();
	cout << info_pop_num << " populations found, population ID, number of the population-specific SNPs (#) and SNP density (#/Mbp) are " << endl;
 	writeLOG << info_pop_num << " populations found, population ID, number of the population-specific SNPs (#) and SNP density (#/Mbp) are " << endl;
	set<string>::iterator it;
	for (it = info_pop_set.begin(); it != info_pop_set.end(); ++ it)
	{
		info_diff_pop.push_back(*it);
		info_pop_snp_num.push_back(info_pop_multiset.count(*it)); 
	}
	for (int i = 0; i < info_pop_num; ++ i)
	{
		cout << info_diff_pop[i] << ": " << info_pop_snp_num[i] << ", " << double(info_pop_snp_num[i] * 100 / ((info_position[info_snp_num -1] - info_position[0]) / 1000000)) / 100.0 << "\t";
		writeLOG << info_diff_pop[i] << ": " << info_pop_snp_num[i] << ", " << double(info_pop_snp_num[i] * 100 / ((info_position[info_snp_num -1] - info_position[0]) / 1000000)) / 100.0 << "\t";
	}
 	cout << endl;
	writeLOG << endl;

	for (int i = 0; i < info_snp_num; ++ i)
	{
		record_num = info_pop_snp.count(info_snp[i]);
		if (record_num > 1)
		{
			cout << "ERROR: " << record_num << " records of SNP " << info_snp[i] << " were found in [" << ifile << "]." << endl;
			writeLOG << "ERROR: " << record_num << " records of SNP " << info_snp[i] << " were found in [" << ifile << "]." << endl;
			exit(0);
		}
		if (info_pop_maf[i] < 0 || info_pop_maf[i] > 1)
		{
			cout << "ERROR: MAF value " << info_pop_maf[i] << " was found for SNP " << info_snp[i]  << "." << endl;
			writeLOG << "ERROR: MAF value " << info_pop_maf[i] << " was found for SNP " << info_snp[i] << "." << endl;
			exit(0);
		}
	}
    
	info_in_bim.resize(info_snp_num, 0);
	not_in_info_num = 0;
	remain_snp_num = 0;
	int diff_allele_num = 0, chr_pos_mismatch_num = 0;
	for (int i = 0; i < snp_num; ++ i)
	{
		if (info_pop_snp.count(snp[i]) == 0)
		{
			not_in_info_num ++;
			pop_maf.push_back(0);
			pop_MiA.push_back("NA");
			pop_MaA.push_back("NA");
			pop.push_back("NA");
			pop_chrom.push_back(0);
			pop_position.push_back(0);
			snp_to_delete[i] = 1;
			continue;
		}
		int idx = info_pop_snp.find(snp[i])->second;
		pop_maf.push_back(info_pop_maf[idx]);
		pop_MiA.push_back(info_pop_MiA[idx]);
		pop_MaA.push_back(info_pop_MaA[idx]);
		pop_chrom.push_back(info_chrom[idx]);
		pop_position.push_back(info_position[idx]);
		info_in_bim[idx] = 1;
		pop.push_back(info_pop[idx]);
		if (info_chrom[idx] != chrom[i] || info_position[idx] != position[i]) {chr_pos_mismatch_num++;}
		if (match_alleles_flag)
		{
			if (allele1[i] == pop_MiA[i])
			{
				if (allele2[i] != pop_MaA[i])
				{
					snp_to_delete[i] = 1;
					diff_allele_num++;
				}
			}
			else if (allele2[i] == pop_MiA[i])
			{
				if (allele1[i] != pop_MaA[i])
				{
					snp_to_delete[i] = 1;
					diff_allele_num++;
				}
			}
			else
			{
				snp_to_delete[i] = 1;
				diff_allele_num++;
			}
		}
		else
		{
			if (allele1[i] != pop_MiA[i] && allele1[i] != pop_MaA[i] && allele2[i] != pop_MiA[i] && allele2[i] != pop_MaA[i])
			{
				snp_to_delete[i] = 1;
				diff_allele_num++;
			}
		}
		if (!snp_to_delete[i])
		{
			pop_set.insert(info_pop[idx]);
			remain_snp_num++;
		}
	}

	if (not_in_info_num > 0)
	{
		cout << not_in_info_num << " SNPs in the haplotype file were excluded, which were not in the SNP panel." << endl;
		writeLOG << not_in_info_num << " SNPs in the haplotype file were excluded, which were not in the SNP panel." << endl;
	}
	if (remain_snp_num == 0)
	{
		cout << "ERROR: 0 SNP available for analysis." << endl;
		writeLOG << "ERROR: 0 SNP available for analysis." << endl;
		exit(0);
	}
	if (diff_allele_num > 0)
	{
		cout << diff_allele_num << " SNPs removed due to discordant alleles with those in reference populations." << endl;
		writeLOG << diff_allele_num << " SNPs removed due to discordant alleles with those in reference populations." << endl;
	}
	cout << "Haplotypes of " << remain_snp_num << " SNPs will be used for the analysis." << endl;
	writeLOG << "Haplotypes of " << remain_snp_num << " SNPs will be used for the analysis." << endl;
	if (chr_pos_mismatch_num >0)
	{
		cout << "Warning: " << chr_pos_mismatch_num << "SNPs have discordant chromosome and position coordinates in [" << ifile << "] and [" << vcfName << "], coordinates in [" << vcfName << "] will be used." << endl;
		writeLOG << "Warning: " << chr_pos_mismatch_num << "SNPs have discordant chromosome and position coordinates in [" << ifile << "] and [" << vcfName << "], coordinates in [" << vcfName << "] will be used." << endl;
	}
}

void pop_allele()
{
	ofstream writefile;
	string outputFilename = ofile + ".pop-allele";
	writefile.open(outputFilename.c_str(), ios::out);
	cout << "Saving population-specific alleles carried by [" << pop_allele_id << "] to [" << outputFilename << "]... " << flush;
	writeLOG << "Saving population-specific alleles carried by [" << pop_allele_id << "] to [" << outputFilename << "]... " << flush;
	writefile << "#population codes: ";
	for (int j = 0; j < info_pop_num; j++)
	{
		writefile << "\t" << info_diff_pop[j] << "=" << j + 1;
	}
	writefile << endl;
	writefile << "#SNP" << "\t" << "POP" << "\t" << "CHROM" << "\t" << "POSITION" << "\t" << "HAPLOTYPE1" << "\t" << "HAPLOTYPE2" << endl;

	vector<string> sample_ref, sample_alt;
	string vcfName = (vfile + ".vcf.gz");
	htsFile *fp = hts_open(vcfName.c_str(), "r");
	if (!fp) 
	{
		cout << endl << "ERROR: Failed to open [" << vcfName << "]." << endl;
		writeLOG << endl << "ERROR: Failed to open [" << vcfName << "]." << endl;
		exit(0);
	}
	int is_vcf = 0, is_bcf = 0;
	int miss_allele_num = 0;
	bcf_hdr_t *hdr = bcf_hdr_read(fp); 
	bcf1_t *rec = bcf_init1();
	int sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, pop_allele_id.c_str());
	if (sample_index  == -1)
	{
		cout << endl << "ERROR: Failed to find [" <<  pop_allele_id << "] in the haplotype file [" << vcfName << "]." << endl;
		writeLOG << endl << "ERROR: Failed to find [" <<  pop_allele_id << "] in the haplotype file [" << vcfName << "]." << endl;
		bcf_close(fp);
		exit(0);
	}
	int vcf_index = 0;
	while (bcf_read(fp, hdr, rec) == 0) 
	{
		if (snp_to_delete[vcf_index])
		{
			vcf_index ++;
			continue;
		}
		int n_fmt = 0;
		int* fmt = NULL;
		bcf_get_genotypes(hdr, rec, &fmt, &n_fmt);
		int hap1_allele = bcf_gt_allele(fmt[2*sample_index]);
		int hap2_allele = bcf_gt_allele(fmt[2*sample_index + 1]);
		if (allele1[vcf_index] == pop_MiA[vcf_index] || allele2[vcf_index] == pop_MaA[vcf_index])
		{
			if (hap1_allele >= 0) {hap1_allele = (hap1_allele + 1) % 2;}
			if (hap2_allele >= 0) {hap2_allele = (hap2_allele + 1) % 2;}
		}
		if (hap1_allele == 1 || hap2_allele == 1)
		{
			writefile << snp[vcf_index] << "\t" << find_position(info_diff_pop, pop[vcf_index]) + 1 << "\t" << pop_chrom[vcf_index] << "\t" << pop_position[vcf_index] << "\t";
			writefile << hap1_allele << "\t" << hap2_allele << endl;
		}

 		free(fmt);
		vcf_index ++;
	}

	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	writefile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}

void lai()
{
	cout << "Calculating local ancestral information vectors by the method of moment... " << flush;
	writeLOG << "Calculating local ancestral information vectors by the method of moment... " << flush;
	unordered_map<string, int> snp_map_lai;
	for (int i = 0; i < snp.size(); ++i) 
	{
		snp_map_lai[snp[i]] = i;
	}

	string vcfName = (vfile + ".vcf.gz");
	string idxName = (vfile + ".vcf.gz.csi");
	htsFile *fp = hts_open(vcfName.c_str(), "r");
	tbx_t *idx = tbx_index_load2(vcfName.c_str(), idxName.c_str());
	if (!idx) 
	{
		idxName = (vfile + ".vcf.gz.tbi");
		idx = tbx_index_load2(vcfName.c_str(), idxName.c_str());
		if (!idx) 
		{
			cout << endl << "ERROR: Failed to open index file of [" << vcfName << "]." << endl;
			writeLOG << endl << "ERROR: Failed to open index file of [" << vcfName << "]." << endl;
			bcf_close(fp);
			exit(0);
		}
	}
	bcf_hdr_t *hdr = bcf_hdr_read(fp); 
	bcf1_t *rec = bcf_init1();
	hts_itr_t *iter;

	ofstream writefile;
	writefile.open((ofile + ".lai").c_str(), ios::out);
	writefile << "#population codes: " << "\t" << "NO_CALL=0" ;
	for (int j = 0; j < info_pop_num; j++)
	{
		writefile << "\t" << info_diff_pop[j] << "=" << j + 1;
	}
	writefile << endl;
//	writefile << "#chr: " << chr << "\t" << "window-size: " << window_size << "\t" << "gap-max: " << gap_max << "\t" << "laiv-min: " << laiv_min << "allele-min: " << allele_min << endl;
	writefile << "#chr" << "\t" << "spos" << "\t" << "epos" << "\t" << "sgpos" << "\t" << "egpos" << "\t" << "n snps";
	for (int j = 0; j < indi_num; ++ j)
	{
		writefile << "\t" << subid[j] << ".0" << "\t" << subid[j] << ".1";
	}
	writefile << endl;

	ostringstream writefile2;
	writefile2 << setprecision(4);
	gzFile gzWrite;
	if (laiv_flag)
	{
		int level = 1; 
		int strategy = Z_DEFAULT_STRATEGY;
		gzWrite = gzopen((ofile + ".laiv.gz").c_str(), "wb");
		gzsetparams(gzWrite, level, strategy);

		writefile2 << "#reference-panel-population: ";
		for (int j = 0; j < info_pop_num; j++)
		{
			writefile2 << "\t" << info_diff_pop[j];
		}
		writefile2 << endl;
		writefile2 << "#chr: " << chr << "\t" << "window-size: " << window_size << endl;
		writefile2 << "#chromosome" << "\t" << "physical_position" << "\t" << "snp_id" << "\t" << "snp_index";
		for (int l = 0; l < indi_num; l ++)
		{
			for (int j = 0; j <info_pop_num; j ++)
			{
				writefile2 << "\t" << subid[l] << ":::hap1:::" << info_diff_pop[j];
			}
			for (int j = 0; j <info_pop_num; j ++)
			{
				writefile2 << "\t" << subid[l] << ":::hap2:::" << info_diff_pop[j];
			}
		}
		writefile2 << endl;
		gzputs(gzWrite, writefile2.str().c_str());
		writefile2.str("");
	}

	ofstream writefile3;
	writefile3.open((ofile + ".Q").c_str(), ios::out);
	writefile3 << "#LAI diploid global ancestry .Q format output" << endl;
	writefile3 << "#sample";
	for (int j = 0; j < info_pop_num; j ++)
	{
		writefile3 << "\t" << info_diff_pop[j];
	}
	writefile3 << endl;

	long long int bim_pos = 0LL;
	int index1, index2;
	double max1, max2;
	int  cache_size = 0, call_start = 0, call_end = 0, call_count = 0;
	bool change_flag;
	int snp_gap = 0;
	eigenArray hap1_d, hap1_n, hap1_laiv, hap2_d, hap2_n, hap2_laiv, total_d, total_n, total_laiv; 
	vector<int> hap1_lai, hap2_lai, hap1_prev_lai, hap2_prev_lai;
	hap1_d.setZero(indi_num, info_pop_num);
	hap1_n.setZero(indi_num, info_pop_num);
	hap1_laiv.setZero(indi_num, info_pop_num);
	hap2_d.setZero(indi_num, info_pop_num);
	hap2_n.setZero(indi_num, info_pop_num);
	hap2_laiv.setZero(indi_num, info_pop_num);
	total_d.setZero(indi_num, info_pop_num);
	total_n.setZero(indi_num, info_pop_num);
	total_laiv.setZero(indi_num, info_pop_num);
	hap1_lai.resize(indi_num, 0);
	hap2_lai.resize(indi_num, 0);
	hap1_prev_lai.resize(indi_num, 0);
	hap2_prev_lai.resize(indi_num, 0);

	Window_info window_info;
	window_info.window_ini(info_chrom,info_position, info_pop, info_diff_pop, window_size, info_snp_num);
	for (int i = 0; i < info_snp_num; i++)
	{
		window_info.window_move(info_chrom, info_position, info_pop, info_diff_pop, i, info_snp_num);
		for (int j = 0; j < info_diff_pop.size(); j++) 
		{
			if (window_info.snp_cnt[j] > cache_size) cache_size = window_info.snp_cnt[j];
		}
	}

	vector<vector<vector<int>>> hap1(info_pop_num, vector<vector<int>>(indi_num,  vector<int>(cache_size, 0)));
	vector<vector<vector<int>>> hap2(info_pop_num, vector<vector<int>>(indi_num,  vector<int>(cache_size, 0)));
	vector<vector<vector<bool>>> missing_bool1(info_pop_num, vector<vector<bool>>(indi_num, vector<bool>(cache_size, true)));
	vector<vector<vector<bool>>> missing_bool2(info_pop_num, vector<vector<bool>>(indi_num, vector<bool>(cache_size, true)));
	vector<vector<double>> maf(info_pop_num, vector<double>(cache_size, 0));
	vector<int> start_ptr(info_pop_num, 0);
	vector<int> end_ptr(info_pop_num, 0);

	for (int i = 0; i < info_snp_num; i++)
	{
		if (i == 0) window_info.window_ini(info_chrom,info_position, info_pop, info_diff_pop, window_size, info_snp_num);
		else window_info.window_move(info_chrom, info_position, info_pop, info_diff_pop, i, info_snp_num);

		for (int j = 0; j < info_pop_num; j++) 
		{
			for (int k = window_info.prev_start; k < window_info.snp_start; k++)
			{
				if (info_pop[k] == info_diff_pop[j])
				{
					for (int l = 0; l < indi_num; ++ l)
					{
						if (!missing_bool1[j][l][start_ptr[j]])
						{
							hap1_n(l, j) = hap1_n(l, j) - hap1[j][l][start_ptr[j]];
							hap1_d(l, j) = hap1_d(l, j) - maf[j][start_ptr[j]];
						}
						if (!missing_bool2[j][l][start_ptr[j]])
						{
							hap2_n(l, j) = hap2_n(l, j) - hap2[j][l][start_ptr[j]];
							hap2_d(l, j) = hap2_d(l, j) - maf[j][start_ptr[j]];
						}
					}
					start_ptr[j] ++;
					if (start_ptr[j] == cache_size) start_ptr[j] = 0;
				}
			}
		}

		for (int j = 0; j < info_pop_num; j++) 
		{
			for (int k = window_info.prev_end; k <= window_info.snp_end; k++)
			{
				if (info_pop[k] == info_diff_pop[j])
				{
					if (k == window_info.prev_end && window_info.prev_end != 0) continue;
					double this_maf = info_pop_maf[k];
					if (info_in_bim[k] != 1)
					{
						for (int l = 0; l < indi_num; ++ l)
						{
							hap1[j][l][end_ptr[j]] = 0;
							hap2[j][l][end_ptr[j]] = 0;
							missing_bool1[j][l][end_ptr[j]] = false;
							missing_bool2[j][l][end_ptr[j]] = false;
							hap1_d(l, j) = hap1_d(l, j) + this_maf;
							hap2_d(l, j) = hap2_d(l, j) + this_maf;
							total_d(l, j) = total_d(l, j) + 2 * this_maf;
						}
						maf[j][end_ptr[j]] = this_maf;
						end_ptr[j] ++;
						if (end_ptr[j] == cache_size) end_ptr[j] = 0;
						continue;
					}

					bim_pos = snp_map_lai[info_snp[k]];
					if (snp_to_delete[bim_pos])
					{
						for (int l = 0; l < indi_num; ++ l)
						{
							hap1[j][l][end_ptr[j]] = 0;
							hap2[j][l][end_ptr[j]] = 0;
							missing_bool1[j][l][end_ptr[j]] = false;
							missing_bool2[j][l][end_ptr[j]] = false;
							hap1_d(l, j) = hap1_d(l, j) + this_maf;
							hap2_d(l, j) = hap2_d(l, j) + this_maf;
							total_d(l, j) = total_d(l, j) + 2 * this_maf;
						}
						maf[j][end_ptr[j]] = this_maf;
						end_ptr[j] ++;
						if (end_ptr[j] == cache_size) end_ptr[j] = 0;
						continue;
					}

					bool alleles_recode = (allele1[bim_pos] == pop_MiA[bim_pos] || allele2[bim_pos] == pop_MaA[bim_pos]);							
					string snp_pos = to_string(chrom[bim_pos]) + ":" + to_string(position[bim_pos]) + "-"  + to_string(position[bim_pos]);
					iter = tbx_itr_querys(idx, snp_pos.c_str());
					int ret = 0;
					kstring_t ksbuf = {0, 0, 0};

					ret = tbx_itr_next(fp, idx, iter, &ksbuf);
					ret = vcf_parse1(&ksbuf, hdr, rec);
					bcf_unpack(rec, BCF_UN_ALL);
					int *gt_array = NULL;
					int ngt = 0;
					bcf_get_genotypes(hdr, rec, &gt_array, &ngt);

					for (int l = 0; l < indi_num; l ++)
					{
						int hap1_allele = bcf_gt_allele(gt_array[l * 2]);
						if (hap1_allele < 0)
						{
							missing_bool1[j][l][end_ptr[j]] = true;
							hap1[j][l][end_ptr[j]] = 0;
						}
						else if (hap1_allele == 0 || hap1_allele == 1)
						{
							if (alleles_recode) hap1_allele = (hap1_allele + 1) % 2; 
							missing_bool1[j][l][end_ptr[j]] = false;
							hap1[j][l][end_ptr[j]] = hap1_allele;
							hap1_n(l, j) = hap1_n(l, j) + hap1_allele;
							hap1_d(l, j) = hap1_d(l, j) + this_maf;
							total_n(l, j) = total_n(l, j) + hap1_allele;
							total_d(l, j) = total_d(l, j) + this_maf;
						}

						int hap2_allele = bcf_gt_allele(gt_array[l * 2 + 1]);
						if (hap2_allele < 0)
						{
							missing_bool2[j][l][end_ptr[j]] = true;
							hap2[j][l][end_ptr[j]] = 0;
						}
						else  if (hap2_allele == 0 || hap2_allele == 1)
						{
							if (alleles_recode) hap2_allele = (hap2_allele + 1) % 2;
							missing_bool2[j][l][end_ptr[j]] = false;
							hap2[j][l][end_ptr[j]] = hap2_allele;
							hap2_n(l, j) = hap2_n(l, j) + hap2_allele;
							hap2_d(l, j) = hap2_d(l, j) + this_maf;
							total_n(l, j) = total_n(l, j) + hap2_allele;
							total_d(l, j) = total_d(l, j) + this_maf;
						}
					}
					maf[j][end_ptr[j]] = this_maf;
					free(gt_array);
					end_ptr[j] ++;
					if (end_ptr[j] == cache_size) end_ptr[j] = 0;
				}
			}
		}

		if (laiv_flag) writefile2 << info_chrom[i] << "\t" << info_position[i] << "\t" <<  info_snp[i] << "\t" << i;

//		hap1_laiv = hap1_n / hap1_d;
//		hap2_laiv = hap2_n / hap2_d;
		for (int l = 0; l < indi_num; l++)
		{
			for (int j = 0; j < info_pop_num; j++) 
			{
				if (window_info.snp_cnt[j] < 1 || (hap1_n(l, j) < allele_min && hap1_n.row(l).sum() > hap1_n(l, j)))
				{
					hap1_laiv(l, j) = 0;
				}
				else
				{
					hap1_laiv(l, j) = hap1_n(l, j) / hap1_d(l, j);
				}
				if (window_info.snp_cnt[j] < 1 || (hap2_n(l, j) < allele_min && hap2_n.row(l).sum() > hap2_n(l, j)))
				{
					hap2_laiv(l, j) = 0;
				}
				else
				{
					hap2_laiv(l, j) = hap2_n(l, j) / hap2_d(l, j);
				}
			}
			if (laiv_flag)
			{
				for (int j = 0; j < info_pop_num; ++j) writefile2 << "\t" << hap1_laiv(l, j);
				for (int j = 0; j < info_pop_num; ++j) writefile2 << "\t" << hap2_laiv(l, j);
			}
			max1 = hap1_laiv.row(l).maxCoeff(&index1);
			max2 = hap2_laiv.row(l).maxCoeff(&index2);
			hap1_lai[l] = 0;
			hap2_lai[l] = 0;
			if (max1 > laiv_min) hap1_lai[l] = index1 + 1;
			if (max2 > laiv_min) hap2_lai[l] = index2 + 1;			
		}
		if (laiv_flag)
		{
			writefile2 << "\n";
			gzputs(gzWrite, writefile2.str().c_str());
			writefile2.str("");
		}

		change_flag = false;
		if (i != 0) 
		{
			snp_gap = info_position[i] - info_position[i - 1];
//			if (snp_gap > gap_max || info_chrom[i] != info_chrom[i - 1]) change_flag = true;
			if (hap1_lai != hap1_prev_lai || hap2_lai != hap2_prev_lai) change_flag = true;
		}
		else call_start = 0;
		if (change_flag)
		{
			call_end = i - 1;
			call_count = call_end - call_start + 1;
			writefile << info_chrom[call_end] << "\t" << info_position[call_start] << "\t" << info_position[call_end] << "\t" << "." << "\t" << "." << "\t" << call_count;
			for (int l = 0; l < indi_num; ++ l) writefile << "\t" << hap1_prev_lai[l] << "\t" << hap2_prev_lai[l];
			writefile << "\n";
			call_start = i;        
		}
		hap1_prev_lai = hap1_lai;
		hap2_prev_lai = hap2_lai;
	}
	call_end = info_snp_num - 1;
	call_count = call_end - call_start + 1;
	writefile << info_chrom[call_end] << "\t" << info_position[call_start] << "\t" << info_position[call_end] << "\t" << "." << "\t" << "." << "\t" << call_count;
	for (int l = 0; l < indi_num; ++ l) writefile << "\t" << hap1_lai[l] << "\t" << hap2_lai[l];
	writefile << "\n";

	total_laiv =  total_n / total_d;
	for (int l = 0; l < indi_num; l ++)
	{
		writefile3 << subid[l];
		for (int j = 0; j < info_pop_num; j++) writefile3 << "\t" << total_laiv(l, j);
		for (int j = 0; j < info_pop_num; j++) writefile3 << "\t" << total_n(l, j) << "\t" << total_d(l, j);
		writefile3 << endl;
	}

	bcf_destroy(rec);
	tbx_itr_destroy(iter);
	bcf_hdr_destroy(hdr);
	tbx_destroy(idx);
	vcf_close(fp);
	writefile.close();
	writefile3.close();
	if (laiv_flag) gzclose(gzWrite);
	cout << "done." << endl;
	writeLOG << "done." << endl;
}

void laiv2lai()
{
	string gzFileName = lgzfile + ".laiv.gz";
	gzFile gzfp = gzopen(gzFileName.c_str(), "rb");
	ofstream writefile;
	writefile.open((ofile + ".lai").c_str(), ios::out);
    	if (gzfp == NULL) 
	{
		cout << "ERROR: Failed to open [" << gzFileName << "]." << endl;
		writeLOG << "ERROR: Failed to open [" << gzFileName << "]." << endl;
		exit(0);
	}

	char buffer[102400];
	int rowIndex = 0;
	int colIndex = 0;
	int gz_pop_num, gz_indi_num;
	vector<string> indi_id;
	vector<int> gz_chrom;
	vector<string> gz_pop;
	vector<int> gz_position, hap1_lai, hap2_lai, hap1_prev_lai, hap2_prev_lai;
	eigenArray hap1_laiv, hap2_laiv;

	while (gzgets(gzfp, buffer, sizeof(buffer)) != NULL) 
	{
		stringstream ss(buffer);
		string token;
		colIndex = 0;
		if (rowIndex == 0)
		{
			getline(ss, token, '\n');
			stringstream tmp1(token);
			string tmp2;
			while (getline(tmp1, tmp2, '\t'))
			{
				if (colIndex == 0) writefile << "#population codes: " << "\t" << "NO_POP=0";
				else
				{
					gz_pop.push_back(tmp2);
					writefile << "\t" << tmp2 << "=" << colIndex;
				}
				colIndex++;
			}
			writefile << endl;
			gz_pop_num  = gz_pop.size();
		}

//		if (rowIndex == 1)
//		{
//			getline(ss, token, '\n');
//			stringstream tmp1(token);
//			string tmp2;
//			while (getline(tmp1, tmp2, '\t')) 
//			{
//				if (colIndex == 0 || colIndex == 1)
//				{
//					writefile << tmp2 << "\t";
//				}
//				colIndex++;
//			}
//			writefile << "gap-max: " << gap_max << "\t" << "laiv-min: " << laiv_min << "allele-min: " << allele_min << endl;
//		}

		if (rowIndex == 1)
		{
			string prev_name;
			getline(ss, token, '\n');
			stringstream tmp1(token);
			string tmp2;
			while (getline(tmp1, tmp2, '\t')) 
			{
				if (colIndex > 3 && ((colIndex - 4) % (gz_pop_num *2) == 0))
				{
					string temp_str = tmp2;
					size_t temp_pos = temp_str.find(":::");
					string indi_name = temp_str.substr(0, temp_pos);
					indi_id.push_back(indi_name);
				}
				colIndex++;
			}
		}
		rowIndex++;
		if (rowIndex > 1 ) break;
	}

	gz_indi_num = indi_id.size();
	writefile << "#chr" << "\t" << "spos" << "\t" << "epos" << "\t" << "sgpos" << "\t" << "egpos" << "\t" << "n snps";
	for (int l = 0; l < gz_indi_num; l++) writefile << "\t" << indi_id[l] << ".0" << "\t" << indi_id[l] << ".1";
	writefile << endl;
	cout << gz_indi_num << " individuals were found in [" << gzFileName << "]. " << endl;
	writeLOG << gz_indi_num << " individuals were found in [" << gzFileName << "]. " << endl;
	cout << "Calling local ancestries with the local ancestral information vectors in [" << gzFileName << "]... " << flush;
	writeLOG << "Calling local ancestries with the local ancestral information vectors in [" << gzFileName << "]... " << flush;

	hap1_laiv.setZero(gz_indi_num, gz_pop_num);
	hap2_laiv.setZero(gz_indi_num, gz_pop_num);
	hap1_lai.resize(gz_indi_num, 0); 
	hap2_lai.resize(gz_indi_num, 0); 
	hap1_prev_lai.resize(gz_indi_num, 0); 
	hap2_prev_lai.resize(gz_indi_num, 0); 
	int i = 0, snp_gap = 0, call_start = 0, call_end = 0, call_count = 0;
	bool change_flag = false;
	double max1, max2;
	int index1, index2;

	while (gzgets(gzfp, buffer, sizeof(buffer)) != NULL) 
	{
		stringstream ss(buffer);
		string token;

		colIndex = 0;
		int count1 = 0, count2 = 0, indi_index = 0;
		bool hap1_flag = true, hap2_flag = false, ana_1 = true;
		while (getline(ss, token, '\t')) 
		{
			if (colIndex == 0) gz_chrom.push_back(stoi(token));
			else if (colIndex == 1) gz_position.push_back(stoi(token));
			else if (colIndex > 3)
			{
				if (hap1_flag && ana_1)
				{
					hap1_laiv(indi_index, count1) = stof(token);
					count1 ++;
					if (count1 == gz_pop_num)
					{
						hap1_flag = false;
						hap2_flag = true;
						count1 = 0;
					}
					ana_1 = false;
				}
				if (hap2_flag && ana_1)
				{
					hap2_laiv(indi_index, count2) = stof(token);
					count2 ++;
					if (count2 == gz_pop_num)
					{
						hap1_flag = true;
						hap2_flag = false;
						count2 = 0;
						indi_index ++;
					}
					ana_1 = false;
				}
			}
			ana_1 = true;
			colIndex++;
		}
		for (int l = 0; l < gz_indi_num; l++)
		{
			max1 = hap1_laiv.row(l).maxCoeff(&index1);
			max2 = hap2_laiv.row(l).maxCoeff(&index2);
			hap1_lai[l] = 0;
			hap2_lai[l] = 0;
			if (max1 > laiv_min) hap1_lai[l] = index1 + 1;
			if (max2 > laiv_min) hap2_lai[l] = index2 + 1;			
		}

		change_flag = false;
		if (i != 0) 
		{
			snp_gap = gz_position[i] - gz_position[i - 1];
//			if (snp_gap > gap_max || gz_chrom[i] != gz_chrom[i - 1]) change_flag = true;
			if (hap1_lai != hap1_prev_lai || hap2_lai != hap2_prev_lai) change_flag = true;
		}
		else call_start = 0;
		if (change_flag)
		{
			call_end = i - 1;
			call_count = call_end - call_start + 1;
			writefile << gz_chrom[call_end] << "\t" << gz_position[call_start] << "\t" << gz_position[call_end] << "\t" << "." << "\t" << "." << "\t" << call_count;
			for (int l = 0; l < gz_indi_num; ++ l) writefile << "\t" << hap1_prev_lai[l] << "\t" << hap2_prev_lai[l];
			writefile << "\n";
			call_start = i;        
		}
		hap1_prev_lai = hap1_lai;
		hap2_prev_lai = hap2_lai;
		i++;
	}
	call_end = i - 1;
	call_count = call_end - call_start + 1;
	writefile << gz_chrom[call_end] << "\t" << gz_position[call_start] << "\t" << gz_position[call_end] << "\t" << "." << "\t" << "." << "\t" << call_count;
	for (int l = 0; l < gz_indi_num; ++ l) writefile << "\t" << hap1_lai[l] << "\t" << hap2_lai[l];
	writefile << "\n";

	gzclose(gzfp);
	writefile.close();
	cout << "done." << endl;
	cout << "Local ancestry calls at " << i << " SNPs were saved to [" << ofile + ".lai" << "]. " << endl;
	writeLOG << "done." << endl;
	writeLOG << "Local ancestry calls at " << i << " SNPs were saved to [" << ofile + ".lai" << "]. " << endl;
}
};  

int main(int argc, char ** argv)
{
	time_t start_time, current_time;
	char * current_time_str;
	LAI lai(argc, argv);
	start_time = time(0);
	current_time_str = ctime(&start_time);
	cout << "Start time: " << current_time_str;
	lai.writeLOG << "Start time: " << current_time_str;
    	if (!lai.laiv2lai_flag)
    	{
		lai.vcf_read();
		lai.snp_info_input();
	}
	else lai.laiv2lai();
	if (lai.pop_allele_flag) lai.pop_allele();
	if (lai.lai_flag) lai.lai();
	current_time = time(0);
	current_time_str = ctime(&current_time);
	cout << "End time: " << current_time_str;
	lai.writeLOG << "End time: " << current_time_str;
	return 0;
}
