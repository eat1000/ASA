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
using namespace std;


class LAI{
private:
const char* delimiter = "\t";

public:
string ofile = "lai", ifile, suffix = ".lai";
bool outfile_flag = false, ifile_flag = false, fill_gaps = false, rfmix = false, palette_flag=false;
int from_chr = 1, to_chr = 22;
ofstream writeLOG;
vector<string> palette{"#FFFFFF", "#0072B2", "#F0E442", "#009E73", "#56B4E9", "#D55E00", "#CC79A7"};

void cover()
{
	cout << "+======================================+" << endl;
	cout << "|                                      |" << endl;
	cout << "|    karyoploteR and Tagore Interface  |" << endl;
	cout << "|    ASA version 1.2.0                 |" << endl;
	cout << "|                                      |" << endl;
	cout << "+======================================+" << endl;
}

void help()
{
	cover();
	cout << "--lai-files\t" << "Prefix of the LAI files." << endl;
	cout << "--from-chr\t" << "Starting chromosome of the LAI files [default: 1]." << endl;
	cout << "--to-chr\t" << "Ending chromosome of the LAI files [default: 22]." << endl;
	cout << "--out\t\t" << "Prefix of the output files for plotting with karyoploteR and Tagore [default: lai]." << endl;
//	cout << "--fill-gaps\t\t" << "Fill the uncalled gaps in the output files [default: false]." << endl;
	cout << "--palette\t" << "Palette for plotting with Tagore." << endl;
	cout << "--suffix\t" << "Suffix of the LAI files [default: .lai]." << endl;
	cout << "--rfmix\t\t" << "RFMix msp file format that uses 0 for population coding." << endl;
}

LAI(int argc, char *argv[]) 
{
	int opt, longindex;
	const char *optstring = "";
	struct option long_options[] =
	{
		{ "lai-files", required_argument, NULL, 0},
		{ "from-chr", required_argument, NULL, 1},
		{ "to-chr", required_argument, NULL, 2},
		{ "fill-gaps", no_argument, NULL, 3},
		{ "out", required_argument, NULL, 4},		
		{ "palette", required_argument, NULL, 5},
		{ "suffix", required_argument, NULL, 6},
		{ "rfmix", no_argument, NULL, 7},
		{ "help", no_argument, NULL, 8},
		{0, 0, 0, 0}
    	};

	while ((opt = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1)
	{
		switch (opt)
		{
		case 0: ifile = optarg; ifile_flag = true; break;
		case 1: from_chr = atoi(optarg); break;
		case 2: to_chr = atoi(optarg); break;
		case 3: fill_gaps = true; break;
		case 4: ofile = optarg; outfile_flag = true; break;
		case 5: 		
		{
			palette_flag = true;
			palette.clear();
			palette.push_back(optarg); 
			for (int i = optind; i < argc; i++) 
			{
				if (argv[i][0] == '-' || argv[i][0] == '\n') break;
				else palette.push_back(argv[i]);
			}
			break;
		}
		case 6: suffix = optarg; break;
		case 7: rfmix = true; break;
		case 8: help(); exit(0); break;
		default: break;
		}
	}
	cout << endl;
	writeLOG << endl;
	cover();

    writeLOG.open((ofile + ".log").c_str(), ios::out);
	cout << "Options specified:" << endl;
	writeLOG << "Options specified:" << endl;
	if (ifile_flag)
	{
		cout << "--lai-files " << ifile << endl;
		writeLOG << "--lai-files " << ifile << endl;
		cout << "--suffix " << suffix << endl;
		writeLOG << "--suffix " << suffix << endl;
	}
   	else
   	{
		cout << "ERROR: use --lai-files to specify the LAI file." << endl;
		writeLOG << "ERROR: use --lai-files to specify the LAI file." << endl;
		exit(0);
	}

	cout << "--from-chr " << from_chr << endl;
	writeLOG << "--from-chr " << from_chr << endl;
	cout << "--to-chr " << to_chr << endl;
	writeLOG << "--to-chr " << to_chr << endl;
//	cout << "--fill-gaps " << fill_gaps << endl;
//	writeLOG << "--fill-gaps " << fill_gaps << endl;
	cout << "--out " << ofile << endl;
	writeLOG << "--out " << ofile << endl;
	if (palette_flag)
	{
		cout << "--palette";
		writeLOG << "--palette";
		for (int i = 0; i < palette.size(); i++) 
		{
			palette[i] = "#" + palette[i];
			cout << " " << palette[i];
			writeLOG << " " << palette[i];
		}
		cout << endl;
		writeLOG << endl;
	}
	if (rfmix)
	{
		cout << "--rfmix" << endl;
		writeLOG << "--rfmix" << endl;
	}
	cout << endl;
	writeLOG << endl;
}
~LAI() {writeLOG.close();}



void LAI_input()
{
	int id_num, pop_num, region_num, chr_snp_num, chr_region_num;
	vector<int> chrom, spos, epos, tmp1, tmp2, nsnps, csidx, ceidx;
	vector<vector<int>> hap1, hap2; 
	vector<string> id;
	string thisline, line1, line2, line3;
	vector <string> parsedline, parsedline1, parsedline2, parsedline3;
	string LAIName = (ifile + "1" + suffix);
	ifstream LAIfile;
	LAIfile.open(LAIName.c_str(),ios::in);
	if (!LAIfile)
	{
		cout << "ERROR: " << LAIName << " was not found!" << endl;
		writeLOG << "ERROR: " << LAIName << " was not found!" << endl;
		exit (0);
	}
	getline(LAIfile, line1, '\n');
	cout << "Population labels found are:" << endl;
	writeLOG << "Population labels found are:" << endl;
	cout << line1 << endl;
	writeLOG << line1 << endl;
//	getline(LAIfile, line2, '\n');
	getline(LAIfile, line3, '\n');
	split(parsedline3, line3, is_any_of(delimiter));
	id_num = (parsedline3.size() - 6) / 2;
	for (int i = 6; i < parsedline3.size(); i++) 
	{
		string str_tmp;
		int pos_tmp;
		str_tmp = parsedline3[i];
		pos_tmp = str_tmp.find(".");
		id.push_back(str_tmp.substr(0, pos_tmp));
		i++;
	}
	LAIfile.close();
	cout << id_num << " individuals found in [" << LAIName << "]." << endl;
	writeLOG << id_num << " individuals found in [" << LAIName << "]." << endl;

	region_num = 0;
	for (int chr = 1; chr <=22; chr++)
	{
		chr_region_num = 0;
		chr_snp_num = 0;
		string LAIName = (ifile + to_string(chr) + suffix);
		ifstream LAIfile;
		LAIfile.open(LAIName.c_str(),ios::in);
		if (!LAIfile)
		{
			cout << "ERROR: " << LAIName << " was not found!" << endl;
			writeLOG << "ERROR: " << LAIName << " was not found!" << endl;
			exit (0);
		}
		csidx.push_back(region_num);
		getline(LAIfile, thisline, '\n');
//		getline(LAIfile, thisline, '\n');
		getline(LAIfile, thisline, '\n');
		while (getline(LAIfile, thisline, '\n'))
		{
			tmp1.clear();
			tmp2.clear();
			split(parsedline, thisline, is_any_of(delimiter));
			chrom.push_back(atoi(parsedline[0].c_str()));
			spos.push_back(atoi(parsedline[1].c_str()));
			epos.push_back(atoi(parsedline[2].c_str()));
			nsnps.push_back(atoi(parsedline[5].c_str()));
			for (int i = 0; i < id_num; i++)
			{
				tmp1.push_back(atoi(parsedline[ 2 * i + 6].c_str()));
				tmp2.push_back(atoi(parsedline[ 2 * i + 7].c_str()));
			}
			hap1.push_back(tmp1);
			hap2.push_back(tmp2);
			chr_region_num++;
			chr_snp_num = chr_snp_num + atoi(parsedline[5].c_str());
			region_num++;
		}
		LAIfile.close();
		ceidx.push_back(region_num - 1);
		cout << chr_region_num << " regions and " << chr_snp_num << " SNPs found in [" << LAIName << "]." << endl;
		writeLOG << chr_region_num << " regions and " << chr_snp_num << " SNPs found in [" << LAIName << "]." << endl;
	}
//	cout << "Total " << region_num << " regions found in the 22 LAI files." << endl;
//	writeLOG << "Total " << region_num << " regions found in the 22 LAI files." << endl;
	
	int this_chrom, this_start, this_end, this_hap1, this_hap2, snp_num;
	vector<int> one_chrom, one_start, one_end, one_hap1, one_hap2, one_snp_num;
	vector<int> out_chrom, out_start, out_end, out_hap1, out_hap2, out_snp_num;
	cout << "Writing individual LAI files for plotting with karyoploteR and Tagore..." << flush;
	writeLOG << "Writing individual LAI files for plotting with karyoploteR and Tagore..." << flush;

	for (int i= 0; i < id_num; i++)
	{
		ofstream tagfile;
		tagfile.open((ofile + "." + id[i] + ".bed").c_str(), ios::out);
		ofstream laifile;
		laifile.open((ofile + "." + id[i] + ".lai").c_str(), ios::out);
		laifile << line1 << endl;
//		laifile << line2 << endl;
		laifile << "#chr" << "\t" << "spos" << "\t" << "epos" << "\t" << "sgpos" << "\t" << "egpos" << "\t" << "n snps" << "\t" << id[i] << ".0" << "\t" << id[i] << ".1" << endl;

		for (int ch = 0; ch < 22; ch++)
		{
			one_chrom.clear();
			one_start.clear();
			one_end.clear();
			one_hap1.clear();
			one_hap2.clear();
			one_snp_num.clear();
			this_chrom = chrom[csidx[ch]];
			this_start = spos[csidx[ch]];
			this_end = epos[csidx[ch]];	
			this_hap1 = hap1[csidx[ch]][i];
			this_hap2 = hap2[csidx[ch]][i];	
			snp_num = nsnps[csidx[ch]];	
			chr_snp_num = nsnps[csidx[ch]];	
			for (int j = csidx[ch] + 1; j <= ceidx[ch]; j++)
			{
				if (hap1[j][i] == this_hap1 && hap2[j][i] == this_hap2)
				{
					snp_num = snp_num + nsnps[j];
					chr_snp_num = chr_snp_num + nsnps[j];
					this_end = epos[j];
					continue;
				}
				else
				{
					one_chrom.push_back(this_chrom);
					one_start.push_back(this_start);
					one_end.push_back(this_end);
					one_hap1.push_back(this_hap1);
					one_hap2.push_back(this_hap2);
					one_snp_num.push_back(snp_num);
					this_chrom = chrom[j];
					this_start = spos[j];
					this_end = epos[j];
					this_hap1 = hap1[j][i];
					this_hap2 = hap2[j][i];	
					snp_num = nsnps[j];
					chr_snp_num = chr_snp_num + nsnps[j];
				}
			}
			one_chrom.push_back(this_chrom);
			one_start.push_back(this_start);
			one_end.push_back(this_end);
			one_hap1.push_back(this_hap1);
			one_hap2.push_back(this_hap2);
			one_snp_num.push_back(snp_num);

			if (rfmix)
			{
				for (int j = 0; j < one_chrom.size(); j++)
				{
					laifile << one_chrom[j] << "\t" << one_start[j] << "\t" << one_end[j]  << "\t" << "." << "\t" << "." << "\t" << one_snp_num[j] << "\t" << one_hap1[j] << "\t" << one_hap2[j] << endl;
					tagfile << one_chrom[j] << "\t" << one_start[j] << "\t" << one_end[j]  << "\t" << "0" << "\t" << "1" << "\t" << palette[one_hap1[j]] << "\t" << "1" << endl;
					tagfile << one_chrom[j] << "\t" << one_start[j] << "\t" << one_end[j]  << "\t" << "0" << "\t" << "1" << "\t" << palette[one_hap2[j]] << "\t" << "2" << endl;
				}
				continue;
			}

			int one_region_num  = one_chrom.size();
			int idx, hap_up, hap_down;
			if (one_hap1[0] == 0)
			{
				idx = 0;
				hap_down = 0;
				for (int j = 1; j < one_region_num; j++)
				{
					if (one_hap1[j] != 0) {	idx = j; hap_down = one_hap1[j]; break;}
				}
				for (int j = 0; j < idx; j++) one_hap1[j] = hap_down;
			}
			if (one_hap2[0] == 0)
			{
				idx = 0;
				hap_down = 0;
				for (int j = 1; j < one_region_num; j++)
				{
					if (one_hap2[j] != 0) {	idx = j; hap_down = one_hap2[j]; break;}
				}
				for (int j = 0; j < idx; j++) one_hap2[j] = hap_down;
			}
			if (one_hap1[one_region_num - 1] == 0)
			{
				idx = one_region_num - 1;
				hap_up = 0;
				for (int j = one_region_num - 2; j >= 0; j--)
				{
					if (one_hap1[j] != 0) {	idx = j; hap_up = one_hap1[j]; break;}
				}
				for (int j = one_region_num - 1; j > idx; j--) one_hap1[j] = hap_up;
			}
			if (one_hap2[one_region_num - 1] == 0)
			{
				idx = one_region_num - 1;
				hap_up = 0;
				for (int j = one_region_num - 2; j >= 0; j--)
				{
					if (one_hap2[j] != 0) {	idx = j; hap_up = one_hap2[j]; break;}
				}
				for (int j = one_region_num - 1; j > idx; j--) one_hap2[j] = hap_up;
			}

			for (int j = 1; j < one_region_num - 1; j++)
			{
				if (one_hap1[j] == 0)
				{
					idx = j;
					hap_up = one_hap1[j - 1];
					hap_down = 0;
					for (int k = j + 1; k < one_region_num; k++)
					{
						if (one_hap1[k] != 0) {idx = k; hap_down = one_hap1[k]; break;}
					}
					if (hap_up == hap_down)
					{
						for (int k = j; k < idx; k++) one_hap1[k] = hap_up;
					}
				}
			}
			for (int j = 1; j < one_region_num - 1; j++)
			{
				if (one_hap2[j] == 0)
				{
					idx = j;
					hap_up = one_hap2[j - 1];
					hap_down = 0;
					for (int k = j + 1; k < one_region_num; k++)
					{
						if (one_hap2[k] != 0) {idx = k; hap_down = one_hap2[k]; break;}
					}
					if (hap_up == hap_down)
					{
						for (int k = j; k < idx; k++) one_hap2[k] = hap_up;
					}
				}
			}

			out_chrom.clear();
			out_start.clear();
			out_end.clear();
			out_hap1.clear();
			out_hap2.clear();
			out_snp_num.clear();
			this_chrom = one_chrom[0];
			this_start = one_start[0];
			this_end = one_end[0];	
			this_hap1 = one_hap1[0];
			this_hap2 = one_hap2[0];	
			snp_num = one_snp_num[0];	
			chr_snp_num = one_snp_num[0];
			for (int j = 1; j < one_region_num; j++)
			{
				if (one_hap1[j] == this_hap1 && one_hap2[j] == this_hap2)
				{
					snp_num = snp_num + one_snp_num[j];
					chr_snp_num = chr_snp_num + one_snp_num[j];
					this_end = one_end[j];
					continue;
				}
				else
				{
					out_chrom.push_back(this_chrom);
					out_start.push_back(this_start);
					out_end.push_back(this_end);
					out_hap1.push_back(this_hap1);
					out_hap2.push_back(this_hap2);
					out_snp_num.push_back(snp_num);
					this_chrom = one_chrom[j];
					this_start = one_start[j];
					this_end = one_end[j];
					this_hap1 = one_hap1[j];
					this_hap2 = one_hap2[j];	
					snp_num = one_snp_num[j];
					chr_snp_num = chr_snp_num + one_snp_num[j];
				}
			}
			out_chrom.push_back(this_chrom);
			out_start.push_back(this_start);
			out_end.push_back(this_end);
			out_hap1.push_back(this_hap1);
			out_hap2.push_back(this_hap2);
			out_snp_num.push_back(snp_num);

			if (fill_gaps)
			{

			}

			for (int j = 0; j < out_chrom.size(); j++)
			{
				laifile << out_chrom[j] << "\t" << out_start[j] << "\t" << out_end[j]  << "\t" << "." << "\t" << "." << "\t" << out_snp_num[j] << "\t" << out_hap1[j] << "\t" << out_hap2[j] << endl;
				tagfile << out_chrom[j] << "\t" << out_start[j] << "\t" << out_end[j]  << "\t" << "0" << "\t" << "1" << "\t" << palette[out_hap1[j]] << "\t" << "1" << endl;
				tagfile << out_chrom[j] << "\t" << out_start[j] << "\t" << out_end[j]  << "\t" << "0" << "\t" << "1" << "\t" << palette[out_hap2[j]] << "\t" << "2" << endl;
			}
		}
		tagfile.close();
		laifile.close();
	}
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


};  

int main(int argc, char ** argv)
{
	time_t start_time, current_time;
	char * current_time_str;
	LAI lai(argc, argv);
	start_time = time(0);
	current_time_str = ctime(&start_time);
	cout << "Start time: " << current_time_str << flush;
	lai.writeLOG << "Start time: " << current_time_str;
	lai.LAI_input();
	current_time = time(0);
	current_time_str = ctime(&current_time);
	double time_difference = difftime(current_time, start_time);
	cout << "End time: " << current_time_str << "Running time: " << time_difference << " s" << endl << endl;
	lai.writeLOG << "End time: " << current_time_str << "Running time: " << time_difference << " s" << endl << endl;
	return 0;
}
