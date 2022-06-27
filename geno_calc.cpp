#include "geno_calc.h"

void ASA::readFAM(string binfilename)
{
	ifstream famfile;
	string thisline;
	vector<string> parsedline;
	famfile.open((binfilename + ".fam").c_str(), ios::in);
	if(!famfile)
	{
		cout << binfilename + ".fam" << " was not found!" << endl;
		writeLOG << binfilename + ".fam" << " was not found!" << endl;
		exit (0);
	}

	while (getline(famfile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		famid.push_back(parsedline[0]);
		subid.push_back(parsedline[1]);
	}

	famfile.close();
	indi_num = subid.size();
	cout << indi_num << " people loaded from FAM file [" + binfilename + ".fam]." << endl;
	writeLOG << indi_num << " people loaded from FAM file [" + binfilename + ".fam]." << endl;
}


void ASA::readBIM(string binfilename)
{
	ifstream bimfile;
	string thisline;
	vector<string> parsedline;
	bimfile.open((binfilename + ".bim").c_str(), ios::in);
	if(!bimfile)
	{
		cout << binfilename + ".bim" << " was not found!" << endl;
		writeLOG << binfilename + ".bim" << " was not found!" << endl;
		exit (0);
	}

	int chrom_this, pos_this;
	int i = 0, chrom_last = -1, pos_last = -1, sex_snp_num = 0, pos_del_num = 0, miss_allele_num = 0;
	while (getline(bimfile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		chrom_this = atoi(parsedline[0].c_str());
		chrom.push_back(chrom_this);
		snp.push_back(parsedline[1]);
		pos_this = atoi(parsedline[3].c_str());
		position.push_back(pos_this);
		allele1.push_back(parsedline[4]);
		allele2.push_back(parsedline[5]);
		snp_to_delete.push_back(0);
		if (allele1[i] == "*" || allele1[i] == "." || allele2[i] == "*" || allele2[i] == ".")
		{
			snp_to_delete[i] = 1;
			miss_allele_num++;
		}
		if (chrom_this < 1 || chrom_this > 22)
		{
			snp_to_delete[i] = 1;
			sex_snp_num++;
		}
		if (chrom_this == chrom_last && pos_this == pos_last)
		{
			snp_to_delete[i] = 1;
			pos_del_num++;
			if (snp_to_delete[i - 1] = 0)
			{
				snp_to_delete[i - 1] = 1;
				pos_del_num++;
			}
		}
		chrom_last = chrom_this;
		pos_last = pos_this;
		i ++;
	}

	bimfile.close();	
	snp_num = snp.size();
	cout << snp_num << " SNPs loaded from BIM file [" + binfilename + ".bim]." << endl;
	writeLOG << snp_num << " SNPs loaded from BIM file [" + binfilename + ".bim]." << endl;
	if (pos_del_num > 0)
	{
		cout << pos_del_num << " SNPs removed due to same physical positions." << endl;
		writeLOG << pos_del_num << " SNPs removed due to same physical positions." << endl;
	}
	if (sex_snp_num > 0)
	{
		cout << sex_snp_num << " SNPs on sex chromosomes removed." << endl;
		writeLOG << sex_snp_num << "SNPs on sex chromosomes removed." << endl;
	}
	if (miss_allele_num > 0)
	{
		cout << miss_allele_num << " SNPs removed due to missing allele information." << endl;
		writeLOG << miss_allele_num << " SNPs removed due to missing allele information." << endl;
	}
}


void ASA::snp_info_check(string ifile)
{
	string thisline;
	vector<string> parsedline;
	infofile.open(ifile.c_str(), ios::in);
	if(!infofile)
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
	if (info_col_num < 5)
	{
		cout << "ERROR: " << ifile << " has " << info_col_num << " columns, 5 columns are expeted." << endl;
		writeLOG << "ERROR: " << ifile << " has " << info_col_num << " columns, 5 columns are expeted." << endl;
		exit(0);
	}
}


void ASA::snp_info_input(string ifile, bool match_alleles)
{
	int record_num;
	string thisline, this_snp;
	vector <string> parsedline, info_snp, info_diff_pop;
	vector <int> info_pop_snp_num;
	multimap<string, int> info_pop_snp;
	set<string> pop_set, info_pop_set;
	unordered_multiset<string> info_pop_multiset;
	infofile.open(ifile.c_str(), ios::in);
	info_snp_num = 0;
	while (getline(infofile, thisline, '\n'))
	{
		split(parsedline, thisline, is_any_of(delimiter));
		this_snp = parsedline[0];
		info_pop_snp.insert(multimap<string, int>::value_type(this_snp, info_snp_num));
		info_snp.push_back(this_snp);
		info_pop_MiA.push_back(parsedline[3]);
		info_pop_MaA.push_back(parsedline[4]);
		info_pop.push_back(parsedline[1]);
		info_pop_set.insert(parsedline[1]);
		info_pop_multiset.insert(parsedline[1]);
		info_pop_maf.push_back(atof(parsedline[2].c_str()));
		info_snp_num++;
	}
	infofile.close(); 	

	cout << info_snp_num << " SNPs loaded from the file [" << ifile << "]." << endl;
	writeLOG << info_snp_num << " SNPs loaded from the file [" << ifile << "]." << endl;
	int info_pop_num = info_pop_set.size();
	cout << info_pop_num << " populations found, population ID and number of the population-specific SNPs are " << endl;
 	writeLOG << info_pop_num << " populations found, population ID and number of the population-specific SNPs are " << endl;

	set<string>::iterator it;
	for (it = info_pop_set.begin(); it != info_pop_set.end(); ++ it)
	{
		info_diff_pop.push_back(*it);
		info_pop_snp_num.push_back(info_pop_multiset.count(*it)); 
	}
	for(int i = 0; i < info_pop_num; ++ i)
	{
		cout << info_diff_pop[i] << ": " << info_pop_snp_num[i] << "\t";
		writeLOG << info_diff_pop[i] << ": " << info_pop_snp_num[i] << "\t";
	}
 	cout << endl;
	writeLOG << endl;

	for(int i = 0; i < info_snp_num; ++ i)
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
	int diff_allele_num = 0;
	for(int i = 0; i < snp_num; ++ i)
	{
		if (info_pop_snp.count(snp[i]) == 0)
		{
			not_in_info_num ++;
			pop_maf.push_back(0);
			pop_MiA.push_back("NA");
			pop_MaA.push_back("NA");
			pop.push_back("NA");
			snp_to_delete[i] = 1;
			continue;
		}
		int idx = info_pop_snp.find(snp[i])->second;
		pop_maf.push_back(info_pop_maf[idx]);
		pop_MiA.push_back(info_pop_MiA[idx]);
		pop_MaA.push_back(info_pop_MaA[idx]);
		info_in_bim[idx] = 1;
		pop.push_back(info_pop[idx]);
		if (match_alleles)
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

	pop_num = pop_set.size();
	for (it = pop_set.begin(); it != pop_set.end(); ++ it) diff_pop.push_back(*it);
	if (not_in_info_num > 0)
	{
		cout << not_in_info_num << " SNPs in the genotype file were excluded, which were not in the SNP panel." << endl;
		writeLOG << not_in_info_num << " SNPs in the genotype file were excluded, which were not in the SNP panel." << endl;
	}
	if(remain_snp_num == 0)
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
	cout << "Genotypes of " << remain_snp_num << " SNPs will be used for the analysis." << endl;
	writeLOG << "Genotypes of " << remain_snp_num << " SNPs will be used for the analysis." << endl;
}


void ASA::make_X_Z(string binfilename, int size, bool centralize, bool standardize, bool ifile_flag)
{
	cout << "Calculating genotype relationship matrix (GRM)... " << flush;
	writeLOG << "Calculating genetic relationship matrix (GRM)... ";

	bedfile.open((binfilename + ".bed").c_str(), ios::in | ios::binary);
	char tmp1char[1];
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);

	int line_byte_size = ceil(double(indi_num) / 4);
	int batch_snp_size = size;
	int batch_byte_size = batch_snp_size * line_byte_size;
	int batch_num = ceil(double(snp_num) / batch_snp_size);
	char *batch_char = new char[batch_byte_size];
	string batch_str;
	eigenMatrix batch_X;

	Z_N = vector<vector<unsigned int>> (indi_num, vector<unsigned int>(indi_num, 0));    
	Z.setZero(indi_num, indi_num);
	vector< vector<int> > batch_miss_pos;
	
	for (int batch_idx = 1; batch_idx <= batch_num; ++ batch_idx)
	{
		int snp_start = (batch_idx - 1) * batch_snp_size;
		int snp_end = batch_idx * batch_snp_size;
		if (batch_idx == batch_num)
		{
			snp_end = snp_num;
			batch_snp_size = snp_end - snp_start;
			batch_byte_size = batch_snp_size * line_byte_size;
		}

		batch_X.resize(indi_num, batch_snp_size);
		vector< vector<int> > X_bool(indi_num, vector<int>(batch_snp_size));
		batch_miss_pos.resize(indi_num);
		bedfile.read(batch_char, batch_byte_size);
		batch_str.assign(batch_char, batch_byte_size);
		dynamic_bitset<unsigned char> batch_bit(batch_str.begin(), batch_str.end());

		#pragma omp parallel for
		for (int i = snp_start; i < snp_end; ++ i)
		{
			int k = i - snp_start;
			if (snp_to_delete[i])
			{
				batch_X.col(k).setZero();
				for (int j = 0; j < indi_num; ++ j) X_bool[j][k] = 1;
				continue;
			}
			 
			eigenVar mu, inv_sd_snp, batch_X_0, batch_X_1, batch_X_2;
			batch_X_0 = 0;
			batch_X_1 = 1;
			batch_X_2 = 2;
			if (centralize)
			{
				mu = 2 * pop_maf[i];
				batch_X_0 = 0 - mu;
				batch_X_1 = 1 - mu;
				batch_X_2 = 2 - mu;
			}
			if (standardize)
			{
				mu = 2 * pop_maf[i];
				if(pop_maf[i] < 1.0e-50) inv_sd_snp = 0.0;
				else inv_sd_snp = sqrt(1 / (2 * pop_maf[i] * (1 - pop_maf[i])));
				batch_X_0 = (0 - mu) * inv_sd_snp;
				batch_X_1 = (1 - mu) * inv_sd_snp;
				batch_X_2 = (2 - mu) * inv_sd_snp;
			}

			int pos = (line_byte_size) * (i - snp_start) * 8;
			for (int j = 0; j < indi_num; ++ j)
			{
				bool tmp1bool = batch_bit[pos ++ ];
				bool tmp2bool = batch_bit[pos ++ ];

				if (!tmp1bool && !tmp2bool)
				{
					if(allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i])
					{
						batch_X(j, k) = batch_X_2;
						X_bool[j][k] = 1;
					}
					else
					{
						batch_X(j, k) = batch_X_0;
						X_bool[j][k] = 1;
					}
					continue;
				}
				if (!tmp1bool && tmp2bool)
				{
					batch_X(j, k) = batch_X_1;
					X_bool[j][k] = 1;
					continue;
				}
				if (tmp1bool && tmp2bool)
				{
					if(allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i])
					{
						batch_X(j, k) = batch_X_0;
						X_bool[j][k] = 1;
					}
					else
					{
						batch_X(j, k) = batch_X_2;
						X_bool[j][k] = 1;
					}
					continue;
				}
				batch_X(j, k) = 0;
				X_bool[j][k] = 0;
			}
		}
		
		int filtered_num = 0;
		for (int i = snp_start; i < snp_end; ++ i)
		{
			if (snp_to_delete[i]) filtered_num = filtered_num + 1;
		}
		
		#pragma omp parallel for
		for (int i = 0; i < indi_num; ++ i)
		{
			for (int j = 0; j < batch_snp_size; ++ j)
			{
				if(X_bool[i][j] == 0) batch_miss_pos[i].push_back(j);
			}
		}
		
		#pragma omp parallel for
		for (int i = 0; i < indi_num; ++ i)
		{
			for (int j = 0; j <= i; ++ j)
			{
				int miss_j = 0;
				for (int k = 0; k < batch_miss_pos[j].size(); ++ k) miss_j = miss_j + X_bool[i][batch_miss_pos[j][k]];
				Z_N[i][j] = Z_N[i][j] + batch_snp_size - batch_miss_pos[i].size() - miss_j - filtered_num;
				Z_N[j][i] = Z_N[i][j];
			}
		}
		
		X_bool.clear();
		batch_miss_pos.clear();
		Z = Z + batch_X * batch_X.transpose();
	}
	
	#pragma omp parallel for
	for (int i = 0; i < indi_num; ++ i)
	{
		for (int j = 0; j <= i; ++ j)
		{
			Z(i, j) = Z(i, j) / Z_N[i][j];
			Z(j, i) = Z(i, j);
		}
	}
	
	delete [] batch_char;
	bedfile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::grm_id_output(string ofile)
{	
	ofstream writefile;
	writefile.open((ofile + ".grm.id").c_str(), ios::out);

	cout << "Saving IDs of the GRM [" + ofile + ".grm.gz] to [" + ofile + ".grm.id]... " << flush;
	writeLOG << "Saving IDs of the GRM [" + ofile + ".grm.gz] to [" + ofile + ".grm.id]... ";

	for (int i = 0; i<indi_num; ++ i) writefile << famid[i] << "\t" << subid[i] << endl;
	cout << "done." << endl;
	writeLOG << "done." << endl;
	writefile.close();
}


void ASA::grm_output(string ofile)
{
	ostringstream writefile;
	gzFile gzWrite;
	gzWrite = gzopen((ofile + ".grm.gz").c_str(), "wb");

	cout << "Saving the GRM to [" + ofile + ".grm.gz]... " << flush;
	writeLOG << "Saving the GRM to [" + ofile + ".grm.gz]... ";
	
	for (int i = 0; i < indi_num; ++ i)
	{
		for (int j = 0; j <= i; ++ j)
		{
			writefile << i + 1 << "\t" << j + 1 << "\t" << Z_N[i][j] << "\t" << Z(i, j) << endl;
			gzputs(gzWrite, writefile.str().c_str());
			writefile.str("");
		}
	}

	gzclose(gzWrite);
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::grm_eigen()
{
	cout << "Performing principal component analysis... " << flush;
	writeLOG << "Performing principal component analysis... ";
	SelfAdjointEigenSolver < eigenMatrix> es(Z);
	Z_V = es.eigenvectors();
	Z_D = es.eigenvalues();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}
	

void ASA::grm_eigvec_pc_eigval_output(string eigval_ofile, string eigvec_ofile, string pc_ofile, int N)
{
	
	ofstream writefile1;
	writefile1.open((eigval_ofile + ".eigenval").c_str(), ios::out);
	ofstream writefile2;
	writefile2.open((eigvec_ofile + ".eigenvec").c_str(), ios::out);
	ofstream writefile3;
	writefile3.open((pc_ofile + ".pc").c_str(), ios::out);

	if (N > indi_num) N = indi_num;
	cout << "Saving top " << N << " eigenvalues of " << indi_num << " individuals to [" + eigval_ofile + ".eigenval] file... " << flush;
	writeLOG << "Saving top " << N << " eigenvalues of " << indi_num << " individuals to [" + eigval_ofile + ".eigenval] file... ";
	double sum_eigenval = 0.0;
	
	for (int i = indi_num - 1; i  >= 0; --i) sum_eigenval = sum_eigenval + Z_D(i);
	for(int j = indi_num - 1; j >= indi_num - N; --j) writefile1 << Z_D(j) << "\t" << Z_D(j) / sum_eigenval << endl;
	
	writefile1.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
	cout << "Saving top " << N << " eigenvectors and PCs of " << indi_num << " individuals to [" + eigvec_ofile + ".eigenvec] and [" + pc_ofile + ".pc] files... " << flush;
	writeLOG << "Saving top " << N << " eigenvectors and PCs of " << indi_num << " individuals to [" + eigvec_ofile + ".eigenvec] and [" + pc_ofile + ".pc] files... ";
	
	for (int i = 0; i < indi_num; ++ i)
	{
		writefile2 << famid[i] << "\t" << subid[i] << "\t";
		writefile3 << famid[i] << "\t" << subid[i] << "\t";
		for(int j = indi_num - 1; j >= indi_num - N; --j)
		{
			writefile2 << Z_V(i, j) << "\t";
			writefile3 << sqrt(remain_snp_num * Z_D(j)) * Z_V(i, j) << "\t";
		}
		writefile2 << endl;
		writefile3 << endl;
	}
	writefile2.close();
	writefile3.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::make_X_psv(string binfilename, int size)
{
	cout << "Calculating principal score vectors... " << flush;
	writeLOG << "Calculating principal score vectors... ";
	eigenMatrix b, c;
	b.setZero(snp_num, pop_num);
	c.setZero(indi_num, pop_num);
	psv.setZero(indi_num, pop_num);

	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < snp_num; ++ i)
		{
			if (!snp_to_delete[i])
			{
				if (pop[i] == diff_pop[j])
				{
					b(i, j) = 2 * pop_maf[i];
					c(0, j) = c(0, j) + 4 * pop_maf[i] * pop_maf[i];
				}
			}
		}
		for (int i = 0; i < info_snp_num; ++ i)
		{
			if (!info_in_bim[i])
			{
				if (info_pop[i] == diff_pop[j]) c(0, j) = c(0, j) + 4 * info_pop_maf[i] * info_pop_maf[i];
			}
		}
		for (int k = 1; k < indi_num; ++ k) c(k, j) = c(0, j);
	}
	
	bedfile.open((binfilename + ".bed").c_str(), ios::in | ios::binary);
	char tmp1char[1];
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	int line_byte_size = ceil(double(indi_num) / 4);
	int batch_snp_size = size;
	int batch_byte_size = batch_snp_size * line_byte_size;
	int batch_num = ceil(double(snp_num) / batch_snp_size);
	char *batch_char = new char[batch_byte_size];
	string batch_str;
	eigenMatrix batch_X;

	for (int batch_idx = 1; batch_idx <= batch_num; ++ batch_idx)
	{
		int snp_start = (batch_idx - 1) * batch_snp_size;
		int snp_end = batch_idx * batch_snp_size;
		if (batch_idx == batch_num)
		{
			snp_end = snp_num;
			batch_snp_size = snp_end - snp_start;
			batch_byte_size = batch_snp_size * line_byte_size;
		}
		
		vector<vector<int>> missing_bool(indi_num, vector<int>(batch_snp_size, 0));
		batch_X.resize(indi_num, batch_snp_size);
		bedfile.read(batch_char, batch_byte_size);
		batch_str.assign(batch_char, batch_byte_size);
		dynamic_bitset<unsigned char> batch_bit(batch_str.begin(), batch_str.end());
		
		#pragma omp parallel for
		for (int i = snp_start; i < snp_end; ++ i)
		{
			int k = i - snp_start;
			if (snp_to_delete[i])
			{
				batch_X.col(k).setZero();
				for (int j = 0; j < indi_num; ++ j) missing_bool[j][k] = 0;
				continue;
			}
			 
			int pos = (line_byte_size) * (i - snp_start) * 8;
			for (int j = 0; j < indi_num; ++ j)
			{
				bool tmp1bool = batch_bit[pos ++ ];
				bool tmp2bool = batch_bit[pos ++ ];
				if (!tmp1bool && !tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i]) batch_X(j, k) = 2;
					else batch_X(j, k) = 0;
					continue;
				}
				if (!tmp1bool && tmp2bool)
				{
					batch_X(j, k) = 1;
					continue;
				}
				if (tmp1bool && tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i]) batch_X(j, k) = 0;
					else batch_X(j, k) = 2;
					continue;
				}
				batch_X(j, k) = 0;
				missing_bool[j][k] = 1;
			}
		}
	
		psv = psv + batch_X * (b.middleRows(snp_start, batch_snp_size));

		#pragma omp parallel for
		for (int k = 0; k < pop_num; ++ k)
		{
			for (int i = 0; i < indi_num; ++ i)
			{
				for (int j = 0; j < batch_snp_size; ++ j)
				{
					if (missing_bool[i][j])
					{
						if (b(snp_start + j, k) != 0) c(i, k) = c(i, k) - b(snp_start + j, k) * b(snp_start + j, k);
					}
				}
			}
		}
	}
	
	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < indi_num; ++ i) psv(i, j) = psv(i, j) / c(i, j);
	}

	delete [] batch_char;
	bedfile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::psv_output(string ofile)
{
	ofstream writefile;
	writefile.open((ofile + ".psv").c_str(), ios::out);

	cout << "Saving the principal score vectors to [" + ofile + ".psv]... " << flush;
	writeLOG << "Saving the principal score vectors to [" + ofile + ".psv]... ";
	writefile << "ID_1" << "\t" << "ID_2";
	for (int j = 0; j < pop_num; ++ j) writefile << "\t" << diff_pop[j];
	writefile << endl;
	
	for (int i = 0; i < indi_num; ++ i)
	{
		writefile << famid[i] << "\t" << subid[i];
		for (int j = 0; j < pop_num; ++ j) writefile << "\t" << psv(i, j);
		writefile << endl;
	}
	
	writefile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::make_X_aiv(string binfilename, int size)
{
	cout << "Calculating ancestral information vectors by the method of moment... " << flush;
	writeLOG << "Calculating ancestral information vectors by the method of moment... ";
	eigenMatrix b, c;
	b.setZero(snp_num, pop_num);
	c.setZero(indi_num, pop_num);
	aiv.setZero(indi_num, pop_num);

	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < snp_num; ++ i)
		{
			if (!snp_to_delete[i])
			{
				if (pop[i] == diff_pop[j])
				{
					b(i, j) = 1;
					c(0, j) = c(0, j) + 2 * pop_maf[i];
				}
			}
		}
		for (int i = 0; i < info_snp_num; ++ i)
		{
			if (!info_in_bim[i])
			{
				if (info_pop[i] == diff_pop[j]) c(0, j) = c(0, j) + 2 * info_pop_maf[i];
			}
		}
		for (int k = 1; k < indi_num; ++ k) c(k, j) = c(0, j);
	}

	bedfile.open((binfilename + ".bed").c_str(), ios::in | ios::binary);
	char tmp1char[1];
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	int line_byte_size = ceil(double(indi_num) / 4);
	int batch_snp_size = size;
	int batch_byte_size = batch_snp_size * line_byte_size;
	int batch_num = ceil(double(snp_num) / batch_snp_size);
	char *batch_char = new char[batch_byte_size];
	string batch_str;
	eigenMatrix batch_X;

	for (int batch_idx = 1; batch_idx <= batch_num; ++ batch_idx)
	{
		int snp_start = (batch_idx - 1) * batch_snp_size;
		int snp_end = batch_idx * batch_snp_size;
		if (batch_idx == batch_num)
		{
			snp_end = snp_num;
			batch_snp_size = snp_end - snp_start;
			batch_byte_size = batch_snp_size * line_byte_size;
		}
		
		vector<vector<int>> missing_bool(indi_num, vector<int>(batch_snp_size, 0));
		batch_X.resize(indi_num, batch_snp_size);
		bedfile.read(batch_char, batch_byte_size);
		batch_str.assign(batch_char, batch_byte_size);
		dynamic_bitset<unsigned char> batch_bit(batch_str.begin(), batch_str.end());
		
		#pragma omp parallel for
		for (int i = snp_start; i < snp_end; ++ i)
		{
			int k = i - snp_start;
			if (snp_to_delete[i])
			{
				batch_X.col(k).setZero();
				for (int j = 0; j < indi_num; ++ j) missing_bool[j][k] = 0;
				continue;
			}
				 
			int pos = (line_byte_size) * (i - snp_start) * 8;
			for (int j = 0; j < indi_num; ++ j)
			{
				bool tmp1bool = batch_bit[pos ++ ];
				bool tmp2bool = batch_bit[pos ++ ];
				if (!tmp1bool && !tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i]) batch_X(j, k) = 2;
					else batch_X(j, k) = 0;
					continue;
				}
				if (!tmp1bool && tmp2bool)
				{
					batch_X(j, k) = 1;
					continue;
				}
				if (tmp1bool && tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i]) batch_X(j, k) = 0;
					else batch_X(j, k) = 2;
					continue;
				}
				batch_X(j, k) = 0;
				missing_bool[j][k] = 1;
			}
		}

		aiv = aiv + batch_X * (b.middleRows(snp_start, batch_snp_size));

		#pragma omp parallel for
		for (int k = 0; k < pop_num; ++ k)
		{
			for (int i = 0; i < indi_num; ++ i)
			{
				for (int j = 0; j < batch_snp_size; ++ j)
				{
					if (missing_bool[i][j])
					{
						if (b(snp_start + j, k) != 0) c(i, k) = c(i, k) - 2 * pop_maf[snp_start + j];
					}
				}
			}
		}
	}
	
	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < indi_num; ++ i)  aiv(i, j) = aiv(i, j) / c(i, j);
	}

	delete [] batch_char;
	bedfile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}

void ASA::aiv_output(string ofile)
{
	ofstream writefile;
	writefile.open((ofile + ".aiv").c_str(), ios::out);
	cout << "Saving the ancestral information vectors to [" + ofile + ".aiv]... " << flush;
	writeLOG << "Saving the ancestral information vectors to [" + ofile + ".aiv]... ";
	writefile << "ID_1" << "\t" << "ID_2";

	for (int j = 0; j < pop_num; ++ j) writefile << "\t" << diff_pop[j];
	writefile << endl;
	
	for (int i = 0; i < indi_num; ++ i)
	{
		writefile << famid[i] << "\t" << subid[i];
		for (int j = 0; j < pop_num; ++ j)   writefile << "\t" << aiv(i, j);
		writefile << endl;
	}
	
	writefile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


void ASA::make_X_mle(string binfilename, int size)
{
	cout << "Calculating ancestral information vectors by the approximate maximum likelihood emtimate... " << flush;
	writeLOG << "Calculating ancestral information vectors by the approximate maximum likelihood emtimate... ";
	eigenMatrix a, b, c, aiv_a, aiv_b, aiv_c;
	a.setZero(snp_num, pop_num);
	b.setZero(snp_num, pop_num);
	c.setZero(snp_num, pop_num);
	aiv.setZero(indi_num, pop_num);
	aiv_a.setZero(indi_num, pop_num);
	aiv_b.setZero(indi_num, pop_num);
	aiv_c.setZero(indi_num, pop_num);

	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < snp_num; ++ i)
		{
			if (!snp_to_delete[i])
			{
				if (pop[i] == diff_pop[j])
				{
					a(i, j) = 4 * pop_maf[i] * pop_maf[i];
					b(i, j) = 2 * pop_maf[i];
					c(i, j) = 1;
				}
			}
		}
		for (int i = 0; i < info_snp_num; ++ i)
		{
			if (!info_in_bim[i])
			{
				if (info_pop[i] == diff_pop[j])
				{
					aiv_a(0, j) = aiv_a(0, j) + 4 * info_pop_maf[i] * info_pop_maf[i];
					aiv_b(0, j) = aiv_b(0, j) + 2 * info_pop_maf[i];
				}
			}
		}
		for (int k = 1; k < indi_num; ++ k)
		{
			aiv_a(k, j) = aiv_a(0, j);
			aiv_b(k, j) = aiv_b(0, j);
		}
	}
	
	bedfile.open((binfilename + ".bed").c_str(), ios::in | ios::binary);
	char tmp1char[1];
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	bedfile.read(tmp1char, 1);
	int line_byte_size = ceil(double(indi_num) / 4);
	int batch_snp_size = size;
	int batch_byte_size = batch_snp_size * line_byte_size;
	int batch_num = ceil(double(snp_num) / batch_snp_size);
	char *batch_char = new char[batch_byte_size];
	string batch_str;
	eigenMatrix batch_X, batch_1_X;

	for (int batch_idx = 1; batch_idx <= batch_num; ++ batch_idx)
	{
		int snp_start = (batch_idx - 1) * batch_snp_size;
		int snp_end = batch_idx * batch_snp_size;
		if (batch_idx == batch_num)
		{
			snp_end = snp_num;
			batch_snp_size = snp_end - snp_start;
			batch_byte_size = batch_snp_size * line_byte_size;
		}
		
		batch_X.resize(indi_num, batch_snp_size);
		batch_1_X.resize(indi_num, batch_snp_size);

		bedfile.read(batch_char, batch_byte_size);
		batch_str.assign(batch_char, batch_byte_size);
		dynamic_bitset<unsigned char> batch_bit(batch_str.begin(), batch_str.end());
		
		#pragma omp parallel for
		for (int i = snp_start; i < snp_end; ++ i)
		{
			int k = i - snp_start;
			if (snp_to_delete[i])
			{
				batch_X.col(k).setZero();
				batch_1_X.col(k).setZero();
				continue;
			}
			 
			int pos = (line_byte_size) * (i - snp_start) * 8;
			for (int j = 0; j < indi_num; ++ j)
			{
				bool tmp1bool = batch_bit[pos ++ ];
				bool tmp2bool = batch_bit[pos ++ ];
				if (!tmp1bool && !tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i])
					{
						batch_X(j, k) = 2;
						batch_1_X(j, k) = -1;
					}
					else
					{
						batch_X(j, k) = 0;
						batch_1_X(j, k) = 1;
					}
					continue;
				}
				if (!tmp1bool && tmp2bool)
				{
					batch_X(j, k) = 1;
					batch_1_X(j, k) = 0;
					continue;
				}
				if (tmp1bool && tmp2bool)
				{
					if (allele1[i] == pop_MiA[i] || allele2[i] == pop_MaA[i])
					{
						batch_X(j, k) = 0;
						batch_1_X(j, k) = 1;
					}
					else
					{
						batch_X(j, k) = 2;
						batch_1_X(j, k) = -1;
					}
					continue;
				}
				batch_X(j, k) = 0;
				batch_1_X(j, k) = 0;
			}
		}
	
		aiv_a = aiv_a + batch_1_X * (a.middleRows(snp_start, batch_snp_size));
		aiv_b = aiv_b + batch_1_X * (b.middleRows(snp_start, batch_snp_size));
		aiv_c = aiv_c + batch_X * (c.middleRows(snp_start, batch_snp_size));
	}

	#pragma omp parallel for
	for (int j = 0; j < pop_num; ++ j)
	{
		for (int i = 0; i < indi_num; ++ i) aiv(i, j) = (sqrt(aiv_b(i, j) * aiv_b(i, j) + 4 * aiv_a(i, j) * aiv_c(i, j)) - aiv_b(i, j)) / aiv_a(i, j) / 2;
	}

	delete [] batch_char;
	bedfile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
	if (not_in_info_num > 0)
	{
		cout << not_in_info_num << " SNPs were excluded, which are not in the SNP panel defined by the information file." << endl;
		writeLOG << not_in_info_num << " SNPs were excluded, which are not in the SNP panel defined by the information file." << endl;
	}
	cout << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
	writeLOG << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
}


void ASA::mle_output(string ofile)
{
	ofstream writefile;
	writefile.open((ofile + ".mle").c_str(), ios::out);
	cout << "Saving the ancestral information vectors to [" + ofile + ".mle]... " << flush;
	writeLOG << "Saving the ancestral information vectors to [" + ofile + ".mle]... ";
	writefile << "ID_1" << "\t" << "ID_2";

	for (int j = 0; j < pop_num; ++ j)  writefile << "\t" << diff_pop[j];
	writefile << endl;
	
	for (int i = 0; i < indi_num; ++ i)
	{
		writefile << famid[i] << "\t" << subid[i];
		for (int j = 0; j < pop_num; ++ j)  writefile << "\t" << aiv(i, j);
		writefile << endl;
	}
	
	writefile.close();
	cout << "done." << endl;
	writeLOG << "done." << endl;
}


