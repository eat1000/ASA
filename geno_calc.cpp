#include "geno_calc.h"

void ASA::allele_info_stat(string binfilename, double thre_num, double low_MAF, double high_MAF, double miss_max, int size)
{
    ifstream famfile, bimfile;
    string thisline;
    vector<string> parsedline;

    bedfile.open((binfilename + ".bed").c_str(), ios::in | ios::binary);
    famfile.open((binfilename + ".fam").c_str(), ios::in);
    bimfile.open((binfilename + ".bim").c_str(), ios::in);
	
    if(!bedfile)
    {
        cout << binfilename + ".bed" << " was not found!" << endl;
        writeLOG << binfilename + ".bed" << " was not found!" << endl;
        exit (0);
    }
	
    if(!famfile)
    {
        cout << binfilename + ".fam" << " was not found!" << endl;
        writeLOG << binfilename + ".fam" << " was not found!" << endl;
        exit (0);
    }
	
    if(!bimfile)
    {
        cout << binfilename + ".bim" << " was not found!" << endl;
        writeLOG << binfilename + ".bim" << " was not found!" << endl;
        exit (0);
    }
	
    while (getline(famfile, thisline, '\n'))
    {
        split(parsedline, thisline, is_any_of(delimiter));
        famid.push_back(parsedline[0]);
        subid.push_back(parsedline[1]);
    }

    indi_num = subid.size();
    famfile.close();
    cout << indi_num << " people loaded from FAM file [" + binfilename + ".fam]." << endl;
    writeLOG << indi_num << " people loaded from FAM file [" + binfilename + ".fam]." << endl;
	
    while (getline(bimfile, thisline, '\n'))
    {
        split(parsedline, thisline, is_any_of(delimiter));
        chrom.push_back(atoi(parsedline[0].c_str()));
        snp.push_back(parsedline[1]);
        position.push_back(atoi(parsedline[3].c_str()));
        allele1.push_back(parsedline[4]);
        allele2.push_back(parsedline[5]);
    }
	
    snp_num = snp.size();
    bimfile.close();
    cout << snp_num << " SNPs loaded from BIM file [" + binfilename + ".bim]." << endl;
    writeLOG << snp_num << " SNPs loaded from BIM file [" + binfilename + ".bim]." << endl;

    cout << "Calculating allele frequencies from BED file [" + binfilename + ".bed]... " << flush;
    writeLOG << "Calculating allele frequencies from BED file [" + binfilename + ".bed]... ";

    allele1_num.resize(snp_num);
    allele2_num.resize(snp_num);
    miss_num.resize(snp_num);
    miss_fre.resize(snp_num);
    allele1_fre.resize(snp_num);
    allele2_fre.resize(snp_num);
    maf.resize(snp_num);
    minor_allele.resize(snp_num);
    snp_to_delete.resize(snp_num);
    min_alle_eqt_alle2.resize(snp_num);
	
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
            allele1_num[i] = 0;
            allele2_num[i] = 0;
            miss_num[i] = 0;
            snp_to_delete[i] = 0;
            int pos = (line_byte_size) * (i - snp_start) * 8;
            for (int j = 0; j < indi_num; ++ j)
            {
                bool tmp1bool = batch_bit[pos ++ ];
                bool tmp2bool = batch_bit[pos ++ ];
                if (!tmp1bool && !tmp2bool)
                {
                    allele1_num[i] = allele1_num[i] + 2;
                    continue;
                }
                if (!tmp1bool && tmp2bool)
                {
                    ++ allele1_num[i];
                    ++ allele2_num[i];
                    continue;
                }
                if (tmp1bool && tmp2bool)
                {
                    allele2_num[i] = allele2_num[i] + 2;
                    continue;
                }
                ++ miss_num[i];
				
            }
			
            if(miss_num[i] == indi_num)
            {
                allele1_fre[i] = 0.0;
                allele2_fre[i] = 0.0;
                miss_fre[i] = 1;
            }
            else
            {
                miss_fre[i] = double(miss_num[i]) / indi_num;
                allele1_fre[i] = double(allele1_num[i]) / (indi_num - miss_num[i]) / 2;
                allele2_fre[i] = double(allele2_num[i]) / (indi_num - miss_num[i]) / 2;
            }
            if(allele1_num[i] > allele2_num[i])
            {
                minor_allele[i] = allele2[i];
                maf[i] = allele2_fre[i];
                min_alle_eqt_alle2[i] = 1;
            }
            else
            {
                minor_allele[i] = allele1[i];
                maf[i] = allele1_fre[i];
                min_alle_eqt_alle2[i] = 0;
            }
        }
		
    }
    delete [] batch_char;
    bedfile.close();

    cout << "done." << endl;
    writeLOG << "done." << endl;
	
    //sorting SNPs
    for(int i = 0; i < snp_num; ++ i)
    {
        chr_idx.push_back(make_pair(chrom[i], make_pair(position[i], i)));
    }
	
    vector<int> chr_reflect;
    vector<int> one_chr_position_reflect;
    vector<string> one_chr_position;

    sort(chr_idx.begin(), chr_idx.end(), cmp_3arr_1<int, int, int>);

    int one_chr_num_begin = 0;
    int one_chr_num_end;
    int current_chr_num = chr_idx[0].first;
	
    for(int i = 0; i < snp_num; ++ i)
    {
        if(chr_idx[i].first != current_chr_num)
        {
            one_chr_num_end = i - 1;
            sort(chr_idx.begin() + one_chr_num_begin, chr_idx.begin() + one_chr_num_end + 1, cmp_3arr_2<int, int, int>);
            one_chr_num_begin = i;
            current_chr_num = chr_idx[i].first;
        }
    }
	
    sort(chr_idx.begin() + one_chr_num_begin, chr_idx.end(), cmp_3arr_2<int, int, int>);
	
// SNPs distance filter
    int thr_pos_del_num = 0;
    int t = 0;
	
    if(thre_num == 0)
    {
        for(int i = 0; i < snp_num - 1; ++ i)
        {
            if(chr_idx[t].first == chr_idx[i + 1].first && (chr_idx[i + 1].second.first - chr_idx[t].second.first) == thre_num)
            {
                thr_pos_del_num = thr_pos_del_num + 1;
                snp_to_delete[chr_idx[i + 1].second.second] = 1;
            }

            else
            {
                t = i + 1;
            }
        }
		
        cout << thr_pos_del_num << " SNPs removed due to distances from previous ones = " << thre_num << ".\n";
        writeLOG << thr_pos_del_num << " SNPs removed due to distances from previous ones = " << thre_num << ".\n";
    }
    else if(thre_num > 0)
    {
        for(int i = 0; i < snp_num - 1; ++ i)
        {
            if(chr_idx[t].first == chr_idx[i + 1].first && (chr_idx[i + 1].second.first - chr_idx[t].second.first) < thre_num)
            {
                thr_pos_del_num = thr_pos_del_num + 1;
                snp_to_delete[chr_idx[i + 1].second.second] = 1;
            }
			
            else
            {
                t = i + 1;
            }
        }
        cout << thr_pos_del_num << " SNPs removed due to distances from previous ones < " << thre_num << ".\n";
        writeLOG << thr_pos_del_num << " SNPs removed due to distances from previous ones < " << thre_num << ".\n";
    }
	
    int low_MAF_del_num = 0;
    int high_MAF_del_num = 0;
    int sex_chrom_snps_num = 0;
    int miss_fre_del_num = 0;
    remain_snp_num = 0;
    for(int i = 0; i < snp_num; ++ i)
    {
// MAF filter
        if(maf[i] < low_MAF)
        {
                low_MAF_del_num = low_MAF_del_num + 1;
                snp_to_delete[i] = 1;
        }
        if(maf[i] > high_MAF)
        {
                high_MAF_del_num = high_MAF_del_num + 1;
                snp_to_delete[i] = 1;
        }

// Sex chrom filter
        if(chrom[i] > 22 || chrom[i] < 1)
        {
            snp_to_delete[i] = 1;
            ++ sex_chrom_snps_num;
        }

// Missing rate filter
        if(miss_fre[i] > miss_max)
        {
            snp_to_delete[i] = 1;
            ++ miss_fre_del_num;
        }

        if(!snp_to_delete[i])
        {
            ++ remain_snp_num;
        }
    }

    cout << low_MAF_del_num << " SNPs removed due to MAFs < " << low_MAF << ".\n";
    writeLOG << low_MAF_del_num << " SNPs removed due to MAFs < " << low_MAF << ".\n";
    cout << high_MAF_del_num << " SNPs removed due to MAFs > " << high_MAF << ".\n";
    writeLOG << high_MAF_del_num << " SNPs removed due to MAFs > " << high_MAF << ".\n";	
    cout << sex_chrom_snps_num << " sex chromosome SNPs removed. " << endl;
    writeLOG << sex_chrom_snps_num << " sex chromosome SNPs removed. " << endl;
    cout << miss_fre_del_num << " SNPs removed due to missing rate > " << miss_max << ".\n";
    writeLOG << miss_fre_del_num << " SNPs removed due to missing rate > " << miss_max << ".\n";
	
    if(remain_snp_num == 0)
    {
        cout << "ERROR: 0 SNPs pass filters." << endl;
        writeLOG << "ERROR: 0 SNPs pass filters." << endl;
        exit(0);
    }

    cout << remain_snp_num << " SNPs passed filters." << endl;
    writeLOG << remain_snp_num << " SNPs passed filters." << endl;
} 


void ASA::allele_info_output(string ofile)
{
    ofstream writefile;
    writefile.open((ofile + ".frq").c_str(), ios::out);

    cout << "Saving allele information to [" + ofile + ".frq] file... " << flush;
    writeLOG << "Saving allele information to [" + ofile + ".frq] file... ";

    writefile << "chrom\t" << "pos\t" << "SNP\t" << "A1\t" << "A2\t" << "A1_num\t" << "A2_num\t" << "A1_freq\t" << "A2_freq\t" << "minor_allele\t" << "MAF\t" << "missing_rate\t" << "filter_flag" << endl;
	
    for(int i = 0; i<snp_num; ++ i)
    {
        writefile << chrom[i] << "\t";
        writefile << position[i] << "\t";
        writefile << snp[i] << "\t";
        writefile << allele1[i] << "\t";
        writefile << allele2[i] << "\t";
        writefile << allele1_num[i] << "\t";
        writefile << allele2_num[i] << "\t";
        writefile << allele1_fre[i] << "\t";
        writefile << allele2_fre[i] << "\t";
        writefile << minor_allele[i] << "\t";
        writefile << maf[i] << "\t";
	writefile << miss_fre[i] << "\t";
        writefile << snp_to_delete[i] << endl;
    }
	 
    cout << "done." << endl;
    writeLOG << "done." << endl;
    writefile.close();
}


template <typename T> void ASA::sort_indexes(const vector<T> &v, vector<int> &idx)
{
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
}


template <typename T1, typename T2, typename T3>
bool ASA::cmp_3arr_1(const pair<T1, pair<T2, T3> > &a, const pair<T1, pair<T2, T3> > &b) 
{
    return a.first < b.first;
}


template <typename T1, typename T2, typename T3>
bool ASA::cmp_3arr_2(const pair<T1, pair<T2, T3> > &a, const pair<T1, pair<T2, T3> > &b) 
{
    return a.second.first < b.second.first;
}


template <typename T1, typename T2>
bool ASA::cmp_2arr_1(const pair<T1, T2> &a, const pair<T1, T2>  &b)
{
    return a.first < b.first;
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

                for (int j = 0; j < indi_num; ++ j)
                {
                    X_bool[j][k] = 1;
                }
                continue;
            }
			 
            eigenVar mu, inv_sd_snp, batch_X_0, batch_X_1, batch_X_2;
            batch_X_0 = 0;
            batch_X_1 = 1;
            batch_X_2 = 2;
            if (centralize)
            {
                mu = 2 * maf[i];
                batch_X_0 = 0 - mu;
                batch_X_1 = 1 - mu;
                batch_X_2 = 2 - mu;
            }
            if (standardize)
            {
                mu = 2 * maf[i];
                if(maf[i] < 1.0e-50)
                {
                    inv_sd_snp = 0.0;
                }
                else
                {
                    inv_sd_snp = sqrt(1 / (2 * maf[i] * (1 - maf[i])));
                }
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
                    if(!min_alle_eqt_alle2[i])
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
                    if(!min_alle_eqt_alle2[i])
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
            if (snp_to_delete[i])
            {
                filtered_num = filtered_num + 1;
            }
        }
		
        #pragma omp parallel for
        for (int i = 0; i < indi_num; ++ i)
        {
            for (int j = 0; j < batch_snp_size; ++ j)
            {
                if(X_bool[i][j] == 0)
                {
                    batch_miss_pos[i].push_back(j);
                }
            }
        }
		
        #pragma omp parallel for
        for (int i = 0; i < indi_num; ++ i)
        {
            for (int j = 0; j <= i; ++ j)
            {
                int miss_j = 0;
                for (int k = 0; k < batch_miss_pos[j].size(); ++ k)
                {
                    miss_j = miss_j + X_bool[i][batch_miss_pos[j][k]];
                }
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
    if (ifile_flag && not_in_info_num > 0)
    {
	cout << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
	writeLOG << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
    }
    else
    {
	cout << "Genotype data of " << remain_snp_num << " SNPs were used." << endl;
	writeLOG << "Genotype data of " << remain_snp_num << " SNPs were used." << endl;
    }
}


void ASA::grm_id_output(string ofile)
{	
    ofstream writefile;
    writefile.open((ofile + ".grm.id").c_str(), ios::out);

    cout << "Saving IDs of the GRM [" + ofile + ".grm.gz] to [" + ofile + ".grm.id]... " << flush;
    writeLOG << "Saving IDs of the GRM [" + ofile + ".grm.gz] to [" + ofile + ".grm.id]... ";

    for(int i = 0; i<indi_num; ++ i)
    {
        writefile << famid[i] << "\t" << subid[i] << endl;
    }

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
	
    for(int i = 0; i < indi_num; ++ i)
    {
        for(int j = 0; j <= i; ++ j)
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
	
    for (int i = indi_num - 1; i  >= 0; --i)
    {
        sum_eigenval = sum_eigenval + Z_D(i);
    }

    for(int j = indi_num - 1; j >= indi_num - N; --j)
    {
        writefile1 << Z_D(j) << "\t" << Z_D(j) / sum_eigenval << endl;
    }
	
    writefile1.close();
    cout << "done." << endl;
    writeLOG << "done." << endl;

    cout << "Saving top " << N << " eigenvectors and PCs of " << indi_num << " individuals to [" + eigvec_ofile + ".eigenvec] and [" + pc_ofile + ".pc] files... " << flush;
    writeLOG << "Saving top " << N << " eigenvectors and PCs of " << indi_num << " individuals to [" + eigvec_ofile + ".eigenvec] and [" + pc_ofile + ".pc] files... ";
	
    for(int i = 0; i < indi_num; ++ i)
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

    if (info_col_num < 3)
    {
	cout << "ERROR: " << ifile << " has " << info_col_num << "columns. Three columns are expeted: the first is SNP ID, the second is population that the SNP is polymorphic in, and the third is MAF of the SNP in the population." << endl;
	writeLOG << "ERROR: " << ifile << " has " << info_col_num << "columns. Three columns are expeted: the first is SNP ID, the second is population that the SNP is polymorphic in, and the third is MAF of the SNP in the population." << endl;
	exit(0);
    }
}


void ASA::snp_info_input(string ifile)
{
    cout << "Reading SNP information file [" << ifile << "]... " << flush;
    writeLOG << "Reading SNP information file [" << ifile << "]... " ;
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
	info_pop.push_back(parsedline[1]);
	info_pop_set.insert(parsedline[1]);
	info_pop_multiset.insert(parsedline[1]);
	info_pop_maf.push_back(atof(parsedline[2].c_str()));
        ++ info_snp_num;
    }
    infofile.close(); 	

    cout << "done." << endl;
    writeLOG << "done." << endl;
    cout << info_snp_num << " SNPs loaded from the information file [" << ifile << "]." << endl;
    writeLOG << info_snp_num << " SNPs loaded from the information file [" << ifile << "]." << endl;
    int info_pop_num = info_pop_set.size();
    cout << info_pop_num << " populations found in the information file, population ID and number of the population-specific SNPs are " << endl;
    writeLOG << info_pop_num << " populations found in the information file, population ID and number of the population-specific SNPs are " << endl;

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
	if (info_pop_maf[i] <= 0 || info_pop_maf[i] > 0.5)
	{
	    cout << "ERROR: MAF value " << info_pop_maf[i] << " was found for SNP " << info_snp[i]  << "." << endl;
	    writeLOG << "ERROR: MAF value " << info_pop_maf[i] << " was found for SNP " << info_snp[i] << "." << endl;
	    exit(0);
	}
    }
    
    info_in_bim.resize(info_snp_num, 0);
    not_in_info_num = 0;
    for(int i = 0; i < snp_num; ++ i)
    {
	if (info_pop_snp.count(snp[i]) == 0)
	{
	    not_in_info_num ++;
	    pop_maf.push_back(0);
	    pop.push_back("NA");
	    if (!snp_to_delete[i])
	    {
		snp_to_delete[i] = 1;
		remain_snp_num --;
	    }
	    continue;
	}
	int idx = info_pop_snp.find(snp[i])->second;
	pop_maf.push_back(info_pop_maf[idx]);
	info_in_bim[idx] = 1;
	pop.push_back(info_pop[idx]);
	if(!snp_to_delete[i])
	{
	    pop_set.insert(info_pop[idx]);
	}
    }

    pop_num = pop_set.size();
    for (it = pop_set.begin(); it != pop_set.end(); ++ it)
    {
	diff_pop.push_back(*it);
    }
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
    for(int j = 0; j < pop_num; ++ j)
    {
        for(int i = 0; i < snp_num; ++ i)
        {
            if(!snp_to_delete[i])
            {
                if(pop[i] == diff_pop[j])
                {
                    b(i, j) = 2 * pop_maf[i];
                    c(0, j) = c(0, j) + 4 * pop_maf[i] * pop_maf[i];
                }
            }
        }
	for(int i = 0; i < info_snp_num; ++ i)
	{
	    if(!info_in_bim[i])
	    {
		if(info_pop[i] == diff_pop[j])
	        {
		    c(0, j) = c(0, j) + 4 * info_pop_maf[i] * info_pop_maf[i];
		}
	    }
	}
        for(int k = 1; k < indi_num; ++ k)
        {
            c(k, j) = c(0, j);
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
                for (int j = 0; j < indi_num; ++ j)
                {
                    missing_bool[j][k] = 0;
                }
                continue;
            }
			 
            int pos = (line_byte_size) * (i - snp_start) * 8;
            for (int j = 0; j < indi_num; ++ j)
            {
                bool tmp1bool = batch_bit[pos ++ ];
                bool tmp2bool = batch_bit[pos ++ ];
                if (!tmp1bool && !tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
                    {
                        batch_X(j, k) = 0;
                    }
		    else
                    {
                        batch_X(j, k) = 2;
                    }
		    continue;
                }
				
                if (!tmp1bool && tmp2bool)
                {
                    batch_X(j, k) = 1;
		    continue;
                }
				
                if (tmp1bool && tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
                    {
                        batch_X(j, k) = 2;
                    }
                    else
                    {
                        batch_X(j, k) = 0;
                    }
		    continue;
                }
				
                batch_X(j, k) = 0;
                missing_bool[j][k] = 1;
            }
        }
	
        psv = psv + batch_X * (b.middleRows(snp_start, batch_snp_size));

        #pragma omp parallel for
        for(int k = 0; k < pop_num; ++ k)
        {
            for(int i = 0; i < indi_num; ++ i)
            {
                for(int j = 0; j < batch_snp_size; ++ j)
                {
                    if(missing_bool[i][j])
                    {
                        if(b(snp_start + j, k) != 0)
                        {
                            c(i, k) = c(i, k) - b(snp_start + j, k) * b(snp_start + j, k);
                        }
                    }
                }
            }
        }
		
    }
	
    #pragma omp parallel for
    for (int j = 0; j < pop_num; ++ j)
    {
        for (int i = 0; i < indi_num; ++ i)
        {
            psv(i, j) = psv(i, j) / c(i, j);
        }
    }

    delete [] batch_char;
    bedfile.close();
    cout << "done." << endl;
    writeLOG << "done." << endl;
    if (not_in_info_num > 0)
    {
	cout << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
	writeLOG << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
    }
    cout << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
    writeLOG << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
}


void ASA::psv_output(string ofile)
{
    ofstream writefile;
    writefile.open((ofile + ".psv").c_str(), ios::out);

    cout << "Saving the principal score vectors to [" + ofile + ".psv]... " << flush;
    writeLOG << "Saving the principal score vectors to [" + ofile + ".psv]... ";

    writefile << "ID_1" << "\t" << "ID_2";

    for (int j = 0; j < pop_num; ++ j)
    {
        writefile << "\t" << diff_pop[j];
    }
    writefile << endl;
	
    for(int i = 0; i < indi_num; ++ i)
    {
        writefile << famid[i] << "\t" << subid[i];
        for (int j = 0; j < pop_num; ++ j)
        {
            writefile << "\t" << psv(i, j);
        }
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
    for(int j = 0; j < pop_num; ++ j)
    {
        for(int i = 0; i < snp_num; ++ i)
        {
            if(!snp_to_delete[i])
            {
                if(pop[i] == diff_pop[j])
                {
                    b(i, j) = 1;
                    c(0, j) = c(0, j) + 2 * pop_maf[i];
                }
            }
        }
	for(int i = 0; i < info_snp_num; ++ i)
	{
	    if(!info_in_bim[i])
	    {
		if(info_pop[i] == diff_pop[j])
		{
		    c(0, j) = c(0, j) + 2 * info_pop_maf[i];
		}
	    }
	}
        for(int k = 1; k < indi_num; ++ k)
        {
            c(k, j) = c(0, j);
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
                for (int j = 0; j < indi_num; ++ j)
                {
                    missing_bool[j][k] = 0;
                }
                continue;
            }
			 
            int pos = (line_byte_size) * (i - snp_start) * 8;
            for (int j = 0; j < indi_num; ++ j)
            {
                bool tmp1bool = batch_bit[pos ++ ];
                bool tmp2bool = batch_bit[pos ++ ];
                if (!tmp1bool && !tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
                    {
                        batch_X(j, k) = 0;
                    }
		    else
                    {
                        batch_X(j, k) = 2;
                    }
		    continue;
                }
				
                if (!tmp1bool && tmp2bool)
                {
                    batch_X(j, k) = 1;
		    continue;
                }
				
                if (tmp1bool && tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
                    {
                        batch_X(j, k) = 2;
                    }
                    else
                    {
                        batch_X(j, k) = 0;
                    }
		    continue;
                }
				
                batch_X(j, k) = 0;
                missing_bool[j][k] = 1;
            }
        }
	
        aiv = aiv + batch_X * (b.middleRows(snp_start, batch_snp_size));

        #pragma omp parallel for
        for(int k = 0; k < pop_num; ++ k)
        {
            for(int i = 0; i < indi_num; ++ i)
            {
                for(int j = 0; j < batch_snp_size; ++ j)
                {
                    if(missing_bool[i][j])
                    {
			if(b(snp_start + j, k) != 0)
			{
			    c(i, k) = c(i, k) - 2 * pop_maf[snp_start + j];
			}
                    }
                }
            }
        }
		
    }
	
    #pragma omp parallel for
    for (int j = 0; j < pop_num; ++ j)
    {
        for (int i = 0; i < indi_num; ++ i)
        {
            aiv(i, j) = aiv(i, j) / c(i, j);
        }
    }

    delete [] batch_char;
    bedfile.close();
    cout << "done." << endl;
    writeLOG << "done." << endl;
    if (not_in_info_num > 0)
    {
	cout << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
	writeLOG << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
    }
    cout << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
    writeLOG << "Genotype data of " << remain_snp_num << " SNPs in the BED file and information of " << info_snp_num << " SNPs in the information file were used." << endl;
}



void ASA::aiv_output(string ofile)
{
    ofstream writefile;
    writefile.open((ofile + ".aiv").c_str(), ios::out);

    cout << "Saving the ancestral information vectors to [" + ofile + ".aiv]... " << flush;
    writeLOG << "Saving the ancestral information vectors to [" + ofile + ".aiv]... ";

    writefile << "ID_1" << "\t" << "ID_2";

    for (int j = 0; j < pop_num; ++ j)
    {
        writefile << "\t" << diff_pop[j];
    }
    writefile << endl;
	
    for(int i = 0; i < indi_num; ++ i)
    {
        writefile << famid[i] << "\t" << subid[i];
        for (int j = 0; j < pop_num; ++ j)
        {
            writefile << "\t" << aiv(i, j);
        }
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
    for(int j = 0; j < pop_num; ++ j)
    {
        for(int i = 0; i < snp_num; ++ i)
        {
            if(!snp_to_delete[i])
            {
                if(pop[i] == diff_pop[j])
                {
		    a(i, j) = 4 * pop_maf[i] * pop_maf[i];
		    b(i, j) = 2 * pop_maf[i];
		    c(i, j) = 1;
                }
            }
        }
	for(int i = 0; i < info_snp_num; ++ i)
	{
	    if(!info_in_bim[i])
	    {
		if(info_pop[i] == diff_pop[j])
		{
		    aiv_a(0, j) = aiv_a(0, j) + 4 * info_pop_maf[i] * info_pop_maf[i];
		    aiv_b(0, j) = aiv_b(0, j) + 2 * info_pop_maf[i];
		}
	    }
	}
        for(int k = 1; k < indi_num; ++ k)
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
		
		
//        vector<vector<int>> missing_bool(indi_num, vector<int>(batch_snp_size, 0));
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
//                for (int j = 0; j < indi_num; ++ j)
//                {
//                    missing_bool[j][k] = 0;
//                }
                continue;
            }
			 
            int pos = (line_byte_size) * (i - snp_start) * 8;
            for (int j = 0; j < indi_num; ++ j)
            {
                bool tmp1bool = batch_bit[pos ++ ];
                bool tmp2bool = batch_bit[pos ++ ];
                if (!tmp1bool && !tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
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
				
                if (!tmp1bool && tmp2bool)
                {
                    batch_X(j, k) = 1;
                    batch_1_X(j, k) = 0;
		    continue;
                }
				
                if (tmp1bool && tmp2bool)
                {
                    if(min_alle_eqt_alle2[i])
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
				
                batch_X(j, k) = 0;
                batch_1_X(j, k) = 0;
//                missing_bool[j][k] = 1;
            }
        }
	
        aiv_a = aiv_a + batch_1_X * (a.middleRows(snp_start, batch_snp_size));
        aiv_b = aiv_b + batch_1_X * (b.middleRows(snp_start, batch_snp_size));
	aiv_c = aiv_c + batch_X * (c.middleRows(snp_start, batch_snp_size));
    }

    #pragma omp parallel for
    for (int j = 0; j < pop_num; ++ j)
    {
        for (int i = 0; i < indi_num; ++ i)
        {
            aiv(i, j) = (sqrt(aiv_b(i, j) * aiv_b(i, j) + 4 * aiv_a(i, j) * aiv_c(i, j)) - aiv_b(i, j)) / aiv_a(i, j) / 2;
        }
    }

    delete [] batch_char;
    bedfile.close();
    cout << "done." << endl;
    writeLOG << "done." << endl;
    if (not_in_info_num > 0)
    {
	cout << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
	writeLOG << not_in_info_num << " SNPs were excluded, which were not in the SNP panel defined by the information file." << endl;
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

    for (int j = 0; j < pop_num; ++ j)
    {
        writefile << "\t" << diff_pop[j];
    }
    writefile << endl;
	
    for(int i = 0; i < indi_num; ++ i)
    {
        writefile << famid[i] << "\t" << subid[i];
        for (int j = 0; j < pop_num; ++ j)
        {
            writefile << "\t" << aiv(i, j);
        }
        writefile << endl;
    }
	
    writefile.close();
    cout << "done." << endl;
    writeLOG << "done." << endl;
}

