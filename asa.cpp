#include "geno_calc.h"

void cover()
{
	cout << "+======================================+" << endl;
	cout << "|                                      |" << endl;
	cout << "|   Ancestral Spectrum Analyzer (ASA)  |" << endl;
	cout << "|   version 1.0.0                      |" << endl;
	cout << "|                                      |" << endl;
	cout << "+======================================+" << endl;
}

void help()
{
	cover();
	cout << "--bfile\t\t" << "Input genotype file in plink binary format." << endl;
	cout << "--ifile\t\t" << "Input information file of a panel of population-specific SNPs. The file is expected to have three columns without headers:" << endl;
	cout << "\t\t" << "the first is SNP ID, the second is population that the SNP is specific to, and the third is MAF of the SNP in the population." << endl;
	cout << "--out\t\t" << "Output file name [default: asa]." << endl;
	cout << "--aiv\t\t" << "Calculate and output ancestral information vectors by the method of moment estimate." << endl;
	cout << "--mle\t\t" << "Calculate and output ancestral information vectors by the approximate maximum likelihood estimate." << endl;
	cout << "--psv\t\t" << "Calculate and output principal score vectors." << endl;
	cout << "--pca\t\t" << "Perform PCA and output top n eigenvectors, eigenvalues and PCs [default: n = 20]." << endl;
	cout << "--freq\t\t" << "Calculate and output allele frequencies." << endl;
	cout << "--grm\t\t" << "Calculate and output genetic relationship matrix." << endl;
	cout << "--maf-min\t" << "Exclude SNPs with MAFs smaller than the specified value [default: 0]." << endl;
	cout << "--maf-max\t" << "Exclude SNPs with MAFs larger than the specified value [default: 0.5]." << endl;
	cout << "--dist-min\t" << "Exclude SNPs with distances from previous ones less or equal to the specified value [default: 0]." << endl;
	cout << "--miss-max\t" << "Exclude SNPs with missing rates larger than the specified value [default: 1]." << endl;
//	cout << "--centralize\t" << "Centralize genotype matrix by SNPs." << endl;
//	cout << "--standardize\t" << "Standardize genotype matrix by SNPs." << endl;
	cout << "--batch-size\t" << "Number of SNPs to be processed in a batch [default: 10000]." << endl;
	cout << "--thread-num\t" << "Number of threads on which the program will be running [default: thread number in your machine - 1]." << endl;
}


int main(int argc, char* argv[])
{
    double main_start_time;
    double main_end_time;
    double main_total_time;
    main_start_time = omp_get_wtime();
    time_t main_current_time;
    char* main_time_int_to_string;
    main_current_time = time(0);
    main_time_int_to_string = ctime(&main_current_time);
    
    int thread_num, pca = 20, max_thread_num = std::thread::hardware_concurrency(), batch_size = 10000;
    double dist_min = 0.0, maf_min = 0.0, maf_max = 0.5, miss_max = 1;
    string bfile, out_file = "asa", ifile;
    bool grm_flag = false, freq_flag = false, std_flag = false, ctr_flag = false, pca_flag = false, bfile_flag = false, maf_min_flag = false, maf_max_flag = false, batch_size_flag = false; 
    bool thread_num_flag = false, dist_min_flag = false, miss_max_flag = false, outfile_flag = false, psv_flag = false, aiv_flag = false, mle_flag = false, ifile_flag = false, analysis_flag = false;
    thread_num = max_thread_num - 1;

    ASA data;

    int opt, longindex;
    const char *optstring = "";
    struct option long_options[] =
    {
	{ "bfile", required_argument,  NULL, 0},
	{ "ifile", required_argument,  NULL, 1},
	{ "out", required_argument, NULL, 2},
	{ "aiv", no_argument, NULL, 3},
	{ "mle", no_argument , NULL, 4},
	{ "psv", no_argument, NULL, 5},
	{ "pca", required_argument, NULL, 6},
	{ "freq", no_argument, NULL, 7},
	{ "grm", no_argument, NULL, 8},
	{ "maf-min", required_argument, NULL, 9},
	{ "maf-max", required_argument, NULL, 10},
	{ "dist-min", required_argument, NULL, 11},
	{ "miss-max", required_argument, NULL, 12},
	{ "centralize", no_argument, NULL, 13},
	{ "standardize", no_argument, NULL, 14},
	{ "batch-size", required_argument , NULL, 15},
	{ "thread-num", required_argument , NULL, 16},
	{ "help", no_argument, NULL, 17},
	{0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1)
    {
	switch (opt)
	{
	case 0: bfile = optarg; bfile_flag = true; break;
	case 1: ifile = optarg; ifile_flag = true; break;
	case 2: out_file = optarg; outfile_flag = true; break;
	case 3: aiv_flag = true; analysis_flag = true; break;
	case 4: mle_flag = true; analysis_flag = true; break;
	case 5: psv_flag = true; analysis_flag = true; break;
	case 6: pca = atoi(optarg); pca_flag = true; analysis_flag = true; break;
	case 7: freq_flag = true; analysis_flag = true; break;
	case 8: grm_flag = true; analysis_flag = true; break;
	case 9: maf_min = atof(optarg); maf_min_flag = true; break;
	case 10: maf_max = atof(optarg); maf_max_flag = true; break;
	case 11: dist_min = atof(optarg); dist_min_flag = true; break;
	case 12: miss_max = atof(optarg); miss_max_flag = true; break;
	case 13: ctr_flag = true; break;
	case 14: std_flag = true; break;
	case 15: batch_size = atoi(optarg); batch_size_flag = true; break;
	case 16: thread_num = atoi(optarg); thread_num_flag = true; break;
	case 17:  help(); exit(0); break;
	default: break;
	}
    }

    cover();
    data.writeLOG.open((out_file + ".log").c_str(), ios::out);
    cout << "Options specified:" << endl;
    data.writeLOG << "Options specified:" << endl;
	
    if (bfile_flag)
    {
	cout << "--bfile " << bfile << endl;
	data.writeLOG << "--bfile " << bfile << endl;
    }
    else
    {
	cout << "ERROR: input genotype file was not specified." << endl;
	data.writeLOG << "ERROR: input genotype file not specified." << endl;
	exit(0);
    }
    if(ifile_flag)
    {
	cout << "--ifile " << ifile << endl;
	data.writeLOG << "--ifile " << ifile << endl;
    }
    cout << "--out " << out_file << endl;
    data.writeLOG << "--out " << out_file << endl;
    if (pca_flag)
    {
	if (pca <= 0)
	{
	    cout << "ERROR: --pca should be larger than 0." << endl;
	    data.writeLOG << "ERROR: --pca should be larger than 0." << endl;
	    exit(0);
	}
	else
	{
	    cout << "--pca " << pca << endl;
	    data.writeLOG << "--pca " << pca << endl;
	}
    }
    if (psv_flag)
    {
	cout << "--psv " << endl;
	data.writeLOG << "--psv " << endl;
    }
    if (aiv_flag)
    {
	cout << "--aiv " << endl;
	data.writeLOG << "--aiv " << endl;
    }
    if (mle_flag)
    {
	cout << "--mle " << endl;
	data.writeLOG << "--mle " << endl;
    }
    if (freq_flag)
    {
	cout << "--freq " << endl;
	data.writeLOG << "--freq " << endl;
    }
    if (grm_flag)
    {
	cout << "--grm " << endl;
	data.writeLOG << "--grm " << endl;
    }
    if (maf_min_flag)
    {
	if (maf_min > 0 && maf_min < 0.5)
	{
	    cout << "--maf-min " << maf_min << endl;
	    data.writeLOG << "--maf-min " << maf_min << endl;
	}
	else
	{
	    cout << "Error: --maf-min should be between 0 and 0.5." << endl;
	    data.writeLOG << "Error: --maf-min should be between 0 and 0.5." << endl;
	    exit(0);
	}
    }
    if (maf_max_flag)
    {
	if (maf_max > 0 && maf_max < 0.5)
	{
	    cout << "--maf-max " << maf_max << endl;
	    data.writeLOG << "--maf-max " << maf_max << endl;
	}
	else
	{
	    cout << "Error: --maf-max should be between 0 and 0.5." << endl;
	    data.writeLOG << "Error: --maf-max should be between 0 and 0.5." << endl;
	    exit(0);
	}
    }
    if (dist_min_flag)
    {
	if (dist_min > 0)
	{
	    cout << "--dist_min " << dist_min << endl;
	    data.writeLOG << "--dist_min " << dist_min << endl;
	}
	else
	{
	    cout << "Error: --dist_min should be larger than 0." << endl;
	    data.writeLOG << "Error: --dist_min should be larger than 0." << endl;
	    exit(0);
	}
    }
    if (miss_max_flag)
    {
	if (miss_max > 0 && miss_max <= 1)
	{
	    cout << "--miss_max " << miss_max << endl;
	    data.writeLOG << "--miss_max " << miss_max << endl;
	}
	else
	{
	    cout << "ERROR: --miss-max should be between 0 and 1." << endl;
	    data.writeLOG << "ERROR: --miss-max should be between 0 and 1." << endl;
	    exit(0);
	}
    }			
    if (ctr_flag)
    {
	    cout << "--centralize " << endl;
	    data.writeLOG << "--centralize " << endl;
    }
    if (std_flag)
    {
	    cout << "--standardize " << endl;
	    data.writeLOG << "--standardize " << endl;
    }		
    if (batch_size_flag)
    {
	if (batch_size > 1)
	{
	    cout << "--batch-size " << batch_size << endl;
	    data.writeLOG << "--batch-size " << batch_size << endl;
	}
	else
	{
	    cout << "Error: --batch-size should be larger than 1." << endl;
	    data.writeLOG << "Error: --batch-size should be larger than 1." << endl;
	    exit(0);
	}
    }
    if (thread_num_flag)
    {
	if (thread_num > 0 && thread_num <= max_thread_num)
	{
	    cout << "--thread-num " << thread_num << endl;
	    data.writeLOG << "--thread-num " << thread_num << endl;
	}
	else
	{
	    cout << "Error: --thread-num should be between 1 and " << max_thread_num << endl;
	    data.writeLOG << "Error: --thread-num should be between 1 and " << max_thread_num << endl;
	    exit(0);
	}
    }

    if(!analysis_flag)
    {
        cout << "ERROR: no analysis has been launched by the options." << endl;
        data.writeLOG << "ERROR: no analysis has been launched by the options." << endl;
        exit(0);
    }
	
    if(maf_max <= maf_min)
    {
        cout << "ERROR: --maf-max should be larger than --maf-min. " << endl;
        data.writeLOG << "ERROR: --maf-max should be larger than --maf-min. " << endl;
        exit(0);
    }
	
    if(std_flag && ctr_flag)
    {
        cout << "ERROR: choose --centralize or --standardize in one analysis." << endl;
        data.writeLOG << "ERROR: choose --centralize or --standardize in one analysis." << endl;
        exit(0);
    }
	
    if(psv_flag && !ifile_flag)
    {
        cout << "ERROR: information file of a panel of population-specific SNPs is necessary for calculating the population score vectors, use --ifile to specify it." << endl;
        data.writeLOG << "ERROR: information file of a panel of population-specific SNPs is necessary for calculating the population score vectors, use --ifile to specify it." << endl;
        exit(0);
    }
	
    if((aiv_flag || mle_flag) && !ifile_flag)
    {
	cout << "ERROR: information file of a panel of population-specific SNPs is necessary for calculating the ancestral information vectors, use --ifile to specify it." << endl;
	data.writeLOG << "ERROR: information file of a panel of population-specific SNPs is necessary for calculating the ancestral information vectors, use --ifile to specify it.." << endl;
	exit(0);
    }


    cout << "Start time: " << main_time_int_to_string;
    data.writeLOG << "Start time: " << main_time_int_to_string;

    stringstream ss;
    ss << thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
    omp_set_num_threads(thread_num);
    cout << "The program will be running on " << thread_num << " threads." << endl;
    data.writeLOG << "The program will be running on " << thread_num << " threads." << endl;

    data.allele_info_stat(bfile, dist_min, maf_min, maf_max, miss_max, batch_size);    
    
    if(ifile_flag)
    {
	data.snp_info_check(ifile);
	data.snp_info_input(ifile);
    }

    if(freq_flag)
    {
	data.allele_info_output(out_file);
    }

    if(psv_flag)
    {
        data.make_X_psv(bfile, batch_size);
        data.psv_output(out_file);
    }

    if(aiv_flag)
    {
	data.make_X_aiv(bfile, batch_size);
	data.aiv_output(out_file);
    }

    if(mle_flag)
    {
	data.make_X_mle(bfile, batch_size);
	data.mle_output(out_file);
    }
	
    if(pca_flag)
    {

        data.make_X_Z(bfile, batch_size, ctr_flag, std_flag, ifile_flag);

	if(grm_flag)
        {
            data.grm_output(out_file);
            data.grm_id_output(out_file);
        }

        data.grm_eigen();
        data.grm_eigvec_pc_eigval_output(out_file, out_file, out_file, pca);

    }
	
	
    if(grm_flag && !pca_flag && !psv_flag && !aiv_flag && !mle_flag)
    {
        data.make_X_Z(bfile, batch_size, ctr_flag, std_flag, ifile_flag);
        data.grm_output(out_file);
        data.grm_id_output(out_file);
    }

    main_current_time = time(0);
    main_time_int_to_string = ctime(&main_current_time);
    main_end_time = omp_get_wtime( );
    main_total_time = main_end_time - main_start_time;
    cout << "End time: " << main_time_int_to_string;
    data.writeLOG << "End time: " << main_time_int_to_string;
    cout << "Program running time: " << main_total_time << endl;
    return 0;
}
