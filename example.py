#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from olivar import build, tiling, save, validate


if __name__ == "__main__":
    build(
        fasta_path = 'example_input/EPI_ISL_402124.fasta', # Path to the fasta reference sequence.

        var_path = 'example_input/delta_omicron_loc.csv', # Optional, path to the csv file of SNP coordinates and frequencies.
        # Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.

        BLAST_db = 'example_input/Human/GRCh38_primary', # Optional, path to the BLAST database. 
        # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").

        out_path = './', # Output directory.
        title = 'olivar-ref', # Name of the Olivar reference file [olivar-ref].
        threads = 1 # Number of threads.
    )


    tiling(
        ref_path = 'example_output/olivar-ref.olvr', # Path to the Olivar reference file (.olvr).
        out_path = './', # Output directory.
        title = 'olivar-design', # Name of this design [olivar-design].
        
        max_amp_len = 420, # Maximum amplicon length [420].
        min_amp_len = 252, # Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        w_egc = 1, # Weight for extreme GC content [1.0].
        w_lc = 1, # Weight for low sequence complexity [1.0].
        w_ns = 1, # Weight for non-specificity [1.0].
        w_var = 1, # Weight for variations [1.0].
        
        temperature = 60, # PCR annealing temperature [60.0].
        salinity = 0.18, # Concentration of monovalent ions in units of molar [0.18].
        dG_max = -11.8, # Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC = 0.2, # Minimum GC content of a primer [0.2].
        max_GC = 0.75, # Maximum GC content of a primer [0.75].
        min_complexity = 0.4, # Minimum sequence complexity of a primer [0.4].
        max_len = 36, # Maximum length of a primer [36].
        check_var = True, # Filter out primer candidates with variations within 5nt of 3' end [False]. 
        # Setting check_var=True is not recommended when a lot of variations are provided, 
        # since this would significantly reduce the number of primer candidates. 

        fP_prefix = '', # Prefix of forward primer. Empty string '' by default.
        rP_prefix = '', # Prefix of reverse primer. Empty string '' by default.
        
        seed = 10, # Random seed for optimizing primer design regions and primer dimer [10].
        threads = 1, # Number of threads.
    )


    # save(
    #     design_out = 'example_output/olivar-design.olvd', # Path to the Olivar design file (.olvd)
    #     out_path = './', # Output directory.
    # )


    validate(
        primer_pool = 'example_output/olivar-design.csv', # Path to the csv file of a primer pool. 
        # Required columns: "amp_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer), "pool" (pool number, e.g., 1).

        pool = 1, # Primer pool number [1]. 
        
        BLAST_db = 'example_input/Human/GRCh38_primary', # Optional, path to the BLAST database. 
        # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        
        out_path = './', # Output directory. 
        title = 'olivar-val', # Name of validation.

        max_amp_len = 1500, # Maximum length of predicted non-specific amplicon [1500].
        temperature = 60, # PCR annealing temperature [60.0].
        threads = 1 # Number of threads.
    )


    validate(
        primer_pool = 'example_output/olivar-design.csv', # Path to the csv file of a primer pool. 
        # Required columns: "amp_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer), "pool" (pool number, e.g., 1).

        pool = 2, # Primer pool number [1]. 
        
        BLAST_db = 'example_input/Human/GRCh38_primary', # Optional, path to the BLAST database. 
        # Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        
        out_path = './', # Output directory. 
        title = 'olivar-val', # Name of validation.

        max_amp_len = 1500, # Maximum length of predicted non-specific amplicon [1500].
        temperature = 60, # PCR annealing temperature [60.0].
        threads = 1 # Number of threads.
    )
