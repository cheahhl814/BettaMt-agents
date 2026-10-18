usage: 
###  Embryophyta plant plastome, 2*(1G raw data, 150 bp) reads
get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_seed.fasta -o plastome_output  -R 15 -k 21,45,65,85,105 -F embplant_pt
###  Embryophyta plant mitogenome
get_organelle_from_reads.py -1 sample_1.fq -2 sample_2.fq -s mt_seed.fasta -o mitogenome_output  -R 30 -k 21,45,65,85,105 -F embplant_mt

GetOrganelle v1.7.7.1 get_organelle_from_reads.py assembles organelle genomes
from genome skimming data. Find updates in
https://github.com/Kinggerm/GetOrganelle and see README.md for more
information.

options:
  -1 FQ_FILE_1          Input file with forward paired-end reads
                        (*.fq/.gz/.tar.gz).
  -2 FQ_FILE_2          Input file with reverse paired-end reads
                        (*.fq/.gz/.tar.gz).
  -u UNPAIRED_FQ_FILES  Input file(s) with unpaired (single-end) reads.
  -o OUTPUT_BASE        Output directory.
  -s SEED_FILE          Input fasta format file as initial seed. Default:
                        /home/cheahhl814/.GetOrganelle/SeedDatabase/*.fasta
  -w WORD_SIZE          Word size (W) for extension. Default: auto-estimated
  -R MAX_ROUNDS         Maximum extension rounds (suggested: >=2). Default: 15
                        (embplant_pt)
  -F ORGANELLE_TYPE     Target organelle genome type(s): embplant_pt/other_pt/
                        embplant_mt/embplant_nr/animal_mt/fungus_mt/fungus_nr/
                        anonym/embplant_pt,embplant_mt/other_pt,embplant_mt,fu
                        ngus_mt
  --max-reads MAX_READS
                        Maximum number of reads to be used per file. Default:
                        1.5E7 (-F
                        embplant_pt/embplant_nr/fungus_mt/fungus_nr); 7.5E7
                        (-F embplant_mt/other_pt/anonym); 3E8 (-F animal_mt)
  --fast FAST_STRATEGY  ="-R 10 -t 4 -J 5 -M 7 --max-n-words 3E7 --larger-
                        auto-ws --disentangle-time-limit 360"
  -k SPADES_KMER        SPAdes kmer settings. Default: 21,55,85,115
  -t THREADS            Maximum threads to use. Default: 1
  -P PRE_GROUPED        Pre-grouping value. Default: 200000
  -v, --version         show program's version number and exit
  -h                    print brief introduction for frequently-used options.
  --help                print verbose introduction for all options.

