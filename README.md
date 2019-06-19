# bestbamhit
Utility for finding the best STAR rna-seq mapping among mappings to several genomes.

    usage: bestbamhit [options] a.bam b.bam ...
      -edit-penalty float
        	multiple for how to penalize edit distance (default 2)
      -keep string
        	file where to write the names of reads matching the first bam file
      -labels string
        	comma-separated list of labels for the BAMs (required)
      -limit int
        	limit the number of sample reads considered (0 = no limit)
      -log string
        	write parameters and stats to a log file
      -max-dist int
        	max edit distance for an alignment (default 5)
      -min-len int
        	min length for an alignment (default 60)
