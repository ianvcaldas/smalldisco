
conda run -n smalldisco python smalldisco.py sirna -o example/sirna.bed -a example/c_elegans-subset.gtf example/bam
conda run -n smalldisco python smalldisco.py tail -o example/tails.tsv -g example/c_elegans-chrom2.fa example/sirna.bed example/bam
