# Take Muscle as an example

# Step 1. Converting bams to juncs
for i in $(seq 1 10)
do
	sh leafcutter/scripts/bam2junc.sh muscle${i}.bam muscle${i}.junc
done


# Step 2. Intron clustering
ls *.junc > muscle_juncfiles.txt
python leafcutter/clustering/leafcutter_cluster.py -j muscle_juncfiles.txt -m 50 -o muscle -l 500000


# Step 3. Calculate the PCA 
python leafcutter/scripts/prepare_phenotype_table.py muscle_perind.counts.gz -p 10 # note: -p 10 specific you want to calculate for sQTL to use as covariates
sh muscle_perind.counts.gz_prepare.sh
head -6 muscle_perind.counts.gz.PCs > muscle_perind.counts.gz.PCs.PC5


# Step 4. sQTL detection
for j in $(seq 1 29)
do
	fastQTL.static --vcf muscle.filtered.vcf.gz --bed muscle_perind.counts.gz.qqnorm_chr${j}.gz --cov muscle_perind.counts.gz.PCs.PC5 --permute 1000 10000 --normal --out muscle_perind.permutation.chr${j} --chunk 1 1
done
