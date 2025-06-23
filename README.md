# Neopestalotiopsis

# Create a working directory and move all FASTQ files there.
mkdir fungal_project
cd fungal_project
mkdir raw_data trimmed_data qc_reports
mv *_R1_001.fastq.gz *_R2_001.fastq.gz raw_data/

# Quality Control and Trimming (fastp)
conda install -c bioconda fastp
cd raw_data
for sample in *_R1_001.fastq.gz
do
  base=$(basename ${sample} _R1_001.fastq.gz)
  fastp -i ${base}_R1_001.fastq.gz -I ${base}_R2_001.fastq.gz \
        -o ../trimmed_data/${base}_R1_trimmed.fastq.gz \
        -O ../trimmed_data/${base}_R2_trimmed.fastq.gz \
        -h ../qc_reports/${base}_fastp.html -j ../qc_reports/${base}_fastp.json
done
fastp -i *_R1.fastq.gz -I *_R2.fastq.gz -o clean_R1.fastq.gz -O clean_R2.fastq.gz -h report.html -j report.json

# Mapping to Reference Genome (BWA + Samtools)
# Index Reference
bwa index ref_genome.fasta

# Align Reads

cd ../trimmed_data
mkdir ../aligned_bam
for sample in 73_S419 80_S420 86_S421 93_S422
do
  echo "Aligning $sample"

  bwa mem -t 4 GCNrosiarefgenome.fna \
    ${sample}_clean_R1_001.fastq.gz \
    ${sample}_clean_R2_001.fastq.gz | \
  samtools view -bS - | \
  samtools sort -o ${sample}_aligned_sorted.bam

  samtools index ${sample}_aligned_sorted.bam
done


bwa mem -t 4 ../ref_genome.fasta \
    ${base}_R1_trimmed.fastq.gz ${base}_R2_trimmed.fastq.gz | \
    samtools view -Sb - > ../aligned_bam/${base}.bam

# For singles 
bwa mem GCNrosiarefgenome.fna 73_S419_clean_R1_001.fastq.gz 73_S419_clean_R2_001.fastq.gz > 73_S419_alignment.sam
# Convert SAM to BAM
samtools view -bS 73_S419_alignment.sam > 73_S419_alignment.bam

# Sort the BAM file
samtools sort 73_S419_alignment.bam -o 73_S419_alignment_sorted.bam

# Index bam files
samtools index 73_S419_alignment_sorted.bam


    
# Sort and Index

cd ../aligned_bam
samtools sort -o sorted_${bam} $bam
samtools index sorted_${bam}

# Variant Calling (bcftools)
conda install -c bioconda bcftools

# Create combined pileup
samtools mpileup -uf reference/ref_genome.fasta aligned_bam/*_sorted.bam | \
    bcftools call -mv -Oz -o variants.vcf.gz

# Index VCF
bcftools index variants.vcf.gz

# Call variants (from previous alignment)
bcftools mpileup -f ref_genome.fasta aligned_bam/*.bam | bcftools call -mv -Oz -o all_samples.vcf.gz

# Filter variants
bcftools filter -e 'QUAL<30 || DP<10' all_samples.vcf.gz -Oz -o filtered.vcf.gz

# Convert to PHYLIP
vcftools --vcf filtered.vcf.gz --out snp_matrix --phylip

# Build SNP tree
fasttree -nt -gtr snp_matrix.phylip > snp_tree.nwk

# Detect Housekeeping Genes & Resistance Genes
# Create BLAST database
makeblastdb -in reference/ref_genome.fasta -dbtype nucl -out blast_results/ref_db

# Create query files
cat > blast_resistance_genes.fa << 'EOF'
>CYP51A
ATGGCTGTATCAGTTCGTTTCTGCTGCTCG...
>SdhB
ATGGACGTGAAGTGGAGGTCAGGATTTTGAG...
EOF
>
# Extract sequences

/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/makeblastdb -in 100contigs.fasta -dbtype nucl
/home/muzhinjin/tikafinal/ncbi-blast-2.15.0+/bin/blastn -query NeopestalotroseITS.fasta -db 100contigs.fasta -outfmt 6 -out 100contigITS.tsv
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx 100contigs.fasta
/home/muzhinjin/tikafinal/samtools-1.19.2/samtools faidx 100contigs.fasta NODE_499_length_3025_cov_4728.818860:450-995 > 100cintigITS_raw.fasta

# Run BLAST searches
blastn -query blast_resistance_genes.fa -db blast_results/ref_db \
       -out blast_results/resistance_blast.txt -evalue 1e-5 -outfmt 6

# Use blastn to search genes
makeblastdb -in ref_genome.fasta -dbtype nucl -out ref_db
blastn -query housekeeping.fasta -db ref_db -out hk_blast.txt -evalue 1e-5 -outfmt 6
blastn -query fungicide_genes.fasta -db ref_db -out resistance_blast.txt -evalue 1e-5 -outfmt 6
# House keeping genes
>ACT1_Fungal_Actin
ATGGTCGGTATGGGTCAAAAAGAAATGGTAAGGTTATGTCACCAACTGGGACGATATGGAGAAGATCTGGCATCACACCTTCTACAACGAGCTGGAA

>EF1A_Translation_Elongation_Factor_1_Alpha
ATGGGTAAGGAGGACAAGACTGGAAGTTGAAGAGGAGCTTGATTTCAAGAACATGACGCTGAAGGTTGACGATATTGATGGAATCGTCGGAGAGGACA

>GAPDH_Glyceraldehyde-3-phosphate_dehydrogenase
ATGGTTGTTGAGGTTGTTGTTGACGTTGGTGTTGGTGATGCTGTTGTTGCTGATGTTGTTGTTGTTGGTTCTGATGCTGTTGGTGTTGGTGGTCGAG

>TUB2_Beta_Tubulin
ATGCGTGAGTGCATCAAGAAGGACTACCTGGGCACCCACCTCGGAGCGTGTGTACAGTCCACCGCCGATGTCTACGTCGTACAGTTGGGAGCCG

>RPB2_RNA_polymerase_II_largest_subunit
ATGGGCAAGGGTCAGATTGATGTCTGGGTCGATTGTTGAGGCGAGTTCGGTGAGTTCGAGGCTGACTTCGAGGAGTTGATGACGAGCTTGAGG


# Fungicide ressistance genes 

>CYP51A_Sterol_14α-demethylase_Azole_target
ATGGCTGTATCAGTTCGTTTCTGCTGCTCGATGTGTCTCGTGCTGGTGCTGTTGGCGGCTCGTTGTCTCGTTGGAGCTGCTGTTGTTGCTGGT

>SdhB_Succinate_Dehydrogenase_Subunit_B_SDHI_target
ATGGACGTGAAGTGGAGGTCAGGATTTTGAGGTCGTTGAGCTGGTTCTGTTCTGGTGGTGGTGATGTGGTGTCTGCTGATGGTCAGGATGTGGAC

>FKS1_1,3-beta-D-glucan_synthase_Echinocandin_target
ATGGTTGGTGCTGGTTGCTGGTGGTTGCTGATGCTGCTGTTGCTGCTGCTGGTCAGGGTTCTGGTGTTGTTGTTGGTTCTGTTGCTGATGTT

>ABC1_ABC_Transporter_Multidrug_Resistance
ATGGCGCTGGTTGTTGGTTCTGATGTTGTTGATGCTGGTGTTGCTGCTGTTGCTGGTGCTGCTGGTGCTGGTGGTGCTGGTGCTGCTGGTTCT

>MDR1_Major_Facilitator_Superfamily_MDR
ATGGTGCTGGTTGTTGCTGCTGCTGGTGTTGTTGTTGCTGGTCAGGTTGGTGCTGGTCAGGGTGTTGGTGGTCAGGGTCAGGGTGGTCAG


  
