#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=viannam@ufl.edu
#SBATCH --account=mresende
#SBATCH --qos=mresende-b
#SBATCH --job-name=hvfc_2_gvcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256gb
#SBATCH --time=96:00:00
#SBATCH --output=phg_find_paths_%j.log
#SBATCH --error=phg_find_paths_%j.err
##SBATCH --array=2-694  # Uncomment if running mapping in array mode

# ---------------------------
# Environment setup
# ---------------------------
module load conda
module load java

cd /blue/mresende/share/viannam/PHG/
conda activate phgv2-conda
export JAVA_OPTS="-Xmx256G"

# ---------------------------
# Step 1: Export hVCF
# ---------------------------
./phg/bin/phg export-vcf \
    --db-path ./tiledb_maize_phgv2 \
    --dataset-type hvcf \
    --sample-names ./tiledb_maize_phgv2/output/updated_assemblies/Ia453.fa \
    --outputDir ./tiledb_maize_phgv2/hvcf_dataset

# ---------------------------
# Step 2: Build K-mer Index
# ---------------------------
./phg/bin/phg build-kmer-index \
    --db-path ./phgv2_B73_ref/ \
    --hvcf-dir ./phgv2_B73_ref/output/vcf_files

# ---------------------------
# Step 3: Map Short Reads (via keyfile)
# ---------------------------
./phg/bin/phg map-kmers \
    --hvcf-dir ./tiledb_maize_phgv2/hvcf_files \
    --kmer-index ./tiledb_maize_phgv2/hvcf_files/kmerIndex.txt \
    --key-file ./mapping_keyfile.txt \
    --output-dir ./tiledb_maize_phgv2/mapping/

# ---------------------------
# Step 4: Find Paths
# ---------------------------
./phg/bin/phg find-paths \
    --path-keyfile path_keyfile.txt \
    --hvcf-dir ./tiledb_maize_phgv2/hvcf_files \
    --reference-genome ./tiledb_maize_phgv2/output/updated_assemblies/B73v5.fa \
    --path-type haploid \
    --threads 24 \
    --out-parents-dir ./tiledb_maize_phgv2/parents/ \
    --output-dir ./tiledb_maize_phgv2/reads_hvcfs

# ---------------------------
# Step 5: Convert hVCF to gVCF
# ---------------------------
./phg/bin/phg hvcf2gvcf \
    --hvcf-dir ./tiledb_maize_phgv2/reads_hvcfs \
    --db-path ./tiledb_maize_phgv2 \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa \
    --output-dir ./tiledb_maize_phgv2/reads_gvcfs

# ---------------------------
# Cleanup
# ---------------------------
conda deactivate
