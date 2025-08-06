#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=viannam@ufl.edu
#SBATCH --account=mresende
#SBATCH --qos=mresende-b
#SBATCH --job-name=phg_B73_load
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH --time=96:00:00
#SBATCH --output=phg_B73_load_%j.log
#SBATCH --error=phg_B73_load_%j.err

# ---------------------------
# Environment setup
# ---------------------------
module load conda
module load java

cd /blue/mresende/share/viannam/PHG/
conda activate phgv2-conda
export JAVA_OPTS="-Xmx320G"

# ---------------------------
# Step 1: Prepare Assemblies
# ---------------------------
./phg/bin/phg prepare-assemblies \
    --keyfile key_file.txt \
    --threads 10 \
    --output-dir ./DataSets/output/updated_assemblies

# ---------------------------
# Step 2: Create BED from GFF
# ---------------------------
./phg/bin/phg create-ranges \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa \
    --gff ./data/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
    --boundary gene \
    --pad 500 \
    --range-min-size 500 \
    -o ./DataSets/output/annotation_B73.bed

# ---------------------------
# Step 3: Align Assemblies
# ---------------------------
./phg/bin/phg align-assemblies \
    --gff ./data/B73v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa \
    --assembly-file-list ./assemblies_list.txt \
    -o ./DataSets/output/B73_aligment/

# ---------------------------
# Step 4: Compress Assemblies
# ---------------------------
./phg/bin/phg agc-compress \
    --db-path ./phgv2_B73_ref \
    --fasta-list ./assemblies_list.txt \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa

# ---------------------------
# Step 5: Create Reference VCF
# ---------------------------
./phg/bin/phg create-ref-vcf \
    --bed ./DataSets/output/annotation_B73.bed \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa \
    --reference-name B73v5 \
    --db-path ./phgv2_B73_ref

# ---------------------------
# Step 6: Create MAF-based VCFs
# ---------------------------
./phg/bin/phg create-maf-vcf \
    --db-path ./phgv2_B73_ref \
    --bed ./DataSets/output/annotation_B73.bed \
    --reference-file ./DataSets/output/updated_assemblies/B73v5.fa \
    --maf-dir ./DataSets/output/B73_aligment \
    -o ./output/vcf_files

# ---------------------------
# Step 7: Load VCFs into PHG DB
# ---------------------------
./phg/bin/phg load-vcf \
    --vcf-dir ./output/vcf_files \
    --db-path ./phgv2_B73_ref \
    --threads 10

# ---------------------------
# Clean up
# ---------------------------
conda deactivate
