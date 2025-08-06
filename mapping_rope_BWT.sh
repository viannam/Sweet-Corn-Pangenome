#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=viannam@ufl.edu
#SBATCH --account=mresende
#SBATCH --qos=mresende-b
#SBATCH --job-name=find_paths
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=160Gb
#SBATCH --time=96:00:00
#SBATCH --output=find_paths_%j.log
#SBATCH --error=find_paths_%j.err
#SBATCH --array=1-693

# ---------------------------
# Environment setup
# ---------------------------
module load conda
module load java

cd /red/mresende/viannam/PHG/
conda activate phgv2-ropebwt-conda
export JAVA_OPTS="-Xmx256G"

# ---------------------------
# Step 1: Export hVCFs (optional, only if new samples were added)
# ---------------------------
./phg/bin/phg export-vcf \
    --db-path ./phgv2_B73_ref \
    --dataset-type hvcf \
    --sample-names ./assemblies_list \
    --output-dir ./phgv2_B73_ref/hvcf_dataset

# ---------------------------
# Step 2: Build ropeBWT index
# ---------------------------
./phg/bin/phg rope-bwt-index \
    --db-path ./phgv2_B73_ref \
    --hvcf-dir ./phgv2_B73_ref/hvcf_files \
    --output-dir ./phgv2_B73_ref/imputation_rope_BWT \
    --index-file-prefix rope_BWT_index

# ---------------------------
# Step 3: Map reads (array job per keyfile)
# ---------------------------
keyfile_path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./key_files/key_file_list.txt)
echo "Using keyfile: ./key_files/${keyfile_path}"

if [[ ! -f "./key_files/${keyfile_path}" ]]; then
    echo "Error: Keyfile does not exist: ./key_files/${keyfile_path}"
    exit 1
fi

./phg/bin/phg map-reads \
    --hvcf-dir ./phgv2_B73_ref/hvcf_files \
    --index ./phgv2_B73_ref/imputation_rope_BWT/rope_BWT_index.fmd \
    --key-file ./key_files/${keyfile_path} \
    --output-dir ./phgv2_B73_ref/read_mappings_rope_BWT

# ---------------------------
# Step 4: Find paths
# ---------------------------
keyfile_path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./key_files/find_paths_keyfile_list.txt)
echo "Using keyfile: ./key_files/${keyfile_path}"

./phg/bin/phg find-paths \
    --path-keyfile /red/mresende/viannam/PHG/key_files/${keyfile_path} \
    --hvcf-dir /red/mresende/viannam/PHG/phgv2_B73_ref/hvcf_files \
    --reference-genome /blue/mresende/share/viannam/PHG/data/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa \
    --path-type haploid \
    --output-dir /red/mresende/viannam/PHG/phgv2_B73_ref/find_paths_rope_BWT/ \
    --out-parents-dir /red/mresende/viannam/PHG/phgv2_B73_ref/find_paths_rope_BWT/parents
