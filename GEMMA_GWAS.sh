#!/bin/bash
#SBATCH --job-name=GWAS_haplo
#SBATCH --mail-user=viannam@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --mem=120GB
#SBATCH --qos=mresende-b
#SBATCH --account=mresende
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=./GWAS_haplo_%a_%A.out
#SBATCH --error=./GWAS_haplo_%a_%A.err
#SBATCH --array=1-20

## Let's run GWAS for sweet corn traits using GEMMA based on haplotypes from the PHG (maize_phg B73 ref)
##the haplotype matrix was generated through the rPHG (rPHGv2-copy.R)
##PHENO FILE ONLY THE BLUES COLUMN!! Replace NA by -9
##Use the Gmatrix regenerated from X*X`\ncol X
#module load R
#Rscript G_matrix_haplo.R

module load gemma

trait=$(cat pheno_list.txt|head -n ${SLURM_ARRAY_TASK_ID}| tail -n 1)

#gemma -g haplo_table_bimbam_format.txt -p $trait -gk 1 -o kinship_matrix_GEMMA_haplo.txt

#gemma -g haplo_table_bimbam_format.txt -p $trait -k ./output/kinship_matrix_GEMMA_haplo.txt.cXX.txt -lmm 2 -o $trait

#gemma -g haplo_table_bimbam_format.txt -p $trait -k ./haplo_kinship_matrix.txt -lmm 2 -o $trait

gemma -g haplo_table_bimbam_format_filetered.txt -p $trait -k ./haplo_kinship_matrix_filtered.txt -lmm 2 -o $trait

#module load R
#Rscript haplo_GWAS_plot.R
