## Clear the environment
rm(list=ls())

## Set working directory
setwd("/red/mresende/viannam/PHG/")
Sys.setenv(JAVA_HOME = "/apps/java/jdk-20.0.1")
Sys.setenv(LD_LIBRARY_PATH = "/apps/java/jdk-20.0.1/lib/server:/opt/slurm/lib64:/usr/local/lib")
.libPaths(c("/home/viannam/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))
# Sys.getenv("JAVA_HOME")
# Sys.getenv("LD_LIBRARY_PATH")
## Install rPHG:
#install.packages("pak")
#pak::pak("maize-genetics/rPHG2")
#remove.packages("rJava")
#install.packages("rJava", type = "source")

## Load package:
library(rPHG2)
library(rJava)
options(java.parameters = "-Xmx500g")
#Initiciate database
initPhg("./phg/lib/")

# Create a connection
locCon <- list.files("./phgv2_B73_ref/find_paths_rope_BWT/", pattern = ".h.vcf", full.names = TRUE )|>
  PHGLocalCon()
locCon@hVcfFiles

# Build a graph
graph <- locCon |> buildHaplotypeGraph()

graph

graph |> javaRefObj()

graph |> readSamples()

graph |> readRefRanges()

#Metadata
graph |> readHapIdMetaData()

#Pos
graph |> readHapIdPosMetaData()

##Haplotype IDs
#To return all haplotype IDs as a “sample × reference range” matrix object, we can use the readHapIds() function

m <- graph |> readHapIds()

write.table(m, file = "haplotype_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

p <- graph |> readHapIdMetaData()
write.table(p, file = "haplotype_ID_assembly.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

phgDs <- graph |> readPhgDataSet()
phgDs
phgDs |> numberOfChromosomes()
phgDs |> numberOfRefRanges()
phgDs |> numberOfHaplotypes()

#Number of unique haplotypes IDs
phgDs |> numberOfHaplotypes(byRefRange = TRUE)

g <- phgDs |> numberOfHaplotypes(byRefRange = TRUE)

write.table(g, file = "number_of_haploypes_by_ref_range.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Removing Scaffolds
# Access the seqnames from the GRanges object
seqnames_phgDs <- as.character(seqnames(phgDs@refRanges))

# Filter out scaffolds (assuming scaffolds contain the word "Scaffold")
filtered_seqnames <- seqnames_phgDs[!grepl("^scaff*", seqnames_phgDs)]

filtered_seqnames <- filtered_seqnames[!is.na(filtered_seqnames)]
filtered_seqnames <- as.vector(filtered_seqnames)
filtered_seqnames <- as.character(filtered_seqnames)
ref_seqnames <- as.character(seqnames(phgDs@refRanges))
ref_seqnames <- ref_seqnames[!is.na(ref_seqnames)]
# Subset the refRanges slot to only include non-scaffold entries
filtered_refRanges <- phgDs@refRanges[ref_seqnames %in% filtered_seqnames]

# Replace the refRanges in the original object with the filtered refRanges
phgDs_filtered <- phgDs
phgDs_filtered@refRanges <- filtered_refRanges

# Now apply numberOfChromosomes to the filtered PHGDataSet
phgDs_filtered |> numberOfChromosomes()

####################################

library(data.table)
library(reshape2)
haplotype_dt <- as.data.table(m, keep.rownames = "Sample")
haplotype_dt[, Sample := sub("\\.R1.*", "", Sample)]
haplotype_long <- reshape2::melt(haplotype_dt, id.vars = "Sample", variable.name = "Refrange", value.name = "HaplotypeID")
haplotype_counts <- aggregate(. ~ HaplotypeID + Sample, data = haplotype_long, FUN = length)
haplotype_long <- as.data.table(haplotype_long)
haplotype_long[1:10,1:3]
haplotype_counts <- haplotype_long[, .N, by = .(HaplotypeID, Sample)]
haplotype_counts[1:10,1:3]
haplotype_summary <- dcast(haplotype_counts, Sample ~ HaplotypeID, value.var = "N", fill = 0)
haplotype_summary <- unique(haplotype_long[, .(Sample, HaplotypeID)])
haplotype_summary[, Presence := 1]

pav_matrix <- dcast(haplotype_summary, Sample ~ HaplotypeID, value.var = "Presence", fill = 0)

write.table(haplotype_summary, file = "haplotype_matrix_phgv2_B73_ref_ropeBWT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(pav_matrix, file = "haplotype_PAV_matrix.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#check unique haplotype IDs
colnames(haplotype_summary)[duplicated(colnames(haplotype_summary))]

any(duplicated(colnames(haplotype_summary)))

######
#Getting ref ranges on the sh2 region chr3:21800000-22800000
ref_ranges <- readRefRanges(phgDs)
ref_ranges_df <- as.data.frame(ref_ranges)

# Filter to chr3:21800000–22800000
##I Increase the range to match the imputation file/plots after
sh2_ranges <- ref_ranges_df[
  ref_ranges_df$seqnames == "chr3" &
    ref_ranges_df$start <= 22800000 &
    ref_ranges_df$end >= 21800000, ]

# Optional: get the reference range IDs
sh2_rr_ids <- sh2_ranges$rr_id

# Find which reference ranges overlap the sh2 region
hap_meta <- readHapIdMetaData(phgDs)

sh2_haps <- hap_meta[hap_meta$ref_range_hash %in% sh2_rr_ids, ]
sh2_haps[, c("hap_id", "sample_name")]
unique(sh2_haps$sample_name)
write.table(sh2_haps, file = "sh2_haplotypes_assemblies_info.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##Getting PAV haplotype matrix from the sh2 haplotypes
sh2_hap_ids <- sh2_haps$hap_id
haplo_annotated <- merge(
  haplotype_summary,
  hap_meta[, c("hap_id", "ref_range_hash"), with = FALSE],  # just needed columns
  by.x = "HaplotypeID", by.y = "hap_id",
  all.x = FALSE, all.y = FALSE
)
haplo_sh2 <- haplo_annotated[ref_range_hash %in% sh2_rr_ids]
sh2_pav_matrix <- dcast(haplo_sh2, Sample ~ HaplotypeID, value.var = "Presence", fill = 0)

write.table(sh2_pav_matrix, file = "pav_matrix_sh2_gene_chr3_21800000_22800000.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###########################################################
#Getting ref ranges on the su1 region chr4:43400000-43800000
ref_ranges <- readRefRanges(phgDs)
ref_ranges_df <- as.data.frame(ref_ranges)

# Filter to chr4:43400000–43800000
##I Increase the range to match the imputation file/plots after
su1_ranges <- ref_ranges_df[
  ref_ranges_df$seqnames == "chr4" &
    ref_ranges_df$start <= 44000000 &
    ref_ranges_df$end >= 43000000, ]

# Optional: get the reference range IDs
su1_rr_ids <- su1_ranges$rr_id

# Find which reference ranges overlap the su1 region
hap_meta <- readHapIdMetaData(phgDs)

su1_haps <- hap_meta[hap_meta$ref_range_hash %in% su1_rr_ids, ]
su1_haps[, c("hap_id", "sample_name")]
unique(su1_haps$sample_name)
write.table(su1_haps, file = "su1_haplotypes_assemblies_info_sweetcap.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##Getting PAV haplotype matrix from the su1 haplotypes

haplo_su1 <- haplo_annotated[ref_range_hash %in% su1_rr_ids]

su1_pav_matrix <- dcast(haplo_su1, Sample ~ HaplotypeID, value.var = "Presence", fill = 0)

# Save SU1 PAV matrix
write.table(su1_pav_matrix, file = "pav_matrix_su1_gene_chr4_40000000_45000000.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Getting ref ranges on the su1 region chr4:40000000-45000000
ref_ranges <- readRefRanges(phgDs)
ref_ranges_df <- as.data.frame(ref_ranges)

# Filter to chr4:40000000–45000000 (5Mb)
##I Increase the range to match the imputation file/plots after
su1_ranges <- ref_ranges_df[
  ref_ranges_df$seqnames == "chr4" &
    ref_ranges_df$start <= 45000000 &
    ref_ranges_df$end >= 40000000, ]

# Optional: get the reference range IDs
su1_rr_ids <- su1_ranges$rr_id

# Find which reference ranges overlap the su1 region
hap_meta <- readHapIdMetaData(phgDs)

su1_haps <- hap_meta[hap_meta$ref_range_hash %in% su1_rr_ids, ]
su1_haps[, c("hap_id", "sample_name")]
unique(su1_haps$sample_name)
write.table(su1_haps, file = "su1_haplotypes_assemblies_info_5M_sweetcap.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##Getting PAV haplotype matrix from the su1 haplotypes

su1_hap_ids <- su1_haps$hap_id
haplo_su1 <- haplo_annotated[ref_range_hash %in% su1_rr_ids]
su1_pav_matrix <- dcast(haplo_su1, Sample ~ HaplotypeID, value.var = "Presence", fill = 0)
write.table(su1_pav_matrix, file = "pav_matrix_su1_gene_chr4_40000000_45000000.txt", sep = "\t", row.names = FALSE, quote = FALSE)
