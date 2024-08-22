suppressMessages(suppressWarnings(library("edgeR", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
ANN_FILE<-args[2]
ORIG_ANN_COL<-args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
OUTPUT_FILE_BASE <- 'edgeR_results'
OUTPUT_NORMALIZED_COUNTS_BASE <- 'edgeR_normalized_counts'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
ANN_COL <- make.names(ORIG_ANN_COL)

# create a string to identify the contrast:
contrast_str = paste0(ANN_COL, '_', CONDITION_B, '_vs_', CONDITION_A)

# read the annotation file:
annotations <- read.table(ANN_FILE, sep='\t', header = T, row.names = 1)

if (!(ANN_COL %in% colnames(annotations))) {
    message(sprintf('The column "%s" was not found in your annotation file. Please check your inputs. Note that this input is case-sensitive and must match exactly.', ORIG_ANN_COL))
    quit(status=1)
}

# filter to only those samples which match the conditions
base_samples_filter <- annotations[, ANN_COL] == CONDITION_A
exp_samples_filter <- annotations[, ANN_COL] == CONDITION_B
orig_base_samples <- rownames(annotations)[base_samples_filter]
orig_exp_samples <- rownames(annotations)[exp_samples_filter]
base_samples <- make.names(orig_base_samples)
exp_samples <- make.names(orig_exp_samples)

if(
    (length(base_samples) == 0)
    ||
    (length(exp_samples) == 0)
){
    message(sprintf('One or both of your sample sets contained zero samples. Check that both %s and %s are valid values in column "%s".', CONDITION_A, CONDITION_B, ORIG_ANN_COL))
    quit(status=1)
}

if(
    (length(base_samples) < 2)
    ||
    (length(exp_samples) < 2)
){
    message(sprintf('One or both of your sample sets contained fewer than 2 samples. To perform differential expression analysis, replicates are required. Check that both conditions %s and %s have 2 or more samples each.', CONDITION_A, CONDITION_B))
    quit(status=1)
}

intersection_list = intersect(base_samples, exp_samples)

if (length(intersection_list) > 0){
    sample_list = paste0(intersection_list, collapse=',')
    message(paste(
       'The following samples were in both contrast groups. Fix this and try again: ',
       sample_list
    ))
    quit(status=1)
}
all_samples <- c(base_samples, exp_samples)
original_sample_names <- c(orig_base_samples, orig_exp_samples)

# this allows us to map back to the original names after
colname_mapping = data.frame(
    orig_names = original_sample_names,
    row.names=all_samples,
    stringsAsFactors=F)

condition_a_list <- rep(CONDITION_A, length(base_samples))
condition_b_list <- rep(CONDITION_B, length(exp_samples))
all_conditions <- c(condition_a_list, condition_b_list)

# full annotation matrix:
annotations <- as.data.frame(cbind(all_samples, all_conditions), stringsAsFactors = F)
colnames(annotations) <- c('sample', 'condition')

# read the raw count matrix, genes as row names:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F)

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols = colnames(count_data)
annotations <- annotations[annotations$sample %in% count_mtx_cols,]

# Check that after subsetting samples, we actually have more than one 'condition' represented.
# This also covers the edge case where only a single sample name was valid. In that case, the 
# column subset below produces an array and `dim(count_data)[2]` throws an error.
fl <- length(levels(as.factor(annotations$condition)))
if(fl < 2){
    message(sprintf('After subsetting the matrix for the samples of interest (%d found), only one cohort of samples was present. Please double-check your inputs or sample names.',  dim(annotations)[1]))
    quit(status=1)
}

# subset to only keep samples corresponding to the current groups in the count_data dataframe
# Also sorts the columns to match the ordering of the annotation dataframe
count_data <- count_data[,annotations$sample]
if (dim(count_data)[2] == 0){
    message('After subsetting the matrix for the samples of interest, the matrix was empty. Please check the input samples and matrix')
    quit(status=1)
}

# create a factor for the conditions:
condition_factor = factor(annotations$condition, levels=c(CONDITION_A,CONDITION_B))

d <- DGEList(count_data, lib.size = as.vector(colSums(count_data, na.rm=T)))
d <- calcNormFactors(d,method='TMM')

design = model.matrix(~condition_factor)
rownames(design) <- colnames(count_data)
d <- estimateDisp(d, design)
fit <- glmFit(d, design)
lrt <- glmLRT(fit)
resOrdered <- topTags(lrt, n=Inf, sort.by='PValue')

# extract and output the normalized counts:
norm_mtx = cpm(d, normalized.lib.size=TRUE)
nc_cols = colnames(norm_mtx)
remapped_cols = colname_mapping[nc_cols, 'orig_names']
colnames(norm_mtx) = remapped_cols
fout2 <- paste(OUTPUT_NORMALIZED_COUNTS_BASE, contrast_str, 'tsv', sep='.')
fout2 <- paste(working_dir, fout2, sep='/')
write.table(norm_mtx, fout2, sep='\t', quote=F, row.names=T)

# merge to create a single table, which makes frontend work easier
m <- merge(resOrdered, norm_mtx, by="row.names")
rownames(m) <- m[,'Row.names']
drops <- c("Row.names", "Gene", "gene")
m <- m[, !(names(m) %in% drops)]

# change column name for the 'stat' column which will match other dge-type analyses
cols <- colnames(m)
cols[which(names(m) == 'LR')] = 'statistic'
cols[which(names(m) == 'logFC')] = 'log2FoldChange'
cols[which(names(m) == 'PValue')] = 'pvalue'
cols[which(names(m) == 'FDR')] = 'padj'
colnames(m) <- cols

output_filename <- paste(OUTPUT_FILE_BASE, contrast_str, 'tsv', sep='.')
output_filename <- paste(working_dir, output_filename, sep='/')
write.table(m, output_filename, sep='\t', quote=F)

json_str = paste0(
       '{"dge_results":"', output_filename, '",',
       '"normalized_counts":"', fout2, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
