suppressMessages(suppressWarnings(library("edgeR", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
BASE_CONDITION_SAMPLES <- args[2]
EXPERIMENTAL_CONDITION_SAMPLES <- args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
OUTPUT_FILE_BASE <- 'edgeR_results'
OUTPUT_NORMALIZED_COUNTS_BASE <- 'edgeR_normalized_counts'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
CONDITION_A <- make.names(CONDITION_A)
CONDITION_B <- make.names(CONDITION_B)

# create a string to identify the contrast:
contrast_str = paste0(CONDITION_B, '_vs_', CONDITION_A)

# the sample names are given as a comma-delimited string. Split them
base_samples <- make.names(strsplit(BASE_CONDITION_SAMPLES, ',')[[1]])
exp_samples <- make.names(strsplit(EXPERIMENTAL_CONDITION_SAMPLES, ',')[[1]])
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
