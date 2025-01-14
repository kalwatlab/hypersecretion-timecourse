# edgeR time course analysis for thapsigargin-treated samples

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))

cnt_type = "DMSO"
trt_type = "Tg"

x <- read.csv(paste0("data/sample_matrix_rawcount_",trt_type,".csv"), header=T, stringsAsFactors=FALSE,row.names=1)
newcols = unlist(lapply(colnames(x),function(x) if(substr(x,1,8)=="Baseline"){return(c(sub("Baseline",cnt_type,x),sub("Baseline",trt_type,x)))}else{return(x)}))
x = x[,unlist(lapply(colnames(x),function(x) if(substr(x,1,8)=="Baseline"){return(c(x,x))}else{return(x)}))]
colnames(x) = newcols
x = as.data.frame(x)

## Setting grouping variables
time<-factor(sub("[_].*","",sub(".*[.]","",colnames(x))),levels=c("0h","1h","2h","6h","24h"))
treat<-factor(sub("[.].*","",colnames(x)),levels=c(cnt_type,trt_type))
levels(time)
levels(treat)

## Creating DGEList object
y<-DGEList(counts=x, genes = row.names(x), group=paste0(as.character(treat),"_",as.character(time)))
y <- calcNormFactors(y)
y

## Filter out lowly expressed genes
keep <- filterByExpr(y)
table(keep)
y <- y[keep,,keep.lib.sizes=FALSE]

## TMM Normalization
y <- calcNormFactors(y)
y$samples

## Extract TMM values for later use
tmm_cpm <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE)

tmm_cpm_df <- as.data.frame(tmm_cpm)
tmm_cpm_df$gene <- rownames(tmm_cpm_df)

# Create a new column for the sample group (e.g., Sample1, Sample2)
tmm_long <- tmm_cpm_df %>%
  pivot_longer(cols = -gene, names_to = "sample_rep", values_to = "tmm_cpm") %>%
  
  # Extract condition (DMSO, SW, Tg) and time point (0h, 1h, etc.)
  mutate(
    condition = gsub("\\..*", "", sample_rep),      # Extract the condition
    time = gsub(".*\\.(\\d+h)_.*", "\\1", sample_rep)  # Extract the time point
  )

# Summarise to get the average for each condition-time point group
tmm_avg_dplyr <- tmm_long %>%
  group_by(gene, condition, time) %>%
  summarise(tmm_avg = mean(tmm_cpm)) %>%
  ungroup()

# Reshape back to wide format for easy viewing
tmm_avg_wide <- tmm_avg_dplyr %>%
  pivot_wider(names_from = c(condition, time), values_from = tmm_avg)

# View the averaged TMM-normalized counts
head(tmm_avg_wide)

# Export the avg log2(TMMs) table
write.table(tmm_avg_wide, file = "output/dmso_tg_log2_avg_tmms.txt", sep = "\t", row.names = F)

## Design for time course analysis
design <- model.matrix(~time+treat+time*treat, data = y$samples) 
design

## Estimate dispersion:
z <- estimateDisp(y,design)
sqrt(z$common.dispersion)
#plotBCV(z)
#plotMDS(z, top = 500, method = "logFC", labels = time, col=ifelse(sub("[_].*","",y$samples$group)==trt_type,"red","blue"))

#To perform quasi-likelihood F-tests:
fit <- glmQLFit(z,design)
#plotQLDisp(fit)

## Time course trend analysis
## Using the design matrix we are choosing genes that are both significantly
## altered over time and based on treatment. i.e. the interaction between drug treatment and time is significant.
test_coefficients = colnames(design)[7:10] # For future analyses this range may need to be adjusted
test_coefficients

qlf <- glmQLFTest(fit,coef=test_coefficients)
tab1 <- as.data.frame(topTags(qlf, n=Inf)) # this table contains all genes regardless of signficance or FC
# use a  lfc = 1 and BH-FDR < 0.05 for downstream DPGP clustering. use FC = 0 and p.value=Inf to export entire gene list
keep_lfc1 = decideTests(qlf,lfc=1,p.value=0.05) # lfc is log2fc, so lfc=1 means a cutoff of a linear fold-change > 2
summary(keep_lfc1)

hist(-log10(tab1$FDR), breaks = 50)
head(-log10(tab1$FDR))

tab2 = tab1[row.names(keep_lfc1)[keep_lfc1==1],]
write.csv(tab2, paste0("output/QLFTest_lfc1_pval0.05_",trt_type,".csv"), row.names = FALSE)

# Need to write a table specifically for use as DPGP input
tab3 = tab2[,2:5]
colnames(tab3) = c("1","2","6","24")
write.table(tab3,file=paste0("output/final_data_matrix_lfc1_pval0.05_",trt_type,".txt"),sep='\t',row.names=FALSE,quote=FALSE) # this is for DPGP

#to get ALL the results, write tab1
write.table(tab1,file=paste0("output/final_data_matrix_nocutoffs_",trt_type,".txt"),sep='\t',row.names=FALSE,quote=FALSE)
