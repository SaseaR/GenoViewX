#' Mplot for three type of QTLs with GWAS association
#' MarkerName : eg. 1:10086:A/C or 1:10086:A:C
#' Usage eg. Mplot_for_multiple_QTL_GWAS_intersect('/A/B/C','a.TBL','QTL1.gz','QTL2.gz,'QTL3.gz,'/a/b/c','plotA',5e-04,10,1,5.5,20)
#' @param file_path the pathway of meta file
#' @param file_name the name of of meta file with suffix
#' @param qtl_file1 the file of type 1 QTL
#' @param qtl_file2 the file of type 2 QTL
#' @param qtl_file3 the file of type 3 QTL
#' @param figure_pathy the pathway to save Mplot
#' @param figure_name the name of of Mplot
#' @param p_col column number of p-value
#' @param snp_col colnum number of MarkerName
#' @param p_threshold P-value threshold of GWAS association
#' @param height height of Mplot
#' @param width width of Mplot
#' @export

Mplot_for_multiple_QTL_GWAS_intersect <- function(file_path,file_name,qtl_file1,qtl_file2,qtl_file3,figure_pathy,figure_name,p_col,snp_col,p_threshold,height,width){

  requireNamespace(library(CMplot))
  requireNamespace(library(data.table))
  dat1 <- data.table::fread(paste0(file_path,'/',file_name))
  names(dat1)[p_col] <- 'pvalue'
  names(dat1)[snp_col] <- 'MarkerName'
  dat1$CHR <- sapply(strsplit(dat1$MarkerName,':'),'[',1)
  dat1$BP <- sapply(strsplit(dat1$MarkerName,':'),'[',2)
  dat_plot1 <- dat1[,c('MarkerName','CHR','BP','pvalue')]
  names(dat_plot1) <- c('SNP','Chromosome','Position',figure_name)
  dat_plot1 <- dat_plot1[which(dat_plot1[,4]<p_threshold),]
  dat_plot1 <- subset(dat_plot1,Chromosome%in%c(1:22))

  qtl1 <- data.table::fread(qtl_file1) # eqtl
  qtl2 <- data.table::fread(qtl_file2) # aqtl
  qtl3 <- data.table::fread(qtl_file3) # sqtl

  group_snps1 <-  dat_plot1$SNP[dat_plot1$SNP%in%qtl1$variant_id] # e
  group_snps2 <-  dat_plot1$SNP[dat_plot1$SNP%in%qtl2$variant_id] # a
  group_snps3 <-  dat_plot1$SNP[dat_plot1$SNP%in%qtl3$variant_id] # s
  group_snps4 <-  dat_plot1$SNP[dat_plot1$SNP%in%intersect(qtl1$variant_id,qtl2$variant_id)] # ae
  group_snps5 <-  dat_plot1$SNP[dat_plot1$SNP%in%intersect(qtl1$variant_id,qtl3$variant_id)] # es
  group_snps6 <-  dat_plot1$SNP[dat_plot1$SNP%in%intersect(qtl3$variant_id,qtl2$variant_id)] # as
  group_snps7 <-  dat_plot1$SNP[dat_plot1$SNP%in%intersect(intersect(qtl3$variant_id,qtl2$variant_id),qtl1$variant_id)] # aes
  highlight_snps <- c(group_snps1,group_snps2,group_snps3)

  highlight_col <- c(rep("#0072B580", length(group_snps1)),   # e
                     rep("#E1872780", length(group_snps3)),   # s
                     rep("#F46E6080", length(group_snps2))   # a
  )

  highlight_pch <- c(rep(15, length(group_snps1)),   # e
                     rep(17, length(group_snps2)),   # a
                     rep(19, length(group_snps3))   # s
  )

  cols=c("#A3A3A3","#C4C4C4","#E2E2E2")
  setwd(figure_pathy)

  CMplot::CMplot(dat_plot1,plot.type="m",LOG10=TRUE,threshold.lwd=c(1,1), threshold.col=c('black','#BC3C29'),main.font=1,ylab=expression(-log[10](italic(P))),
                 threshold=c(5e-04,5e-08),threshold.lty=c(1,2),chr.den.col=NULL,col = cols,cex =c(0.5,0.5,0.5),lab.font=1,main="Chromosome",height = height,width = width,
                 amplify= TRUE, signal.cex=0.5, file="pdf",file.output=TRUE,verbose=TRUE,legend.pos='middle',  highlight = highlight_snps,
                 highlight.col = highlight_col,highlight.pch = highlight_pch)
}
