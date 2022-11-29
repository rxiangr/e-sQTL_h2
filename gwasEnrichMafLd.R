#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("4 argument must be supplied", call.=FALSE)
}

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
library(ggplot2)
options(bitmapType='cairo')

ciseqtlpath <- args[1] #'/group/dairy/variantid/for_Ruidong/cgtexdt/gwas/Blood.topeQTL'
ciseqtlpatt <- args[2] #'p5e-08top0.cis.topsel.txt.gz'
transeqtlpath <- args[3] #'/group/dairy/variantid/for_Ruidong/cgtexdt/gwas/Blood.topeQTL'
transeqtlpatt <- args[4] #'p5e-08top3.trans.topsel.txt.gz'
gwasfn <- args[5] #'/group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gctawtgwas/cbGwasRes/cow/raw/cb.tr01.mlma'
ldmafpath <- args[6] #'/group/dairy/variantid/for_Ruidong/ldscore/ldscout'
ldmafpatt <- args[7] #'score.ld'
chipdtfn <- args[8] #'/group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/funcgenom/chipsqSNP/ChipseqSNP.ovlp.txt.gz'
Npkcut <- as.numeric(args[9])
Nrep <- as.numeric(args[10]) #3
outpath <- args[11] #'/group/dairy/variantid/for_Ruidong/gwasEnrichMafLd/gwasenrich'
outpref <- args[12] #'Blood.tr01'

#---readin chip-seq data
cat(paste0('Reading chIP-seq data at ',Sys.time()),'\n')
chipdt <- fread(chipdtfn)
#---read in eqtl data
cat(paste0('Reading eqtl data at ',Sys.time()),'\n')
ciseqtllist <- list.files(ciseqtlpath,ciseqtlpatt)
#ciseqtldt <- setDT(bind_rows(lapply(paste0(ciseqtlpath,'/',ciseqtllist),fread)))
ciseqtldt <- setDT(do.call(rbind,lapply(paste0(ciseqtlpath,'/',ciseqtllist),fread)))
ciseqtlSNP <- setDT(data.frame(SNP=ciseqtldt[,unique(SNP)]))
ciseqtlSNP <- ciseqtlSNP[SNP %in% chipdt[Npk>=Npkcut]$SNP]
#--trans
transeqtllist <- list.files(transeqtlpath,transeqtlpatt)
#transeqtldt <- setDT(bind_rows(lapply(paste0(transeqtlpath,'/',transeqtllist),fread)))
transeqtldt <- setDT(do.call(rbind,lapply(paste0(transeqtlpath,'/',transeqtllist),fread)))
transeqtlSNP <- setDT(data.frame(SNP=transeqtldt[(!SNP %in% unlist(ciseqtlSNP)),unique(SNP)]))
transeqtlSNP <- transeqtlSNP[SNP %in% chipdt[Npk>=Npkcut]$SNP]

#---read in gwas data
cat(paste0('Reading GWAS data at ',Sys.time()),'\n')
gwasdt <- fread(gwasfn)
nonzero.pmin <- gwasdt[p>0,min(p)]
gwasdt[,p:=ifelse(p==0,nonzero.pmin,p)]
gwasdt[,mlogp:=-log10(p)]
gwasdt[,tval:=abs(b/se)]
#---read in mafld data
cat(paste0('Reading MAF-LD data at ',Sys.time()),'\n')
ldmaflist <- list.files(ldmafpath,ldmafpatt)
ldmafdt <- setDT(bind_rows(lapply(paste0(ldmafpath,'/',ldmaflist),fread)))
setorder(ldmafdt,MAF,ldscore)
#ldmafdt[,MAFbin:=paste0('m',ntile(MAF, 20))]
ldmafdt[,MAFbin:=paste0('chr',chr,'m',ntile(MAF, 20)),by=chr]
ldmafdt[,MAFLDbin:=paste0(MAFbin,'l',as.numeric(ntile(ldscore, 20),by=MAFbin))]
#---get bin data for cis and trans eqtl
#--function to produce test data
cat(paste0('Generating test data at ',Sys.time()),'\n')
dtconvfunc <- function(ldmafdt,gwasdt,taglist1,taglist2,snpcol,bincol,Ncol){
bindt <- ldmafdt[get(snpcol) %in% unlist(taglist1)]
bindt[,N:=length(unique(get(snpcol))),by=bincol]
bin.ct <- bindt[,.N,by=bincol]
ldmafdt.bin <- merge(ldmafdt,bin.ct,by=bincol,sort=F)
ldmafdt.bin1 <- merge(ldmafdt.bin[,grep(paste0('^',snpcol,'$|','^',bincol,'$|','^',Ncol,'$'),colnames(ldmafdt.bin)),with=F],gwasdt[(!SNP %in% unlist(taglist1))&(!SNP %in% unlist(taglist2))&(!SNP %in% chipdt$SNP),c(2,9:11)],by=snpcol,sort=F)
return(ldmafdt.bin1)
}
testdt.cis <- dtconvfunc(ldmafdt,gwasdt,ciseqtlSNP,transeqtlSNP,'SNP','MAFLDbin','N')
testdt.trans <- dtconvfunc(ldmafdt,gwasdt,transeqtlSNP,ciseqtlSNP,'SNP','MAFLDbin','N')

#----function to sample random variants
cat(paste0('Sampling random variants with matching bins at ',Sys.time()),'\n')
set.seed(123)
samfunc <- function(ldmafdt.tagbindt,pcol) {
randsmpdt <- ldmafdt.tagbindt[,list(SNP=sample(SNP,N,replace=T)),by=MAFLDbin]
ldmafdt.tagbindt_rowsub <- ldmafdt.tagbindt[SNP %in% randsmpdt$SNP]
set2dt <- ldmafdt.tagbindt_rowsub[,list(m.est=mean(get(pcol),na.rm=T),se=sd(get(pcol),na.rm=T)/sqrt(.N),nSNP=.N)]
return(set2dt)
}
#Nrep <- 3
#---cis random
smplist.cis <- list()
cat(paste0('Sampling and analysis of cis started at ',Sys.time()),'\n')
#smplist.cis <- replicate(3,data.frame(samfunc(testdt.cis,'mlogp')))
smplist.cis <- replicate(Nrep,data.frame(samfunc(testdt.cis,'mlogp')))
cisresdim <- dim(smplist.cis)[-1]
cbres.cis <- data.frame(t(matrix(unlist(smplist.cis),cisresdim)))
colnames(cbres.cis) <- c('m.est','se','nSNP')
cat(paste0('Sampling and analysis of cis finished at ',Sys.time()),'\n')
#---trans random
smplist.trans <- list()
cat(paste0('Sampling and analysis of trans started at ',Sys.time()),'\n')
#smplist.trans <- replicate(3,data.frame(samfunc(testdt.trans,'mlogp')))
smplist.trans <- replicate(Nrep,data.frame(samfunc(testdt.trans,'mlogp')))
transresdim <- dim(smplist.trans)[-1]
cbres.trans <- data.frame(t(matrix(unlist(smplist.trans),transresdim)))
colnames(cbres.trans) <- c('m.est','se','nSNP')
cat(paste0('Sampling and analysis of trans finished at ',Sys.time()),'\n')

#---get a random set of cis and trans
randsmpdt.cis <- unique(testdt.cis[,list(SNP=sample(SNP,N,replace=T)),by=MAFLDbin])
randsmpdt.trans <- unique(testdt.trans[!SNP %in% unlist(randsmpdt.cis$SNP),list(SNP=sample(SNP,N,replace=T)),by=MAFLDbin])

#---mean effects of cis
cat(paste0('Analysis targeted eQTL data started at ',Sys.time()),'\n')
cat(paste0(' at ',Sys.time()),'\n')
setfunc <- function(taglist,gwasdt,snpcol,pcol) {
set1dt <- gwasdt[get(snpcol) %in% unlist(taglist),list(m.est=mean(get(pcol),na.rm=T),se=sd(get(pcol),na.rm=T)/sqrt(.N),nSNP=.N)]
return(set1dt)
}
cistagdt <- setfunc(ciseqtlSNP,gwasdt,'SNP','mlogp')
#---mean effects of trans
transtagdt <- setfunc(transeqtlSNP,gwasdt,'SNP','mlogp')

#---save analysis results
cat(paste0('Generating analysis results at ',Sys.time()),'\n')
resdt <- setDT(rbind(data.frame(type='cis',cistagdt),data.frame(type='cis.rand',cbres.cis),data.frame(type='trans',transtagdt),data.frame(type='trans.rand',cbres.trans)))
write.table(data.frame(fn=outpref,resdt),paste0(outpath,'/',outpref,'.mest.txt'),row.names=F,quote=F,sep=' ')
cat(paste0('individual m.est results saved to ',outpath,'/',outpref,'.mest.txt'),sep='\n')

#----generate overall results
tag.overall <- resdt[!grepl('rand',type),list(mean(m.est),sqrt(sum(se^2))),by=type]
#rand.overall <- resdt[grepl('rand',type),list(mean(m.est),sqrt(sum(se^2))),by=type]
rand.overall <- resdt[grepl('rand',type),list(mean(m.est),sd(m.est)/sqrt(.N)),by=type]
tdt <- abs(tag.overall[,2]-rand.overall[,2])/sqrt(tag.overall[,3]^2+rand.overall[,3]^2)
colnames(tdt) <- 'enrich.t'
tval.overall <- cbind(type=paste0(unlist(tag.overall[,1]),'VSrand'),tdt)
tval.overall[,enrich.p:=2*pnorm(abs(enrich.t), lower.tail = F)]
write.table(data.frame(fn=outpref,tval.overall),paste0(outpath,'/',outpref,'.overallres.txt'),row.names=F,quote=F,sep=' ')
cat(paste0('Overall enrichment results saved to ',outpath,'/',outpref,'.overallres.txt'),sep='\n')

#----output a class file
ciseqtlSNP[, c("chr", "bp") := tstrsplit(SNP, ":", fixed=TRUE)]
ciseqtlSNP[,chr:=gsub('Chr','',chr)]
ciseqtlSNP[,class:=1]
ciseqtlSNP[,type:='cis']
transeqtlSNP[, c("chr", "bp") := tstrsplit(SNP, ":", fixed=TRUE)]
transeqtlSNP[,chr:=gsub('Chr','',chr)]
transeqtlSNP[,class:=2]
transeqtlSNP[,type:='trans']
randsmpdt.cis[, c("chr", "bp") := tstrsplit(SNP, ":", fixed=TRUE)]
randsmpdt.cis[,chr:=gsub('Chr','',chr)]
randsmpdt.cis[,class:=3]
randsmpdt.cis[,type:='cis.rand']
randsmpdt.trans[, c("chr", "bp") := tstrsplit(SNP, ":", fixed=TRUE)]
randsmpdt.trans[,chr:=gsub('Chr','',chr)]
randsmpdt.trans[,class:=3]
randsmpdt.trans[,type:='trans.rand']
classdt <- rbind(ciseqtlSNP,transeqtlSNP,randsmpdt.cis[,-1],randsmpdt.trans[,-1])
classdt1 <- merge(classdt,gwasdt[,c(2,9:ncol(gwasdt)),with=F],by='SNP',sort=F)
colnames(classdt1)[ncol(classdt1)-2] <- 'GWA.p'
fwrite(classdt1,paste0(outpath,'/',outpref,'.classdt.gz'),row.names=F,quote=F,compress='gzip',sep=' ',na=NA)
cat(paste0('SNP class data saved to ',outpath,'/',outpref,'.classdt.gz'),sep='\n')


