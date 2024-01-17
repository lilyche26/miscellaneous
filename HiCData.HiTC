library(HiTC)
library(Matrix)


#####----------1.generate HiTC object---#####
##day00.HERV1.KO.HiC.Rep1
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3263nnn/GSM3263087/suppl/GSM3263087_RH_350.hic
# ##day00.HERV2.KO.HiC.Rep1
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3263nnn/GSM3263089/suppl/GSM3263089_RH_352.hic
# ##day00.ctrl.HiC.Rep1
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3263nnn/GSM3263085/suppl/GSM3263085_RH_348.hic


system('/volume2/volume2/tools2/straw/C++/straw KR /GSM3263085_RH_348.hic chr4 chr4 BP 10000 | /HiCProcessing/col2mat.awk -v chr=chr4 -v bin_size=10000 -v chr_size=191154276  > GSM3263085_RH_348.norm.asc')
system('/volume2/volume2/tools2/straw/C++/straw KR /GSM3263087_RH_350.hic chr4 chr4 BP 10000 | /HiCProcessing/col2mat.awk -v chr=chr4 -v bin_size=10000 -v chr_size=191154276  > GSM3263087_RH_350.chr4.norm.asc')
system('/volume2/volume2/tools2/straw/C++/straw KR /GSM3263089_RH_352.hic chr4 chr4 BP 10000 | /HiCProcessing/col2mat.awk -v chr=chr4 -v bin_size=10000 -v chr_size=191154276  > GGSM3263089_RH_352.chr4.norm.asc')


mat_ctrl <- read.table('GSM3263085_RH_348.chr4.norm.asc',sep = '\t', header = T, row.names = 1)
mat_ctrl <- as.matrix(mat_ctrl)
head(mat_ctrl)
grng <- lapply(rownames(mat_ctrl), function(x){
    splt <- strsplit(x, '\\|')
    return(unlist(splt)[3])
})
grng <- gsub('chrchr', 'chr', grng)
grng <- StringToGRanges(grng, sep = c(":", "-"))
gnames <- lapply(rownames(mat_ctrl), function(x){
    s <- gsub('\\|chrchr.+', '', x)
    return(s)
})
names(grng) <- as.character(gnames)
rownames(mat_ctrl) <- gnames
colnames(mat_ctrl) <- gnames
hic_ctrl <- HTCexp(mat_ctrl, xgi=grng, ygi=grng)
detail(hic_ctrl)
head(summary(hic_ctrl))




#####----------2.QC---#####
pdf(file="QC.Control.pdf")
par(mfrow=c(2,2))
CQC(hic_ctrl, winsize = 1e+06, dev.new=FALSE, hist.dist=FALSE)
dev.off()



#####---------3.Data Tranformation---#####
# hic_ctrl_hervh1 <- extractRegion(hic_ctrl, c(1,2), chr="chr4", from=89000000, to=97000000)
hic_ctrl_hervh1 <- extractRegion(hic_ctrl, c(1,2), chr="chr4", from=91000000, to=95000000)
# hic_ctrl_hervh1 <- extractRegion(hic_ctrl, c(1,2), chr="chr4", from=54500000, to=58500000)
hic_ctrl_hervh1.binned <-  binningC(hic_ctrl_hervh1, binsize=10000, method="median", step = 1)
# mapC(hic_ctrl_hervh1.binned)
hic_ctrl_hervh1.exp <- getExpectedCounts(hic_ctrl_hervh1.binned, method="loess", stdev=TRUE, plot=TRUE)
hic_ctrl_hervh1.norm.binned <- normPerExpected(hic_ctrl_hervh1.binned, method="mean")
max(hic_ctrl_hervh1.norm.binned@intdata, na.rm = T)
min(hic_ctrl_hervh1.norm.binned@intdata, na.rm = T)
# mapC(hic_ctrl_hervh1.norm.binned)


pdf(file="Control.Norm.91-95m.pdf")
mapC(hic_ctrl_hervh1.norm.binned)
dev.off()


#####---------4.Comparing Experiments---##### 
mat_ko_D00_HERV1_Rep1 <- read.table('GSM3263087_RH_350.chr4.norm.asc', 
                             sep = '\t', header = T, row.names = 1)
mat_ko_D00_HERV1_Rep1 <- as.matrix(mat_ko_D00_HERV1_Rep1)

grng <- lapply(rownames(mat_ko_D00_HERV1_Rep1), function(x){
    splt <- strsplit(x, '\\|')
    return(unlist(splt)[3])
})
grng <- gsub('chrchr', 'chr', grng)
grng <- StringToGRanges(grng, sep = c(":", "-"))
gnames <- lapply(rownames(mat_ko_D00_HERV1_Rep1), function(x){
    s <- gsub('\\|chrchr.+', '', x)
    return(s)
})

names(grng) <- gnames
rownames(mat_ko_D00_HERV1_Rep1) <- gnames
colnames(mat_ko_D00_HERV1_Rep1) <- gnames



hic_ko_D00_HERV1_Rep1 <- HTCexp(mat_ko_D00_HERV1_Rep1, xgi=grng, ygi=grng)
hic_ko_hervh1 <- extractRegion(hic_ko_D00_HERV1_Rep1, c(1,2), chr="chr4", from=91000000, to=95000000)
hic_ko_hervh1.binned <-  binningC(hic_ko_hervh1, binsize=10000, method="median", step = 3)
hic_ko_hervh1.norm.binned <- normPerExpected(hic_ko_hervh1.binned, method="mean")


pdf(file="KO.D00.HERV1.Rep1.Norm.91-95m.pdf")
mapC(hic_ko_hervh1.norm.binned)
dev.off()


##compare:
pdf(file="Ctrl.vs.KO.D00.HERV1.Rep1.NotNorm.91-95m.pdf")
mapC(hic_ctrl_hervh1.binned, hic_ko_hervh1.binned)
dev.off()
pdf(file="Ctrl.vs.KO.D00.HERV1.Rep1.Norm.91-95m.pdf")
mapC(hic_ctrl_hervh1.norm.binned, hic_ko_hervh1.norm.binned)
dev.off()

pdf(file="Ctrl.D00.Rep1.Norm.91-95m.pdf")
mapC(hic_ctrl_hervh1.norm.binned)
dev.off()
pdf(file="KO.D00.HERV1.Rep1.Norm.91-95m.pdf")
mapC(hic_ko_hervh1.norm.binned)
dev.off()


#####---------5. Call TADs---#####
# Detection of TADs
# Directionality Index
# Using HiTC
##calculate DI:
di_ctrl <- directionalityIndex(hic_ctrl_hervh1)
##plot DI:
pdf(file="Ctrl.D00.Rep1.DI.TAD.pdf")
barplot(di_ctrl, col=ifelse(di_ctrl>0,"red","blue"), main = 'WT',border=NA)
barplot(di_ctrl, col=ifelse(di_ctrl>0,"darkred","darkgreen"), main = 'WT')
mapC(hic_ctrl_hervh1.norm.binned)
dev.off()

# di_ctrl <- directionalityIndex(hic_ctrl_hervh1.binned, winup = 1e+09, windown = 1e+09)
# barplot(di_ctrl, col=ifelse(di_ctrl>0,"darkred","darkgreen"), main = 'WT')

di_ko <- directionalityIndex(hic_ko_hervh1)
pdf(file="KO.D00.HERV1.Rep1.DI.TAD.pdf")
barplot(di_ko, col=ifelse(di_ko>0,"red","blue"), main = 'KO',border=NA)
barplot(di_ko, col=ifelse(di_ko>0,"darkred","darkgreen"), main = 'KO')
mapC(hic_ko_hervh1.norm.binned)
dev.off()

sum(is.na(test))

