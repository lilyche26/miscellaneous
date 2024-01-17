### Step 1: Building the STAR index.*
gtf=/volume2/volume2/GencodeData/gencode.v40.annotation.gtf
genome=/volume2/volume2/GencodeData/GRCh38.primary_assembly.genome.fa
index=/volume2/volume2/STARIndex.HumanGRCh38
STAR --runThreadN 13 --runMode genomeGenerate --sjdbGTFfile  $gtf --genomeFastaFiles $genome --genomeDir $index --sjdbOverhang 100

###alternative index: Ensembl
gtf=/volume2/volume2/EnsemblData/Homo_sapiens.GRCh38.106.chr.gtf
genome=/volume2/volume2/EnsemblData/Homo_sapiens.GRCh38.dna.primary_assembly.fa
index=/volume2/volume2/STARIndex.HumanGRCh38.Ensembl
STAR --runThreadN 13 --runMode genomeGenerate --sjdbGTFfile $gtf --genomeFastaFiles $genome --genomeDir $index --sjdbOverhang 100



###------- Step 2 Align for loop-----############
#!/bin/bash
# define variables
index=/volume2/volume2/STARIndex.HumanGRCh38.Ensembl
# get our data files
FILES=/volume2/volume2/LiTing.RNA-seq/igm-storage2.ucsd.edu/220426_A01535_0115_AHTHJWDSX3/*_L003_R1_001.fastq.gz
out_path=/volume2/volume2/LiTing.RNA-seq/STAR_result

for f in $FILES
do
    echo $f
    #base=$(basename $f .fastq.gz)
    base=$(basename $f _L003_R1_001.fastq.gz)
    pathbase=$(dirname $f)
    read2=$pathbase"/"$base"_L003_R2_001.fastq.gz"
    echo $base
    STAR --runThreadN 13 --genomeDir $index --readFilesIn $f $read2 --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix $out_path/$base"_"
done
echo "done!"



###convert bam file to bigwig file:
###Chip-seq
bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2150570000
    --ignoreForNormalization chrX
    --extendReads

##RNA-seq
bamCoverage -b reads.bam -o coverage.bw




###------- Step 3: STAR results extract-----############
FILES=*_Log.final.out
for f in $FILES
do
    #echo $f
    grep "Uniquely mapped reads % |" $f
done



###------- Step 4: Covert STAR output for loop : *ReadsPerGene.out.tab to Deseq2 count table;
###R
ff <- list.files( path = "./STAR_result", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame(sapply( counts.files, function(x) x[ , 4] ) )
####trascript stand-specifit V4; could see in N_noFeature
ff <- gsub( "[_]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/STAR_result/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
write.table(counts,file="STAR.count_matrix.18samples.WT.KO.61498genes.xls",quote=F,sep="\t",row.names=T,col.names=T)

