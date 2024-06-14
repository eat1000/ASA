library(karyoploteR)

studylanc <- "/data/lai/ASA."
outdir <- "/data/pdf/ASA."
id_list <- read.delim("/data/admix_id.lst", header = FALSE)

id_num <- dim(id_list)[1]
for (i in c(1:id_num))
{
tmp <- strsplit(readLines(paste(studylanc,id_list[i,1],".lai",sep=""), n=1), " |\t|=")
pop <- tmp[[1]]
npop <- (length(pop)-5)/2
pop_list <- pop[rep(1:npop)*2+4]

DM <- read.delim(paste(studylanc,id_list[i,1],".lai",sep=""), col.names=c("chr", "spos", "epos", "sgpos", "egpos", "nsnps", "hap1", "hap2"), sep = "\t", comment.char = "#", header = FALSE)
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 50
pp$data2outmargin <- 50   
nsub <- (dim(DM)[2] -6) /2 

hap_switches <- toGRanges(data.frame(chr=paste("chr",DM[,1],sep=""), start=DM$spos, end=DM$epos))
colors <- c('#FFFFFF', '#0072B2', '#E69F00', '#F0E442', '#009E73', '#56B4E9', '#D55E00', '#CC79A7', '#000000')

pdf(file=paste(outdir,id_list[i,1],".pdf",sep=""),7,7)
kp <- plotKaryotype(genome="hg38", chromosomes=c("autosomal"),plot.type="2", plot.params=pp)
kpRect(kp, hap_switches, y0=0, y1=1, col=colors[DM[, 7] + 1], data.panel = 1, border=NA, clipping=TRUE) 
kpRect(kp, hap_switches, y0=0, y1=1, col=colors[DM[, 8] + 1], data.panel = 2, border=NA, clipping=TRUE)
legend(x = "bottomright", fill = colors[1:npop + 1], legend = pop_list, bty = "n")
dev.off()
}



