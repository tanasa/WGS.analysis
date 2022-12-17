################################################################################################
################################################################################################
################################################################################################
################################################################################################

library("karyoploteR")
library("plyr")
library("doBy")

################################################ some examples are displayed at the link below :
#################################### https://bernatgel.github.io/karyoploter_tutorial/#Examples
################################################################################################
################################################################################################
################################################################################################
################################################################################################

args <- commandArgs(TRUE)
FILE <- args[1]  

################################################################################################
####################################################### FILE <- "SPCG.all.filtered.SV.AF0.1.AD8"
################################################################################################
       
name <- basename(FILE)

x <- read.delim(FILE, header=T, sep="\t")

head(x)
dim(x)

length(x$SV_ID)
length(unique(x$SV_ID))

### in order to ISOLATE only the UNIQUE SV_ID, as the INVERSIONS occur 2 times :
### http://r.789695.n4.nabble.com/Extracing-only-Unique-Rows-based-on-only-1-Column-td1015815.html

library("doBy")
y <- x[firstobs(x$SV_ID),]
dim(y)

library("plyr")
z <- ddply(x, .(SV_ID), head, n = 1)
dim(z)

###############################################################################################
###############################################################################################

# colnames(z)
# [1] "Sample....names.list.of.files...index.list.of.files.."
# [2] "SV_ID"                                                
# [3] "seqnames.x"                                           
# [4] "start.x"                                              
# [5] "CHR2.x"                                               
# [6] "END.x"                                                
# [7] "SVTYPE.x"                                             
# [8] "CT.x"    
# [72] "AD"                                                   
# [73] "AF"                                                   
# [74] "Sample"
### to print for verification purposes : 

write.table(z, file = paste(name, ".verify.unique.ID.txt",sep=""), 
                 sep ="\t",
                 row.names = FALSE,
                 col.names = TRUE)   

DEL <- subset(z, SVTYPE.x=="DEL")
DUP <- subset(z, SVTYPE.x=="DUP")
INS <- subset(z, SVTYPE.x=="INS")
INV <- subset(z, SVTYPE.x=="INV")
TRA <- subset(z, SVTYPE.x=="TRA")

dim(DEL)
dim(DUP)
dim(INS)
dim(INV)
dim(TRA)

### the number of records are :

dim(DEL)[1]
dim(DUP)[1]
dim(INS)[1]
dim(INV)[1]
dim(TRA)[1]

################################################################################################
################################################################################################
################################################################################################
################################################################################################
# [3] "seqnames.x"                                           
# [4] "start.x"                                              
# [5] "CHR2.x"                                               
# [6] "END.x"  
# [74] "Sample"

### to extract the coordinates for all SV :

coord.SV.partner1 <- makeGRangesFromDataFrame(data.frame(chr=z$seqnames.x, 
                                                         start=z$start.x,
                                                         end=z$start.x ))

coord.SV.partner2 <- makeGRangesFromDataFrame(data.frame(chr=z$CHR2.x, 
                                                         start=z$END.x,
                                                         end=z$END.x ))

#### in order to merge these GRANGES objects :

coord.SV.partners <- c(coord.SV.partner1, coord.SV.partner2)
length(coord.SV.partners)

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(DEL)[1] > 0)
{
coord.DEL.partner1 <- makeGRangesFromDataFrame(data.frame(chr=DEL$seqnames.x, 
                                                         start=DEL$start.x,
                                                         end=DEL$start.x ))

coord.DEL.partner2 <- makeGRangesFromDataFrame(data.frame(chr=DEL$CHR2.x, 
                                                         start=DEL$END.x,
                                                         end=DEL$END.x ))

#### in order to merge these GRANGES objects :

coord.DEL.partners <- c(coord.DEL.partner1, coord.DEL.partner2)
}

# length(coord.DEL.partners)

##### if there are no deletions, introducing a DUMMY variable in order to avoid error messages:

if (dim(DEL)[1] == 0)
{
coord.DEL.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(DUP)[1] > 0)
{
coord.DUP.partner1 <- makeGRangesFromDataFrame(data.frame(chr=DUP$seqnames.x, 
                                                         start=DUP$start.x,
                                                         end=DUP$start.x ))

coord.DUP.partner2 <- makeGRangesFromDataFrame(data.frame(chr=DUP$CHR2.x, 
                                                         start=DUP$END.x,
                                                         end=DUP$END.x ))

#### in order to merge these GRANGES objects:

coord.DUP.partners <- c(coord.DUP.partner1, coord.DUP.partner2)
}

# length(coord.DUP.partners)

##### if there are no duplications, introducing a DUMMY variable in order to avoid error messages:

if (dim(DUP)[1] == 0)
{
coord.DUP.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(INS)[1] > 0)
{
coord.INS.partner1 <- makeGRangesFromDataFrame(data.frame(chr=INS$seqnames.x, 
                                                         start=INS$start.x,
                                                         end=INS$start.x ))

coord.INS.partner2 <- makeGRangesFromDataFrame(data.frame(chr=INS$CHR2.x, 
                                                         start=INS$END.x,
                                                         end=INS$END.x ))

#### in order to merge these GRANGES objects :

coord.INS.partners <- c(coord.INS.partner1, coord.INS.partner2)
}

# length(coord.INS.partners)

##### if there are no insertions, introducing a DUMMY variable in order to avoid error messages:

if (dim(INS)[1] == 0)
{
coord.INS.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(INV)[1] > 0)
{

coord.INV.partner1 <- makeGRangesFromDataFrame(data.frame(chr=INV$seqnames.x, 
                                                         start=INV$start.x,
                                                         end=INV$start.x ))

coord.INV.partner2 <- makeGRangesFromDataFrame(data.frame(chr=INV$CHR2.x, 
                                                         start=INV$END.x,
                                                         end=INV$END.x ))

### in order to merge these GRANGES objects :

coord.INV.partners <- c(coord.INV.partner1, coord.INV.partner2)
}

# length(coord.INV.partners)

##### if there are no inversions, introducing a DUMMY variable in order to avoid error messages:

if (dim(INV)[1] == 0)
{
coord.INV.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(TRA)[1] > 0)
{
coord.TRA.partner1 <- makeGRangesFromDataFrame(data.frame(chr=TRA$seqnames.x, 
                                                         start=TRA$start.x,
                                                         end=TRA$start.x ))

coord.TRA.partner2 <- makeGRangesFromDataFrame(data.frame(chr=TRA$CHR2.x, 
                                                         start=TRA$END.x,
                                                         end=TRA$END.x ))

#### in order to merge these GRANGES objects :

coord.TRA.partners <- c(coord.TRA.partner1, coord.TRA.partner2)
}

# length(coord.TRA.partners)

##### if there are no translocations, introducing a DUMMY variable in order to avoid error messages:

if (dim(TRA)[1] == 0)
{
coord.TRA.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
### in this version we print only SV and CNV

pdf(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.pdf", sep=""), width=10, height=12)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="green", data.panel=1, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )
# kpPlotRegions(kp, coord.SV.partners, col="magenta", data.panel=1, r0=0.21, r1=0.4 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=2, r0=0.21, r1=0.4 ) 

dev.off()

### PNG :

png(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.png", sep=""), 
                          width=10, height=12, units="cm", res = 300, pointsize = 10)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="green", data.panel=1, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )
# kpPlotRegions(kp, coord.SV.partners, col="magenta", data.panel=1, r0=0.21, r1=0.4 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=2, r0=0.21, r1=0.4 ) 

dev.off()

###############################################################################################################################################
###############################################################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# [3] "seqnames.x"                                           
# [4] "start.x"                                              
# [5] "CHR2.x"                                               
# [6] "END.x"  
# [74] "Sample"
################################################################################################
######################## in these displays we show the LENGTH of DEL < DUP < INS < INV

if (dim(DEL)[1] > 0)
{
coord.DEL.length <- makeGRangesFromDataFrame(data.frame(chr=DEL$seqnames.x, 
                                                         start=DEL$start.x,
                                                         end=DEL$END.x ))
} else { coord.DEL.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1)) }

####

if (dim(DUP)[1] > 0)
{
coord.DUP.length <- makeGRangesFromDataFrame(data.frame(chr=DUP$seqnames.x, 
                                                         start=DUP$start.x,
                                                         end=DUP$END.x ))
} else {coord.DUP.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

####

if (dim(INS)[1] > 0)
{
coord.INS.length <- makeGRangesFromDataFrame(data.frame(chr=INS$seqnames.x, 
                                                         start=INS$start.x,
                                                         end=INS$END.x ))
} else {coord.INS.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

####

if (dim(INV)[1] > 0)
{
coord.INV.length <- makeGRangesFromDataFrame(data.frame(chr=INV$seqnames.x, 
                                                         start=INV$start.x,
                                                         end=INV$END.x ))
} else {coord.INV.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

############################################################################################
############################################################################################
################## for TRA we do keep the display of the breakpoints :
################## as we have defined those breakpoints above
############################################################################################
############################################################################################

if (dim(TRA)[1] > 0)
{
coord.TRA.partner1 <- makeGRangesFromDataFrame(data.frame(chr=TRA$seqnames.x, 
                                                         start=TRA$start.x,
                                                         end=TRA$start.x ))

coord.TRA.partner2 <- makeGRangesFromDataFrame(data.frame(chr=TRA$CHR2.x, 
                                                         start=TRA$END.x,
                                                         end=TRA$END.x ))

#### in order to merge these GRANGES :

coord.TRA.partners <- c(coord.TRA.partner1, coord.TRA.partner2)
}

# length(coord.TRA.partners)

##### if there are no translocations, introducing a DUMMY variable in order to avoid error messages:

if (dim(TRA)[1] == 0)
{
coord.TRA.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

pdf(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.pdf", sep=""), 
                                                                  width=10, height=12)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=2, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=1, r0=0.21, r1=0.4 )
kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=1, r0=0.21, r1=0.4 ) 

dev.off()

#################################################################################################
#################################################################################################

png(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.png", sep=""), 
                 width=10, height=12, units="cm", res = 300, pointsize = 10)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=2, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=1, r0=0.21, r1=0.4 )
kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                     data.panel=1, r0=0.21, r1=0.4 )

dev.off()

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
### to display the TRANSLOCATIONS as LINKS :

# coord.TRA.partner1 
# coord.TRA.partner2 

# kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

# kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="red", data.panel=2, r0=0, r1=0.2 )

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
### to display the TRANSLOCATIONS as ARCS, the chromosomes are represented linearly  

pdf(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.pdf", sep=""), 
                                                                  width=30, height=20)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg38", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=1, r0=0.61, r1=0.8 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0, r1=0.8) 

dev.off()

### PNG :

png(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.png", 
          sep=""), width=30, height=20, units="cm", res = 300, pointsize = 10)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg38", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=1, r0=0.61, r1=0.8 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0, r1=0.8) 

dev.off()

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
################################################################################################# CN.MOPS
### we are going to add the information on CNV (AMP and DEL) from cn.mops

FILE2 <- args[2]  
 
## FILE2 <- "SPCG.list.of.CNV.segments.by.cn.MOPS.LIST.before.filtering.info"
      
name2 <- basename(FILE2)

x2 <- read.delim(FILE2, header=T, sep="\t")

### x2 <- read.delim("ACC303_6M.list.of.CNV.segments.by.cn.MOPS.LIST.before.filtering.info",
###                                                         header=T, sep="\t") 

head(x2)
dim(x2)

###############################################################################################
###############################################################################################
### to assign DEL and AMP :

x2$CNV <- ifelse( ((x2$CN=="CN0") | (x2$CN=="CN1")), "DELETION", "AMPLIFICATION")  

#### and writing for verifications :

write.table(x2, file = paste(name2, ".assign.CNV.txt",sep=""), 
                 sep ="\t",
                 row.names = FALSE,
                 col.names = TRUE) 

###############################################################################################
### to separate into 2 dataframes :

DELETION <- subset(x2, CNV=="DELETION")
head(DELETION)
dim(DELETION)

AMPLIFICATION <- subset(x2, CNV=="AMPLIFICATION")
head(AMPLIFICATION)
dim(AMPLIFICATION)

###############################################################################################
###############################################################################################
###############################################################################################
#### to write for verification purposes :

write.table(DELETION, file = paste(name2, ".assign.CNV.DELETION.txt",sep=""), 
                 sep ="\t",
                 row.names = FALSE,
                 col.names = TRUE) 

write.table(AMPLIFICATION, file = paste(name2, ".assign.CNV.AMPLIFICATION.txt",sep=""), 
                 sep ="\t",
                 row.names = FALSE,
                 col.names = TRUE) 

################################################################################################
################################################################################################
################################################################################################
 
DELETION_r <- makeGRangesFromDataFrame(data.frame(chr=DELETION$seqnames, 
                                                  start=DELETION$start,
                                                  end=DELETION$end ))  


AMPLIFICATION_r <- makeGRangesFromDataFrame(data.frame(chr=AMPLIFICATION$seqnames, 
                                                       start=AMPLIFICATION$start,
                                                       end=AMPLIFICATION$end )) 

############################################################################################
############################################################################################
############################################################################################
### saving as PNG file

png(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.and.CNV.png", sep=""), 
    width=10, height=12, units="cm", 
    res = 300, pointsize = 10)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="brown", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="green", data.panel=2, r0=0.61, r1=0.8 )
# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )

kpPlotRegions(kp, DELETION_r, col="red", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, AMPLIFICATION_r, col="black", data.panel=1, r0=0.41, r1=0.6 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=1, r0=0, r1=0.21 ) 

dev.off()

################################################################################################
############################################################################################
### saving as PDF file

pdf(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.and.CNV.pdf", sep=""), 
    width=10, height=12)

kp <- plotKaryotype(genome="hg38", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="brown", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="green", data.panel=2, r0=0.61, r1=0.8 )
# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )

kpPlotRegions(kp, DELETION_r, col="red", data.panel=1, r0=0.21, r1=0.40 )
kpPlotRegions(kp, AMPLIFICATION_r, col="black", data.panel=1, r0=0.41, r1=0.6 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=1, r0=0, r1=0.21 ) 

dev.off()

############################################################################################
############################################################################################
############################################################################################
############################################################################################
### to make the same plots on the HORIZONTAL AXIS, by adding AMP and DEL 
### PNG 

png(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.and.CNV.png", 
          sep=""), width=30, height=20, units="cm", res = 300, pointsize = 10)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg38", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=1, r0=0.61, r1=0.8 )

kpPlotRegions(kp, DELETION_r, col="red", data.panel=2, r0=0, r1=0.1 )
kpPlotRegions(kp, AMPLIFICATION_r, col="black", data.panel=2, r0=0.11, r1=0.2 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0.21, r1=1) 

dev.off()

############################################################################################
############################################################################################
########## to make the same plots on the HORIZONTAL AXIS, by adding AMP and DEL 
########## save as PDF
############################################################################################

pdf(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.and.CNV.pdf", sep=""), 
                                                                  width=30, height=20)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg38", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="brown", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="green", data.panel=1, r0=0.61, r1=0.8 )

kpPlotRegions(kp, DELETION_r, col="red", data.panel=2, r0=0, r1=0.1 )
kpPlotRegions(kp, AMPLIFICATION_r, col="black", data.panel=2, r0=0.11, r1=0.2 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0.21, r1=1) 

dev.off()

#################################################################################################################
#################################################################################################################
#################################################################################################################
