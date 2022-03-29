##This script produces plots useful to evaluate the distribution of ref and alt allelic depthsin a vcf. Must be used together with vcf_sites_refalt_counts.py which preprocesses vcf files.


library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

for(a in args){
refalt <- read.table(a)
names(refalt) <- c("cont", "pos", "ref", "alt")

png(filename=sub("\\.[^\\.]*$", ".png", a), width=300, height=500)

scatter <- ggplot(refalt) + 
	geom_point(aes(x=ref, y=alt)) +
	geom_abline(slope=1) +
	theme_bw()

histo <- ggplot(refalt) +
	geom_histogram(aes(x=ref/(ref+alt))) +
	theme_bw()

grid.arrange(scatter, histo)
dev.off()
}

