
library(reshape2)
library(vioplot)
data_all <- read.table('FRiP_data_file.tsv', header=TRUE, sep="\t")
c1 <- subset(data_all, group == 'real')$value
c2 <- subset(data_all, group == 'random')$value
#pdf('FRiP.pdf')
vioplot(c1, c2, col = "#91D1C2B2")
#dev.off()


