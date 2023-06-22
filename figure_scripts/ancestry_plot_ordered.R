## SUPPLEMENTAL FIG 1 ### 

library(ggplot2)
library(dplyr)
library(reshape2)

#cols <- c(EUR='#CA4136',AFR='#1F619E',EAS='#496849', SAS='#819D86',AMR='#B7E0EB')
#cols <- c(EUR='#CA4136',AFR='#1F619E',EAS='#819D86', SAS='#8B4E85',AMR='#FDC652')
cols <- c(EUR='#F9C5C2',AFR='#B7E0EB',EAS='#819D86', SAS='#8B4E85',AMR='#FDC652') # Using color palette colors but not sticking to our ancestry colors

dat <- read.csv('mmc2_modified.csv', header=T, row.names=1, check.names=F)
dat$Sample <- rownames(dat)


# Pick a cancer to plot
#df <- melt(dat[ dat$tumor_type=='BRCA',c(1:2,8:13)]) # Subset to cancer type of interest
#df <- melt(dat[ dat$tumor_type=='UCEC',c(1:2,8:13)]) # Subset to cancer type of interest
df <- melt(dat[ dat$tumor_type=='THCA',c(1:2,8:13)]) # Subset to cancer type of interest

# Remove 'Admixture_percent_' from the ancestry labels
df$variable <- substring(df$variable,19,22)
df$variable <- factor(df$variable, levels=c('EUR','AFR','EAS','SAS','AMR'))


# TODO: need to replace color palette with paper colors
df <- df[with(df, order(consensus_ancestry,variable,-value)), ]
df <- df[complete.cases(df),]
df$Sample <- factor( as.character(df$Sample), levels = unique(df$Sample))

ggplot(df, aes(fill=variable, y=value, x=Sample)) + geom_bar(position="fill", stat="identity") +   scale_fill_manual(values=cols)
#ggsave('ancestry_barplot_BRCA.png', width=10, height=4)
#ggsave('ancestry_barplot_UCEC.png', width=10, height=4)
ggsave('ancestry_barplot_THCA.png', width=10, height=4)
# don't forget to save - wide and long!

