suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(qqman)))
suppressWarnings(suppressMessages(library(tidyverse)))

toptrait <- read_tsv("GWAS_resfile_1.txt")
bottomtrait <- read_tsv("GWAS_resfile_2.txt")

#colnames(toptrait)[which(colnames(toptrait)=="mtag_pval")]="P"
#colnames(bottomtrait)[which(colnames(bottomtrait)=="mtag_pval")]="P"

toptrait<-subset(toptrait, !is.na(P))
bottomtrait<-subset(bottomtrait, !is.na(P))

nsitestop=nrow(toptrait)
nsitesbottom=nrow(bottomtrait)

toptrait<-subset(toptrait,P <= 1e-2)
bottomtrait<-subset(bottomtrait,P <= 1e-2)

toptrait$CHR=gsub("chrX|X","chr23",toptrait$CHR)
toptrait$CHR=as.numeric(gsub("chr","",toptrait$CHR))

bottomtrait$CHR=gsub("chrX|X","chr23",bottomtrait$CHR)
bottomtrait$CHR=as.numeric(gsub("chr","",bottomtrait$CHR))

toptrait=toptrait[toptrait[,do.call(order, .SD), .SDcols = c("CHR","BP")]]
bottomtrait=bottomtrait[bottomtrait[,do.call(order, .SD), .SDcols = c("CHR","BP")]]

toptrait$BP <- as.numeric(toptrait$BP)
bottomtrait$BP <- as.numeric(bottomtrait$BP)

toptrait <- toptrait %>%
    select(SNP,CHR,BP,P) %>%
    mutate(P = -log(P, 10))

combineddat <- bottomtrait %>%
    select(SNP,CHR,BP,P) %>%
    mutate(P = log(P, 10)) %>%
    rbind(toptrait)

library(grid)

combineddat$CHR <- as.factor(combineddat$CHR)

png("Miami_plot.png",type='cairo',width=12,height=6,units='in',res=300)

ggplot(combineddat, aes(x=BP, y=P, color=CHR)) +
    geom_point(alpha=0.8, size=0.5) +
    facet_grid(~CHR, scales= "free_x", space="free_x",switch="x") +
    scale_color_manual(values = rep(c("black", "darkgray"), 23 )) +
    theme_bw() +
    scale_y_continuous(breaks=c(-50,-20,-8,-4,0,4,8,20,50), label=c(50,20,8,4,0,4,8,20,50)) +
    geom_hline(yintercept=7.30103, linetype="dashed") +
    geom_hline(yintercept=-7.30103, linetype="dashed") +
    labs(x = "Chromosome") +
    scale_y_log10() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text.x = element_text(size=8),
      strip.background = element_rect(colour="white", fill="white"),
      panel.margin = unit(0, "lines"),
      axis.ticks.x=element_blank()
    )

dev.off()
