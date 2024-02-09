setwd("E:/CHUG_reduction/MA/")
library(ggplot2)
require(ggpubr)
library(dplyr)


##new approach using dplyr package
##For CHUG
library(dplyr)
df <- read.table('CHUG/CHUG_mut_summary_refNCBI.txt',sep='\t',header = T,stringsAsFactors = F)
cds <- df[which(df$genetype =='CDS'),]
#cds_cnt <- count(cds, ~genepos+muttype)
cds_cnt <- cds %>%
  group_by(genepos, muttype) %>%
  summarise(freq = n())

cds_cnt$genepos <- as.factor(cds_cnt$genepos)#run this step, otherwise the bar would be too thick to be seen.

ncds <- df[which(df$genetype !='CDS'),]
#ncds_cnt <- count(ncds, ~genepos+muttype)
ncds_cnt <- ncds %>% group_by(genepos,muttype) %>% summarise(freq = n())

p1 <- ggplot() +
  geom_bar(data=cds_cnt, aes(x=as.numeric(as.character(genepos)), y=freq, fill=muttype), stat="identity") +
  scale_fill_manual(values=c('bps'="#333333",'del'="#FFD000",'ins'="#FF6670")) +
  geom_point(data=ncds_cnt,aes(x=genepos,y=0,shape=muttype,color=muttype,fill=muttype),alpha=0.65,size=2) +
  scale_shape_manual(values=c('bps'=23,'del' = 24,'ins'= 25))+
  scale_colour_manual(values=c('bps'="#333333",'del'='#FFD000','ins'="#FF6670")) +
  coord_cartesian(ylim=c(0,41),xlim = (c(0,2663031))) + 
  labs(x = "Genome position (bases)",y = "No. of mutation per CDS") +
  theme_classic()


#hotspots <- c(53992.5, 144618.5, 211390, 656508, 854157, 900612, 954924.5, 1219769.5, 1581725.5, 1616488,1859572.5,1886915,1927163.5,1949103,2383738,2475118.5,2625995) #p<0.05, obsolete ref
hotspots <- c(53992.5, 144618.5, 656508, 854157, 900612, 954924.5, 1219769.5, 1581725.5, 1587338, 1616488, 1886915, 1927163.5, 1859572.5, 1949103, 2383738, 2475118.5, 2625980) #p<0.05, NCBI ref
hotspots_mut <- df[which(df$genepos %in% hotspots),]
#hot_mut_cnt <- count(hotspots_mut, ~genepos)
hot_mut_cnt <- hotspots_mut %>% group_by(genepos,locus) %>% summarise(freq = n())
## add annotation of hotspot genes
for (i in 1:length(hotspots)) {
  p <- p1 + annotate("segment", x = as.numeric(hot_mut_cnt[i,1]), xend = as.numeric(hot_mut_cnt[i,1]), y = as.numeric(hot_mut_cnt[i,3])+2, yend = as.numeric(hot_mut_cnt[i,3]), colour = "purple",size=0.3, arrow=arrow()) +
    annotate("text",x=as.numeric(hot_mut_cnt[i,1]),y=as.numeric(hot_mut_cnt[i,3])+5,label= hot_mut_cnt[i,2],size=2.5)
  print(as.numeric(hot_mut_cnt[i,1]))
  p1 <- p
}

pdf(file = "CHUG/CHUG_mut_dist_refNCBI.pdf",width = 10,height = 4)
p1
dev.off()

############Sulfitobacter ######
chr_len = data.frame(chr=c('tig00000001', 'tig00000020'),
                     len = c(2988286,349750))
bound <- chr_len
colnames(bound) <- c('chr','end')
for (i in 2:dim(chr_len)[1]) {
  bound[i,2] <- sum(chr_len[1:i,2])
}

##plot mutations CDS by CDS, using dots to represent mut in non-CDS
df.ori <- read.table('Sulfito/Sulfito_mut_summary.txt',sep='\t',header = T,stringsAsFactors = F)

chr <- df.ori %>% filter(replicon == chr_len[1,1]) 
df <- chr
for (i in 2:dim(chr_len)[1]) {
  pls <- df.ori %>% filter(replicon == chr_len[i,1]) %>% 
    mutate(genepos = genepos + sum(chr_len[(1:(i-1)),2]) + 1)
  df <- rbind(df,pls)
}


cds <- df[which(df$genetype =='CDS'),]
cds_cnt <- cds %>% group_by(genepos,muttype) %>% count()
cds_cnt$genepos <- as.factor(cds_cnt$genepos)#run this step, otherwise the bar would be too thick to be seen.

ncds <- df[which(df$genetype !='CDS'),]
ncds_cnt <- ncds %>% group_by(genepos,muttype) %>% count()

p1 <- ggplot() +
  geom_bar(data=cds_cnt, aes(x=as.numeric(as.character(genepos)), y=n, fill=muttype), stat="identity") +
  scale_fill_manual(values=c('bps'="#333333",'del'="#FFD000",'ins'="#FF6670")) +
  geom_point(data=ncds_cnt,aes(x=genepos,y=0,shape=muttype,color=muttype,fill=muttype),alpha=0.65,size=2) +
  scale_shape_manual(values=c('bps'=23,'del' = 24,'ins'= 25))+
  scale_colour_manual(values=c('bps'="#333333",'del'='#FFD000','ins'="#FF6670")) +
  coord_cartesian(ylim=c(0,41),xlim = (c(0,3338036))) + 
  geom_vline(xintercept = chr_len[1,2],linetype="dotted") +
  labs(x = "Genome position (bases)",y = "No. of mutation per CDS") +
  theme_classic()


hotspots <- c(2088497.5, 110400, 769196) #p<0.05 #bug in check 2088498
hotspots_mut <- df[which(df$genepos %in% hotspots),]
hot_mut_cnt <- count(hotspots_mut, genepos)

## add annotation of hotspot genes
for (i in 1:length(hotspots)) {
  p <- p1 + annotate("segment", x = as.numeric(hot_mut_cnt[i,1]), xend = as.numeric(hot_mut_cnt[i,1]), y = hot_mut_cnt[i,2]+1, yend = hot_mut_cnt[i,2], colour = "purple",size=0.3, arrow=arrow())
  p1 <- p
}

pdf(file = "Sulfito/Sulfito_mut_dist_uniscale.pdf",width = 10,height = 4)
p1
dev.off()



#for DFL12
chr_len = data.frame(chr=c('NC_009952', 'NC_009955', 'NC_009956', 'NC_009957', 'NC_009958', 'NC_009959'),
                     len = c(3789584,190506,152970,126304,86208,72296))
bound <- chr_len
colnames(bound) <- c('chr','end')
for (i in 2:dim(chr_len)[1]) {
  bound[i,2] <- sum(chr_len[1:i,2])
}

##plot mutations CDS by CDS, using dots to represent mut in non-CDS
df.ori <- read.table('DFL12/DFL12_mut_summary.txt',sep='\t',header = T,stringsAsFactors = F)

chr <- df.ori %>% filter(replicon == chr_len[1,1]) 
df <- chr
for (i in 2:dim(chr_len)[1]) {
  pls <- df.ori %>% filter(replicon == chr_len[i,1]) %>% 
    mutate(genepos = genepos + sum(chr_len[(1:(i-1)),2]) + 1)
  df <- rbind(df,pls)
}


cds <- df[which(df$genetype =='CDS'),]
cds_cnt <- cds %>% group_by(genepos,muttype,locus) %>% count()
cds_cnt$genepos <- as.factor(cds_cnt$genepos)#run this step, otherwise the bar would be too thick to be seen.

ncds <- df[which(df$genetype !='CDS'),]
ncds_cnt <- ncds %>% group_by(genepos,muttype) %>% count()

## plot the lower part of histogram due to the hotspot RS15215

p1 <- ggplot() +
  geom_bar(data=cds_cnt, aes(x=as.numeric(as.character(genepos)), y=n, fill=muttype),stat="identity") +
  scale_fill_manual('CDS',values=c('bps'="#333333",'del'="#FFD000",'ins'="#FF6670")) +
  geom_point(data=ncds_cnt,aes(x=genepos,y=0,shape=muttype,color=muttype,fill=muttype),
             alpha=0.65,size=2) +
  scale_shape_manual(values=c('bps'=23,'del' = 24,'ins'= 25))+
  scale_colour_manual(values=c('bps'="#333333",'del'='#FFD000','ins'="#FF6670")) +
  coord_cartesian(ylim=c(0,41),xlim = (c(0,4417868))) + 
  geom_vline(data = bound, mapping = aes(xintercept = end),linetype="dotted") + 
  labs(x = "Genome position (bases)",y = "No. of mutation per CDS") +
  theme_classic()

## plot the upper part if necesary
#p2 <- ggplot(data=cds_cnt, aes(x=as.numeric(as.character(genepos)),y=n)) +
#  geom_bar(stat="identity",aes(fill=muttype)) +
#  scale_fill_manual('mut.type',values=c('bps'="#333333",'del'='#FFD000','ins'="#FF6670")) +
#  coord_cartesian(ylim=c(25,41),xlim = (c(0,4417868))) + 
#  labs(x=NULL,y=NULL,fill=NULL) +
#  theme_classic() +
#  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y=element_text(size=10))

## add annotation of hotspot genes  
sorted_cnt <- cds_cnt %>% arrange(desc(n)) #p<0.05 
hotspots_mut <- df[which(df$genepos %in% head(sorted_cnt$genepos,8)),]
hot_mut_cnt <- hotspots_mut %>% group_by(genepos,locus) %>% count()

hot_mut_cnt <- as.data.frame(hot_mut_cnt)
for (i in 1:dim(hot_mut_cnt)[1]) {
  p <- p1 + annotate("segment", x = hot_mut_cnt[i,1], xend = hot_mut_cnt[i,1], y = hot_mut_cnt[i,3]+5, yend = hot_mut_cnt[i,3], colour = "purple",size=0.3, arrow=arrow()) +
    annotate("text",x=hot_mut_cnt[i,1],y=hot_mut_cnt[i,3]+5,label= hot_mut_cnt[i,2],size=2.5)
  p1<- p
}


#p <- ggarrange(p2,p1,heights=c(1/4, 3/4),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")

pdf("DFL12/DFL12_mut_dist.pdf",height = 4,width = 10)
p
dev.off()


#######plot mutations by genomic windows (10Kbp) ##########
df <- read.table('CHUG/CHUG_mut_summary.txt',sep='\t',header = T,stringsAsFactors = F)

## plot the last part of histogram due to the hotspot
p_by_window <- ggplot() + 
  geom_histogram(data = df,aes(x=mutpos,fill=muttype),breaks=seq(1,2663031,by=10000),alpha=.6,position = 'identity') +
  scale_fill_manual(values=c("#333333","#FFD000","#FF6670")) +
  coord_cartesian(ylim=c(0,35),xlim = (c(1,2663031))) + 
  labs(x = "Genome position (bases)",y = "Number of mutations  per 10 kilobases") +
  theme_classic()

cds <- df[which(df$genetype =='CDS'),]
cds_cnt <- cds %>% group_by(genepos,muttype) %>% count()
cds_cnt$genepos <- as.factor(cds_cnt$genepos)#run this step, otherwise the bar would be too thick to be seen.

ncds <- df[which(df$genetype !='CDS'),]
ncds_cnt <- ncds %>% group_by(genepos,muttype) %>% count()
ncds_cnt$genepos <- as.factor(ncds_cnt$genepos)

p_by_gene <- ggplot() +
  geom_bar(data=cds_cnt, aes(x=as.numeric(as.character(genepos)), y=n, fill=muttype), stat="identity") +
  scale_fill_manual(values=c('bps'="#333333",'del'="#FFD000",'ins'="#FF6670")) +
  geom_point(data=ncds_cnt,aes(x=as.numeric(as.character(genepos)),y=0,shape=muttype,color=muttype,fill=muttype),alpha=0.65,size=2) +
  scale_shape_manual(values=c('bps'=23,'del' = 24,'ins'= 25))+
  scale_colour_manual(values=c('bps'="#333333",'del'='#FFD000','ins'="#FF6670")) +
  coord_cartesian(ylim=c(0,27),xlim = (c(0,2663031))) + 
  labs(x = "Genome position (bases)",y = "No. of mutation per CDS") +
  theme_classic()


hotspots <- c(53992.5, 144618.5, 211390, 656508, 854157, 900612, 954924.5, 1219769.5, 1581725.5, 1616488,1859572.5,1886915,1927163.5,1949103,2383738,2475118.5,2625995) #p<0.05
hotspots_mut <- df[which(df$genepos %in% hotspots),]
hot_mut_cnt <- hotspots_mut %>% count(genepos)

## add annotation of hotspot genes
for (i in 1:length(hotspots)) {
  p1 <- p_by_gene + annotate("segment", x = hot_mut_cnt[i,1], xend = hot_mut_cnt[i,1], y = hot_mut_cnt[i,2]+2, yend = hot_mut_cnt[i,2], colour = "purple",size=0.3, arrow=arrow())
  p_by_gene <- p1
}

p <- ggarrange(p_by_window,p_by_gene,ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v") 

pdf(file = "CHUG/CHUG_mut_dist_combined.pdf",width = 12,height = 7.5)
p
dev.off()


####comparison between real mut rate and expected mut rate
##CHUG
m <- 505    #num of BPS
T <- 180 * 1288 # #MA lines * #generations
s <- 2663031    # genome size
rate <- poisson.test(m)
rate$estimate/(T*s)
rate$conf.int/(T*s)
#7.481570e-10 8.925012e-10

m2 <- 259
T2 <- 94 * 1288
s <- 2663031
rate2 <- poisson.test(m2)
rate2$estimate/(T2*s)
rate2$conf.int/(T2*s)
#7.084424e-10 9.073292e-10

##Sulfito
m <- 243    #num of BPS
T <- 51 * 4320 # #MA lines * #generations
s <- 3338036    # genome size
rate <- poisson.test(m)
rate$estimate/(T*s)
rate$conf.int/(T*s)
#2.901749e-10 3.746777e-10

m2 <- 280
T2 <- 58 * 4320
s <- 3338036
rate2 <- poisson.test(m2)
rate2$estimate/(T2*s)
rate2$conf.int/(T2*s)
#2.967087e-10 3.763747e-10


##DFL12
m <- 80    #num of BPS
T <- 62 * 1701 # #MA lines * #generations
s <- 4417868    # genome size
rate <- poisson.test(m)
rate$estimate/(T*s)
rate$conf.int/(T*s)
#1.361509e-10 2.137009e-10

m2 <- 218
T2 <- 149 * 1701
s <- 4417868
rate2 <- poisson.test(m2)
rate2$estimate/(T2*s)
rate2$conf.int/(T2*s)
#1.697055e-10 2.223266e-10


col <- c("#A52A2A","#800000","#bca052")
