setwd("E:/CHUG_reduction/MA/synthesis/PGLS")
library(ape)
library(nlme)
library(phytools)
library(picante)
library(car)
library(MuMIn)
library(ggplot2)
library(scales)
library(caper)
library(ggrepel)
library(rsq)
library(dplyr)

###########for Roseo Only #######
df0 <-read.table("gnm_feaures_26_species.txt",sep="\t",header=T,stringsAsFactors = F)
df <- df0[c(13,24,25,26),]
tree <- read.newick("Roseo_GTDB_tree_rooted.nwk")
write.nexus(tree,file="caper_roseo.nex",translate = T)  ## to get the translate nexus
tree <- read.nexus("caper_roseo.nex")
tree$node.label <- NULL ##otherwise error was printed
cmp <- comparative.data(phy = tree,data = df,
                        names.col = "Species",warn.dropped=TRUE)

#########Mut rate vs Gsize ############
pgls.mutload.size<-pgls(log10(Mutation.rate)~ log10(Genome.size), data = cmp, lambda ='ML')
summary(pgls.mutload.size)
a <- summary(pgls.mutload.size)  
intercept <- a$coefficients[1,1]
slope <- a$coefficients[2,1]

res <- residuals(pgls.mutload.size,phylo = T)
rownames(res) <- rownames(pgls.mutload.size$residuals)
res<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res)[(abs(res)>3)]

##GLM regression
BB <- glm(log10(Mutation.rate)~ log10(Genome.size), data = df)
rsqr_glm <- rsq::rsq(BB) 
b <- summary(BB)
#log10(Mut_load) -0.14756    0.09534  -1.548    0.132
outlierTest(BB)


p <- ggplot(data = df, aes(x = Genome.size, y = Mutation.rate)) +
  geom_point(color='red', size =3) + 
  geom_text_repel(aes(label=Species),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, linewidth =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, linewidth =1.2) +
  labs(x="Genome size (megabases)",y="Base-substitution mutation rate \n per nucleotide site per generation (μ)") +
  annotation_logticks(sides = 'b') +# Show log tick marks at the bottom
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "GenomeSize_Mutrate_Roseo.pdf",width = 9.5,height = 6)
p
dev.off()

##### Gsize vs Ne #####
pgls.mut.ne<-pgls(log10(Genome.size) ~ log10(Ne_median), data = cmp, lambda='ML')
a <- summary(pgls.mut.ne)  
intercept <- a$coefficients[1,1]
slope <- a$coefficients[2,1]
res <- residuals(pgls.mut.ne,phylo = T)
rownames(res) <- rownames(pgls.mut.ne$residuals)
res<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res)[(abs(res)>3)]#gives the names of the outliers

A <- glm(log10(Genome.size)~ log10(Ne_median), data = df)
Rsqr_glm <- rsq::rsq(A) #0.796755
b <- summary(A)
outlierTest(A)


p <- ggplot(data = df, aes(x = Ne_median, y = Genome.size)) + 
  geom_point(color='red', size =3) + 
  geom_text_repel(aes(label=Species),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, linewidth =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  labs(x="Effective population size (Ne)",y="Genome size (megabases)") +
  annotation_logticks(sides = 'l') +# Show log tick marks in the left side 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "Gsize_Ne_Roseo.pdf",width = 9.5,height = 6)
p
dev.off()

##### Gsize vs indel rate #####
df0 <-read.table("gnm_feaures_35_species.txt",sep="\t",header=T,stringsAsFactors = F)
df <- df0[c(13,33:35),]
tree <- read.newick("Roseo_GTDB_tree_rooted.nwk")
write.nexus(tree,file="caper_roseo.nex",translate = T)  ## to get the translate nexus
tree <- read.nexus("caper_roseo.nex")
tree$node.label <- NULL ##otherwise error was printed
cmp <- comparative.data(phy = tree,data = df,
                        names.col = "Species",warn.dropped=TRUE)

pgls.indel.ne<-pgls(log10(Genome.size) ~ log10(Indel_rate), data = cmp, lambda='ML')
a <- summary(pgls.indel.ne)  
intercept <- a$coefficients[1,1]
slope <- a$coefficients[2,1]
res <- residuals(pgls.indel.ne,phylo = T)
rownames(res) <- rownames(pgls.indel.ne$residuals)
res<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res)[(abs(res)>3)]#gives the names of the outliers

A <- glm(log10(Genome.size)~ log10(Indel_rate), data = df)
Rsqr_glm <- rsq::rsq(A) #0.796755
b <- summary(A)
outlierTest(A)


p <- ggplot(data = df, aes(x = Indel_rate, y = Genome.size)) + 
  geom_point(color='red', size =3) + 
  geom_text_repel(aes(label=Species),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, size =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  labs(x="Indel rate",y="Genome size (megabases)") +
  annotation_logticks(sides = 'l') +# Show log tick marks in the left side 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "Gsize_indelrate_Roseo.pdf",width = 9.5,height = 6)
p
dev.off()


#########Mutrate vs Ne ######
pgls.mut.ne<-pgls(log10(Mutation.rate) ~ log10(Ne_median), data = cmp, lambda='ML')
a <- summary(pgls.mut.ne)  
intercept <- a$coefficients[1,1]
slope <- a$coefficients[2,1]


#res1 <- pgls.mut.ne$residuals ##non-phylogenetic residual
res <- residuals(pgls.mut.ne,phylo = T)
rownames(res) <- rownames(pgls.mut.ne$residuals)
res<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res)[(abs(res)>3)]#gives the names of the outliers

##GLM regression
A <- glm(log10(Mutation.rate)~ log10(Ne_median), data = df)
rsqr_glm <- rsq::rsq(A) 
b <- summary(A)
#log10(Ne_median) -0.78538    0.08656  -9.073 1.03e-08 ***
outlierTest(A)
rsq::rsq(A)



p <- ggplot(data = df, aes(x = Ne_median, y = Mutation.rate)) + 
  geom_point(color='red', size =3) + 
  geom_text_repel(aes(label=Species),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, linewidth =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  labs(x="Effective population size (Ne)",y="Base-substitution mutation rate \n per nucleotide site per generation (μ)") + theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
                                                                                                                                 axis.line = element_line(colour = 'black', linewidth = 1),
                                                                                                                                 axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "Mut_Ne_Roseo.pdf",width = 9.5,height = 6)
p
dev.off()
