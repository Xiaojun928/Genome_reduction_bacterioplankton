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

df <-read.table("gnm_feaures_35_species.txt",sep="\t",header=T,stringsAsFactors = F)
tree <- read.newick("35spe_GTDB_tree_rooted.nwk")
write.nexus(tree,file="caper_35.nex",translate = T)  ## to get the translate nexus
tree <- read.nexus("caper_35.nex")
tree$node.label <- NULL ##otherwise error was printed
cmp <- comparative.data(phy = tree,data = df,names.col = "Species", 
                        warn.dropped=TRUE)

##### mutrate vs Gsize  ######
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


df.Roseo <- df[which(df[,1] %in% c("CHUG_sp_HKCCA1288","Ruegeria_pomeroyi_DSS_3","Dinoroseobacter_shibae_DFL_12","Sulfitobacter_sp_EE-36")),c(3,2)]

#prepare R^2 and slope of PGLS and GLM to show it in the plot  
##PGLS
xpos_pgls <- rep(6,3)
ypos_pgls <- c(1E-8,10^(-8.1),10^(-8.2))
Rsqr <- paste("R\u00B2"," = ",round(a$adj.r.squared,3),sep="")
lamda <- paste("lambda = ",round(a$param.CI$lambda$opt,3),sep="")
lamda_n_R_chr <- paste(lamda,Rsqr,sep=", ")
if(a$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(a$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",round(slope,3), slope_p, sep=" ")
labels_pgls <- c("PGLS:", lamda_n_R_chr ,slope_chr) 

#GLM
xpos_glm <- rep(6,3)
ypos_glm <- c(10^(-8.42),10^(-8.5),10^(-8.6))
Rsqr_glm <- paste("R\u00B2"," = ",round(rsqr_glm,3),sep="")
slope_glm <- round(b$coefficients[2,1],3)
if(b$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(b$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",slope_glm,slope_p,sep=" ")
labels_glm <- c("GLM:", Rsqr_glm ,slope_chr) 


p <- ggplot(data = df, aes(x = Genome.size, y = Mutation.rate)) +
  geom_point(color='blue', size =3) + 
  geom_text_repel(aes(label=id),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, size =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  geom_point(data = df.Roseo, col = c('red','red','red','red'),size =3) + labs(x="Genome size (megabases)",y="Base-substitution mutation rate \n per nucleotide site per generation (μ)") +
  annotation_logticks(sides = 'b') +# Show log tick marks at the bottom
  annotate("text", x = xpos_pgls, y = ypos_pgls,label = labels_pgls, hjust = 0) +
  annotate("rect", xmin = 5.9, xmax = 11.5, ymin = 10^(-7.95), ymax = 10^(-8.25),fill =NA, color = "blue",linewidth = 1.5) +
  annotate("text", x = xpos_glm, y = ypos_glm,label = labels_glm, hjust = 0) +
  annotate("rect", xmin = 5.9, xmax = 11.5, ymin = 10^(-8.385), ymax = 10^(-8.65),fill =NA, color = "grey",linetype = "dashed",linewidth = 1.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "GenomeSize_Mutrate.pdf",width = 9.5,height = 6)
p
dev.off()




######### with Ne ############
df <-read.table("gnm_feaures_26_species.txt",sep="\t",header=T,stringsAsFactors = F)
tree <- read.newick("26spe_GTDB_tree_rooted.nwk")
write.nexus(tree,file="caper_26.nex",translate = T)  ## to get the translate nexus
tree <- read.nexus("caper_26.nex")
tree$node.label <- NULL ##otherwise error was printed
cmp <- comparative.data(phy = tree,data = df,
                        names.col = "Species",warn.dropped=TRUE)

###### Mutrate vs Ne ####
pgls.mut.ne<-pgls(log10(Mutation.rate) ~ log10(Ne_median), data = cmp, lambda='ML')
#pgls.mut.ne<-pgls(log10(Mutation.rate) ~ log10(Ne_median), data = cmp, lambda=0.968)
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
rsqr_glm <- rsq::rsq(A) #0.796755
b <- summary(A)
#log10(Ne_median) -0.78538    0.08656  -9.073 1.03e-08 ***
outlierTest(A)
rsq::rsq(A)

df.Roseo <- df[which(df[,1] %in% c("CHUG_sp_HKCCA1288","Ruegeria_pomeroyi_DSS_3","Dinoroseobacter_shibae_DFL_12","Sulfitobacter_sp_EE-36")),c('Ne_median','Mutation.rate')]

#prepare R^2 and slope of PGLS and GLM to show it in the plot  
##PGLS
xpos_pgls <- rep(10^(7.7),3)
ypos_pgls <- c(1E-8,10^(-8.1),10^(-8.2))
Rsqr <- paste("R\u00B2"," = ",round(a$adj.r.squared,3),sep="")
#lamda <- paste("lambda = ",round(a$param.CI$lambda$opt,3),sep="")
lamda <- paste("lambda = ","0.968",sep="")
lamda_n_R_chr <- paste(lamda,Rsqr,sep=", ")
if(a$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(a$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",round(slope,3), slope_p, sep=" ")
labels_pgls <- c("PGLS:", lamda_n_R_chr ,slope_chr) 

#GLM
xpos_glm <- rep(10^(7.7),3)
ypos_glm <- c(10^(-8.42),10^(-8.5),10^(-8.6))
Rsqr_glm <- paste("R\u00B2"," = ",round(rsqr_glm,3),sep="")
slope_glm <- round(b$coefficients[2,1],3)
if(b$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(b$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",slope_glm,slope_p,sep=" ")
labels_glm <- c("GLM:", Rsqr_glm ,slope_chr) 

p <- ggplot(data = df, aes(x = Ne_median, y = Mutation.rate)) + 
  geom_point(color='blue', size =3) + 
  geom_text_repel(aes(label=id),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, linewidth =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotate("text", x = xpos_pgls, y = ypos_pgls,label = labels_pgls, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.4), ymin = 10^(-7.95), ymax = 10^(-8.25),fill =NA, color = "blue",linewidth = 1.5) +
  annotate("text", x = xpos_glm, y = ypos_glm,label = labels_glm, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.4), ymin = 10^(-8.385), ymax = 10^(-8.65),fill =NA, color = "grey",linetype = "dashed",linewidth = 1.5) +
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  geom_point(data = df.Roseo, col = c('red','red','red','red'),size =3)+
  labs(x="Effective population size (Ne)",y="Base-substitution mutation rate \n per nucleotide site per generation (μ)") + theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
                                                                                                                                 axis.line = element_line(colour = 'black', linewidth = 1),
                                                                                                                                 axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "Mut_Ne.pdf",width = 9.5,height = 6)
p
dev.off()



####### MutLoad vs Ne ##########
pgls.mutload.ne <- pgls(log10(Mut_load) ~ log10(Ne_median), 
                        data = cmp, lambda='ML')
a <- summary(pgls.mutload.ne)  
intercept <- a$coefficients[1,1]
slope <- a$coefficients[2,1]
res <- residuals(pgls.mutload.ne,phylo = T)
rownames(res) <- rownames(pgls.mutload.ne$residuals)
res<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res)[(abs(res)>3)]#gives the names of the outliers

B <- glm(log10(Mut_load)~ log10(Ne_median), data = df)
b <- summary(B)
rsqr_glm <-  rsq::rsq(B)
outlierTest(B)

#prepare R^2 and slope of PGLS and GLM to show it in the plot  
##PGLS
xpos_pgls <- rep(10^(7.7),3)
ypos_pgls <- c(1E-2,10^(-2.1),10^(-2.2))

Rsqr <- paste("R\u00B2"," = ",round(a$adj.r.squared,3),sep="")
lamda <- paste("lambda = ",round(a$param.CI$lambda$opt,3),sep="")
lamda_n_R_chr <- paste(lamda,Rsqr,sep=", ")
if(a$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(a$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",round(slope,3), slope_p, sep=" ")
labels_pgls <- c("PGLS:", lamda_n_R_chr ,slope_chr) 

#GLM
xpos_glm <- rep(10^(7.7),3)
ypos_glm <- c(10^(-2.42),10^(-2.5),10^(-2.6))
Rsqr_glm <- paste("R\u00B2"," = ",round(rsqr_glm,3),sep="")
slope_glm <- round(b$coefficients[2,1],3)
if(b$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(b$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",slope_glm,slope_p,sep=" ")
labels_glm <- c("GLM:", Rsqr_glm ,slope_chr) 


df.Roseo <- df[which(df[,1] %in% c("CHUG_sp_HKCCA1288","Ruegeria_pomeroyi_DSS_3","Dinoroseobacter_shibae_DFL_12","Sulfitobacter_sp_EE-36")),c('Ne_median','Mut_load')]
p <- ggplot(data = df, aes(x = Ne_median, y = Mut_load)) + 
  geom_point(color='blue', size =3) + 
  geom_text_repel(aes(label=id),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, size =1) + ##regression line
  annotate("text", x = xpos_pgls, y = ypos_pgls,label = labels_pgls, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.355), ymin = 10^(-1.95), ymax = 10^(-2.25),fill =NA, color = "blue",linewidth = 1.5) +
  annotate("text", x = xpos_glm, y = ypos_glm,label = labels_glm, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.355), ymin = 10^(-2.4), ymax = 10^(-2.65),fill =NA, color = "grey",linetype = "dashed",linewidth = 1.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  geom_point(data = df.Roseo, col = c('red','red','red','red'),size =3)+
  labs(x="Effective population size (Ne)",
       y="Genome-wide mutation rate in protein-coding \n DNA per geenration (Up)") +theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
                                                                                          axis.line = element_line(colour = 'black', linewidth = 1),
                                                                                          axis.ticks = element_line(colour = "black", linewidth = 1))


pdf(file = "Mutload_Ne.pdf",width = 9.5,height = 6)
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

#prepare R^2 and slope of PGLS and GLM to show it in the plot  
##PGLS
xpos_pgls <- rep(10^(7.7),3)
ypos_pgls <- c(10^(0.3),10^(0.26),10^(0.22))
Rsqr <- paste("R\u00B2"," = ",round(a$adj.r.squared,3),sep="")
lamda <- paste("lambda = ",round(a$param.CI$lambda$opt,3),sep="")
lamda_n_R_chr <- paste(lamda,Rsqr,sep=", ")
if(a$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(a$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",round(slope,3), slope_p, sep=" ")
labels_pgls <- c("PGLS:", lamda_n_R_chr ,slope_chr) 

#GLM
xpos_glm <- rep(10^(7.7),3)
ypos_glm <- c(10^(0.16),10^(0.12),10^(0.08))
Rsqr_glm <- paste("R\u00B2"," = ",round(rsqr_glm,3),sep="")
slope_glm <- round(b$coefficients[2,1],3)
if(b$coefficients[2,4] < 0.001)
{
  slope_p <- "(p < 0.001)"
}else{
  slope_p <- paste("(p = ",round(b$coefficients[2,4],3),")",sep="")
}
slope_chr <- paste("slope =",slope_glm,slope_p,sep=" ")
labels_glm <- c("GLM:", Rsqr_glm ,slope_chr) 

df.Roseo <- df[which(df[,1] %in% c("CHUG_sp_HKCCA1288","Ruegeria_pomeroyi_DSS_3","Dinoroseobacter_shibae_DFL_12","Sulfitobacter_sp_EE-36")),c("Genome.size","Ne_median")]
p <- ggplot(data = df, aes(x = Ne_median, y = Genome.size)) + 
  geom_point(color='blue', size =3) + 
  geom_text_repel(aes(label=id),size=5,max.overlaps = Inf) +
  geom_abline(color='blue',intercept = intercept, slope = slope, size =1) + ##regression line
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  annotate("text", x = xpos_pgls, y = ypos_pgls,label = labels_pgls, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.3), ymin = 10^(0.188), ymax = 10^(0.34),fill =NA, color = "blue",linewidth = 1.5) +
  annotate("text", x = xpos_glm, y = ypos_glm,label = labels_glm, hjust = 0) +
  annotate("rect", xmin = 10^(7.65), xmax = 10^(8.3), ymin = 10^(0.06), ymax = 10^(0.18),fill =NA, color = "grey",linetype = "dashed",linewidth = 1.5) +
  theme_classic() + geom_smooth(method='glm',formula=y~x,se=FALSE,colour="grey",linetype=2, size =1.2) +
  geom_point(data = df.Roseo, col = c('red','red','red','red'),size =3)+
  labs(x="Effective population size (Ne)",y="Genome size (megabases)") +
  annotation_logticks(sides = 'l') +# Show log tick marks in the left side 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 1))

pdf(file = "Gsize_Ne.pdf",width = 9.5,height = 6)
p
dev.off()


############test by leaving one out ##########
df0 <-read.table("gnm_feaures_26_species.txt",sep="\t",header=T,stringsAsFactors = F)

setwd("E:/CHUG_reduction/MA/synthesis/PGLS/PGLS_base_GTDB_conct_tree")
tree_list <- list.files("E:/CHUG_reduction/MA/synthesis/PGLS/PGLS_base_GTDB_conct_tree",pattern = "prune.*nwk")
mtx <- data.frame(species = df0$Species,R2 = rep(0,nrow(df0)),
                  slope = rep(0,nrow(df0)), p = rep(0,nrow(df0)))
rownames(mtx) <- df0$Species

tree_list <- tree_list[-11] #Mesoplasma failed to run PGLS when lambda was setting as 'ML'

for (trefile in tree_list) {
  tree <- read.tree(trefile)
  print(trefile)
  species <- sub("prune_","",sub(".nwk","",basename(trefile)))
  nex_filename <- paste(species,"pruned","nex",sep=".")
  write.nexus(tree,file=nex_filename,translate = T)  ## to get the translate nexus
  tree <- read.nexus(nex_filename)
  tree$node.label <- NULL ##otherwise error was printed
  df  <- df0[which(df0$Species != species),]
  cmp <- comparative.data(phy = tree,data = df,names.col = "Species", warn.dropped=TRUE)
  pgls.mut.ne<-pgls(log10(Genome.size) ~ log10(Ne_median), data = cmp, lambda='ML')
  a <- summary(pgls.mut.ne)  
  mtx[species,"slope"] <- a$coefficients[2,1]
  mtx[species,"p"] <- a$coefficients[2,4]
  mtx[species,"R2"] <- a$adj.r.squared[1,1]
}

mtx[,"gnm_size"] <- df0$Genome.size
#if mtx$p < 0.05, suggesting that the species support the negative relationship between Ne adn genome size
mtx %>% filter(p > 0.05)




