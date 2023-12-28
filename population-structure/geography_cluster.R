library('ggtree')
library('tidyverse')
library(ggpubr)
library("RColorBrewer")
library(pophelper)
library(scatterpie)
library(sf)
library(raster)
library(tidyr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
world=st_read('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/shp/States_shapefile_shp/States_shapefile.shp')%>%
  subset(State_Code!='AK'&State_Code!='HI')
pop=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/allKMap.csv')%>%
  dplyr::select('Population.ID','Latitude','Longitude')
setwd('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/samples')
EV=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/environmental_factors.csv')
sp_name='pet_subs_pet'
file=paste('./all_',sp_name,'_samples',sep='')
sp_data=read.csv(file,sep=' ',header = F)
colnames(sp_data)=c('ind','Population.ID','group','sp')
#####nj
#tree <- groupClade(tree, .node=c(295,342,333,251,203))
nj_tree=paste('nj/all_',sp_name,'.nwk',sep='')
tree=read.tree(nj_tree)
tree=split(dna_pca$ind, dna_pca$group_new) %>% groupOTU(tree, .)
nj_p=ggtree(tree, aes(color=group),layout = "circular",branch.length='none')+
  scale_color_manual(values=c(lzq_colors,'grey'))+
  #geom_tiplab(size=1)
  theme(
    legend.position="none",
    panel.grid= element_blank()
  )
ggsave(nj_p,file=paste('./FIGURE/',sp_name,'_nj.pdf',sep=''),width=4,height=4)

####cv to difine best K
cv=read.csv(paste("./admixture/all_",sp_name,'_individuals.CV.txt',sep=''),header=FALSE)
cv$K=1:dim(cv)[1]
cv_p=ggplot(cv,aes(x = K, y = V1)) + geom_line() + geom_point(size=3)+
  theme_bw()
ggsave(cv_p,filename=paste('./FIGURE/',sp_name,'_cv.pdf',sep=''),height=3,width=4)
######barplot for admixture
bestk=2
lzq_colors1=brewer.pal(bestk, 'Paired')
alist <- readQ(files=list.files(path=paste('./admixture/all_',sp_name,'_individuals',sep=''), full.names=T),indlabfromfile=F)
inds <- read.delim(paste('./admixture/all_',sp_name,'_individuals.individuals',sep=''),header=FALSE,stringsAsFactors=F)
ind_order_from_nj=tree$tip.label
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,function(x) {x%>%mutate(ind=inds$V1,name=inds$V1)%>%column_to_rownames('name')})
alist<-lapply(alist,function(x) {x[ind_order_from_nj, ]%>% dplyr::select(-ind)})
plotQ(alignK(alist[1:bestk]),imgoutput="join",barsize=1,height=2.6,width=10, 
      splab=paste0("K=",sapply(alist[1:bestk],ncol)),showlegend=T,clustercol = lzq_colors1,
      showindlab=F,useindlab=T,indlabsize=1,indlabspacer=-1,sharedindlab=T,
      barbordersize=0,outputfilename=paste(sp_name,'.admixture_new',sep=''),imgtype="pdf",exportpath="./FIGURE")

######map
pie_plot<-alist[[bestk]]%>%
  rownames_to_column('ind')%>%
  merge(.,sp_data,by='ind')%>%
  merge(.,pop)%>%
  dplyr::select(-ind,-group,-sp)%>%
  group_by(.,Population.ID)%>%
  summarize_each(funs(mean))
temp_n=bestk+1
map_k=ggplot()+geom_sf(data=world,fill='transparent',position = "identity")+
  scale_x_continuous(labels = ~ .x) +
  scale_y_continuous(labels = ~ .x) +
  coord_sf(ylim = c(35,43),xlim = c(-110,-95))+
  geom_scatterpie(aes(x=Longitude, y=Latitude,r=0.2), data=pie_plot,
                  cols=colnames(pie_plot)[2:temp_n],alpha=0.9,size=0.15)+
  scale_fill_manual(values=lzq_colors1)+
  theme_bw()+
  theme(
    legend.position="none",
    panel.grid= element_blank()
  )
ggsave(map_k,file=paste('./FIGURE/',sp_name,'_map.pdf',sep=''),width=8,height=4)
################

##########new map based on group 
pie=alist[[bestk]]%>%
  rownames_to_column('ind')%>%
  merge(.,sp_data,by='ind')%>%
  merge(.,pop)%>%
  dplyr::select(-ind,-group,-sp)%>%
  group_by(.,Population.ID)%>%
  summarize_each(funs(mean))
pie$group_new=group_new=apply(pie[,grep('Cluster',colnames(pie))],1,function(x) {if (max(x)>0.7){which.max(x)}else{return(99)}})
  
map_group=ggplot()+geom_sf(data=world,fill='transparent',position = "identity")+
  scale_x_continuous(labels = ~ .x) +
  scale_y_continuous(labels = ~ .x) +
  coord_sf(ylim = c(35,43),xlim = c(-110,-95))+
  geom_point(pie,mapping=aes(x=Longitude, y=Latitude,color=factor(group_new)),size=4)+
    scale_color_manual(values=c(lzq_colors1,'grey'))+
  theme_bw()+
  theme(
    legend.position="none",
    panel.grid= element_blank()
  )
 
ggsave(map_group,file=paste('./FIGURE/',sp_name,'_map_group.pdf',sep=''),width=8,height=4.5)
######################pca based on EVs
data=EV %>% 
  merge(.,pie[c('Population.ID','group_new')],by='Population.ID') %>%
  column_to_rownames('Population.ID')%>%
  na.omit()
pca=prcomp(x=data[c(12,19:53)],center=T,scale=TRUE)%>% summary()
pca_plot <- pca$x %>% data.frame()%>%
  mutate(ecotype=data$group_new)
plot_rotation<- pca$rotation%>% data.frame() %>%
  rownames_to_column('EVs')
pc1=data.frame(pca$importance)['Proportion of Variance','PC1']*100
pc2=data.frame(pca$importance)['Proportion of Variance','PC2']*100
evs_p=ggplot()+
  geom_segment(data=plot_rotation,mapping = aes(x=0,y=0,xend=PC1*10,yend=PC2*10),arrow=arrow(length = unit(0.2,"cm")))+
  geom_text(data=plot_rotation,aes(x=PC1*10,y=PC2*10,label=EVs),fontface='italic',vjust=-4,hjust=0.2,size=2.2,check_overlap=T,alpha=0.8,color='red')+
  geom_point(data=pca_plot,aes(x=PC1,y=PC2,color=factor(ecotype)),size=3,shape=19)+
 # geom_text(data=pca_plot,aes(x=PC1,y=PC2,label=row.names(pca_plot)),vjust=-1,hjust=0.2,size=3,color="black")+
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  scale_color_manual(values = c(lzq_colors,'grey'))+
  labs(x = paste('PC1(',pc1,'%)'),y =  paste('PC2(',pc2,'%)'),size = 12)+
  theme_bw()+  
  theme(
    legend.position="none",
    panel.grid= element_blank()
  )
ggsave(evs_p,file=paste('./FIGURE/',sp_name,'_env_pca.pdf',sep=''),width=4,height=3)
################dna pca
dna_file=paste('./pca/all_',sp_name,'_individuals.evec',sep='')
explain_por=paste('./pca/all_',sp_name,'_individuals.eval',sep='')
ex_por=read.csv(explain_por,header=F)
dna_pca<-read.table(dna_file,header = F,skip=1)%>% 
  mutate(ind=str_split(V1,':',simplify = T)[,1])%>%
  dplyr::select(ind,V2,V3)%>%
  rename(PC1=V2,PC2=V3)%>%
  merge(sp_data,.,by='ind')%>%
  merge(.,pie[c('Population.ID','group_new')],by='Population.ID')
gene_p=ggplot()+
  geom_point(data=dna_pca,aes(x=PC1,y=PC2,color=factor(group_new)),shape=19,size=2,alpha=0.7)+
  labs(x = paste('PC1(',ex_por[1,1],'%)'),y =  paste('PC2(',ex_por[2,1],'%)'),size = 12)+
  scale_color_manual(values = c(lzq_colors,'grey'))+
  theme_bw()+  
  theme(
    legend.position="none",
    panel.grid= element_blank()
  )
ggsave(gene_p,file=paste('./FIGURE/',sp_name,'_gene_pca.pdf',sep=''),width=4,height=3)
