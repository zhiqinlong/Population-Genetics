squash_axis <- function(from, to, factor) { 
  trans <- function(x) {    
    isq <- x > from & x < to
    ito <- x >= to
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)}
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}


combine_chrs<-function (dir_path){
  out_sweep=data.frame()
  for (i in 1:17){
    if (i < 10){i=paste('0',i,sep='')}
    path=paste(dir_path,'Ha412HOChr',i,'.sweep.output.2kb',sep='')
    sweep=read.table(path,header=T)
    sweep$CHR=paste('Chr',i,sep = '')
    out_sweep=rbind(out_sweep,sweep)
  }
  threshold=quantile(out_sweep$LR,0.95)
  out_sweep= out_sweep %>%  mutate(sig=case_when(LR>threshold  ~ 'sig',
                                                 TRUE~ 'non_sig'),
                                   location=round(location,0))
  chr_len <- out_sweep%>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(location))
  chr_pos <- chr_len  %>%
    mutate(total = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len)
  window_pos <- chr_pos %>%
    left_join(out_sweep, ., by="CHR") %>%
    arrange(CHR, location) %>%
    mutate( BPcum = location + total)%>%
    mutate(end=location+2000)%>%
    select(CHR,location,end,LR,alpha,BPcum,sig)
  X_axis <-  window_pos  %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  group=strsplit(dir_path,'/')
  out_path=paste('./sweepfinder2/sweep_out_sig/',group[[1]][4],'_sig_sweepout.txt',sep='')
  out_sweep_sig=window_pos %>% filter(sig=='sig')%>%
    select(CHR,location,end,LR,alpha,BPcum,sig)
  write.table(out_sweep_sig, out_path,quote=FALSE,row.names =FALSE,sep="\t",col.names =FALSE)
  
  return(list(window_pos,X_axis,threshold,chr_pos))
}
##########inversion positions


####GEA location
inv_gea_loc<-function(inv_path,gea_path,species,bpcum,x_text){
  inversions=read.csv(inv_path)
  inversion_sp=inversions%>%
    filter(Species==species)%>%
    mutate(ystart=0,yend=1,START=START*1000000,END=END*1000000)%>%
    mutate(CHR=case_when(CHR<10~paste('Chr0',CHR,sep=''),
                         CHR>=10~paste('Chr',CHR,sep='')))
  out_path=gsub(' ','',paste('./Marco_nature/inv_',species,'.txt',sep=''))
  write.table(inversion_sp%>%select(CHR,START,END,Haploblock_ID),out_path,quote=FALSE,row.names =FALSE,sep="\t",col.names =FALSE)
  inversion_sp_new<-bpcum %>%
    left_join(inversion_sp, ., by="CHR") %>%
    arrange(CHR, START) %>%
    mutate( start_cum = START + total,
            end_cum = END + total)
  inv_plot<-ggplot(bpcum,aes(x=total))+
    geom_rect(inversion_sp_new,mapping=aes(xmin = start_cum,xmax = end_cum , ymin = ystart , ymax = yend ),color='darkred',fill='darkblue',size=0.1) +
    scale_x_continuous(label = x_text$CHR, breaks= x_text$center )+
    theme_bw()+facet_grid(Species~., scales = "free_y")+
    theme(
      legend.position="none",
      axis.title.x =element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y= element_blank(),
      panel.grid= element_blank()
    )
  GEA=read.table(gea_path,sep=',',skip=1,header=T,na.strings='',fill=T)
  GEA_new=GEA%>% filter(Species==species) %>% 
    mutate(across(Chromosome,~gsub('Ha412HO','',.)))%>%
    rename(CHR=Chromosome)%>%
    select(c(Species,CHR,Gene.start,Gene.end))%>%
    left_join(.,bpcum, by="CHR")%>%
    arrange(CHR, Gene.start) %>%
    mutate( start_cum = Gene.start + total,
            end_cum = Gene.end + total,
            ystart=0,yend=1)
  out_path_gea=gsub(' ','',paste('./Marco_nature/gea_',species,'.txt',sep=''))
  write.table(GEA_new%>%select(c(CHR,Gene.start,Gene.end)),out_path_gea,quote=FALSE,row.names =FALSE,sep="\t",col.names =FALSE)
  gea_plot=ggplot(bpcum,aes(x=total))+
    geom_rect(GEA_new,mapping=aes(xmin = start_cum,xmax = end_cum , ymin = ystart , ymax = yend),color='darkred',fill='darkblue',size=0.1) +
    scale_x_continuous( label = x_text$CHR, breaks= x_text$center )+
    theme_bw()+facet_grid(Species~., scales = "free_y")+
    theme(
      legend.position="none",
      axis.title.x =element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y= element_blank(),
      panel.grid= element_blank()
    )
  return(list(inv_plot,gea_plot))
}
