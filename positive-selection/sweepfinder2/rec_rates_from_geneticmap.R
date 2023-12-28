library(splines)
library(data.table)
library(dplyr)


genetic_map=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt',sep='\t')
files=list.files("/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/",pattern='.*vcf',full.names = F)
for (file in files){
  chrnum=strsplit(file,split='\\.')[[1]][1]
  genetic_map_chr=subset(genetic_map,chr==chrnum)
  fitted_spline<-smooth.spline(genetic_map_chr$pos,genetic_map_chr$cM)

  plot(fitted_spline)
  vcf_path=paste('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/',file,sep='')
  vcf=read.table(vcf_path,header=F)
  vcf=vcf[c(1,2)]
  predicted_rates<-predict(fitted_spline,vcf$V2,deriv = 1)
  predicted_rates=as.data.frame(predicted_rates)%>% 
    mutate( rec_rates= ifelse(y < 0, 0, y))%>%
    rename(position=x)%>% 
    select(rec_rates,position)

  write.table(predicted_rates,file=paste('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/sweepfinder2/recombination_rates/',chr,'_recombination_rates',sep=''),row.names=F)
}
?write.table

a=lm(genetic_map_chr$pos~genetic_map_chr$cM)
plot(a)
rates<-predict(a,x=vcf$V2)
?predict()
