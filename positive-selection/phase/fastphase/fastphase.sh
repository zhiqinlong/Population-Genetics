workdir="/SundanceScratch/zhiqin/phd_work/selective_sweeps"
for chr in Ha412HOChr{01..17}
do
#       python vcf2fastphase.py /TEMP_SCRATCH/zhiqin/vcf/${chr}.vqsr90.pass.highmappability.biallelic.four_sps.recode.vcf.gz ./vcf_impute_miss/${chr} 
echo "/SundanceScratch/zhiqin/software/fastPHASE $workdir/ihs/vcf_impute_miss/${chr} -H-4 -K8 -o./fastphase_out/${chr}.fastphase"
done
