# AtINO80_LightFuntion2020_source_code
The code is built for the paper named "The chromatin remodeling factor AtINO80 represses photomorphogenesis by modulating nucleosome density and H2A.Z incorporation in light-related genes"

Note: the bioinformatics scripts conclude shell and R scripts which are given in directories named "shell_scripts" and "R_scripts" respectively. R custom functions used in R scripts are listed in directory named "R_functions". Please source the related funcitons befor perform R script. In addition to, the provided scripts here may require additional steps and modifications to perform well on different platforms.

1.analysis of RNAseq data for atino80-5 mutant at white light
  1)copy fastq.gz file to the current directory
     cp /mnt/USB1/WL_atino80/*.gz ./
  2)trim adapters, map reads to genome, sort bam files, get bigwig files and align reads to genes
     sh RNAseq_paired_stranded.sh
  3)get differentially expressed genes comparing mutant with wild-type
     Rscript RNAseq_DEG_atino80_WL.r

2.overlap between AtINO80-regulated and light-induced genes for Figure 3A-3C. 
  1)use R to get venn diagram and p value for two sets of genes, analyze GO functions for overlap genes
source("OverlapPlot_function.r")# import function firstly
 
Plot("/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_INO80_5/","/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_regulate/PNAS/up_all.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_diff.xls","light_induced","atino80_diff","atino80-diff and light_induced")
  2)use R to get bubble figure for GO analysis result 
source("GOBubblePlot_function.r")# import function firstly
     
BubblePlot("/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_INO80_5/GO/","/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_INO80_5/GO/GO_choose_term.txt",c(0,50),c(1,18),c(1,3,6,9))
     
3.analysis of RNAseq data for atino80-6 mutant at dark
  1)copy fastq.gz file to the current directory
     cp /mnt/USB1/Dark_atino80/sequence_2/*.gz ./
  2)trim adapters, map reads to genome, sort bam files, get bigwig files and align reads to genes
     sh RNAseq_paired_unstranded.sh
  3)get differentially expressed genes comparing mutant with wild-type
     Rscript RNAseq_DEG_atino80_dark.r

4. overlap between AtINO80-regulated and dark-induced genes for Figure 3D-3F. 
  1)use R to get venn diagram and p value for two sets of genes, analyze GO functions for overlap genes
     source("OverlapPlot_function.r")# import function firstly
     Plot("/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_DarkINO80_5/","/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_regulate/PNAS/down_all.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/Dark_atino80_sequence2/DEG/15bei_pvalue/Darkino80_diff.xls","light_repressed","atino80_diff","atino80-diff and light_repressed")
  2)use R to get bubble figure for GO analysis result 
     source("GOBubblePlot_function.r")# import function firstly
     BubblePlot("/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_DarkINO80_5/GO/","/Volumes/OpheliaData/atino80/atino80_0317/light_regulate_genes/light_DarkINO80_5/GO/GO_choose_term.txt",c(0,20),c(1,14),c(0,3,6))
      
5.analysis of ChIPseq data for atino80-5 mutant at white light
  1)copy fastq.gz file to the current directory
     cp /mnt/USB1/ara_ino80_H2AZ_repeat2/*.gz ./
  2)map reads to genome, sort bam files, remove duplicates, get bigwig files and bed files
     sh ChIPseq_paired_nocut.sh
  3)get fragment size for Figure S8 
     sh FragmentSize.sh 
  4)find enriched regions (peaks) of histone variants (H2A.Z and H3)
     sh PeaksSicer.sh
  5)find overlap peaks between two replicates by R
     source("PeakOverlap_function.r")# import function firstly
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/overlappeaks/") 
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_WT1_H2AZ_repeat1.bed","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_WT1_H2AZ_repeat2.bed","WT1_repeat1","WT1_repeat2","WT1_H2AZ_repeats_overlap_peaks")
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_ino80_H2AZ_repeat1.bed","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_ino80_H2AZ_repeat2.bed","ino80_repeat1","ino80_repeat2","ino80_H2AZ_repeats_overlap_peaks")
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_WT1_H3_repeat1.bed","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_WT1_H3_repeat2.bed","WT1_repeat1","WT1_repeat2","WT1_H3_repeats_overlap_peaks")
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_ino80_H3_repeat1.bed","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/rawpeak/sig_ino80_H3_repeat2.bed","ino80_repeat1","ino80_repeat2","ino80_H3_repeats_overlap_peaks")
  6)remove the prefix "chr" of seqnames in overlap enrichment peaks    
     ls -1 *.bed |while read id
     do
     echo "${id%.bed}" 
     awk -F "hr" '{print  $2}' ${id} >  ${id%.bed}_nochr.bed
     done  
     wc -l *_nochr.bed
  7)get merged bam files, bigwig files and bed files
     sh merge.sh   

6.analysis of published ChIPseq data (GSE108450)
  1)download data from NCBI, map reads to genome, sort bam files, remove duplicates and get bigwig files and bed files
     sh ChIPseq_single_nocut.sh
  2)find enriched regions (peaks) of H2A.Z
     sh PeakSicer_GSE108450.sh
  3)find overlap peaks between two replicates by R
     source("PeakOverlap_function.r")# import function firstly
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/published_ChIPseq/ara/peaks/overlappeaks/") 
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/published_ChIPseq/ara/peaks/repeats/sig_SRR6412321.bed","/Volumes/OpheliaData/atino80/atino80_0317/published_ChIPseq/ara/peaks/repeats/sig_SRR6412322.bed","WT1_repeat1","WT1_repeat2","WT1_H2AZ_repeats_overlap_peaks")
  4)remove the prefix "chr" of seqnames in overlap enrichment peaks    
     ls -1 *.bed |while read id
     do
     echo "${id%.bed}" 
     awk -F "hr" '{print  $2}' ${id} >  ${id%.bed}_nochr.bed
     done  
     wc -l *_nochr.bed 

7. compare H2A.Z peaks in our ChIPseq data with published data for Figure S5C
     source("PeakOverlapTest_function.r")# import function firstly  
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/overlappeaks/compareWithPublised/") 
     Peakover("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/overlappeaks/WT1_H2AZ_repeats_overlap_peaks_nochr.bed","/Volumes/OpheliaData/atino80/atino80_0317/published_ChIPseq/ara/peaks/overlappeaks/WT1_H2AZ_repeats_overlap_peaks_nochr.bed","our","published","WTH2AZ_enriched_overlap_peaks_our_published")   
     
8.fragment size for Figure S5B
     sh FragmentSize.sh 
     
9.normolized reads density of histone variants(H2AZ and H3) across genes for Figure 4A-4B and 5A-5C
  1)get matrix of normolized reads density on gene body, around TSS and TES
     sh matrix.sh  
  2)get matrix for heatmap by R
     Rscript HeatmapMatrix.r
  3)plot heatmap
     sh HeatmapPlot.sh
  4)plot average density profile by R
     Rscript PlotProfile.r     
  
10.differentially enriched regions of histone modifcations (H2A.Z and H3) for Figure 4C-4E and Figure S6.
  1)get fold change of histone variants(H2AZ and H3) on enrichment regions (peaks) in mutant and wild-type
     sh PeaksSicerDiff.sh
  2)find differentially enriched regions by R
     source("DiffPeak_function.r")# import function firstly  
     DiffPeak1("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/","H2AZ","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/FCH2AZ_merge_nochr.summary",1.2,0.05)
     DiffPeak1("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/","H3","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/FCH3_merge_nochr.summary",1.2,0.05)   
  3)annotate differentially enriched regions by R
     source("PeakAnnotate_function.r")# import function firstly 
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/annotation/")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/H2AZ_sigPeakdown.xls","H2AZ_downpeak")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/H2AZ_sigPeakup.xls","H2AZ_uppeak")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/H3_sigPeakdown.xls","H3_downpeak")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/H3_sigPeakup.xls","H3_uppeak")     
  4)plot venn diagram and perform GO analysis for the overlap between H2A.Z, H3 diff-peaks and AtINO80-regulated genes
     source("VennPlot1_function.r")# import function firstly 
     VennPlot1("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/venn_histiondiff_with_expressiondiff/atino80_diff/","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/annotation/H2AZ_diffpeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_diff.xls","H2AZ_diff","atino80-diff","H2AZdiff_and_atino80diff")
     VennPlot1("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/venn_histiondiff_with_expressiondiff/atino80_diff/","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/annotation/H3_diffpeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_diff.xls","H3_diff","atino80-diff","H3diff_and_atino80diff")
  5)plot venn diagram for the overlap between H2A.Z diff-peaks and H3 diff-peaks 
     source("VennPlot2_function.r")# import function firstly 
     VennPlot2("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/venn_histiondiff_with_expressiondiff/H2AZdiff_H3diff/","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/annotation/H2AZ_diffpeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/sicer_diff/merge/annotation/H3_diffpeak_annoOvgene_genelist.xls","H2AZ_diff","H3_diff","H2AZdiff_and_H3diff")
  
11.differentially enriched regions of H2A.Z/H3 for Figure 4F
  1)find differentially enriched regions of H2A.Z considering H3 as input
     Rscript DiffBind_H2AZvsH3.r 
  2)plot MAplot by R
     Rscript MAplot.r  
  3)annotate differentially enriched regions
     source("PeakAnnotate_function.r")# import function firstly 
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/annotation/")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/H2AZ_H3_Down_nochr.xls","H2AZ:H3_downpeak")
     peakanno("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/H2AZ_H3_Up_nochr.xls","H2AZ:H3_uppeak")

12.enrichment analysis of H2AZ/H3 differentially enriched regions on DEGs in atino80-5 for Figure 4G and 4H
     source("EnrichmentAnalyze_function.r")# import function firstly 
     setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/relation_with_diffregulatedgenes/")
     EnrichAna("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/annotation/H2AZ:H3_downpeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_up.xls","downH2AZ:H3","up-regulated genes","H2AZ:H3down_genes_with_upregulated")
     EnrichAna("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/annotation/H2AZ:H3_downpeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_down.xls","downH2AZ:H3","down-regulated genes","H2AZ:H3down_genes_with_downregulated")
     EnrichAna("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/annotation/H2AZ:H3_uppeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_up.xls","upH2AZ:H3","up-regulated genes","H2AZ:H3up_genes_with_upregulated")
     EnrichAna("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/peaks_sicer2/diffbind/overlappeak/annotation/H2AZ:H3_uppeak_annoOvgene_genelist.xls","/Volumes/OpheliaData/atino80/atino80_0317/RNAseq/WL_atino80/DEG/15bei_pvalue/WLino80_down.xls","upH2AZ:H3","down-regulated genes","H2AZ:H3up_genes_with_downregulated")



