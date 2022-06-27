#args <- commandArgs()
#print(args)
#datasets_path <- args[6]
#print(datasets_path)

#("Y:/Melissa-abacus/Barbara_CPP/Bulk-RNA-Seq")
setwd("Y:/NOBACKUP/Melissa-nobackup/Bulk_RNA-seq/JPA_GRCh38")
library(dplyr)
library(tidyr)
library(ggplot2)
require(data.table) ## 1.9.4+

#datasets = read.delim(file= "20200214_Pipeline_input.tsv")
datasets = read.delim(file= "GRCh38_RNA_Pipeline_input_JPAs.tsv")
#datasets = datasets %>% filter(Sample=="N10316_S10")

load.calls <- function(x){
  # load STAR-Fusion ####
  print(paste0(x$Sample, " - STAR-Fusion"))
  star.fusion = read.delim(file=paste0("star-fusion/",x$Sample,"/star-fusion.fusion_predictions.abridged.coding_effect.tsv"))
  colnames(star.fusion)[1]="FusionName"
  star.fusion = star.fusion %>% separate(col = LeftGene, into = c("LeftGene", "LeftGene_ID"), sep = "\\^") %>% 
    separate(col = RightGene, into = c("RightGene", "RightGene_ID"), sep = "\\^") %>% 
    separate(col = LeftBreakpoint, into = c("LeftCHROM", "LeftBreak","LeftOrientation"), sep = "\\:") %>% 
    separate(col = RightBreakpoint, into = c("RightCHROM", "RightBreak","RightOrientation"), sep = "\\:") %>%
    mutate(ID=paste0(x$Sample,"_star.fusion_",rownames(star.fusion)), Confidence= NA, Caller="STAR.Fusion")
  
  # load Arriba ####
  print(paste0(x$Sample, " - Arriba"))
  arriba = read.delim(file=paste0("arriba/",x$Sample,"/",x$Sample,".fusions.tsv"))
  colnames(arriba)[1]="gene1"
  arriba = arriba %>% separate(col = breakpoint1, into = c("Chrom1", "Breakpoint1"), sep = "\\:") %>% 
    separate(col = breakpoint2, into = c("Chrom2", "Breakpoint2"), sep = "\\:") %>% 
    mutate(ID=paste0(x$Sample,"_Arriba_",rownames(arriba)), Chrom1=paste0("chr",Chrom1), Chrom2=paste0("chr",Chrom2),
           FusionName=paste0(gene1,"--",gene2), split_reads=as.numeric(split_reads1)+as.numeric(split_reads2), Caller= "Arriba")
  
  # load InFusion ####
  print(paste0(x$Sample, " - InFusion"))
  infusion =  read.delim(file=paste0("InFusion/",x$Sample,"/gene-gene.fusions.detailed.tsv")) #pre-filtered for only gene-gene fusions
  colnames(infusion)[1]="ID"
  infusion["break_on_exon"] = gsub("yes", "ExonBreak", infusion$break_on_exon)
  infusion = infusion %>% mutate(ID=paste0(x$Sample,"_InFusion_",ID), ref1=paste0("chr",ref1), ref2=paste0("chr",ref2), FusionName=paste0(gene_1,"--",gene_2))
  infusion_stringent =  read.delim(file=paste0("InFusion/",x$Sample,"_stringent/fusions.detailed.txt ")) #pre-filtered for only gene-gene fusions
  infusion_stringent = infusion_stringent %>% mutate(ID=paste0("InFusion_",X.id))
  infusion_stringent = unique(infusion_stringent$ID)
  infusion = infusion %>% mutate(stringent= ifelse(ID %in% infusion_stringent, "stringent", NA), Caller="InFusion")
  
  # Give files common headers ####
  #FusionName= , ID= , chr1=  ,break1= , gene1=, strand1= ,chr2= , break2= , gene2= ,strand2=  , 
  # SplitReads= , DiscordantMates= , FusionType= , Feature1= , Feature2= , Confidence=
  common.star.fusion = star.fusion %>% select(c(FusionName=FusionName,  Caller=Caller, ID=ID, chr1=LeftCHROM, break1=LeftBreak, gene1=LeftGene,
                                             strand1=LeftOrientation, chr2=RightCHROM, break2=RightBreak, gene2=RightGene,
                                             strand2=RightOrientation, SplitReads=JunctionReadCount,
                                             DiscordantMates=SpanningFragCount, FusionType=PROT_FUSION_TYPE,
                                             Feature1=CDS_LEFT_ID, Feature2=CDS_RIGHT_ID, Confidence=Confidence))
  common.arriba = arriba %>% select(c(FusionName=FusionName, Caller=Caller, ID=ID, chr1=Chrom1, break1=Breakpoint1, gene1=gene1,
                                    strand1=strand1.gene.fusion., chr2=Chrom2, break2=Breakpoint2, gene2=gene2,
                                    strand2=strand2.gene.fusion., SplitReads=split_reads,
                                    DiscordantMates=discordant_mates, FusionType=type, Feature1=site1,
                                    Feature2=site2 , Confidence=confidence))
  common.infusion = infusion %>% select(c(FusionName=FusionName, Caller=Caller, ID=ID, chr1=ref1, break1=break_pos1, gene1=gene_1,
                                        strand1=gene_1_strand, chr2=ref2, break2=break_pos2, gene2=gene_2,
                                        strand2=gene_2_strand, SplitReads=num_split, DiscordantMates=num_paired, 
                                        FusionType=break_on_exon, Feature1=feature_1, Feature2=feature_2, Confidence=stringent))
  masterlist = as.data.frame(rbind(common.star.fusion, common.arriba, common.infusion))
  masterlist = masterlist %>% mutate(patient=x$Sample, Inter.intra=ifelse(chr1==chr2,"Intra", "Inter")) %>% 
    select(Patient=patient, Caller, FusionName, everything())
  return(masterlist)
  }



#this dataset does not contain normal tissue
AllSamples_AllCalls=NULL
for(i in 1:nrow(datasets)) {
  data = datasets[i,]
  print(data$Sample)
  # put blood and tumor file in the same place
  dir.create("Fusion_filtering_Output", showWarnings = FALSE) #stops warnings if folder already exists
  dir.create(paste0("Fusion_filtering_Output/",data$Sample), showWarnings = FALSE) #stops warnings if folder already exists
  #Call load.calls function ####
  masterlist = load.calls(data)
  fwrite(masterlist, paste0("Fusion_filtering_Output/",data$Sample,"/",data$Sample,"_AllCallers.csv"))
  
  AllSamples_AllCalls = as.data.frame(rbind(AllSamples_AllCalls,masterlist))
  # Sometimes two caller detect the same event but gene order is different 
  # ex. AllSamples_AllCalls %>% filter(gene1=="MDM4" | gene2=="MDM4"))
  AllGenes = as.data.frame(c(unique(AllSamples_AllCalls$gene1),unique(AllSamples_AllCalls$gene2)))
  AllGenes = AllGenes %>% mutate(Duplicated= duplicated(AllGenes[,1])) %>% filter(Duplicated==T)
  AllGenes = as.vector(AllGenes$`c(unique(AllSamples_AllCalls$gene1), unique(AllSamples_AllCalls$gene2))`)
  t = AllSamples_AllCalls %>% filter(gene1 %in% AllGenes | gene2%in% AllGenes)
  t2= AllSamples_AllCalls %>% group_by(FusionName, Caller) %>% summarize(count=n()) %>% 
    spread(key=Caller, value = count)
  t2["Total.Count"]=rowSums(t2[,2:4],na.rm = T)
  t2 = t2 %>% mutate(called.by.Arriba = ifelse(!is.na(Arriba),as.numeric(1), 0),
                     called.by.InFusion = ifelse(!is.na(InFusion),as.numeric(1), 0),
                     called.by.STAR.Fusion = ifelse(!is.na('STAR-Fusion'),as.numeric(1), 0),
                     called.by.N.callers = (as.numeric(called.by.Arriba) + as.numeric(called.by.InFusion) + as.numeric(called.by.STAR.Fusion)))
  t2[is.na(t2)] = 0
}

fwrite(AllSamples_AllCalls, paste0("Fusion_filtering_Output/AllSamples_AllCallers.noFiltering.csv"))
BRAF_fusion = AllSamples_AllCalls %>% filter(gene1 == "BRAF" | gene2 == "BRAF") %>% filter(FusionName != "BRAF--BRAF")
fwrite(BRAF_fusion, paste0("Fusion_filtering_Output/AllSamples_AllCallers.BRAFfusions.csv"))
