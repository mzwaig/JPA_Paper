args <- commandArgs()
print(args)
data_path <- args[6]
#subtypes <- as.vector(unique(as.data.frame(strsplit(args[7], ",")))[,1])
#print(subtypes)
cluster = T
#subtypes = "LGG"

setwd("Y:/")
cluster = F
library(dplyr)
library(tidyr)
#library(plyr) #causes issues with svim.bnd
library(ggplot2)
require(data.table) ## 1.9.4+
library(gridExtra)
library(grid)
library(gsubfn)
library(limma)
library(openxlsx)
library(tibble)
library(stringr)
library(tidyselect)
library(vcfR)
options(digits=10)
path = ifelse(cluster == T, "/lb/project/ioannisr/", "Y:/")
#library(wordcloud2)

#conda ####
#conda install -c r r-stringi r-dplyr r-tidyr r-ggplot2 r-data.table r-gridextra
#conda install -c bioconda grid bioconductor-limma
#conda install -c conda-forge r-gsubfn
#####

#datasets = read.delim(file= data_path, sep=",",header=T, na.strings=c("","NA")) #use this to open files
datasets = read.delim(file= paste0(path,"Melissa-abacus/Pipelines_MZ/GRCh38_JPA_datasets.csv"), sep=",",header=T, na.strings=c("","NA")) #use this to open files
datasets= datasets %>% filter(!is.na(ONT) | !is.na(WGS) | !is.na(PacBio))
#datasets= datasets %>% filter(Patient=="EPT2061")
#datasets= datasets %>% filter(Patient=="MDT-AP-2673")
#datasets= datasets %>% filter(diagnosis %in% subtypes)
datasets= datasets %>% filter(class == "Tumor") #%>% filter(Patient == "MDT-AP-2673")
print(paste0(datasets$Patient,"_", datasets$class))

#datasets=datasets[8,]
#datasets=datasets[23:38,]
list.of.samples = datasets %>% mutate(Patient=paste0(Patient,".",class)) 
list.of.samples = list.of.samples[order(list.of.samples$Patient),]
list.of.samples = as.vector(unique(list.of.samples$Patient))
list.callers= c("LongRanger","SvABA-10X","LinkedSV","NAIBR","GROC_SV", "svim", "CuteSV","NanoVar", "delly","lumpy", "wham","SvABA-WGS","manta","CNVkit" ,"sniffles")
#x=datasets[1,]
#data=x

# Parameters ####
#chromosomes = paste0("chr", c(1:22,"X","Y"))
chromosomes = c(1:22,"X","Y")
v.small=10
small=30
medium=100
large=1000

v.small.kb= paste(v.small, "kb", sep= "")
small.kb= paste(small, "kb", sep= "")
medium.kb= paste(medium, "kb", sep= "")
large.kb=paste(large, "kb", sep= "")

v.small= as.numeric(paste(v.small, "000", sep= ""))
small= as.numeric(paste(small, "000", sep= ""))
medium=as.numeric(paste(medium, "000", sep= ""))
large=as.numeric(paste(large, "000", sep= ""))

# Both breakpoints of SV must fall within
interval= as.numeric(1000)
#interval= as.numeric(paste(margin, "000", sep= ""))
#####

# File check ####
for(i in 1:nrow(datasets)) {
  x = datasets[i,]
  print(x$Patient)
  file_delly = list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".delly.*.noheader.tsv"))
  delly = paste0(path, x$WGS,"/",file_delly)
  file_lumpy <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".lumpy.*.noheader.tsv"))
  lumpy = paste0(path, x$WGS,"/",file_lumpy)
  file_wham <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".wham.*.noheader.tsv"))
  wham = paste0(path, x$WGS,"/",file_wham)
  file_SvABA_WGS <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".svaba.*.noheader.tsv"))
  SvABA_WGS = paste0(path, x$WGS,"/",file_SvABA_WGS)
  file_manta <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".manta.*.noheader.tsv"))
  manta = paste0(path, x$WGS,"/",file_manta)
  file_cnvkit <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".cnvkit.*.noheader.tsv"))
  cnvkit = paste0(path, x$WGS,"/",file_cnvkit)
  check = rbind.data.frame(cbind.data.frame(caller = "LongRanger", Exists = file.exists(paste0(x$LongRanger))),
                           cbind.data.frame(caller = "LongRanger DELs", Exists = file.exists(paste0(x$file_path,"_custom_files/", x$LR_name, "_LR_dels_parsed.vcf"))),
                           cbind.data.frame(caller = "SvABA-10X", Exists = file.exists(paste0(x$SvABA))),
                           cbind.data.frame(caller = "LinkedSV", Exists = file.exists(paste0(x$LinkedSV, "LSV.filtered_large_svcalls.bedpe"))),
                           cbind.data.frame(caller = "NAIBR", Exists = file.exists(paste0(x$NAIBR))),
                           cbind.data.frame(caller = "GROC-SV", Exists = file.exists(paste0(path,"NOBACKUP/Melissa-nobackup/10X_linked-reads/GROC_SV/",x$LR_name,"/classes.txt"))),
                           cbind.data.frame(caller = "SVIM", Exists = file.exists(paste0(path, x$ONT, "svim/", x$ONT_name, "/", x$ONT_name, ".final_results.noheader.vcf"))),
                           cbind.data.frame(caller = "CuteSV", Exists = file.exists(paste0(path, x$ONT, "CuteSV/", x$ONT_name, "/", x$ONT_name, ".CuteSV.vcf"))),
                           cbind.data.frame(caller = "NanoVar", Exists = file.exists(paste0(path, x$ONT, "NanoVar/", x$ONT_name, "/", x$ONT_name, ".sorted.nanovar.PASS.vcf"))),
                           cbind.data.frame(caller = "Delly", Exists = file.exists(paste0(delly))),
                           cbind.data.frame(caller = "Lumpy", Exists = file.exists(paste0(lumpy))),
                           cbind.data.frame(caller = "Wham", Exists = file.exists(paste0(wham))),
                           cbind.data.frame(caller = "SvABA-WGS", Exists = file.exists(paste0(SvABA_WGS))),
                           cbind.data.frame(caller = "Manta", Exists = ifelse(x$PairedWGS=="Yes",file.exists(paste0(manta)), TRUE)),
                           cbind.data.frame(caller = "CNVkit", Exists = file.exists(paste0(cnvkit))),
                           cbind.data.frame(caller = "Sniffles", Exists = file.exists(paste0(path,x$PacBio_name, "/", x$PacBio ,".vcf"))))
  print(check)
  ifelse(all(check$Exists) == TRUE, "", print(paste0(x$Patient, " is missing files")))
  #stopifnot(all(check$Exists) == TRUE)
}




# Load SV call files (A.1) ####
load.calls.by.sample <- function(x){
  masterlist=NULL
  class = ifelse(x$class=="Tumor", "somatic", "germline")
  # Load LongRanger - bedpe ####
  print(paste(x$LR_name, x$class," - Load LongRanger",sep=" "))
  LR.bedpe = fread(file= paste0(x$LongRanger), stringsAsFactors = FALSE, check.names=T)
  colnames(LR.bedpe)[1:12]= c("chrom1",	"start1",	"stop1",	"chrom2",	"start2",	"stop2",	"name",	"qual",	"strand1",	"strand2",	"filters",	"info")
  LR.bedpe = LR.bedpe%>% filter(!(str_detect(chrom1, "#")))
  LR.bedpe= LR.bedpe %>% mutate(caller="LongRanger", class=x$class, Technology="10X")
  LR.bedpe["COPY"] = LR.bedpe$info %>% strapplyc("COPY=(.*)", simplify = TRUE) %>% substr(1,1)
  LR.bedpe["COPY"] = gsub("c", NA, LR.bedpe$COPY)
  LR.bedpe["SVTYPE"] = ifelse(grepl("TYPE=DEL", LR.bedpe$info), "DEL", 
                               ifelse(grepl("TYPE=DUP", LR.bedpe$info), "DUP",
                                      ifelse(grepl("TYPE=INV", LR.bedpe$info), "INV",
                                             ifelse(grepl("TYPE=UNK", LR.bedpe$info), "UNK",
                                                    ifelse(grepl("TYPE=DISTAL", LR.bedpe$info), "DISTAL", "???")))))
   # select common headers
  LR.bedpe = LR.bedpe %>%
    select(CHROM=chrom1, START=start1, CHR2=chrom2 , END=start2 , ID=name, QUAL=qual, FILTER=filters, 
           SVTYPE=SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO=info) %>% 
    mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
  masterlist=as.data.frame(rbind(masterlist,LR.bedpe)) 
  # Load LongRanger - dels ####
  print(paste(x$LR_name, x$class," - Load LongRanger, dels",sep=" "))
  LR.dels.all = fread(file= paste0(x$file_path,"_custom_files/", x$LR_name, "_LR_dels_parsed.vcf"), stringsAsFactors = FALSE, check.names=T)
  colnames(LR.dels.all)[15] = "GT"
  LR.dels = LR.dels.all %>% filter(SVTYPE=="DEL") %>% 
    mutate(pair_ID=gsub("call", "del", ID), CHR2=CHROM, caller="LongRanger", class=x$class, Technology="10X", COPY= NA,
           INFO = paste0("HAP_ALLELIC_FRAC=",HAP_ALLELIC_FRAC,"ALLELIC_FRAC=",ALLELIC_FRAC, "PS=", PS)) %>% 
    select(CHROM=CHROM, START=POS, CHR2=CHR2, END , ID=pair_ID, QUAL, FILTER, 
           SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO) %>% 
    mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
           SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
  masterlist=as.data.frame(rbind(masterlist,LR.dels)) 
  # BNDs
  print(paste(x$LR_name,x$class," - Load LongRanger, BNDs",sep=" "))
  LR.BNDs = LR.dels.all %>% filter(SVTYPE=="BND") %>% mutate(pair_ID=gsub("_[12]", "", ID)) %>% mutate(pair_ID=gsub("call", "del",pair_ID))
  df1 = LR.BNDs[grep("_1",LR.BNDs$ID),]
  df1 = df1 %>% select(-END, -ID, -REF, -ALT, PS1= PS, -GT)
  colnames(df1)[1:2] = c("CHROM", "START")
  df2 = LR.BNDs[grep("_2",LR.BNDs$ID),]
  df2 = df2 %>% select(-END, - ID, -REF, -ALT, PS2=PS, -GT)
  colnames(df2)[1:2] = c("CHR2", "END")
  LR.BNDs2 = merge(df1, df2, by = c("pair_ID","QUAL", "FILTER","SVTYPE","SVTYPE2","SVLEN", "HAP_ALLELIC_FRAC","ALLELIC_FRAC"))
  LR.BNDs2 = LR.BNDs2  %>%  mutate(caller="LongRanger", CHR2=CHROM, class=x$class, Technology="10X", COPY= NA,
           INFO = paste0("HAP_ALLELIC_FRAC=",HAP_ALLELIC_FRAC,"ALLELIC_FRAC=",ALLELIC_FRAC, "PS1=", PS1, "PS2=", PS2)) %>% 
    select(CHROM=CHROM, START=START, CHR2=CHR2, END , ID=pair_ID, QUAL, FILTER, 
           SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO) %>% 
    mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
           SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
  masterlist=as.data.frame(rbind(masterlist,LR.BNDs2)) 
  # Load 10X SvABA ####
  if (!is.na(x$SvABA)) {
    print(paste(x$LR_name,x$class," - Load SvABA",sep=" "))
    SvABA = fread(file= paste0(x$SvABA), stringsAsFactors = FALSE, check.names=T)
    colnames(SvABA)[1:12]= c("CHROM","POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO", "FORMAT", paste0(x$Patient,".Tumor"), paste0(x$Patient,".Blood"), "SUPPORT")
    SvABA = SvABA %>% mutate(pair_ID=as.numeric(gsub(":[12]", "", ID)))
    df1 = SvABA[grep(":1",SvABA$ID),]
    colnames(df1)[1:2] = c("CHROM", "POS")
    df2 = SvABA[grep(":2",SvABA$ID),]
    colnames(df2)[1:2] = c("CHR2", "END")
    SvABA.parsed = merge(df1, df2, by = c("pair_ID","QUAL","FILTER","FORMAT",paste(x$tumor), paste(x$blood), "SUPPORT"))
    SvABA.parsed = SvABA.parsed %>% mutate(caller="SvABA-10X", class=x$class, Technology="10X",SVTYPE=NA, COPY=NA, ID=paste("SvABA_",pair_ID,sep=""))
    #SvABA.parsed["INFO"] = SvABA.parsed$INFO.x %>% strapplyc(";(.*)", simplify = TRUE) # drop supporting barcodes
    SvABA.parsed = SvABA.parsed %>% mutate(INFO=strapplyc(SvABA.parsed$INFO.x, ";(.*)", simplify = TRUE))
    # select common headers
    SvABA.parsed = SvABA.parsed %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID,QUAL=QUAL,  FILTER=FILTER, 
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO) %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    masterlist=as.data.frame(rbind(masterlist,SvABA.parsed))
    }
  # Load LinkedSV tumor ####
  if (!is.na(x$LinkedSV)) {
    print(paste(x$LR_name,x$class," - Load LinkedSV",sep=" "))
    linkedSV.sv = fread(file= paste0(x$LinkedSV, "LSV.filtered_large_svcalls.bedpe"), header=T, fill= T, stringsAsFactors = FALSE, na.strings=c("","NA"))
    colnames(linkedSV.sv) = c("chrom1","start1", "stop1",	"chrom2",	"start2", "stop2",	"SVTYPE",	"ID", "SVLEN", "QUAL", "FILTER", "INFO")
    if (nrow(linkedSV.sv)>0) {
      linkedSV.sv = linkedSV.sv  %>% drop_na(chrom1) %>% mutate(caller = "LinkedSV", class =  as.character(x$class),
                                                                Technology = "10X",COPY = NA)
      } else {
        linkedSV.sv = data.frame(chrom1=character(0), start1=character(0), stop1=character(0),
                                 chrom2=character(0),	start2=character(0), stop2=character(0),
                                 SVTYPE=character(0),	ID=character(0), SVLEN=character(0),
                                 QUAL=character(0), FILTER=character(0), INFO=character(0),
                                 caller=character(0), class=character(0), Technology=character(0), COPY=character(0))
        }
    #linkedSV.cnv = fread(file= paste0(x$LinkedSV, "phased_possorted_bam.bam.large_cnv.bedpe"), header=F,stringsAsFactors = FALSE)
    linkedSV.cnv = fread(file= paste0(x$LinkedSV, "LSV.large_cnv.bedpe"), header=F,stringsAsFactors = FALSE)
    if (nrow(linkedSV.cnv)>0) {
      colnames(linkedSV.cnv) = c("chrom1","start1", "stop1",	"chrom2",	"start2", "stop2",	"SVTYPE",	"ID", "PASS", "SVLEN", "QUAL", "FILTER", "INFO")
      linkedSV.cnv = linkedSV.cnv %>% mutate(caller="LinkedSV", class=x$class, Technology="10X") %>% select(-PASS)
      #linkedSV.cnv["COPY"] = linkedSV.cnv$INFO %>% strapplyc("COPY_NUMBER=(.*)", simplify = TRUE) %>% substr(1,1)
      linkedSV.cnv = linkedSV.cnv %>% mutate(COPY=strapplyc(linkedSV.cnv$INFO, "COPY_NUMBER=(.*)", simplify = TRUE))
      }
      linkedSV =as.data.frame(rbind(linkedSV.sv, linkedSV.cnv))
      linkedSV["ID"] = paste0("LSV_",rownames(linkedSV))# select common headers
      linkedSV = linkedSV %>% 
        select(CHROM=chrom1, START=start1, CHR2=chrom2 , END=start2 , ID=ID,QUAL=QUAL,FILTER=FILTER, 
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO) %>% 
        mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
               SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
      masterlist=as.data.frame(rbind(masterlist,linkedSV))
      }
  # Load NAIBR tumor ####
  if (!is.na(x$NAIBR)) {
    print(paste(x$LR_name,x$class," - Load NAIBR",sep=" "))
    NAIBR = fread(file=paste0(x$NAIBR), sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names=T)
    NAIBR = NAIBR %>% mutate(SVTYPE=NA, caller="NAIBR", class=x$class, Technology="10X", COPY=NA,ID=paste("NAIBR_",rownames(NAIBR),sep=""),INFO=NA)
    # select common headers
    NAIBR = NAIBR %>% 
      select(CHROM=Chr1, START=Break1, CHR2=Chr2 , END=Break2 , ID=ID,QUAL=Score, FILTER=Pass.filter,
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO) %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    masterlist=as.data.frame(rbind(masterlist,NAIBR))
    }
  # Load GROC-SV ####
  if (!is.na(x$GROC_SV)) {
    print(paste(x$LR_name, " ", x$class," - Load GROC_SV",sep=""))
    grocsv.classes = fread(file=paste0(path,"NOBACKUP/Melissa-nobackup/10X_linked-reads/GROC_SV/",x$LR_name,"/classes.txt"), sep="\t", stringsAsFactors = FALSE, check.names=T, header=TRUE)
    grocsv.classes$classes <- gsub('11', 'Blood', grocsv.classes$classes)
    grocsv.classes$classes <- gsub('10', 'Blood', grocsv.classes$classes)
    grocsv.classes$classes <- gsub('1', 'Tumor', grocsv.classes$classes)
    grocsv.genotypes = fread(file= paste0(path,"NOBACKUP/Melissa-nobackup/10X_linked-reads/GROC_SV/",x$LR_name,"/genotypes.tsv"), sep="\t", stringsAsFactors = FALSE, check.names=T, header=TRUE)
    grocsv.full = data.frame(chromx = grocsv.genotypes$chromx, chromy= grocsv.genotypes$chromy, cluster = grocsv.genotypes$cluster,
                             x= grocsv.genotypes$x, y= grocsv.genotypes$y, class = grocsv.classes$classes, QUAL=grocsv.genotypes$blacklist , 
                             SVLEN=grocsv.genotypes$dist, FILTER=grocsv.genotypes$quality)
    #grocsv.full$cluster <- sub("^", "cluster_", grocsv.full$cluster) # cluster is not unique since GROC-sv can detect complex rearrangements
    grocsv.full = grocsv.full %>% mutate(caller="GROC_SV", SVTYPE=NA, COPY=NA,Technology="10X",ID=paste0("GROCSV_",rownames(grocsv.full)),INFO=NA)
    # If SVLEN comes back negative, SV is Intrachromosomal and the order of the events is reversed 
    grocsv.full["Inter.intra"] = ifelse(grocsv.full$chromx==grocsv.full$chromy, "Intra", "Inter")
    grocsv.full["SVLEN"] = ifelse(grocsv.full$Inter.intra=="Intra", grocsv.full$y - grocsv.full$x, "")
    grocsv.full = grocsv.full %>% mutate(CHROM=ifelse(SVLEN<0,as.character(chromy),as.character(chromx)),
                                         CHR2=ifelse(SVLEN<0,as.character(chromx),as.character(chromy)), 
                                         START=ifelse(SVLEN<0,y,x), END=ifelse(SVLEN<0,x,y))
    grocsv.full["SVLEN2"] = ifelse(grocsv.full$Inter.intra=="Intra", grocsv.full$END - grocsv.full$START, "") #check that all length are pos
    # select common headers
    grocsv.full = grocsv.full %>% 
      select(CHROM=CHROM , START=START , CHR2=CHR2 , END=END , ID=ID,  QUAL=QUAL, FILTER=FILTER, 
             SVTYPE= SVTYPE, caller= caller, class= class , COPY= COPY, Technology=Technology, INFO=INFO) %>%
      filter(class == as.character(x$class)) %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    print(paste0("GROC-SV classes :", unique(grocsv.full$class)))
    masterlist=as.data.frame(rbind(masterlist,grocsv.full))
    }
  
  # Load ONT - SVIM ####
  if (!is.na(x$ONT)) {
  print(paste0("Load SVIM calls - ", x$ONT_name))
  #svim = fread(paste0("/lb/project/ioannisr/", x$ONT, "final_results.noheader.vcf"), stringsAsFactors = FALSE, check.names=T)
  svim = fread(paste0(path, x$ONT, "svim/", x$ONT_name, "/",x$ONT_name, ".final_results.noheader.vcf"), stringsAsFactors = FALSE)
  colnames(svim)=c("#CHROM","POS","ID", "REF", "ALT", "QUAL","FILTER","INFO","FORMAT","Sample")
  svim = svim %>% mutate(SVTYPE = ifelse(grepl("SVTYPE=DEL", svim$INFO), "DEL", 
                                         ifelse(grepl("SVTYPE=INS", svim$INFO), "INS",
                                                ifelse(grepl("SVTYPE=INV", svim$INFO), "INV",
                                                       ifelse(grepl("SVTYPE=DUP:TANDEM", svim$INFO), "DUP:TANDEM",
                                                              ifelse(grepl("\\bSVTYPE=DUP_INT;CUTPASTE\\b", svim$INFO), "DUP_INT:CUTPASTE",
                                                                     ifelse(grepl("\\bSVTYPE=DUP_INT\\b", svim$INFO), "DUP_INT",
                                                                            ifelse(grepl("SVTYPE=BND", svim$INFO), "BND", "???"))))))))
  svim.separated = NULL
  # split DELs
  svim.separated = svim %>% filter(SVTYPE=="DEL" | SVTYPE=="DUP:TANDEM" | SVTYPE=="DUP_INT") %>% separate(INFO, c("TYPE", "END","SVLEN", "SUPPORT", "STD_SPAN", "STD_POS", "READS"), ";") %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = gsub('SVLEN=', '', SVLEN),
           caller= "svim",  Technology="ONT_DNA", CHR2=`#CHROM` ,
           INFO=paste("SVLEN=",SVLEN,SUPPORT,STD_SPAN,STD_POS, sep=";"), Inter.intra="Intra") %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
           SVTYPE=SVTYPE, caller=caller , Technology=Technology, INFO=INFO)
  # split INSs
  svim.ins = svim %>% filter(SVTYPE=="INS")%>% separate(INFO, c("TYPE", "END","SVLEN", "SUPPORT", "STD_SPAN", "STD_POS", "SEQS" ,"READS"), ";") %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = gsub('SVLEN=', '', SVLEN),
           END= as.numeric(as.character(END)) + as.numeric(as.character(SVLEN)),
           caller= "svim",  Technology="ONT_DNA", CHR2=`#CHROM` ,
           INFO=paste("SVLEN=",SVLEN,SUPPORT,STD_SPAN,STD_POS, sep=";"), Inter.intra="Intra") %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
           SVTYPE=SVTYPE, caller=caller , Technology=Technology, INFO=INFO)
  # split INVs
  svim.inv = svim %>% filter(SVTYPE=="INV") %>% separate(INFO, c("TYPE", "END", "SUPPORT", "STD_SPAN", "STD_POS","READS"), ";") %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = as.numeric(END)-as.numeric(POS),
           caller= "svim",  Technology="ONT_DNA", CHR2=`#CHROM` ,
           INFO=paste(SUPPORT,STD_SPAN,STD_POS, sep=";"), Inter.intra="Intra") %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
           SVTYPE=SVTYPE, caller=caller , Technology=Technology, INFO=INFO)
  # split DUP_INT:CUTPASTE
  svim.int = svim %>% filter(SVTYPE=="DUP_INT:CUTPASTE") %>% separate(INFO, c("TYPE", "CUTPASTE", "END", "SVLEN", "SUPPORT", "STD_SPAN", "STD_POS","READS"), ";") %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = gsub('SVLEN=', '', SVLEN),
           caller= "svim",  Technology="ONT_DNA", CHR2=`#CHROM` ,
           INFO=paste("SVLEN=",SVLEN,SUPPORT,STD_SPAN,STD_POS, sep=";"), Inter.intra="Intra") %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
           SVTYPE=SVTYPE, caller=caller , Technology=Technology, INFO=INFO)
  # split BNDs
  svim.bnd = svim %>% filter(SVTYPE=="BND") %>% separate(INFO, c("TYPE", "SUPPORT", "STD_POS1", "STD_POS2","READS"), ";") %>% 
    separate(ALT, c("CHR2", "END"), sep = '[[:]]', remove = F) %>% 
    mutate(CHR2 = gsub(as.character("N"), "", CHR2), CHR2 = gsub(as.character("\\]*"), "", CHR2), CHR2 = gsub(as.character("\\[*"), "", CHR2),
           END = gsub(as.character("N"), "", END), END = gsub(as.character("\\]*"), "", END), END = gsub(as.character("\\[*"), "", END)) %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), caller= "svim", Technology="ONT_DNA", 
           INFO=paste(SUPPORT,STD_POS1,STD_POS2, sep=";"), 
           Inter.intra  = ifelse(as.numeric(as.character(`#CHROM`)) == as.numeric(as.character(CHR2)), "Intra", "Inter"),
           SVLEN = ifelse(Inter.intra  == "Intra",as.numeric(END)-as.numeric(POS), NA )) %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
           SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra )
  # still duplicated! order breakpoints for smallest to big and then deduplicate
  svim.BND2 = svim.bnd %>% mutate(CHROM.NUM=ifelse(CHROM=="X", 23, ifelse(CHROM=="Y", 24, as.numeric(svim.bnd$CHROM))), 
                                  CHR2.NUM=ifelse(CHR2=="X", 23,ifelse(CHR2=="Y", 24, as.numeric(svim.bnd$CHR2))))
  svim.BND2["SWAP"] = as.numeric(svim.BND2$CHR2.NUM - svim.BND2$CHROM.NUM) # need to do this for inter and Intra since the line requires a values
  svim.BND2 = svim.BND2 %>% filter(!is.na(as.character(SWAP))) %>% # removes non-chromosomal events
    mutate(CHROM.fix=ifelse(SWAP<0,CHR2,CHROM),CHR2.fix=ifelse(SWAP<0,CHROM,CHR2),
           START.fix=ifelse(SWAP<0,END,START), END.fix=ifelse(SWAP<0,START,END),
           CHROM.NUM.fix=ifelse(SWAP<0,CHR2.NUM,CHROM.NUM),CHR2.NUM.fix=ifelse(SWAP<0,CHROM.NUM,CHR2.NUM)) # Swap these too to re-test
  svim.BND2["SWAP2"] = as.numeric(svim.BND2$CHR2.NUM.fix - svim.BND2$CHROM.NUM.fix)
  svim.BND2 = svim.BND2 %>% 
    select(CHROM=CHROM.fix , START=START.fix , CHR2=CHR2.fix , END=END.fix,  QUAL=QUAL,
           FILTER=FILTER, SVTYPE= SVTYPE, SVLEN=SVLEN, caller= caller,
           Technology=Technology, INFO=INFO, Inter.intra ) %>% distinct() 
  svim.BND2 = svim.BND2 %>% mutate(ID= paste0("svimBND_", rownames(svim.BND2))) %>% 
    select(CHROM, START, CHR2, END, ID, QUAL=QUAL, FILTER=FILTER,
           SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra ) 
  
  # still duplicated! order breakpoints for smallest to big and then deduplicate
  svim.BND3 = svim.BND2
  svim.BND3["SWAP2"] = ifelse(svim.BND3$SVLEN <= 0, "Yes", NA) # make swap 2 for Intra events to order breakpoints
  svim.BND3 = svim.BND3 %>% 
    mutate(START.fix=ifelse(is.na(SWAP2),START,
                            ifelse(SWAP2=="Yes", END, "?")),
           END.fix=ifelse(is.na(SWAP2),END,
                          ifelse(SWAP2=="Yes",START, "?")),
           ABS.SVLEN =abs(SVLEN)) %>%
    select(CHROM=CHROM, START=START.fix , CHR2=CHR2 , END=END.fix,  QUAL=QUAL,
           FILTER=FILTER, SVTYPE= SVTYPE, SVLEN=SVLEN, caller= caller,
           Technology=Technology, INFO=INFO, ABS.SVLEN, Inter.intra ) %>% distinct(CHROM, START, CHR2, END, .keep_all = TRUE) 
  svim.BND3 = svim.BND3 %>% mutate(ID= paste0("svimBND_", rownames(svim.BND3)),
                                   SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END)-as.numeric(START), "")) %>% 
    select(CHROM, START, CHR2, END, ID, QUAL=QUAL, FILTER=FILTER, #ABS.SVLEN,
           SVTYPE=SVTYPE, caller=caller , Technology=Technology, INFO=INFO)
  
  # put together parsed vcf
  svim.separated = rbind.data.frame(svim.separated, svim.ins, svim.inv, svim.int, svim.BND3)
  svim.separated = svim.separated %>% mutate(COPY=NA,  class=x$class) %>% 
    select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO) %>% 
    mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
           SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
  masterlist=as.data.frame(rbind(masterlist,svim.separated))
  }
  
  # Load ONT - CuteSV ####
  if (!is.na(x$ONT)) {
    print("Load CuteSV calls - DEL, DUP, INS, INV, BND" ) 
    y = paste0(path, x$ONT, "CuteSV/", x$ONT_name, "/",x$ONT_name, ".CuteSV.vcf")
    cuteSV = read.vcfR(file = y)
    cuteSV = as.data.frame(cbind(cuteSV@fix,cuteSV@gt))
    #colnames(cuteSV)=c("#CHROM","POS","ID", "REF", "ALT", "QUAL","FILTER","INFO","FORMAT",sample)
    cuteSV = cuteSV %>% mutate(SVTYPE = ifelse(grepl("SVTYPE=DEL", cuteSV$INFO), "DEL",
                                               ifelse(grepl("SVTYPE=DUP", cuteSV$INFO), "DUP",
                                                      ifelse(grepl("SVTYPE=INS", cuteSV$INFO), "INS",
                                                             ifelse(grepl("SVTYPE=INV", cuteSV$INFO), "INV",
                                                                    ifelse(grepl("SVTYPE=BND", cuteSV$INFO), "BND", "???"))))))
    # split DEL, DUP, INS & INV
    cuteSV.separated = cuteSV %>% filter(SVTYPE=="DEL" | SVTYPE=="DUP" | SVTYPE=="INS" | SVTYPE=="INV") %>% mutate(INFO2 = INFO) %>%
      separate(INFO2, c("PRECISION", "SVTYPE", "SVLEN",	"END"), ";") %>%
      mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = gsub('SVLEN=', '', SVLEN),
             caller= "cuteSV",  Technology="ONT_DNA", CHR2=CHROM, Inter.intra="Intra") %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller, Technology=Technology, INFO=INFO, Inter.intra=Inter.intra)
    
    # split BNDs, already ordered by position and chromosome, each on one line
    cuteSV.bnd = cuteSV %>% filter(SVTYPE=="BND") %>% mutate(INFO2 = INFO) %>% separate(INFO2, c("PRECISION", "SVTYPE"), ";") %>%
      separate(ALT, c("CHR2", "END"), sep = '[[:]]', remove = F) %>% 
      mutate(CHR2 = gsub(as.character("N"), "", CHR2), CHR2 = gsub(as.character("\\]*"), "", CHR2), CHR2 = gsub(as.character("\\[*"), "", CHR2),
             END = gsub(as.character("N"), "", END), END = gsub(as.character("\\]*"), "", END), END = gsub(as.character("\\[*"), "", END)) %>% #Marco PFA aligned to GRCh38!!! need to remove chr then add it again to fix orientation
      mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), caller= "cuteSV",  Technology="ONT_DNA", 
             Inter.intra  = ifelse(CHROM == CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra  == "Intra",as.numeric(END)-as.numeric(POS), NA )) %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra)
    
    # put together parsed vcf
    cuteSV.separated = rbind.data.frame(cuteSV.separated, cuteSV.bnd)
    cuteSV.separated = cuteSV.separated %>%  mutate(COPY=NA,  class=x$class) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO,SVLEN,Inter.intra)
    masterlist=as.data.frame(rbind(masterlist,cuteSV.separated))
  }
  
  # Load ONT - NanoVar ####
  if (!is.na(x$ONT)) {
    print("Load NanoVar calls - DEL, DUP, INS, INV, BND" ) 
    y = paste0(path, x$ONT, "NanoVar/", x$ONT_name, "/",x$ONT_name, ".sorted.nanovar.PASS.vcf")
    nanoVar = read.vcfR(file = y)
    nanoVar = as.data.frame(cbind(nanoVar@fix,nanoVar@gt))
    #colnames(nanoVar)=c("#CHROM","POS","ID", "REF", "ALT", "QUAL","FILTER","INFO","FORMAT",sample)
    nanoVar = nanoVar %>% mutate(SVTYPE = ifelse(grepl("SVTYPE=DEL", nanoVar$INFO), "DEL",
                                                 ifelse(grepl("SVTYPE=DUP", nanoVar$INFO), "DUP",
                                                        ifelse(grepl("SVTYPE=INS", nanoVar$INFO), "INS",
                                                               ifelse(grepl("SVTYPE=INV", nanoVar$INFO), "INV",
                                                                      ifelse(grepl("SVTYPE=BND", nanoVar$INFO), "BND", "???"))))))
    # split DEL, DUP & INV
    nanoVar.separated = nanoVar %>% filter(SVTYPE=="DEL" | SVTYPE=="DUP" | SVTYPE=="INV") %>% 
      mutate(INFO2 = INFO) %>% separate(INFO2, c("SVTYPE",	"END", "SVLEN"), ";") %>%
      mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), SVLEN = gsub('SVLEN=', '', SVLEN),
             caller= "nanoVar",  Technology="ONT_DNA", CHR2=CHROM, Inter.intra="Intra") %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra=Inter.intra)
    # split INSs
    nanoVar.ins = nanoVar %>% filter(SVTYPE=="INS") %>% mutate(INFO2 = INFO) %>% separate(INFO2, c("SVTYPE",	"END", "SVLEN"), ";") %>%
      mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END), 
             SVLEN = gsub('SVLEN=', '', SVLEN),  SVLEN = gsub('>', '', SVLEN), SVLEN = as.numeric(as.character(SVLEN)),
             END= as.numeric(as.character(END)) + as.numeric(as.character(SVLEN)),
             caller= "nanoVar",  Technology="ONT_DNA", CHR2=CHROM, Inter.intra="Intra") %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra=Inter.intra)
    
    # split BNDs, already ordered by position and chromosome, each on one line
    nanoVar.bnd = nanoVar %>% filter(SVTYPE=="BND") %>% mutate(INFO2 = INFO) %>% separate(INFO2, c("SVTYPE",	"END", "SVLEN"), ";") %>%
      separate(ALT, c("CHR2", "END"), sep = '[[:]]', remove = F) %>% 
      mutate(CHR2 = gsub(as.character("N"), "", CHR2), CHR2 = gsub(as.character("\\]*"), "", CHR2), CHR2 = gsub(as.character("\\[*"), "", CHR2),
             END = gsub(as.character("N"), "", END), END = gsub(as.character("\\]*"), "", END), END = gsub(as.character("\\[*"), "", END),
             CHROM = gsub('chr', '', `CHROM`), CHR2 = gsub('chr', '', CHR2)) %>% #need to remove chr then add it again to fix orientation
      mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), caller= "nanoVar",  Technology="ONT_DNA", 
             Inter.intra  = ifelse(CHROM == CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra  == "Intra",as.numeric(END)-as.numeric(POS), NA )) %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID, QUAL=QUAL, FILTER=FILTER, 
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra)
    # duplicated! order chromosomes smallest to big and then deduplicate
    nanoVar.BND2 = nanoVar.bnd %>% mutate(CHROM.NUM=ifelse(CHROM=="X", 23, ifelse(CHROM=="Y", 24, as.numeric(nanoVar.bnd$CHROM))), 
                                          CHR2.NUM=ifelse(CHR2=="X", 23,ifelse(CHR2=="Y", 24, as.numeric(nanoVar.bnd$CHR2))))
    nanoVar.BND2["SWAP"] = as.numeric(nanoVar.BND2$CHR2.NUM - nanoVar.BND2$CHROM.NUM) # need to do this for inter and Intra since the line requires a values
    nanoVar.BND2 = nanoVar.BND2 %>% filter(!is.na(as.character(SWAP))) %>% # removes non-chromosomal events
      mutate(CHROM.fix=ifelse(SWAP<0,CHR2,CHROM),
             CHR2.fix=ifelse(SWAP<0,CHROM,CHR2),
             START.fix=ifelse(SWAP<0,as.numeric(as.character(END)),as.numeric(as.character(START))), 
             END.fix=ifelse(SWAP<0,as.numeric(as.character(START)),as.numeric(as.character(END))),
             CHROM.NUM.fix=ifelse(SWAP<0,CHR2.NUM,CHROM.NUM),
             CHR2.NUM.fix=ifelse(SWAP<0,CHROM.NUM,CHR2.NUM)) # Swap these too to re-test
    nanoVar.BND2["SWAP2"] = as.numeric(nanoVar.BND2$CHR2.NUM.fix - nanoVar.BND2$CHROM.NUM.fix)
    nanoVar.BND2 = nanoVar.BND2 %>% 
      select(CHROM=CHROM.fix , START=START.fix , CHR2=CHR2.fix , END=END.fix,ID=ID,  QUAL=QUAL,
             FILTER=FILTER, SVTYPE= SVTYPE, SVLEN=SVLEN, caller= caller,
             Technology=Technology, INFO=INFO, Inter.intra=Inter.intra) %>% distinct() 
    nanoVar.BND2 = nanoVar.BND2 %>% mutate(CHROM = paste0("chr", CHROM), CHR2 = paste0("chr", CHR2)) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL=QUAL, FILTER=FILTER,
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller , Technology=Technology, INFO=INFO, Inter.intra=Inter.intra) 
    
    # still duplicated! order breakpoints for smallest to big and then deduplicate
    nanoVar.BND3 = nanoVar.BND2 %>% mutate(SVLEN = ifelse(Inter.intra=="Intra", END-START, NA))
    nanoVar.BND3["SWAP2"] = ifelse(nanoVar.BND3$SVLEN <= 0, "Yes", NA) # make swap 2 for Intra events to order breakpoints
    nanoVar.BND3 = nanoVar.BND3 %>% 
      mutate(START.fix=ifelse(is.na(SWAP2),START,
                              ifelse(SWAP2=="Yes", END, "?")),
             END.fix=ifelse(is.na(SWAP2),END,
                            ifelse(SWAP2=="Yes",START, "?")),
             ABS.SVLEN =abs(SVLEN)) %>%
      select(CHROM=CHROM, START=START.fix , CHR2=CHR2 , END=END.fix, ID=ID, QUAL=QUAL,
             FILTER=FILTER, SVTYPE= SVTYPE, caller= caller,
             Technology=Technology, INFO=INFO, Inter.intra=Inter.intra, ABS.SVLEN) %>% distinct(CHROM, START, CHR2, END, .keep_all = TRUE) 
    nanoVar.BND3 = nanoVar.BND3 %>% mutate(SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END)-as.numeric(START), "")) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL=QUAL, FILTER=FILTER, #ABS.SVLEN,
             SVTYPE=SVTYPE, SVLEN=SVLEN, caller=caller, Technology=Technology, INFO=INFO, Inter.intra=Inter.intra)
    
    # put together parsed vcf
    nanoVar.separated = rbind.data.frame(nanoVar.separated, nanoVar.ins, nanoVar.BND3)
    nanoVar.separated = nanoVar.separated %>% mutate(COPY=NA,  class=x$class) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO, Inter.intra, SVLEN)
    masterlist=as.data.frame(rbind(masterlist,nanoVar.separated))
    
  }
  # Load Delly calls ####
  if (!is.na(x$WGS)) {
  print("Load WGS delly calls")
  file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".delly.*.noheader.tsv"))
  delly = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
  delly = delly %>% separate(INFO, c("Precision", "SVTYPE", "SVMETHOD","CHR2", "END","PE", "MAPQ", "CT", "CIPOS", "CIEND"), ";") %>%
    mutate(SVTYPE = gsub('SVTYPE=', '', SVTYPE), CHR2 = gsub('CHR2=', '', CHR2), END = gsub('END=', '', END), ID = paste0("delly_",ID),
           caller= "delly", class=x$class, Technology="WGS", COPY=NA,
           INFO=paste(Precision, PE, MAPQ, CT, CIPOS, CIEND, sep=";")) %>% 
    select(CHROM=`#CHROM`, START=POS, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO) %>% 
    mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
           SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
  masterlist=as.data.frame(rbind(masterlist,delly))
  }
  # Load Lumpy calls ####
  if (!is.na(x$WGS)) {
    print("Load WGS lumpy calls")
    #lumpy = fread(paste0(path, x$WGS, "somatic_calls/", x$Patient, ".lumpy.",class,".noheader.tsv"), stringsAsFactors = FALSE)
    file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".lumpy.*.noheader.tsv"))
    lumpy = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
    lumpy = lumpy %>% mutate(INFO2 = INFO) %>% separate(INFO2, c("SVTYPE", "X2", "X3"), ";") %>%
      mutate(CHR2=`#CHROM`, SVTYPE = gsub('SVTYPE=', '', SVTYPE), ID = paste0("lumpy-",ID)) 
    lumpy2 = lumpy %>% filter(SVTYPE != "BND") %>% 
      mutate(END = gsub('END=', '', X3), caller= "lumpy", class=x$class, Technology="WGS", COPY=NA) %>% 
      select(CHROM=`#CHROM`, START=POS, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO)
    
    lumpy.bnd = lumpy %>% filter(SVTYPE == "BND") %>% mutate(MATEID = ID) %>%
      separate(MATEID, c("pair_ID"), "_") %>%
      select(`#CHROM`, POS, ID, QUAL, FILTER, SVTYPE, INFO, pair_ID)
    df1 = lumpy.bnd[grep("_1",lumpy.bnd$ID),]
    colnames(df1)[1:2] = c("CHROM", "START")
    df2 = lumpy.bnd[grep("_2",lumpy.bnd$ID),]
    colnames(df2)[1:2] = c("CHR2", "END")
    lumpy.bnd.parsed = merge(df1, df2, by = c("pair_ID","QUAL","FILTER", "SVTYPE"))
    lumpy.bnd.parsed = lumpy.bnd.parsed %>% 
      mutate(caller="lumpy", class=x$class, Technology="WGS", COPY=NA, ID=pair_ID) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO=INFO.x)
    lumpy2 = as.data.frame(rbind(lumpy.bnd.parsed, lumpy2))
    lumpy2= lumpy2 %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    
    masterlist=as.data.frame(rbind(masterlist,lumpy2))
  }
  # Load Wham calls ####
  if (!is.na(x$WGS)) {
    print("Load WGS wham calls")
    #wham = fread(paste0(path, x$WGS, "somatic_calls/", x$Patient, ".wham.",class,".noheader.tsv"), stringsAsFactors = FALSE)
    file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".wham.*.noheader.tsv"))
    #print(paste0(path, x$WGS,"/",file_list))
    wham = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
    wham = wham %>% separate(INFO, c("SVTYPE", "SVLEN","ID2", "SUPPORT", "MERGED", "REFINED", "END","POS2", "FIVE","THREE","LID", "RID","CIPOS", "CIEND","COLLAPSED"), ";") %>%
      mutate(CHR2=`#CHROM`,SVTYPE = gsub('SVTYPE=', '', SVTYPE), END = gsub('END=', '', END),
           caller= "wham", class=x$class, Technology="WGS", COPY=NA,
           INFO=paste(SVLEN,SUPPORT, MERGED,REFINED, CIPOS, CIEND,COLLAPSED, sep=";")) %>% 
      select(CHROM=`#CHROM`, START=POS, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO) %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    masterlist=as.data.frame(rbind(masterlist,wham))
  }
  # Load SvABA WGS calls ####
  if (!is.na(x$WGS)) {
    print(paste(x$Patient,x$class," - Load SvABA WGS",sep=" "))
    #SvABA_WGS =  fread(paste0(path, x$WGS, "somatic_calls/", x$Patient, ".svaba.",class,".noheader.tsv"), stringsAsFactors = FALSE)
    file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".svaba.*.noheader.tsv"))
    #print(paste0(path, x$WGS,"/",file_list))
    SvABA_WGS = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
    colnames(SvABA_WGS)[1:12]= c("CHROM","POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO", "FORMAT", paste0(x$Patient,".Tumor"), paste0(x$Patient,".Blood"), "SUPPORT")
    SvABA_WGS = SvABA_WGS %>% mutate(pair_ID=as.numeric(gsub(":[12]", "", ID)))
    df1 = SvABA_WGS[grep(":1",SvABA_WGS$ID),]
    colnames(df1)[1:2] = c("CHROM", "POS")
    df2 = SvABA_WGS[grep(":2",SvABA_WGS$ID),]
    colnames(df2)[1:2] = c("CHR2", "END")
    SvABA_WGS.parsed = merge(df1, df2, by = c("pair_ID","QUAL","FILTER","FORMAT",paste(x$tumor), paste(x$blood), "SUPPORT"))
    SvABA_WGS.parsed = SvABA_WGS.parsed %>% mutate(caller="SvABA-WGS", class=x$class, Technology="WGS",SVTYPE=NA, COPY=NA, ID=paste("SvABA_WGS_",pair_ID,sep=""))
    #SvABA_WGS.parsed["INFO"] = SvABA_WGS.parsed$INFO.x %>% strapplyc(";(.*)", simplify = TRUE) # drop supporting barcodes
    SvABA_WGS.parsed = SvABA_WGS.parsed %>% mutate(INFO=strapplyc(SvABA_WGS.parsed$INFO.x, ";(.*)", simplify = TRUE))
    # select common headers
    SvABA_WGS.parsed = SvABA_WGS.parsed %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID,QUAL=QUAL,  FILTER=FILTER, 
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO) %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    masterlist=as.data.frame(rbind(masterlist,SvABA_WGS.parsed))
  }
  
  # Load Manta calls ####
  if (!is.na(x$WGS) & x$PairedWGS=="Yes") {
    print("Load WGS manta calls")
    #manta = fread(paste0(path, x$WGS, "somatic_calls/", x$Patient, ".manta.",class,".noheader.tsv"), stringsAsFactors = FALSE)
    file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".manta.*.noheader.tsv"))
    manta = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
    manta = manta %>% mutate(INFO2 = INFO) %>% separate(INFO2, c("X1", "X2"), ";") %>%
      mutate(SVTYPE = ifelse(grepl("SVTYPE=DEL", manta$INFO), "DEL",
                             ifelse(grepl("SVTYPE=DUP", manta$INFO), "DUP",
                                    ifelse(grepl("SVTYPE=INS", manta$INFO), "INS",
                                           ifelse(grepl("SVTYPE=INV", manta$INFO), "INV",
                                                  ifelse(grepl("SVTYPE=BND", manta$INFO), "BND", "???")))))) # ADD TRANSLOCATIONS
    manta2 = manta %>% filter(SVTYPE != "BND") %>% 
      mutate(CHR2=`#CHROM`, END = gsub('END=', '', X1), caller= "manta", class=x$class, Technology="WGS", COPY=NA) %>% 
      select(CHROM=`#CHROM`, START=POS, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO)
    
    manta.bnd = manta %>% filter(SVTYPE == "BND") %>% mutate(MATEID = gsub('MATEID=', '', X2)) %>% 
      separate(MATEID, c("v1", "v2", "v3","v4","v5", "v6", "v7","v8"), ":") %>%
      mutate(pair_ID = paste(v1, v2, v3,v4,v5,sep= ":"), ID_fix = paste0(v1,":",v2,":",v3,":",v4,":",v5,"_",v8)) %>%
      select(`#CHROM`, POS, ID_fix, QUAL, FILTER, SVTYPE, INFO, pair_ID, ID)
    df1 = manta.bnd[grep("_0",manta.bnd$ID_fix),]
    colnames(df1)[1:2] = c("CHROM", "START")
    df2 = manta.bnd[grep("_1",manta.bnd$ID_fix),]
    colnames(df2)[1:2] = c("CHR2", "END")
    manta.bnd.parsed = merge(df1, df2, by = c("pair_ID","QUAL","FILTER", "SVTYPE"))
    manta.bnd.parsed = manta.bnd.parsed %>% 
      mutate(caller="manta", class=x$class, Technology="WGS",SVTYPE, COPY=NA, ID=pair_ID) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, caller,class,COPY,Technology, INFO=INFO.x)
    manta2 = as.data.frame(rbind(manta.bnd.parsed, manta2))
    manta2 = manta2 %>% 
      mutate(Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA))
    
    masterlist=as.data.frame(rbind(masterlist, manta2))
  }
  # Load CNVkit calls ####
  if (!is.na(x$WGS) & x$PairedWGS=="Yes") {
    print("Load WGS CNVkit calls")
    file_list <- list.files(path=paste0(path, x$WGS), pattern = paste0(x$WGS_name,".cnvkit.*.noheader.tsv"))
    cnvkit = fread(paste0(path, x$WGS,"/",file_list), stringsAsFactors = FALSE)
    cnvkit = cnvkit %>% 
      mutate(SVTYPE = ifelse(grepl("SVTYPE=DEL", cnvkit$INFO), "DEL",
                             ifelse(grepl("SVTYPE=DUP", cnvkit$INFO), "DUP", NA))) %>%
      separate(INFO, c("X1", "X2", "X3"), ";", remove = FALSE) 
    cnvkit = cnvkit %>% mutate(END = gsub("END=", "", cnvkit$X3)) %>% 
      mutate(CHROM=`#CHROM`, START=POS, CHR2=`#CHROM`, caller="cnvkit", class=x$class, Technology="WGS",COPY=NA, ID=paste0("cnvkit_", rownames(cnvkit)),
             Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA)) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, SVLEN, caller,class,COPY,Technology, INFO, Inter.intra)
    
    masterlist=as.data.frame(rbind(masterlist, cnvkit))
  }
  # Load Sniffles calls ####
  if (!is.na(x$PacBio_name)) {
    print("Load PacBio Sniffles calls")
    sniffles = fread(paste0(path, x$PacBio_name, "/", x$PacBio ,".noheader.vcf"), stringsAsFactors = FALSE)#[-1:-2,]
    colnames(sniffles) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample")
    sniffles2 = sniffles %>% 
      mutate(SVTYPE = gsub("*>", "", ALT), SVTYPE = gsub("*<", "", SVTYPE)) %>% 
      filter(SVTYPE == "DUP" | SVTYPE == "DEL" | SVTYPE == "INS" | SVTYPE == "INVDUP" | SVTYPE == "INV" | SVTYPE == "DUP/INS" | SVTYPE == "DEL/INV" ) %>%
       separate(INFO, c("X1", "X2", "X3", "X4"), ";", remove = FALSE)
    sniffles2 = sniffles2 %>% mutate(CHR2 = gsub("CHR2=", "", sniffles2$X3), END = gsub("END=", "", sniffles2$X4)) %>% 
      mutate(START=POS, caller="sniffles", class=x$class, Technology="PacBio",COPY=NA, ID=paste0("sniffles_", ID),
             Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA)) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, SVLEN, caller,class,COPY,Technology, INFO, Inter.intra)
    
    sniffles.bnd = sniffles %>% mutate(SVTYPE = gsub("*>", "", ALT), SVTYPE = gsub("*<", "", SVTYPE)) %>% 
      filter(SVTYPE != "DUP" & SVTYPE != "DEL" & SVTYPE != "INS" & SVTYPE != "INVDUP" & SVTYPE != "INV" & SVTYPE != "DUP/INS" & SVTYPE != "DEL/INV" ) %>%
      mutate(SVTYPE="BND") %>% separate(ALT, c("CHR2", "END"), sep = '[[:]]', remove = F) %>% 
      mutate(CHR2 = gsub(as.character("N"), "", CHR2), CHR2 = gsub(as.character("\\]*"), "", CHR2), CHR2 = gsub(as.character("\\[*"), "", CHR2),
             END = gsub(as.character("N"), "", END), END = gsub(as.character("\\]*"), "", END), END = gsub(as.character("\\[*"), "", END)) %>%
      filter(CHROM %in% chromosomes & CHR2 %in% chromosomes) %>%
      mutate(START=POS, caller="sniffles", class=x$class, Technology="PacBio",COPY=NA, ID=paste0("sniffles_", ID),
             Inter.intra=ifelse(CHROM==CHR2, "Intra", "Inter"),
             SVLEN = ifelse(Inter.intra=="Intra", as.numeric(END) - as.numeric(START), NA)) %>% 
      select(CHROM, START, CHR2, END, ID, QUAL, FILTER, SVTYPE, SVLEN, caller,class,COPY,Technology, INFO, Inter.intra)
      
    masterlist=as.data.frame(rbind(masterlist, sniffles2, sniffles.bnd))
  }
  #Make master list with common headers ####
  print(paste(x$Patient,x$class," - Make masterlist",sep=" "))
  masterlist = masterlist %>% mutate(caller.class=paste(caller,".",class, sep=""), Patient = x$Patient)
  #masterlist["Inter.intra"] = ifelse(masterlist$CHROM==masterlist$CHR2, "Intra", "Inter")
  masterlist["SVLEN"] = ifelse(masterlist$Inter.intra=="Intra", as.numeric(masterlist$END) - as.numeric(masterlist$START), NA)
  masterlist["ABS.SVLEN"] = ifelse(masterlist$Inter.intra=="Intra", abs(as.numeric(masterlist$END) - as.numeric(masterlist$START)), NA)
  masterlist["Size.class"] = ifelse(masterlist$Inter.intra=="Inter","Interchr",
                                    ifelse(masterlist$ABS.SVLEN<v.small,"Very Small",
                                    ifelse(masterlist$ABS.SVLEN>=v.small & masterlist$ABS.SVLEN< small,"Small",
                                         ifelse(masterlist$ABS.SVLEN>=small & masterlist$ABS.SVLEN< medium, "Medium",
                                                ifelse(masterlist$ABS.SVLEN>=medium & masterlist$ABS.SVLEN< large, "Large",
                                                       ifelse(masterlist$ABS.SVLEN>= large, "Very Large","???"))))))
  masterlist = masterlist %>% filter(FILTER=="PASS" | FILTER=="." | FILTER=="SEG_DUP") %>%  subset(CHROM %in% chromosomes) %>% subset(CHR2 %in% chromosomes)
  # some callers will call y chrom deletions in females, filter these calls in females
  if (x$sex=="F") {
    masterlist = masterlist %>% filter(CHROM!="Y" | CHR2!="Y")
    } else {
      masterlist = masterlist
      }
  orderlist = c("Very Small","Small","Medium","Large","Very Large","Interchr")
  masterlist <- transform(masterlist, Size.class = factor(Size.class, levels = orderlist))
  return(masterlist)
  }


# Overlap SV calls (A.2) ####
overlapping.calls.by.sample <- function(x,y){
  #x=data
  #y = masterlist
  # Orient interchromosomal SVs ####
  #   Process interchromosomal events since order of breakpoints may be different, want the smaller chr first
  #   (need to test then in two orientations since some might write them 1,10,etc), GROC-SV Intra delt with 
  y2 = y %>% mutate(CHROM.NUM=ifelse(CHROM=="X", 23,ifelse(CHROM=="Y", 24, as.numeric(y$CHROM))), 
                   CHR2.NUM=ifelse(CHR2=="X", 23,ifelse(CHR2=="Y", 24, as.numeric(y$CHR2))))
  y2["SWAP"] = as.numeric(y2$CHR2.NUM - y2$CHROM.NUM) # need to do this for inter and Intra since the line requires a values
  y2 = y2 %>% mutate(CHROM.fix=ifelse(SWAP<0,CHR2,CHROM),CHR2.fix=ifelse(SWAP<0,CHROM,CHR2),
                   START.fix=ifelse(SWAP<0,END,START), END.fix=ifelse(SWAP<0,START,END),
                   CHROM.NUM.fix=ifelse(SWAP<0,CHR2.NUM,CHROM.NUM),CHR2.NUM.fix=ifelse(SWAP<0,CHROM.NUM,CHR2.NUM)) # Swap these too to re-test
  y2["SWAP2"] = as.numeric(y2$CHR2.NUM.fix - y2$CHROM.NUM.fix)
  y2 = y2 %>% select(CHROM=CHROM.fix , START=START.fix , CHR2=CHR2.fix , END=END.fix , ID=ID,  QUAL=QUAL,
                   FILTER=FILTER, SVTYPE= SVTYPE, caller= caller, class= class , COPY= COPY, Technology=Technology, 
                   caller.class=caller.class, Patient= Patient, Inter.intra=Inter.intra, SVLEN=SVLEN, Size.class=Size.class, INFO=INFO)
  
  # Make interval for overlap ####
  #   Easier to breakup masterlist so all columns will have the same names 
  interval.small = 50
  masterlist.interval = y2 %>% filter(Size.class!="Very Small") %>% mutate(START.lower=as.numeric(START)-interval, START.upper=as.numeric(START)+interval, END.lower=as.numeric(END)-interval, 
         END.upper=as.numeric(END)+interval) #%>% select(c("CHROM","START.lower","START.upper","CHR2","END.lower", "END.upper","ID", "caller"))
  masterlist.interval.small = y2 %>% filter(Size.class=="Very Small") %>% mutate(START.lower=as.numeric(START)-interval.small, START.upper=as.numeric(START)+interval.small, END.lower=as.numeric(END)-interval.small, 
                                                                           END.upper=as.numeric(END)+interval.small)
  masterlist.interval = as.data.frame(rbind(masterlist.interval, masterlist.interval.small))
  # Compare each caller to the masterlist one at a time ####
  #   drop calls "overlapping" because they are from the same caller, merge results to master list
  #   will have issue that an SV called by LR and NAIBR will also be reported separately as NAIBR and LR, etc.
  #   can make unique ID if you use the ID columns from the merged files so information will be in the same order
  #   doing it this way since I want to know if the breakpoints overlap, NOT if the entire events overlap
  masterlist.overlap=masterlist.interval
  for(i in 1:length(list.callers)) {
    print(paste0("Overlap calls by ", list.callers[i]))
    df = masterlist.interval %>% filter(caller==list.callers[i]) %>% 
    select(c("CHROM","START.lower","START.upper","CHR2","END.lower", "END.upper","START", "END","ID","caller","QUAL","SVTYPE","COPY", "INFO"))
    
    setDT(masterlist.interval)[,START.lower]
    setDT(masterlist.interval)[,START.upper]
    setDT(df)[,START.lower]
    setDT(df)[,START.upper]
    setkey(masterlist.interval, CHROM, CHR2,START.lower, START.upper)
    start.overlap = foverlaps(df, masterlist.interval, by.x=c("CHROM", "CHR2", "START.lower", "START.upper"))
    start.filtered = start.overlap %>% filter(!is.na(START.lower)) %>% 
      select(-c("i.START", "i.END","i.QUAL","i.SVTYPE","i.COPY","i.INFO")) %>% 
      mutate(Pair_ID=ifelse(!is.na(ID),paste0(ID, ".", i.ID), "")) # need ID to merge start and end overlaps # %>% filter(caller!=i.caller), can't filter since you need these to be in unique_ID
    
    setDT(masterlist.interval)[,END.lower]
    setDT(masterlist.interval)[,END.upper]
    setDT(df)[,END.lower]
    setDT(df)[,END.upper]
    setkey(masterlist.interval, CHROM, CHR2, END.lower, END.upper)
    end.overlap = foverlaps(df, masterlist.interval, by.x=c("CHROM", "CHR2", "END.lower", "END.upper"))
    #end.overlap = end.overlap %>% mutate(Pair_ID=paste(ID, ".", i.ID, sep="")) # need ID to merge start and end overlaps
    end.filtered = end.overlap %>% filter(!is.na(START.lower)) %>% 
      select(-c("i.START", "i.END","i.QUAL","i.SVTYPE","i.COPY","i.INFO")) %>%
      mutate(Pair_ID=ifelse(!is.na(ID),paste0(ID, ".", i.ID), "")) # need ID to merge start and end overlaps # %>% filter(caller!=i.caller), can't filter since you need these to be in unique_ID
    
    colnames(start.filtered) = paste(colnames(start.filtered),".S",sep="")
    colnames(end.filtered) = paste(colnames(end.filtered),".E",sep="")
    IDs=merge(start.filtered, end.filtered, by.x="Pair_ID.S", by.y="Pair_ID.E", all=TRUE)
    IDs["both.ends.match"]= ifelse(!is.na(IDs$CHROM.S), ifelse(!is.na(IDs$CHROM.E), "Yes", "No"), "No")
    IDs.filt = IDs %>% filter(both.ends.match=="Yes") %>% select(-c(contains(".E",),contains("Pair_ID"), "i.START.lower.S","i.START.upper.S"))
    colnames(IDs.filt)=gsub(".S", "", colnames(IDs.filt)) # removes S from SVTYPE
    colnames(IDs.filt)[23:24]=c(paste(list.callers[i], ".ID",sep=""),paste(list.callers[i], ".caller",sep=""))
    
    colnames(df)=paste(list.callers[i], ".", colnames(df),sep="")
    IDs.filt.info = merge(IDs.filt, df, by=(c(paste(list.callers[i], ".ID",sep=""),paste(list.callers[i], ".caller",sep="")))) %>% 
      select("CHROM",colnames(masterlist.interval),everything()) %>%
      select(-c(paste0(list.callers[i],".START.upper"), paste0(list.callers[i],".START.lower"),
                                                                 paste0(list.callers[i],".END.upper"),
                                                                 paste0(list.callers[i],".END.lower")))
    masterlist.overlap = merge(masterlist.overlap, IDs.filt.info, by=colnames(masterlist.interval), all=TRUE, allow.cartesian=TRUE)
    }
  # Processing overlapping calls ####
  print("Make unique_ID")
  masterlist.overlap2 = masterlist.overlap %>% unite("unique_ID", paste0(list.callers,".ID"), sep = ".", remove = FALSE)
  masterlist.overlap2$unique_ID = gsub(".NA\\w(*SKIP)(*FAIL)|.NA", "", masterlist.overlap2$unique_ID, perl=TRUE)
  masterlist.overlap2$unique_ID = gsub("NA.", "", masterlist.overlap2$unique_ID, fixed = TRUE)
  print("Count callers which detect SV")
  masterlist_unique_calls= masterlist.overlap2 %>% 
    select(-c("START","END","ID", "QUAL","FILTER", "SVTYPE","caller","COPY", "caller.class","SVLEN","INFO"))  %>%
    mutate(called.by.LongRanger = ifelse(!is.na(LongRanger.ID),as.numeric(1), 0), 
           `called.by.SvABA-10X` = ifelse(!is.na(`SvABA-10X.ID`),as.numeric(1), 0), 
           called.by.LinkedSV = ifelse(!is.na(LinkedSV.ID),as.numeric(1), 0), 
           called.by.NAIBR = ifelse(!is.na(NAIBR.ID),as.numeric(1), 0), 
           called.by.GROC_SV = ifelse(!is.na(GROC_SV.ID),as.numeric(1), 0),
           called.by.svim = ifelse(!is.na(svim.ID),as.numeric(1), 0),
           called.by.CuteSV = ifelse(!is.na(CuteSV.ID),as.numeric(1), 0),
           called.by.NanoVar = ifelse(!is.na(NanoVar.ID),as.numeric(1), 0),
           called.by.delly = ifelse(!is.na(delly.ID),as.numeric(1), 0),
           called.by.lumpy = ifelse(!is.na(lumpy.ID),as.numeric(1), 0),
           called.by.wham = ifelse(!is.na(wham.ID),as.numeric(1), 0),
           `called.by.SvABA-WGS` = ifelse(!is.na(`SvABA-WGS.ID`),as.numeric(1), 0),
           called.by.manta = ifelse(!is.na(manta.ID),as.numeric(1), 0),
           called.by.cnvkit = ifelse(!is.na(CNVkit.ID),as.numeric(1), 0),
           called.by.sniffles = ifelse(!is.na(sniffles.ID),as.numeric(1), 0),
           called.by.N.callers=(as.numeric(called.by.LongRanger) + 
                                  as.numeric(`called.by.SvABA-10X`) + 
                                  as.numeric(called.by.LinkedSV) + 
                                  as.numeric(called.by.NAIBR) + 
                                  as.numeric(called.by.GROC_SV) + 
                                  as.numeric(called.by.svim) + 
                                  as.numeric(called.by.CuteSV) + 
                                  as.numeric(called.by.NanoVar) + 
                                  as.numeric(called.by.delly) +
                                  as.numeric(called.by.lumpy) +
                                  as.numeric(called.by.wham) + 
                                  as.numeric(`called.by.SvABA-WGS`) + 
                                  as.numeric(called.by.manta) +
                                  as.numeric(called.by.cnvkit) +
                                  as.numeric(called.by.sniffles))) %>%
    select("unique_ID" , "Patient","class" ,"Inter.intra", "Size.class", CHROM = CHROM, CHR2= CHR2, contains("called"), everything())
  print("Remove Duplicates")
  masterlist_unique_calls = masterlist_unique_calls[!duplicated(masterlist_unique_calls$unique_ID),]
  masterlist_unique_calls = masterlist_unique_calls %>% filter(unique_ID!="NA")
  print("SVTYPE consensus")
  # only LongRanger and LinkedSV return SVTYPEs so only need to do one comparison 
  masterlist_unique_calls = masterlist_unique_calls %>% rowwise() %>% 
    mutate(SVTYPE = ifelse(grepl(!is.na(LongRanger.SVTYPE),!is.na(LinkedSV.SVTYPE))==TRUE, LongRanger.SVTYPE, 
                           paste0(LongRanger.SVTYPE,"/",LinkedSV.SVTYPE)))
  masterlist_unique_calls$SVTYPE = gsub("/NA", "", masterlist_unique_calls$SVTYPE, fixed=TRUE)
  masterlist_unique_calls$SVTYPE = gsub("NA/", "", masterlist_unique_calls$SVTYPE, fixed = TRUE)
  masterlist_unique_calls = masterlist_unique_calls %>% 
    select("unique_ID" , "Patient","class" , "Inter.intra", "Size.class", "CHROM", "CHR2", "SVTYPE",
           "called.by.N.callers", contains("called"), everything(), -"Technology")
  
  return(masterlist_unique_calls)
  }

# Loop to make masterlist of calls and overlapp within sample ####
for(i in 1:nrow(datasets)) {
  data = datasets[i,]
  print(data$Patient)
  # put blood and tumor file in the same place
  #dir.create(paste0(data$file_path,"_custom_files/Comparing_SVs"), showWarnings = FALSE) #stops warnings if folder already exists
  #Call load.calls function (A.1) ####
  masterlist = NULL
  masterlist = load.calls.by.sample(data)
  # output 1
  fwrite(masterlist, paste0(path,"Melissa-abacus/SV_database_pipeline/Cross_Technology_Consensus/", data$Patient, "_AllCalls_AllTechnologies.csv"))
  #masterlist = fread(file= paste0(path,"Melissa-abacus/SV_database_pipeline/Cross_Technology_Consensus/", data$Patient, "_AllCalls_AllTechnologies.csv"), stringsAsFactors = FALSE, check.names=T)
  #Call overlapping calls on masterlist generate by load.calls function (A.2) ####
  print("Overlapping calls")
  masterlist_unique_calls = NULL
  masterlist_unique_calls = overlapping.calls.by.sample(data,masterlist)
  masterlist_unique_calls = masterlist_unique_calls %>% 
    mutate(POS = paste0("chr",CHROM, ":",START.lower,"-",START.upper,";chr",CHR2, ":",END.lower,"-",END.upper)) %>%
    select(unique_ID, POS, everything())
  fwrite(masterlist_unique_calls, paste0(path,"Melissa-abacus/SV_database_pipeline/Cross_Technology_Consensus/", data$Patient, "_OverlappingCalls_AllCalls_AllTechnologies.csv"))
  #masterlist_unique_calls = fread(file= paste0(path,"Melissa-abacus/SV_database_pipeline/Cross_Technology_Consensus/", data$Patient, "_OverlappingCalls_AllCalls_AllTechnologies.csv"), stringsAsFactors = FALSE, check.names=T)
  
  #Summary tables ####
  #Summary.by.size.class = masterlist %>% group_by(class, caller, Size.class) %>% summarize(count=n())
  print("Make summary")
  Summary.by.size.class = masterlist %>% group_by(class, Technology, caller, Size.class) %>% dplyr::summarise(count = n())
  Summary.by.size.class =  Summary.by.size.class %>% spread(key=Size.class,value = count)
  Summary.by.size.class[is.na(Summary.by.size.class)] = 0
  Summary.by.size.class['Total.Calls'] = rowSums(Summary.by.size.class[4:9])
  Summary.by.size.class['More.Than.10kb'] = rowSums(Summary.by.size.class[5:8])
  # discrepencies with common. files are due to filtering of contigs
  Summary.by.sample.caller = masterlist %>% group_by(class, Technology, caller) %>% 
    dplyr::summarise(count=n(), Mean.SVLEN = trunc(mean(SVLEN, na.rm=TRUE)), 
              Median.SVLEN = trunc(median(SVLEN, na.rm=TRUE)), Min.SVLEN = min(SVLEN, na.rm=TRUE), Max.SVLEN = max(SVLEN, na.rm=TRUE))
  Summary.by.sample.caller[is.na(Summary.by.sample.caller)] = 0
  
  orderlist.chromosomes = c(1:22,"X","Y", "Inter")
  masterlist <- transform(masterlist, CHROM = factor(CHROM, levels = orderlist.chromosomes))
  orderlist = c("Small","Medium","Large","Very Large","Interchr")
  masterlist <- transform(masterlist, Size.class = factor(Size.class, levels = orderlist))
  masterlist <- transform(masterlist, CHROM = factor(caller, levels = list.callers))
  
  
  order.callers = c(1:max(unique(masterlist_unique_calls$called.by.N.callers)))
  masterlist_unique_calls <- transform(masterlist_unique_calls, called.by.N.callers = factor(called.by.N.callers, levels = order.callers))
  
  # Plots ####
  print(paste(data$Patient," - Plot",sep=""))
   p1 =  masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=caller, fill= caller)) + 
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + 
    #facet_grid(caller  ~ . ,scales = "free", space = "free_x") + 
    scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + ylab("Count (log10)") + 
    #facet_grid(caller ~ .) + scale_y_log10(labels = function(x) format(x, scientific = FALSE))  + ylab("Count (log10)") + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + theme_bw() + theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)

   p2 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=CHROM, fill= caller)) + 
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + 
    facet_grid(caller  ~ . ,scales = "free", space = "free_x") + 
    scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + ylab("Count (log10)") + 
    #scale_y_log10(labels = function(x) format(x, scientific = FALSE))  + ylab("Count (log10)") + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls across callers by chromosome", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + theme_bw() + theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)
   
  p3 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=caller, fill= caller)) + facet_grid(SVTYPE ~. ,scales = "free", space = "free_x") +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per SVTYPE across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + scale_x_discrete(drop=FALSE) +
    geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)
  
  p4 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=caller, fill= caller)) + facet_grid(. ~ Size.class) +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per size group across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)
  
  p5 = masterlist_unique_calls %>% ggplot(aes(x=called.by.N.callers)) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    geom_bar(fill = "purple", position = position_dodge2(width = 1.5, preserve = "single"))  + 
    ggtitle(paste(data$Patient," ",data$class, " - SVs called by across callers", sep="")) + ylab("Count (log10)") + 
    geom_text(stat='count', aes(label=..count..), vjust=0, colour = "black")  + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + scale_x_discrete(drop=FALSE) +
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1))
  
  p6 = masterlist_unique_calls %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=SVTYPE, fill= as.factor(called.by.N.callers))) + 
    facet_grid(.~ as.factor(called.by.N.callers), scales = "free", space = "free_x") +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + 
    scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)

  #PDF of plots ####
  size.breakdown = paste0("Very Small < ", v.small.kb, "  | ", v.small.kb, "< Small  <", small.kb, "  | ", small.kb, "< Medium < ",medium.kb,"  | ",medium.kb, "< Large <", large.kb, " | Very Large >",large.kb)
  pdf(file = paste0(path,"Melissa-abacus/SV_database_pipeline/Cross_Technology_Consensus/", data$Patient, ".SVcalling.plots.pdf"), bg = "white", width = 10, height = 7)
  grid.arrange(p1, bottom = paste0(size.breakdown,"\n"))
  grid.arrange(p2, bottom = "Not all callers detect all variant types \n")
  grid.arrange(p3, bottom = "Not all callers detect all variant types \n")
  grid.arrange(p4, bottom = paste0(size.breakdown,"\n"))
  t1= tableGrob(Summary.by.sample.caller, rows=NULL)
  grid.arrange(t1, top = "Summary of calls by caller")
  t2=  tableGrob(Summary.by.size.class, rows=NULL)
  grid.arrange(t2, top = "Summary of calls by size class")
  grid.arrange(p5)
  grid.arrange(p6)
  # Venn diagrams
  #grid.arrange(grobTree(vennDiagram(vennCounts(masterlist_unique_calls[,10:14]), circle.col = 1:length(list.callers), names=list.callers)),top=paste(data$Patient, " Tumor - SVs called by across callers", sep=""))
  dev.off()
  }

