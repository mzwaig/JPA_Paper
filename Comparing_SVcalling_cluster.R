unlink(".RData")
rm(list = ls())
args <- commandArgs()
print(args)
data_path <- args[6]
#subtypes <- as.vector(unique(as.data.frame(strsplit(args[7], ",")))[,1])
#print(subtypes)
cluster = T

#setwd("Y:/Melissa-abacus/")
#cluster = F
#subtypes = "GBM"
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
require(data.table) ## 1.9.4+
library(gridExtra)
library(grid)
library(gsubfn)
library(limma)
library(openxlsx)
library(tibble)
library(stringr)
options(digits=10)
#library(wordcloud2)
path = ifelse(cluster == T, "/lb/project/ioannisr/", "Y:/")


#conda ####
#conda install -c r r-stringi r-dplyr r-tidyr r-ggplot2 r-data.table r-gridextra
#conda install -c bioconda grid bioconductor-limma
#conda install -c conda-forge r-gsubfn
#####

datasets = read.delim(file= paste0(data_path), sep=",",header=T, na.strings=c("","NA")) #use this to open files
#datasets = read.delim(file= paste0(path,"Melissa-abacus/Pipelines_MZ/All_datasets.csv"), sep=",",header=T, na.strings=c("","NA")) #use this to open files
#datasets= datasets %>% filter(diagnosis %in% subtypes) #%>% filter(class == "Tumor")
#datasets= datasets %>% filter(Patient=="MDT-AP-2878") # | Patient=="CPP_3T-N_103_16")

#datasets=datasets[55:56,]
print(paste0(datasets$Patient,"_", datasets$class))

list.of.samples = datasets %>% mutate(Patient=paste0(Patient,".",class)) 
list.of.samples = list.of.samples[order(list.of.samples$Patient),]
list.of.samples = as.vector(unique(list.of.samples$Patient))
list.callers= c("LongRanger","SvABA","LinkedSV","NAIBR","GROC_SV")
#x=datasets[1,]
#data=x

# File check ####
for(i in 1:nrow(datasets)) {
  x = datasets[i,]
  print(x$Patient)
  check = rbind.data.frame(cbind.data.frame(caller = "LongRanger", Exists = file.exists(paste0(x$LongRanger))),
                           cbind.data.frame(caller = "LongRanger DELs", Exists = file.exists(paste0(x$file_path,"_custom_files/", x$LR_name, "_LR_dels_parsed.vcf"))),
                           cbind.data.frame(caller = "SvABA", Exists = file.exists(paste0(x$SvABA))),
                           cbind.data.frame(caller = "LinkedSV", Exists = file.exists(paste0(x$LinkedSV, "LSV.filtered_large_svcalls.bedpe"))),
                           cbind.data.frame(caller = "NAIBR", Exists = file.exists(paste0(x$NAIBR))),
                           cbind.data.frame(caller = "GROC-SV", Exists = file.exists(paste0(path,"NOBACKUP/Melissa-nobackup/10X_linked-reads/GROC_SV/",x$LR_name,"/classes.txt"))))
  print(check)
  ifelse(all(check$Exists) == TRUE, "", print(paste0(x$Patient, " is missing files")))
  #stopifnot(all(check$Exists) == TRUE)
  }

# Paramters ####
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


# Load SV call files (A.1) ####
load.calls.by.sample <- function(x){
  masterlist=NULL
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
           SVTYPE=SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO=info)
  masterlist=as.data.frame(rbind(masterlist,LR.bedpe)) 
  # Load LongRanger - dels ####
  print(paste(x$LR_name, x$class," - Load LongRanger, dels",sep=" "))
  LR.dels.all = fread(file= paste0(x$file_path,"_custom_files/", x$LR_name, "_LR_dels_parsed.vcf"), stringsAsFactors = FALSE, check.names=T)
  colnames(LR.dels.all)[15] = "GT"
  LR.dels = LR.dels.all %>% filter(SVTYPE=="DEL") %>% 
    mutate(pair_ID=gsub("call", "del", ID), CHR2=CHROM, caller="LongRanger", class=x$class, Technology="10X", COPY= NA,
           INFO = paste0("HAP_ALLELIC_FRAC=",HAP_ALLELIC_FRAC,"ALLELIC_FRAC=",ALLELIC_FRAC, "PS=", PS)) %>% 
    select(CHROM=CHROM, START=POS, CHR2=CHR2, END , ID=pair_ID, QUAL, FILTER, 
           SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO)
  masterlist=as.data.frame(rbind(masterlist,LR.dels)) 
  # BNDs
  print(paste(x$Patient,x$class," - Load LongRanger, BNDs",sep=" "))
  LR.BNDs = LR.dels.all %>% filter(SVTYPE=="BND") %>% mutate(pair_ID=gsub("_[12]", "", ID)) %>% mutate(pair_ID=gsub("call", "del",pair_ID))
  if(nrow(LR.BNDs)>0){
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
             SVTYPE, caller=caller , class=class, COPY= COPY, Technology=Technology, INFO)
    masterlist=as.data.frame(rbind(masterlist,LR.BNDs2))
  }
  # ploidy based on LongRanger CNvs ####
  #CNV.10X = LR.bedpe %>% filter(!is.na(COPY)) %>% mutate(SVLEN=as.numeric(END)-as.numeric(START)) %>% 
  #  group_by(CHROM, COPY) %>% summarise(count=n(), bases= sum(SVLEN))
  #t2 = merge( CNV.10X, chr_lengths, by.x= "CHROM", by.y= "Chromosome")
  #t3 = t2 %>% group_by(CHROM, `Total.length.(bp)`) %>% summarise(CNA.bases=sum(bases)) %>% 
  #  mutate(chr2.bases= (`Total.length.(bp)`-CNA.bases), COPY=2)
  #t4=as.data.frame(rbind(cbind(chrom= t2$CHROM, COPY=t2$COPY, bases=t2$bases, Total.len=t2$`Total.length.(bp)`),
  #                       cbind(chrom= t3$CHROM, COPY=t3$COPY, bases=t3$chr2.bases, Total.len=t3$`Total.length.(bp)`)))
  #t4= t4 %>% mutate(chr.portion=(as.numeric(as.character(bases)))/as.numeric(as.character(Total.len))*100) 
  #t5 = t4 %>% group_by(chrom=as.character(chrom), chr.start=1, chr.end=as.numeric(as.character(Total.len))) %>% summarise(ploidy = weighted.mean(as.numeric(as.character(COPY)),
  #                                                                 as.numeric(as.character(bases))))
  #t5  %>% ggplot(aes(y= ploidy)) + geom_rect(aes(xmin = as.numeric(chr.start), xmax= as.numeric(chr.end), 
  #                                           ymin=as.numeric(ploidy)-0.5 , ymax=as.numeric(ploidy)+0.5)) + facet_grid(~ chrom)
  
  # Load SvABA ####
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
    SvABA.parsed = SvABA.parsed %>% mutate(caller="SvABA", class=x$class, Technology="10X",SVTYPE=NA, COPY=NA, ID=paste("SvABA_",pair_ID,sep=""))
    #SvABA.parsed["INFO"] = SvABA.parsed$INFO.x %>% strapplyc(";(.*)", simplify = TRUE) # drop supporting barcodes
    SvABA.parsed = SvABA.parsed %>% mutate(INFO=strapplyc(SvABA.parsed$INFO.x, ";(.*)", simplify = TRUE))
    # select common headers
    SvABA.parsed = SvABA.parsed %>% 
      select(CHROM=CHROM, START=POS, CHR2=CHR2 , END=END , ID=ID,QUAL=QUAL,  FILTER=FILTER, 
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO)
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
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO)
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
             SVTYPE=SVTYPE, caller=caller, class=class, COPY= COPY, Technology=Technology, INFO=INFO)
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
    #grocsv.genotypes = fread(file= paste0("Y:/NOBACKUP/Melissa-nobackup/10X_linked-reads/GROC_SV/",x$Patient,"/genotypes.tsv"), sep="\t", stringsAsFactors = FALSE, check.names=T, header=TRUE)
    grocsv.full = data.frame(chromx = grocsv.genotypes$chromx, chromy= grocsv.genotypes$chromy, cluster = grocsv.genotypes$cluster,
                             x= grocsv.genotypes$x, y= grocsv.genotypes$y, class = grocsv.classes$classes, QUAL=grocsv.genotypes$blacklist , 
                             SVLEN=grocsv.genotypes$dist, FILTER=grocsv.genotypes$quality)
    #grocsv.full$cluster <- sub("^", "cluster_", grocsv.full$cluster) # cluster is not unique since GROC-sv can detect complex rearrangements
    grocsv.full = grocsv.full %>% mutate(caller="GROC_SV", SVTYPE=NA, COPY=NA,Technology="10X",ID=paste0("GROCSV_",rownames(grocsv.full)),INFO=NA)
    # If SVLEN comes back negative, SV is Intrachromosomal and the order of the events is reversed 
    grocsv.full["Inter.intra"] = ifelse(as.character(grocsv.full$chromx)==as.character(grocsv.full$chromy), "Intra", "Inter")
    grocsv.full["SVLEN"] = ifelse(grocsv.full$Inter.intra=="Intra", grocsv.full$y - grocsv.full$x, "")
    grocsv.full = grocsv.full %>% mutate(CHROM=ifelse(SVLEN<0,as.character(chromy),as.character(chromx)),
                                         CHR2=ifelse(SVLEN<0,as.character(chromx),as.character(chromy)), 
                                         START=ifelse(SVLEN<0,y,x), END=ifelse(SVLEN<0,x,y))
    grocsv.full["SVLEN2"] = ifelse(grocsv.full$Inter.intra=="Intra", grocsv.full$END - grocsv.full$START, "") #check that all length are pos
    # select common headers
    grocsv.full = grocsv.full %>% 
      select(CHROM=CHROM , START=START , CHR2=CHR2 , END=END , ID=ID,  QUAL=QUAL, FILTER=FILTER, 
             SVTYPE= SVTYPE, caller= caller, class= class , COPY= COPY, Technology=Technology, INFO=INFO) %>%
      filter(class == as.character(x$class))
    #print(paste0("GROC-SV classes :", unique(grocsv.full$class)))
    masterlist=as.data.frame(rbind(masterlist,grocsv.full))
    }
  #Make master list with common headers ####
  print(paste(x$Patient,x$class," - Make masterlist",sep=" "))
  masterlist = masterlist %>% mutate(caller.class=paste(caller,".",class, sep=""), Patient = x$Patient)
  masterlist["Inter.intra"] = ifelse(masterlist$CHROM==masterlist$CHR2, "Intra", "Inter")
  masterlist["SVLEN"] = ifelse(masterlist$Inter.intra=="Intra", as.numeric(masterlist$END) - as.numeric(masterlist$START), NA)
  masterlist["Size.class"] = ifelse(masterlist$Inter.intra=="Inter","Interchr",
                                    ifelse(masterlist$SVLEN<v.small,"Very Small",
                                    ifelse(masterlist$SVLEN>=v.small & masterlist$SVLEN< small,"Small",
                                         ifelse(masterlist$SVLEN>=small & masterlist$SVLEN< medium, "Medium",
                                                ifelse(masterlist$SVLEN>=medium & masterlist$SVLEN< large, "Large",
                                                       ifelse(masterlist$SVLEN>= large, "Very Large","???"))))))
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
    masterlist.overlap = merge(masterlist.overlap, IDs.filt.info, by=colnames(masterlist.interval), all=TRUE)
    }
  # Processing overlapping calls ####
  print("Make unique_ID")
  masterlist.overlap2 = masterlist.overlap %>% unite("unique_ID", paste0(list.callers,".ID"), sep = ".", remove = FALSE)
  masterlist.overlap2$unique_ID = gsub(".NA\\w(*SKIP)(*FAIL)|.NA", "", masterlist.overlap2$unique_ID, perl=TRUE)
  masterlist.overlap2$unique_ID = gsub("NA.", "", masterlist.overlap2$unique_ID, fixed = TRUE)
  print("Count callers which detect SV")
  masterlist_unique_calls= masterlist.overlap2 %>% 
    select(-c("START","END","ID", "QUAL","FILTER", "SVTYPE","caller","COPY", "caller.class","SVLEN","INFO"))  %>%
    mutate(called.by.LongRanger=ifelse(!is.na(LongRanger.ID),as.numeric(1), 0), called.by.SvABA=ifelse(!is.na(SvABA.ID),as.numeric(1), 0), 
           called.by.LinkedSV=ifelse(!is.na(LinkedSV.ID),as.numeric(1), 0), called.by.NAIBR=ifelse(!is.na(NAIBR.ID),as.numeric(1), 0), 
           called.by.GROC_SV=ifelse(!is.na(GROC_SV.ID),as.numeric(1), 0), 
           called.by.N.callers=(as.numeric(called.by.LongRanger) + as.numeric(called.by.SvABA) + as.numeric(called.by.LinkedSV) + as.numeric(called.by.NAIBR) + as.numeric(called.by.GROC_SV))) %>%
    select("unique_ID" , "Patient","class" , "Technology", "Inter.intra", "Size.class", CHROM = CHROM, CHR2= CHR2, contains("called"), everything())
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
    select("unique_ID" , "Patient","class" , "Technology", "Inter.intra", "Size.class", "CHROM", "CHR2", "SVTYPE",
           contains("called"), everything())
  
  return(masterlist_unique_calls)
  }


#OverlappingCalls_AllCaller_AllSamples=NULL
#AllSamples_AllCalls = NULL
for(i in 1:nrow(datasets)) {
  data = datasets[i,]
  print(data$Patient)
  # put blood and tumor file in the same place
  dir.create(paste0(data$file_path,"_custom_files/Comparing_SVs"), showWarnings = FALSE) #stops warnings if folder already exists
  #Call load.calls function (A.1) ####
  masterlist = load.calls.by.sample(data)
    # output 1
  fwrite(masterlist, paste(data$file_path,"_custom_files/Comparing_SVs/Masterlist_AllCallers.csv",sep=""))
  #masterlist = fread(file= paste0(data$file_path,"_custom_files/Comparing_SVs/Masterlist_AllCallers.csv"), stringsAsFactors = FALSE, check.names=T)
  #Call overlapping calls on masterlist generate by load.calls function (A.2) ####
  print("Overlapping calls")
  masterlist_unique_calls = overlapping.calls.by.sample(data,masterlist)
  fwrite(masterlist_unique_calls, paste(data$file_path,"_custom_files/Comparing_SVs/Overlapping_calls_across_callers.csv",sep=""))
  
  #masterlist_unique_calls = fread(file= paste0(data$file_path,"_custom_files/Comparing_SVs/Overlapping_calls_across_callers.csv"), stringsAsFactors = FALSE, check.names=T)
  #Make across sample masterlist ####
  #OverlappingCalls_AllCaller_AllSamples = as.data.frame(rbind(OverlappingCalls_AllCaller_AllSamples,masterlist_unique_calls))
  #AllSamples_AllCalls = as.data.frame(rbind(AllSamples_AllCalls,masterlist))
  #Summary tables ####
  #Summary.by.size.class = masterlist %>% group_by(class, caller, Size.class) %>% summarize(count=n())
  Summary.by.size.class = masterlist %>% group_by(class, caller, Size.class) %>% dplyr::summarise(count = n())
  Summary.by.size.class =  Summary.by.size.class %>% spread(key=Size.class,value = count)
  Summary.by.size.class[is.na(Summary.by.size.class)] = 0
  Summary.by.size.class['Total.Calls'] = rowSums(Summary.by.size.class[3:8])
  Summary.by.size.class['More.Than.10kb'] = rowSums(Summary.by.size.class[4:7])
  # discrepencies with common. files are due to filtering of contigs
  Summary.by.sample.caller = masterlist %>% group_by(class, caller) %>% 
    dplyr::summarise(count=n(), Mean.SVLEN = trunc(mean(SVLEN, na.rm=TRUE)), 
              Median.SVLEN = trunc(median(SVLEN, na.rm=TRUE)), Min.SVLEN = min(SVLEN, na.rm=TRUE), Max.SVLEN = max(SVLEN, na.rm=TRUE))
  Summary.by.sample.caller[is.na(Summary.by.sample.caller)] = 0
  
  
  # Plots ####
  print(paste(data$Patient," - Plot",sep=""))
  p1 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=Size.class, fill= caller)) + facet_wrap( ~ class) +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE))  + ylab("Count (log10)") + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per size group across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + theme_bw() + theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) 
  
  p2 = masterlist %>% ggplot(aes(x=SVTYPE, fill= caller)) + facet_wrap( ~ class) + 
    scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + ylab("Count (log10)") +  geom_bar(position = position_dodge2(width = 2, preserve = "single")) + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per SVTYPE across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + theme_bw() + theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) 
  
  p3 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=SVTYPE, fill= caller)) + facet_grid(caller ~ class) +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per SVTYPE group across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count', aes(label=..count..), vjust=0.5, colour = "black", size=2)
  
  p4 = masterlist %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=caller, fill= caller)) + facet_grid(class ~ Size.class) +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per size group across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count', aes(label=..count..), vjust=0.5, colour = "black", size=2)
  
  p5 = masterlist_unique_calls %>% ggplot(aes(x=called.by.N.callers)) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    geom_bar(fill = "purple", position = position_dodge2(width = 1.5, preserve = "single"))  + 
    ggtitle(paste(data$Patient," ",data$class, " - SVs called by across callers", sep="")) + ylab("Count (log10)") + 
    geom_text(stat='count', aes(label=..count..), vjust=0, colour = "black")  + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1))
  
  p6 = masterlist_unique_calls %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=SVTYPE, fill= as.factor(called.by.N.callers))) + 
    facet_grid(as.factor(called.by.N.callers) ~ class) +
    geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
    ggtitle(paste(data$Patient," ",data$class, " - Number of calls per SVTYPE group across callers", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
    theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
    geom_text(stat='count', aes(label=..count..), vjust=0.5, colour = "black", size=2)
  
  #PDF of plots ####
  size.breakdown = paste0("Very Small < ", v.small.kb, "  | ", v.small.kb, "< Small  <", small.kb, "  | ", small.kb, "< Medium < ",medium.kb,"  | ",medium.kb, "< Large <", large.kb, " | Very Large >",large.kb)
  pdf(file = paste(data$file_path,"_custom_files/Comparing_SVs/", data$Patient, ".SVcalling.plots.pdf", sep= ""), bg = "white", width = 10, height = 7)
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
  grid.arrange(grobTree(vennDiagram(vennCounts(masterlist_unique_calls[,10:14]), circle.col = 1:length(list.callers), names=list.callers)),top=paste(data$Patient, " Tumor - SVs called by across callers", sep=""))
  dev.off()
   }

# 
# #fwrite(OverlappingCalls_AllCaller_AllSamples, paste("Comparing_SVs/Overlapping_SVs_AllCallers_AllSamples.csv",sep=""))
# 
# orderlist = c("1","2","3","4","5")
# OverlappingCalls_AllCaller_AllSamples <- transform(OverlappingCalls_AllCaller_AllSamples, 
#                                                    called.by.N.callers = factor(called.by.N.callers, levels = orderlist))
# orderlist = c("Very Small","Small","Medium","Large","Very Large","Interchr")
# AllSamples_AllCalls <- transform(AllSamples_AllCalls, Size.class = factor(Size.class, levels = orderlist))
# 
# #### Across sample summary ####
# # Tables ####
# print("Make summary table")
# # calculate SVs and CNVs together since each caller will have different cutoffs
# Summary.by.size.class.ALL= AllSamples_AllCalls %>% group_by(Patient,class, caller, Size.class) %>% dplyr::summarise(count=n())
# Summary.by.size.class.ALL =  Summary.by.size.class.ALL %>% spread(key=Size.class,value = count)
# Summary.by.size.class.ALL[is.na(Summary.by.size.class.ALL)] = 0
# Summary.by.size.class.ALL['Total.Calls'] = rowSums(Summary.by.size.class.ALL[4:9])
# Summary.by.size.class.ALL['More.Than.10kb'] = rowSums(Summary.by.size.class.ALL[5:8])
# Summary.by.size.class.ALL['Less.Than.10kb'] = Summary.by.size.class.ALL$`Very Small`
# 
# # discrepencies with common. files are due to filtering of contigs
# Summary.by.sample.caller.ALL = AllSamples_AllCalls %>% group_by(Patient, class, caller) %>% 
#   dplyr::summarise(Total.Calls=n(), Mean.SVLEN = trunc(mean(SVLEN, na.rm=TRUE)), 
#             Median.SVLEN = trunc(median(SVLEN, na.rm=TRUE)), Min.SVLEN = min(SVLEN, na.rm=TRUE), Max.SVLEN = max(SVLEN, na.rm=TRUE))
# Summary.by.sample.caller.ALL[is.na(Summary.by.sample.caller.ALL)] = 0
# 
# Summary.by.callers = AllSamples_AllCalls  %>% group_by(Patient, class, caller) %>% 
#   dplyr::summarise(count=n()) %>% spread(key=caller,value = count)
# Summary.by.callers[is.na(Summary.by.callers)] = 0
# len.callers=length(Summary.by.callers)
# Summary.by.callers['Total.Calls'] = rowSums(Summary.by.callers[3:len.callers])
# 
# # Summarize SVs called by >1 callers
# Summary.by.N.callers = OverlappingCalls_AllCaller_AllSamples  %>% group_by(class, called.by.N.callers) %>% 
#   dplyr::summarise(count=n()) %>% spread(key=called.by.N.callers,value = count)
# Summary.by.N.callers[is.na(Summary.by.N.callers)] = 0
# len.N.callers=length(Summary.by.N.callers)
# Summary.by.N.callers['Total.Calls'] = rowSums(Summary.by.N.callers[2:len.N.callers])
# Summary.by.N.callers['Made.by.2plus.callers'] = rowSums(Summary.by.N.callers[3:len.N.callers])
# Summary.by.N.callers = Summary.by.N.callers %>% mutate(Frac.Made.by.2plus.callers= (Made.by.2plus.callers / Total.Calls)*100)
# 
# print("Make xlsx")
# list_of_datasets <- list("By.size.class" = Summary.by.size.class.ALL, "Size.by.caller" = Summary.by.sample.caller.ALL,
#                          "Count.by.caller"= Summary.by.callers, "Count.by.number.of.callers"= Summary.by.N.callers)
# write.xlsx(list_of_datasets, file = "SV_database_pipeline/Comparing_SVs_Summary.xlsx")
# 
# # Plot ####
# print("Plot Graphs")
# #  want plot the looks at the total number of calls made by >1 callers, split by Tumor/blood
# p1 = OverlappingCalls_AllCaller_AllSamples  %>% filter(class=="Tumor") %>% ggplot(aes(x=called.by.N.callers)) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
#   geom_bar(fill = "purple", position = position_dodge2(width = 1.5, preserve = "single"))  + 
#   ggtitle(paste("Tumor - SVs called by across callers", sep="")) + ylab("Count (log10)") + 
#   geom_text(stat='count', aes(label=..count..), vjust=0, colour = "black",  size=4)  + theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
#   theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + scale_x_discrete(drop=FALSE)
# 
# p2 = OverlappingCalls_AllCaller_AllSamples %>% filter(class=="Blood") %>% ggplot(aes(x=called.by.N.callers)) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
#   geom_bar(fill = "purple", position = position_dodge2(width = 1.5, preserve = "single"))  + 
#   ggtitle(paste("Blood - SVs called by across callers", sep="")) + ylab("Count (log10)") + 
#   geom_text(stat='count', aes(label=..count..), vjust=0, colour = "black",  size=4)  + theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
#   theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + scale_x_discrete(drop=FALSE)
# 
# p3 =  AllSamples_AllCalls  %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=Patient, fill= caller)) + facet_grid(Size.class ~ class) +
#   geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
#   ggtitle("Number of calls per Sample across callers") + 
#   theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
#   theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
#   geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=2, angle = 90)
# 
# p4 =  AllSamples_AllCalls  %>% filter(!is.na(Size.class)) %>% ggplot(aes(x=caller, fill= caller)) + facet_grid(class ~ Size.class) +
#   geom_bar(position = position_dodge2(width = 1.5, preserve = "single")) + scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
#   ggtitle("Number of calls per size range across callers") + 
#   theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + ylab("Count (log10)") + 
#   theme_bw() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   theme(axis.text.y.left = element_text(angle = 0, hjust = 1)) + 
#   geom_text(stat='count',  position = position_dodge(width = 1), aes(label=..count..), hjust=1, colour = "black", size=3, angle = 90)
# 
# pdf(file = paste("SV_database_pipeline/Comparing_SVs_AllCalls_AllCallers_AllSamples.plots.pdf", sep= ""), bg = "white", width = 10, height = 7)
# grid.arrange(p1, p2)
# grid.arrange(p3, bottom = paste0("Very Small < 10kb  | 10kb< Small  <30kb  | 30kb< Medium < 100kb  | 100kb< Large <1000kb | Very Large >1000kb \n"))
# grid.arrange(p4, bottom = paste0("Very Small < 10kb  | 10kb< Small  <30kb  | 30kb< Medium < 100kb  | 100kb< Large <1000kb | Very Large >1000kb \n"))
# tumor= OverlappingCalls_AllCaller_AllSamples  %>% filter(class=="Tumor")
# grid.arrange(grobTree(vennDiagram(vennCounts(tumor[,10:14]), circle.col = 1:length(list.callers), names=list.callers)),top=paste(" Tumor - SVs called by across callers", sep=""))
# blood= OverlappingCalls_AllCaller_AllSamples  %>% filter(class=="Blood")
# grid.arrange(grobTree(vennDiagram(vennCounts(blood[,10:14]), circle.col = 1:length(list.callers), names=list.callers)),top=paste(" Blood - SVs called by across callers", sep=""))
# dev.off()
# #####
