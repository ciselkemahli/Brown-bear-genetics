<header>

<!--
  <<< Author notes: Course header >>>
  Include a 1280×640 image, course title in sentence case, and a concise description in emphasis.
  In your repository settings: enable template repository, add your 1280×640 social image, auto delete head branches.
  Add your open source license, GitHub uses MIT license.
-->

# Genome analysis of brown bears of Sarıkamış, Kars

This tutorial includes the analysis of RAD-sequencing (A next-generation sequencing method) of brown bears. The samples from brown bears were collected from Sarıkamış, Kars region since 2012 by Kuzey Doğa Derneği (Association). The captured bears were GPS-collared to determine their foraging behavior. The brown bear populations found in these regions with high anthropogenic impact show two basic life history strategies in terms of feeding behavior. The first strategy is ‘sedentary’, in which Eastern Anatolian brown bears use urban garbage dumps as their primary feeding area and rarely forage in natural habitats. The second type of strategy is ‘migratory’, in which Eastern Anatolian brown bears do not normally go to the dump but actively migrate to find food in natural habitats. 

Here, you can see the population genetics analysis for the Eastern Anatolia brown bear population as well as behavioral differences within the population. 

Most of the analysis were performed on bash environment. For this study, we will use some programs:
- [ANGSD - angsd](http://www.popgen.dk/angsd/index.php/ANGSD#Overview)
- [ngsTools/TUTORIAL.md at master · mfumagalli/ngsTools · GitHub](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)
- [PCAngsd - software](http://www.popgen.dk/software/index.php/PCAngsd)
- ngsParalog 

## Reference Genome

I used *Ursus americanus* (American black bear) reference genome as link below. ASM334442v1_HiC.fasta.gz is the file that I used for alignment (Dudchenko et al., 2017; Dudchenko et al., 2018; Srivastava et al., 2019).

https://www.dnazoo.org/assemblies/Ursus_americanus

## World Brown Bears

I added the bear sequence data from different studies. 

1. Barlow et al. (2018) -> 3 *Ursus arctos* from Europe

GEO1 -> Georgia, Great Caucasus

RUS1 -> Balakhtinsky District, Russia

SLO1 -> Hrusica, Slovenia

2. Cahill et al. (2013) -> 7 *Ursus maritimus*, 2 *Ursus arctos*-Alaska

AK1 -> F - Admiralty Island, Alaska

AK2 -> F - Denali National Park, Alaska

UM1 -> M - Lancaster Sound, USA (Smithsonian Natural History Museum)

UM2 -> F - Wrangel Island, Russia (Wrangel Island State Nature Reserve)

UM3 -> M - North Beaufort Sea, Canada (Canada Wildlife Service)

UM4 -> M - Chukchi Sea, USA (US Fish and Wildlife Services)

UM5 -> M - North Beaufort Sea, Canada (US Fish and Wildlife Services)

UM6 -> M - West Hudson Bay, Canada (Canada Wildlife Service)

UM7 -> F - West Hudson Bay, Canada (Canada Wildlife Service)

3. Benazzo et al. (2017) -> 5 Italy Apenine mountains *Ursus arctos* population and 3 from Europe

APN1 -> F - Central Italy

APN2 -> M - Central Italy

APN3 -> F - Central Italy

APN4 -> M - Central Italy

APN5 -> F - Central Italy

GRE1 -> M - Greece

GRE2 -> M - Greece

SLK1 -> M - Slovakia

Download the data from NCBI and aligned to *Ursus americanus* reference genome. All fastq files are formed using SRA toolkit.

In this study, European population includes GRE1, GRE2, SLK1, and SLO1.
Alaskan population includes AK1 and AK2.
RUS1 and GEO1 are stated as two different individuals.
Apennine population includes APN1, APN2, APN3, APN4 and APN5.
*Ursus maritimus* individuals are used for determining ancestral population because the polar bear is the most close relative to brown bear.

## Alignment to *Ursus americanus* ref-seq

```bash
#!/bin/bash
#This script creates different shell scripts for each individuals to decrease the processing time. 

pop=$1 #List includes the name of the individuals stated as in fastq files.
n=$(wc -l ${pop} | awk '{print $1}')
ref=$2 #Fasta reference sequence
        bwa index -a bwtsw ${ref}  #Index reference genome

x=1 
while [ $x -le $n ] 
do
        string="sed -n ${x}p ${pop}"
        str=$($string)
        var=$(echo $str | awk -F"\t" '{print $1}')
        set -- $var
        c1=$1   ### ID ###
                echo "#!/bin/bash" > ${c1}.sh
                echo "" >> ${c1}.sh
	# Align the fastq files to the reference genome using BWA.
                echo "bwa mem ${ref} ../sequences/${c1}_R1.fastq ../sequences/${c1}_R2.fastq > ${c1}.sam" >> ${c1}.sh 
	# Create bam files using Samtools
                echo "~/bin/samtools-1.9/samtools view -bS ${c1}.sam > ${c1}.bam" >> ${c1}.sh
	# Sort and fixmate (add mate scores) using to be used in markdup Samtools
                echo "~/bin/samtools-1.9/samtools sort -n -o ${c1}_namesort.bam ${c1}.bam" >> ${c1}.sh
                echo "~/bin/samtools-1.9/samtools fixmate -m ${c1}_namesort.bam ${c1}_fixmate.bam" >> ${c1}.sh
	# Sort the positions and remove duplicates
                echo "~/bin/samtools-1.9/samtools sort -o ${c1}_positionsort.bam ${c1}_fixmate.bam" >> ${c1}.sh
                echo "~/bin/samtools-1.9/samtools markdup -r ${c1}_positionsort.bam ${c1}_sorted_rmdup.bam" >> ${c1}.sh
	# Index last aligned bam file
                echo "~/bin/samtools-1.9/samtools index ${c1}_sorted_rmdup.bam" >> ${c1}.sh
	# Remove some used but no needed files from the folder
                echo "rm ${c1}.sam ${c1}.bam ${c1}_namesort.bam ${c1}_fixmate.bam" >> ${c1}.sh
	# Send each script to the server
                sbatch -J cisel-align -n 1 -N 1 -p mid -t 720 -o cisel-%j.out ${c1}.sh
        x=$(( $x + 1 ))
done
```

## Filters and Removing Paralogs

```bash
#!/bin/bash
# First to determine paralogs we should have snp position file.
	~/bin/angsd/angsd -GL 1 -ref UrsMar_ASM334442v1_HiC_zoo.fasta -out results_snp/btr.allsite -doMajorMinor 1 -doMaf 1 -minInd 28 -doCounts 1 -minMapQ 10 -minQ 20 -bam btr.bamlist -nThreads 10 -setMinDepth 220 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 2 -doGlf 2

# Get the first two column of btr.allsite.mafs.gz file to obtain all SNPs of the population. -> btr.all.pos (n=21,769,097)

# Determine the paralog sites using ngsParalog
	samtools mpileup -b btr.bamlist -l results_snp/btr.all.pos -f UrsMar_ASM334442v1_HiC_zoo.fasta > results_paralog/btr.all.depth
	~/bin/ngsParalog/ngsParalog calcLR -infile results_paralog/btr.all.depth > results_paralog/btr.all.paralogs

#btr.all.paralogs (n=21,769,097)

# We remove the sites have LRT values bigger than 24 as paralogs. (paralog24.list)

# Remove paralog sites from the all position list and then index it. This list will be used in the later analysis to select all important sites.
	less btr.all.paralogs | awk '$5>24  {print $1, $2}'> paralog24.list (n=23,622)
	grep -Fvxf paralog24_site btr.all.pos > btr.site.list (n=21,745,475)
	~/bin/angsd/angsd sites index btr.site.list
```

## PCA (Principle Component Analysis)

```bash
#!/bin/bash
# Output the genotype likelihoods to obtain beagle file for all brown bear individuals 
~/bin/angsd/angsd -GL 1 -ref UrsMar_ASM334442v1_HiC_zoo.fasta -out arctos -doGlf 2 -doMajorMinor 1 -doMaf 1 -bam ursusarctos.bamlist -nThreads 10 -minInd 22 -minMaf 0.05 -SNP_pval 1e-6 -minMapQ 10 -minQ 20 -setMinDepth 220 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -baq 2 -only_proper_pairs 1 -sites results_paralog/btr.site.list

# Estimate covariance matrix using 10 threads
python ~/bin/pcangsd/pcangsd.py -beagle arctos.beagle.gz -o arctos -threads 10
```

Using R, I plot the PCA plot.

```R
# PLOT PCA using covariance matrix and cluster file.
library(ggplot2)
########## Read input file ##########
covar <- read.table("arctos.cov", stringsAsFact=F)
########## Read annot file ##########
annot <- read.table("arctos.clst", header=T) 
# File includes FID,IID,CLUSTER information for each individual.
 
########## Parse components to analyze ##########
comp <- as.numeric(strsplit("1-2", "-", fixed=TRUE)[[1]])
# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");
# Write eigenvalues
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Behavior <- factor(annot$CLUSTER)
PC$Tra <- factor(annot$IID)
PC$Lab <- factor(annot$FID)

########## Write PC component ##########
title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")
col8=c("#0072B2","#56B4E9","#D55E00","#009E73","#E69F00","#CC79A7","#999999","#000000")

########## Plot ##########
ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Behavior"),size = 3) + ggtitle(title) + theme_classic() + annotate(geom="text", x=0.36, y=-0.2, label="APN", color="black",face="bold") + annotate(geom="text", x=0.1, y=0.23, label="EUR", color="black",face="bold") + annotate(geom="text", x=-0.023, y=-0.08, label="TUR", color="black",face="bold") + geom_point(aes(x=-0.06, y =-0.05), colour = "black", shape = 1, size = 35, stroke=1.5) + geom_point(aes(x=0.385, y=-0.2), colour = "black", shape = 1, size = 14, stroke=1.5) + geom_point(aes(x=0.13, y=0.24), colour = "black", shape = 1, size = 22, stroke=1.5) + scale_color_manual(values=col8, labels = c("Alaska", "Apennine","Turkey","Georgia","Greece","Russia","Slovakia","Slovenia")) + theme(legend.position = c(0.92, 0.85), legend.title = element_blank(), legend.text=element_text(size=15), panel.border = element_rect(colour = "black", fill=NA, size=1)) 
```

## Admixture Analysis

```bash
#!/bin/bash
# Use the beagle file obtained grom genotype likelihood explained in the PCA part.

n=$1    ### Number of replicates ###   30
m=$2    ### Number of ancestral populations ###  8
beag=$3 ### beagle.gz file ###

mkdir results_admixture

for j in $(seq 2 $m)
do
        mkdir results_admixture/admixture_arctos${j}
        for i in $(seq 1 $n)
        do
                ~/bin/NGSadmix -likes ${beag} -K ${j} -P 8 -o results_admixture/admixture_arctos${j}/arctos_admixture${j}_${i} -minMaf 0.05
        done
done
```

Create a log file including best like values from the output files using arctos_admixturej_i.log. File should be in order of each best like values from 2 ancestral populations to 8, each have 30 different values.

Using CLUMPAK (http://clumpak.tau.ac.il/bestK.html) (Kopelman et al., 2015), upload log_prob folder (log prob table file) and find Best K estimation. After Best K was determined, one of qopt file from these folders was dowloaded, and apply below R script to obtain admixture plot.

Using R, Admixture plot is plotted.

```R
########## Packages to be used ########## 
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)  # melt

col6=c("#D55E00","#0072B2","#E69F00","#009E73","#56B4E9","#F0E442","#CC79A7")

########## Admixture plot for K6 ########## 
pop<-read.table("arctos_org2.info",as.is=T) #file contains behavior and name
admix6<-as.matrix(read.table("arctos_admixture6_14.qopt"))

k6=as.data.frame(cbind(pop[,1],pop[,2],admix6))
colnames(k6)=c("Ind","Order","V1","V2","V3","V4","V5","V6")
k6$Order=as.numeric(as.character(factor(k6$Order, levels = unique(k6$Order))))
k6<-k6[order(k6[,2]),]
ad_k6 <- melt(data = k6, id.vars = c("Ind","Order"))
ad_k6$value = as.numeric(as.character(factor(ad_k6$value, levels = unique(ad_k6$value))))
ad_k6$order=factor(ad_k6$Order, levels = unique(ad_k6$Order))
ad_k6$Ind=factor(ad_k6$Ind, levels = unique(ad_k6$Ind))

k6plot = ggplot(ad_k6, aes(x=Ind, y=value, fill = variable, group=Order)) + geom_bar(aes(group=Order), stat="identity",width=1, color = "white", size =0.3, position = position_stack(vjust =1)) +bscale_fill_manual(values=col6) + theme(panel.spacing.x = unit(0.1, "lines"), axis.text.x = element_blank(), panel.grid = element_blank()) + theme_minimal() + labs(x = "", y = "Ancestry") + scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size =8, hjust = 0.5)) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + scale_x_discrete(expand = expansion(add = 1)) + theme(legend.position = "none") + geom_line(aes(x = 5.5), color = "black", size=1) + geom_line(aes(x = 13.5), color = "black", size=1) + geom_line(aes(x = 26.5), color = "black", size=1) 
k6plot
```

## Genetic diversity analysis

```bash
#!/bin/bash
# Calculate saf files and the ML estimate of the sfs using the EM algorithm for each population
# Populations used as bamlist are Europe, Apennine, Turkey, Georgia, Alaska, Russia

# First perform sfs to be used in diversity analysis
~/bin/angsd/angsd -b ${pop}.bamlist -ref UrsMar_ASM334442v1_HiC_zoo.fasta -anc UrsMar_ASM334442v1_HiC_zoo.fasta -out results_sfs/${pop} -GL 1 -doSaf 1 -minMapQ 10 -minQ 20 -nThreads 8 -sites results_paralog/btr.site.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 2
~/bin/angsd/misc/realSFS results_sfs/${pop}.saf.idx -maxIter 100 -P 8 > results_sfs/${pop}.sfs
# Estimate the site allele frequency likelihood
~/bin/angsd/angsd -bam ${pop}.bamlist -out results_diversity/${pop} -ref UrsMar_ASM334442v1_HiC_zoo.fasta -anc UrsMar_ASM334442v1_HiC_zoo.fasta -GL 1 -doThetas 1 -doSaf 1 -pest results_sfs/${pop}.sfs -minMapQ 10 -minQ 20 -nThreads 8 -only_proper_pairs 1 -sites results_paralog/btr.site.list
#Estimate Tajimas D and other statistics for each scaffold
~/bin/angsd/misc/thetaStat do_stat results_diversity/${pop}.thetas.idx -outnames results_diversity/${pop}.thetas
# Calculate diversity statistics with window scan
~/bin/angsd/misc/thetaStat do_stat results_diversity/${pop}.thetas.idx -outnames results_diversity/${pop}_50kb.thetas -win 50000 -step 25000
```

Using python, the obtained file was divided to each scaffolds.

'all_50scaf_tw.txt' file looks like: (tw/nsites)

TURKEY_sc	TURKEY_tw	ALASKA_sc	ALASKA_tw	EUROPE_sc	EUROPE_tw	GEORGIA_sc	GEORGIA_tw	RUSSIA_sc	RUSSIA_tw	APENINE_sc	APENINE_tw

HiC_scaffold_1	0.00146608	HiC_scaffold_1	5.17143e-05	HiC_scaffold_1	0.000279876	HiC_scaffold_1	0.000423312	HiC_scaffold_1	0.000237703	HiC_scaffold_1	1.84962e-05

HiC_scaffold_1	0.00146608	HiC_scaffold_1	5.17143e-05	HiC_scaffold_1	0.000279876	HiC_scaffold_1	0.000423312	HiC_scaffold_1	0.000237703	HiC_scaffold_1	1.84962e-05

Using R, diversity plot for each scaffold was obtained.

```R
########## Function to plot multiple plot at the same time ########## 
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)  
  numPlots = length(plots) 
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y) 
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }  
}
########## Plot with loop and multiplot function ###############
library(ggplot2)
library(scales) ##ggplot scientific
setwd("~/Desktop/Analysis/diversity/scaffolds_20210707/") #Folder includes theta watterson values for the first 37 scaffolds
filelist = list.files(pattern = ".*.txt")
scaflist = lapply(filelist, function(x)read.table(x, header=F)) 
namelist = paste("Scaffold",1:37,sep="")
col6=c("#0072B2","#56B4E9","#E69F00","#009E73","#CC79A7","#D55E00")
myplots <- list()  # new empty list
for (i in 1:37) {
  p1 <- ggplot() + theme_classic() + geom_boxplot(data=scaflist[[i]], aes(x=V1, y=V2, fill=V1), outlier.shape=NA) + scale_y_continuous(limits = seq(0, 0.01, 0.005)) + scale_fill_manual(values=col6,labels = c("AK", "APN", "EUR","GEO", "RUS", "TUR")) + theme(axis.text.x = element_text(angle = 0, vjust = 1, size=5, hjust = 0.5)) + theme(axis.text.y = element_text(angle = 0, vjust = 1, size =5, hjust = 0.5)) + ggtitle(namelist[i]) + ylab(expression(paste(theta,"w"))) + xlab("") + theme(plot.title = element_text(hjust = 0.5, size=7, face="bold"), legend.position = "none") + scale_x_discrete(labels=c("AK", "APN", "EUR","GEO", "RUS", "TUR"))
  myplots[[i]] <- p1  # add each plot into plot list
}
multiplot(plotlist = myplots, cols = 5) #Change cols according to your visualization
```

## Genome-wide Association Study

```bash
#!/bin/bash
# btr3.bamlist includes 30 individuals (13 migratory and 17 sedentary - BTR36 removed due to undetermined behavior)
# btr3.ybin -> A file containing the case control status. 0 meaning migratory and 1 meaning sedentary. Check http://www.popgen.dk/angsd/index.php/Association for detailed information.
~/bin/angsd/angsd -GL 1 -doAsso 4 -out results_association/wr -doMajorMinor 1 -doMaf 1 -bam btr3.bamlist -nThreads 8 -doPost 1 -yBin results_association/btr3.ybin -SNP_pval 1e-6 -minMaf 0.05 -minHigh 10 -doCounts 1 -minCount 3 -minMapQ 10 -minQ 20 -only_proper_pairs 1 -sites results_paralog/btr.site.list
# Outputs are lrt0 and maf files. Our interest is lrt0 file. -999 values should be omitted, indicating no information. 
# Sort lrt0 file based on LRT (log-likelihood ratios). Determine the significance value by chi-squared with one degree of freedom. Select the sites based on the significance value. 
# Create top selected file (wr.top) including scaffold and site names. Perform genotyping for these determined sites with 90% posterior probability.
~/bin/angsd/angsd -ref UrsMar_ASM334442v1_HiC_zoo.fasta -GL 1 -out results_association/wr_top90 -doGlf 2 -doMajorMinor 1 -doMaf 1 -bam btr3.bamlist -nThreads 8 -doGeno 5 -doPost 1 -postCutoff 0.90 -only_proper_pairs 1 -sites results_association/wr.top
```

R
```R
# Table includes genotypes as Anc, Der, Het and NN values for each individuals at each site. 
########## Reading tables ########## 
wr_geno=read.table("wr_90", header = T)
wr <- melt(wr_geno, id.vars=c("Site","Code"))
wr$Code = factor(wr$Code, levels = unique(wr$Code))
wr$value = factor(wr$value, levels = unique(wr$value))
ind_name=unique(wr[,3])
position_name=unique(wr[,1])
col3 = c("#E69F00","#0072B2","#999999","#009E73") #c("Het","Anc","NN","Der")
g <- ggplot(wr, aes(x=Code, y=variable, fill=value))
genplot = g + geom_tile(color=1, size=0.2) + scale_fill_manual(values=col3, labels = c("Het","Anc","NN","Der")) + guides(fill=guide_legend(title="Alleles")) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + theme(axis.text.x = element_text(angle = 60, vjust = 1, size = 10, hjust = 1)) + scale_y_discrete(name ="Individuals", labels=ind_name) + scale_x_discrete(name ="Position", labels=position_name)
genplot
grid.text("Migratory", x = unit(0.94, "npc"), y = unit(0.46, "npc"), rot=270,gp=gpar(fontsize=13))
grid.text("Sedentary", x = unit(0.94, "npc"), y = unit(0.80, "npc"), rot=270,gp=gpar(fontsize=13))
```

## Population branch statistics

```bash
#!/bin/bash
# PBS for European-Apenine-Turkey/Wild/Resident 
~/bin/angsd/misc/realSFS results_sfs/${pop1}.saf.idx results_sfs/${pop2}.saf.idx -P 8 > results_fst/${pop1}_${pop2}.ml
~/bin/angsd/misc/realSFS results_sfs/${pop1}.saf.idx results_sfs/${pop3}.saf.idx -P 8 > results_fst/${pop1}_${pop3}.ml
 ~/bin/angsd/misc/realSFS results_sfs/${pop2}.saf.idx results_sfs/${pop3}.saf.idx -P 8 > results_fst/${pop2}_${pop3}.ml
~/bin/angsd/misc/realSFS fst index results_sfs/${pop1}.saf.idx results_sfs/${pop2}.saf.idx results_sfs/${pop3}.saf.idx -fstout results_fst/${pop1}_${pop2}_${pop3}.pbs -sfs results_fst/${pop1}_${pop2}.ml -sfs results_fst/${pop1}_${pop3}.ml -sfs results_fst/${pop2}_${pop3}.ml
#Window scan
~/bin/angsd/misc/realSFS fst stats2 results_fst/${pop1}_${pop2}_${pop3}.pbs.fst.idx -P 8 -win 50000 -step 25000 > results_fst/${pop1}_${pop2}_${pop3}.pbs.txt
```

Asso-PBS-wrt file includes above values. We just checked the PBS values for the scanned windows includes previously found associated sites.

→ AssociationSites	Nsites	Fst01_Tur	Fst02_Tur	Fst12_Tur	PBS2_Tur	Fst02_Sed	Fst12_Sed	PBS2_Sed	Fst02_Mig	Fst12_Mig	PBS2_Mig

```R
library(ggplot2)
library(reshape2)
########## PBS Plot ##########
asso<-read.table("Asso-PBS-wrt",as.is=T,header=T)
global_res = 0.084006
global_mig = 0.081761
global_tur = 0.055098

asso1=matrix(0,nrow=nrow(asso),ncol=4)
asso1[,1]=asso$AssociationSites
asso1[,2]=asso$PBS2_Sed
asso1[,3]=asso$PBS2_Mig
asso1[,4]=asso$PBS2_Tur
colnames(asso1)=c("Position","Sedentary", "Migratory","Turkey")
asso1=as.data.frame(asso1)

asso2 <- melt(data = asso1, id.vars = "Position")
asso2$Position=factor(asso2$Position, levels = unique(asso2$Position))
asso2$value = as.numeric(as.character(factor(asso2$value, levels = unique(asso2$value))))

ggplot(asso2, aes(y=Position, x=value, color=variable, group=variable)) + theme_bw() + geom_point(aes(shape=variable, color=variable), size = 4, stroke = 1) + scale_shape_manual(values=c(15,16,3)) + theme(axis.text.x = element_text(angle = 0, vjust = 1, size =10, hjust = 0.5)) + theme(axis.text.y = element_text(angle = 0,size =10)) + scale_color_manual(values=col3, labels = c("Sedentary","Migratory","Turkey")) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + xlab("PBS") + ylab("Positions") + geom_vline(xintercept=global_res,linetype="longdash",color="#CC79A7",size=1.5) + geom_vline(xintercept=global_mig,linetype="longdash",color="#0072B2",size=1.5) + geom_vline(xintercept=global_tur,linetype="longdash",color="#D55E00",size=1.5) + theme(legend.position = c(0.9, 0.92), legend.title = element_blank(), legend.background = element_rect(linetype="solid",colour="black"))
```


```bash
library(caret)
library(mlbench)
library(randomForest)
```


## Referanslar

Barlow, A., Cahill, J. A., Hartmann, S., Theunert, C., Xenikoudakis, G., Fortes, G. G., Paijmans, J. L., Rabeder, G., Frischauf, C., & Grandal-d’Anglade, A. (2018). Partial genomic survival of cave bears in living brown bears. Nature ecology & evolution, 2(10), 1563-1570.

Benazzo, A., Trucchi, E., Cahill, J. A., Maisano Delser, P., Mona, S., Fumagalli, M., Bunnefeld, L., Cornetti, L., Ghirotto, S., & Girardi, M. (2017). Survival and divergence in a small group: The extraordinary genomic history of the endangered Apennine brown bear stragglers. Proceedings of the National Academy of Sciences, 114(45), E9589-E9597.

Cahill, J. A., Green, R. E., Fulton, T. L., Stiller, M., Jay, F., Ovsyanikov, N., Salamzade, R., St. John, J., Stirling, I., & Slatkin, M. (2013). Genomic evidence for island population conversion resolves conflicting theories of polar bear evolution. PLoS genetics, 9(3), e1003345.

Dudchenko, O., Batra, S. S., Omer, A. D., Nyquist, S. K., Hoeger, M., Durand, N. C., Shamim, M. S., Machol, I., Lander, E. S., & Aiden, A. P. (2017). De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds. science, 356(6333), 92-95.

Dudchenko, O., Shamim, M. S., Batra, S. S., Durand, N. C., Musial, N. T., Mostofa, R., Pham, M., St Hilaire, B. G., Yao, W., & Stamenova, E. (2018). The Juicebox Assembly Tools module facilitates de novo assembly of mammalian genomes with chromosome-length scaffolds for under $1000. BioRxiv, 254797.

Kopelman, N. M., Mayzel, J., Jakobsson, M., Rosenberg, N. A., & Mayrose, I. (2015). Clumpak: a program for identifying clustering modes and packaging population structure inferences across K. Molecular ecology resources, 15(5), 1179-1191.

Srivastava, A., Kumar Sarsani, V., Fiddes, I., Sheehan, S. M., Seger, R. L., Barter, M. E., Neptune-Bear, S., Lindqvist, C., & Korstanje, R. (2019). Genome assembly and gene expression in the American black bear provides new insights into the renal response to hibernation. DNA Research, 26(1), 37-44.




