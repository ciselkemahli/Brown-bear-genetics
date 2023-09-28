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

### Reference Genome

I used *Ursus americanus* (American black bear) reference genome as link below. ASM334442v1_HiC.fasta.gz is the file that I used for alignment (Dudchenko et al., 2017; Dudchenko et al., 2018; Srivastava et al., 2019).

https://www.dnazoo.org/assemblies/Ursus_americanus

### World Brown Bears

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

```ruby
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

```ruby
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

```ruby
# Output the genotype likelihoods to obtain beagle file for all brown bear individuals 
~/bin/angsd/angsd -GL 1 -ref UrsMar_ASM334442v1_HiC_zoo.fasta -out arctos -doGlf 2 -doMajorMinor 1 -doMaf 1 -bam ursusarctos.bamlist -nThreads 10 -minInd 22 -minMaf 0.05 -SNP_pval 1e-6 -minMapQ 10 -minQ 20 -setMinDepth 220 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -baq 2 -only_proper_pairs 1 -sites results_paralog/btr.site.list

# Estimate covariance matrix using 10 threads
python ~/bin/pcangsd/pcangsd.py -beagle arctos.beagle.gz -o arctos -threads 10
```

Using R, I plot the PCA plot.

```ruby
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

```ruby
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

```ruby
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



```ruby
library(caret)
library(mlbench)
library(randomForest)
```



```ruby
library(caret)
library(mlbench)
library(randomForest)
```



```ruby
library(caret)
library(mlbench)
library(randomForest)
```



```ruby
library(caret)
library(mlbench)
library(randomForest)
```



```ruby
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
