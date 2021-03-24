# script to read and parse data files
rm(list=ls())
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

#############################################################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#chrom_list <- c("chr9","chr13","chr22")
chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
#cov_file <- "./coverage_phased_tenx_"
#cov_file <- "./coverage_phased_tenx_bam1_"
cov_file <- "./coverage_phased_tenx_bam2_"
#cov_file <- "./coverage_phased_illumina_"

#chr_choice <- "chr10"
chr_choice <- "chr17"
binsize <- 100000

cov_tracks = tibble();
for (chr in chrom_list) {
  file_string <- paste( cov_file , chr , ".dat" , sep="" )
  cov_data <- read_delim( file_string , delim="\t", col_names = FALSE )
  cov_data <- dplyr::rename( cov_data ,"chromosome"="X1","pos"="X2","hap"="X3","hapA"="X6","hapB"="X7")
  cov_data <- cov_data %>% transform(pos_bin = as.integer(pos/binsize))
  #cov_data <- cov_data %>% transform(pos_bin = cut(pos,breaks=linear_cut))
  cov_tracks <- bind_rows( cov_tracks , cov_data )
}


#chr_track <- cov_tracks %>% filter(chromosome == "chr3")  #chr_track <- chr_track %>% transform(pos_bin = cut(pos,breaks=linear_cut))
bin_track <- cov_tracks %>% group_by(chromosome,pos_bin) %>% summarise(A = mean(hapA),B = mean(hapB),num = n(),pos=first(pos))
bin_track <- bin_track %>% gather(A,B,key = hap,value = cov)
cov_tracks <- cov_tracks %>% gather(hapA,hapB,key = hap,value = cov)
#bin_track <- bin_track %>% filter(num > 10)
bin_track <- bin_track %>% mutate(cov_norm = cov/6.25)
bin_track <- bin_track %>% filter( num > 10 )

#########################################
unphased_cn_file <- "./K562_cn_binned_100000_tenx_bam1.dat"
unphased_data <- read_delim( unphased_cn_file , delim="\t", col_names = FALSE )
unphased_data <- dplyr::rename( unphased_data ,"chromosome"="X1","pos_bin"="X2","cov"="X3")
unphased_data <- unphased_data %>% mutate(pos = pos_bin*binsize,hap = "U")
unphased_data <- unphased_data %>% mutate(cov_norm = cov/6200)
unphased_data <- unphased_data %>% filter(chromosome %in% chrom_list)


loh_data1 <- unphased_data %>% filter(chromosome == "chr2") %>% filter(pos > 188000000) %>% mutate(hap = "A")
loh_data2 <- unphased_data %>% filter(chromosome == "chr10") %>% filter(pos > 86000000) %>% mutate(hap = "B")
loh_data3 <- unphased_data %>% filter(chromosome == "chr12") %>% filter(pos < 22000000) %>% mutate(hap = "A")
loh_data4 <- unphased_data %>% filter(chromosome == "chr20") %>% filter(pos < 29000000) %>% mutate(hap = "A")
loh_data5 <- unphased_data %>% filter(chromosome == "chr22") %>% filter(pos > 23500000) %>% mutate(hap = "B")
loh_data6 <- unphased_data %>% filter(chromosome == "chr4") %>% filter(pos > 159500000) %>% filter(pos < 162500000) %>% mutate(hap = "A")
loh_data7 <- unphased_data %>% filter(chromosome == "chr17") %>% filter(pos < 27000000) %>% mutate(hap = "A")
loh_data <- bind_rows(loh_data1,loh_data2,loh_data3,loh_data4,loh_data5,loh_data6,loh_data7)

######################################################
bin_track <- bin_track %>% filter(chromosome == chr_choice)
unphased_data <- unphased_data %>% filter(chromosome == chr_choice)
loh_data <- loh_data %>% filter(chromosome == chr_choice)

#############################################################
chromosome1_max = max(unphased_data$pos)
ylimits <- c(0,6)
#ylimits <- c(0,5000)
#xlimits <- c(0,chromosome1_max)
xlimits <- c(chromosome1_max,0)

log_scale_break <- 2
addition_num <- 1
bin_track <- bin_track %>% mutate( ln2_mean = case_when( cov_norm > log_scale_break ~ log2(cov_norm) + addition_num, cov_norm <= log_scale_break ~ cov_norm ) )
unphased_data <- unphased_data %>% mutate( ln2_mean = case_when( cov_norm > log_scale_break ~ log2(cov_norm) + addition_num, cov_norm <= log_scale_break ~ cov_norm ) )
loh_data <- loh_data %>% mutate( ln2_mean = case_when( cov_norm > log_scale_break ~ log2(cov_norm) + addition_num, cov_norm <= log_scale_break ~ cov_norm ) )

xlab <- c( 0, 30, 60, 90, 120, 150, 180, 210, 240)

bin_track$hap = factor(bin_track$hap, levels=c('A','B','U'))
#bin_track$hap = factor(bin_track$hap, levels=c('U','A','B')) 

loh_data$hap = factor(loh_data$hap, levels=c('A','B','U'))
#loh_data$hap = factor(loh_data$hap, levels=c('U','A','B')) 

unphased_data$hap = factor(unphased_data$hap, levels=c('A','B','U'))
#unphased_data$hap = factor(unphased_data$hap, levels=c('U','A','B')) 

#############################################################
p1 <- ggplot(data = bin_track) + 
  geom_point(aes(x = pos,y = ln2_mean,color = hap),size=0.2) +
  geom_point(data = unphased_data,aes(x = pos,y = ln2_mean,color = hap),size = 0.2) +
  geom_point(data = loh_data,aes(x = pos,y = ln2_mean,color = hap),size = 0.2) +   #scale_color_manual(values=c("grey35", "firebrick2", "dodgerblue3")) +
  scale_color_manual(values=c("firebrick2", "dodgerblue3", "grey35")) +
  facet_grid(.~hap) +   #facet_grid(hap~.) +
  scale_y_continuous(labels = c("0","1","2","4","8","16","32"),breaks = c(0,1,2,3,4,5,6),expand = c(0.0,0.0),limits = ylimits) +
  #scale_x_continuous(limits = xlimits,labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0)) +
  scale_x_reverse(limits = xlimits,labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0)) +
  xlab(chr_choice) +
  ylab("Copy Number") +
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"),legend.position = "none",axis.text.x = element_text(angle = 90),panel.grid.major = element_line(colour="grey", size=0.5)) #,panel.grid.major = element_line(colour="grey", size=0.5)

p1


#############################################################
chromosome1_scale = chromosome1_max/5000000
#output_file <- paste("~/2020_2_february_workdir/K562/phased_cn/phased_cn_",chr_choice,".pdf",sep = "")
output_file <- paste("~/2020_2_february_workdir/K562/phased_cn/phased_cn_rotate_",chr_choice,".pdf",sep = "")
ggsave(output_file, width = chromosome1_scale, height = 5, units = "cm")


























######################################################
#p3 <- ggplot(data = bin_track) + geom_point(aes(x = pos,y = cov_norm,color = haplotype),size = 0.2) +
#geom_point(data = unphased_data,aes(x = pos,y = cov_norm,color = haplotype),size = 0.2) +  #haplotype  #num  #,color = haplotype
#  geom_point(data = loh_data,aes(x = pos,y = cov_norm,color = haplotype),size = 0.2) +
#  scale_color_manual("haplotype",values = c("firebrick2", "royalblue3","grey35")) +
#  theme_classic() + 
#geom_vline(xintercept=180000000) +
#  facet_grid(factor(chromosome,levels = chrom_list)~.) +
#  ylab("Copy Number") + 
#  xlab("Genome Coordinate") +
#  scale_x_continuous(labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0),limits = c(0,290000000)) +
#  scale_y_continuous(expand = c(0.0,0.0),limits = c(0,6)) +
#scale_y_continuous(expand = c(0.0,0.0),limits = c(0,10)) +
#  theme(strip.text.y = element_text(size=8, angle=0, face="bold"),strip.background = element_rect(fill="#ffffff"))


#############################################################
#############################################################
#############################################################
#############################################################

#binned <- jdata %>% group_by(bin) %>% dplyr::summarise(A = mean(hapA_cov), B = mean(hapB_cov), U = mean(tot_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
#binned <- binned %>% mutate(middle_pos = ((min_pos + max_pos)/2.0)) %>% filter(n>5)
#binned_gather <- binned %>% gather(A,B,U,key = "hap",value = "cov")
#binned_gather <- binned_gather %>% mutate(copy_number = cov/5.5)
#binned_gather$hap = factor(binned_gather$hap, levels=c('U','A','B'))   
#binned_gather$hap = factor(binned_gather$hap, levels=c('A','B','U'))

#############################################################
#chromosome1_max = max(bin_track$pos)
#ylimits <- c(0,50)
#ylimits <- c(0,6)
#ylimits <- c(0,5000)
#ylimits <- c(0,3.5)   #c(0,50)  #c(0,180)
#xlimits <- c(0,chromosome1_max)
#xlimits <- c(chromosome1_max,0)
#xlimits <- c(38000000,41000000)

#log_scale_break <- 2
#addition_num <- 1

#binned_gather <- binned_gather %>% mutate( ln2_mean = case_when( cov_norm > log_scale_break ~ log2(cov_norm) + addition_num, cov_norm <= log_scale_break ~ cov_norm ) )

#xlab <- c( 0, 30, 60, 90, 120, 150, 180, 210, 240)
#xlab <- c( 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50)
#xlab <- c( 38,38.5, 39, 39.5, 40, 40.5, 41)







#############################################################
#p1 <- ggplot(data=binned_gather) + geom_point(aes(x = middle_pos,y = ln2_mean,color = hap),size=0.5) +
  #scale_color_manual(values=c("grey35", "firebrick2", "dodgerblue3")) +
#  scale_color_manual(values=c("firebrick2", "dodgerblue3", "grey35")) +
#  facet_grid(.~hap) +
  #facet_grid(hap~.) +
#  scale_y_continuous(labels = c("0","1","2","4","8","16","32"),breaks = c(0,1,2,3,4,5,6),expand = c(0.0,0.0),limits = ylimits) +
  #scale_y_continuous(limits = ylimits,expand = c(0.0,0.0)) +
  #scale_x_continuous(limits = xlimits,labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0)) +
#  scale_x_reverse(limits = xlimits,labels = paste0(xlab, "Mb"),breaks = 10^6*xlab,expand = c(0.0,0.0)) +
#  xlab(chromosome) +
#  ylab("Copy Number") +
#  theme_classic() +
#  theme(panel.spacing = unit(1, "lines"),legend.position = "none",axis.text.x = element_text(angle = 90),panel.grid.major = element_line(colour="grey", size=0.5)) #,panel.grid.major = element_line(colour="grey", size=0.5)


#p1

#############################################################
#chromosome1_scale = chromosome1_max/5000000
#output_file <- paste("~/2019_11_november_workdir/HCC1954/phased_hic_dec/phased_cn_",chromosome,".pdf",sep = "")
#output_file <- paste("~/2019_11_november_workdir/HCC1954/phased_hic_dec/phased_cn_rotate_",chromosome,".pdf",sep = "")
#ggsave(output_file, width = chromosome1_scale, height = 5, units = "cm")
