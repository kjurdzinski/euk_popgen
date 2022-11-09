
library(rworldmap)
library("plotrix")
library(cluster)
library(pheatmap)
library(stringi)
library(spaa)
library(ape)
library(phytools)
library(vegan)
library(randomForest)
library(stringi)
library(lubridate)
library(marmap)

infofile1 <- "../tara_metadata/TARA_SAMPLES_CONTEXT_SEQUENCING_20170515_mod.tab"
infofile2 <- "../tara_metadata/TARA_reg_stations_mod.tab"
#mag_info_file = "~/aquatic/mag2.0/mag_summary_tables/Eukaryotic_MAG_summary.tab"
#eggnog_annotations_file = "~/aquatic/popdiv/annotations/NOG.annotations.txt"
#cog_categories_file = "~/aquatic/popdiv/annotations/cog_categories.txt"
subsample = 10 # level of subsmapling/min_cov used in pogenom

genome = c(
  "A_anophagefferens",
  "B_prasinos",
  "E_huxleyi",
  "G_theta",
  "N_gaditana",
  "O_tauri",
  "P_glacialis",
  "P_tricornutum",
  "T_pseudonana"
)

prefix <- paste("../pogenom_output/mf1_mc10_ss10/", genome, "_mf1_mc10_ss10", sep = "")
#prefix <- paste("~/aquatic/popdiv/euka_pogenom/pogenom_output/mf1_mc10/", genome, "_mf1_mc10", sep = "")

totpifile <- paste(prefix,  "intradiv.txt", sep = ".")
locuspifile <- paste(prefix,  "intradiv-per-locus.txt", sep = ".")
totfstfile <- paste(prefix,  "fst.txt", sep = ".")
allelefile <- paste(prefix,  "allele-freqs.txt", sep = ".")
#pifile <- paste(prefix,  "intradiv-per-gene.txt", sep = ".")
#genefstfile <- paste(prefix,  "fst-per-gene.txt", sep = ".")
#aminofreqfile <- paste(prefix,  "aminoacid-freqs.txt", sep = ".")
#aapifile <- paste(prefix,  "intradiv-aminoacid-per-gene.txt", sep = ".")
#pNpSfile <- paste(prefix,  "pNpS-per-gene.txt", sep = ".")
#permgenefstfile <- paste(prefix,  "permuted-fst-per-gene.txt.gz", sep = ".")
#gff_file <- paste("~/aquatic/popdiv/annotations/gff_files/", gff_file, sep = "")
#eggnog_file <- paste("~/aquatic/popdiv/annotations/emapper_files/", paste(bacl, "_proteins.faa.emapper.annotations.txt", sep=""), sep = "")

###

## read genome info 
#tab = read.delim(mag_info_file)
#ix = match(mag, tab[,1])
#bacl = tab[ix,2]
#tax = tab[ix,3]
#a<-strsplit(as.character(tax), ";")
#res<-as.data.frame(t(stri_list2matrix(a)))
#phylum = res[,2]
#genome_size = tab[ix,4]

####

# get all sample names from totpifiles
tot_sample = c()
for (i in 1:length(prefix)) {
  tab = read.delim(totpifile[i])
  sample =  tab[1:(nrow(tab)-1),1]
  sample = gsub("^X", "", sample)
  tot_sample = c(tot_sample, sample)
}
sample = sort(unique(tot_sample))

# read sample info
tab <- read.delim(infofile1)
ix = match(sample, tab[,21])
station = tab[ix,4]
depth = tab[ix,7]
filt = tab[ix,9]
filt = tab[ix,10]
prefilt = tab[ix,11]
method = tab[ix,18]
tab <- read.delim(infofile2)
ix = match(station, tab[,2])
date_time = tab[ix,3]
lat = tab[ix,5]
lon = tab[ix,6]
dist_to_coast = tab[ix,14]
region = as.character(tab[ix,16])
station = as.character(station)
station = gsub("^TARA_", "", station)
ix = grep("P21817", sample)
region[ix] = "Baltic Sea"
#lat[ix] = 
#lon[ix] = 


# sample groups
ms = grep("(MS)", region)
nao = grep("(NAO)", region)
sao = grep("(SAO)", region)
npo = grep("(NPO)", region)
spo = grep("(SPO)", region)
so = grep("(SO)", region)
io = grep("(IO)", region)
bs = grep("Baltic Sea", region)

region[ms] = "MS"
region[nao] = "NAO"
region[sao] = "SAO"
region[npo] = "NPO"
region[spo] = "SPO"
region[so] = "SO"
region[io] = "IO"
region[bs] = "BS"

# colors & shapes
col_region[ms] = "#fdbf6f"
col_region[nao] = "#a6cee3"
col_region[sao] = "#1f78b4"
col_region[npo] = "#b2df8a"
col_region[spo] = "#33a02c"
col_region[so] = "#fb9a99"
col_region[io] = "#e31a1c"
col_region[bs] = "black"

##############################
#### reading pogenom data ####

# read tot nt-pi
totpi_matr = list()
pi_matr = matrix(ncol = length(sample), nrow = length(prefix))
colnames(pi_matr) = sample
rownames(pi_matr) = genome
for (i in 1:length(prefix)) {
  tab = read.delim(totpifile[i])
  tab = tab[1:(nrow(tab)-1),]
  tab[,1] = gsub("^X", "", tab[,1])
  ix = match(tab[,1], sample)
  pi_matr[i,ix] = tab[,3] # for using normalised pi, otherwise tab[,2]
  totpi_matr[[i]] = tab
}

# read per locus nt-pi
locuspi_matr = list()
for (i in 1:length(prefix)) {
  tab = read.delim(locuspifile[i])
  rownames(tab) = paste(tab[,1],tab[,2], sep = "|")
  tab = tab[,3:(ncol(tab) - 2)]
  tab = tab[,seq(from=1, to=ncol(tab), by=2)]
  colnames(tab) = gsub(".Intra_locus_pi", "", colnames(tab))
  locuspi_matr[[i]] = tab
}

## read per gene nt-pi
#pi_gene_matr = pi_gene = pi_contig = pi_gene_num_loci = list()
#for (i in 1:length(prefix)) {
#  tab <- read.delim(pifile[i])
#  pi_gene[[i]] <- tab[,2]
#  pi_gene[[i]] = gsub("ID=PROKKA_MOD_", "", pi_gene[[i]])
#  pi_contig[[i]] <- tab[,1]
#  pi_gene_num_loci[[i]] <- tab[,7]
#  pi_gene_matr[[i]] = as.matrix(tab[,8:(ncol(tab)-1)])
#  colnames(pi_gene_matr[[i]]) = gsub("^X", "", colnames(pi_gene_matr[[i]]))
#  colnames(pi_gene_matr[[i]]) = gsub(".Intra_gene_pi", "", colnames(pi_gene_matr[[i]]))
#}

# reading fst data
fst_matr = list()
for (i in 1:length(prefix)) {
  tab <- read.delim(totfstfile[i], row.names = 1)
  fst_matr[[i]] = tab
  colnames(fst_matr[[i]]) = gsub("X", "", colnames(fst_matr[[i]]))
}

# reading allele data
allele_dist = list()
unique_alleles = list()
tot_alleles = list()
major_allele_counts = list()
minor_allele_counts = list()
major_allele_norm_counts = list()
major_allele = list()
minor_allele = list()
for (i in 1:length(prefix)) {
  tab = read.delim(allelefile[i])
  contig = tab[,1]
  pos = tab[,2]
  locus = paste(contig, pos)
  allele_counts = as.matrix(tab[,3:(ncol(tab)-4)])
  s = colnames(allele_counts)
  s = s[seq(from = 1, to = length(s), by = 4)]
  s = gsub("X", "", s)
  s = gsub(".A", "", s)
  A_counts = allele_counts[,seq(1, ncol(allele_counts), 4)]
  T_counts = allele_counts[,seq(2, ncol(allele_counts), 4)]
  C_counts = allele_counts[,seq(3, ncol(allele_counts), 4)]
  G_counts = allele_counts[,seq(4, ncol(allele_counts), 4)]
  A_norm_counts = A_counts
  T_norm_counts = T_counts
  C_norm_counts = C_counts
  G_norm_counts = G_counts
  for (j in 1:ncol(A_counts)) {
    sum = A_counts[,j] + T_counts[,j] + C_counts[,j] + G_counts[,j]
    A_norm_counts[,j] = A_counts[,j] / sum
    T_norm_counts[,j] = T_counts[,j] / sum
    C_norm_counts[,j] = C_counts[,j] / sum
    G_norm_counts[,j] = G_counts[,j] / sum
  }
  Max_counts = Max_norm_counts = matrix(ncol = ncol(A_counts), nrow = nrow(A_counts))
  Max_allele = Min_allele = rep(NA, nrow(A_counts))
  for (j in 1:nrow(A_counts)) {
    list1 = list(A_norm_counts[j,], T_norm_counts[j,], C_norm_counts[j,], G_norm_counts[j,])
    list2 = list(A_counts[j,], T_counts[j,], C_counts[j,], G_counts[j,])
    ix = sort(unlist(lapply(list1, mean, na.rm=T)), index.return=T, decreasing=T)$ix
    Max_norm_counts[j,] = list1[[ix[1]]]
    Max_counts[j,] = list2[[ix[1]]]
    Max_allele[j] = c("A","T","C","G")[ix[1]]
    Min_allele[j] = c("A","T","C","G")[ix[2]]
  }
  major_allele_counts[[i]] = Max_counts
  major_allele_norm_counts[[i]] = Max_norm_counts
  major_allele[[i]] = Max_allele
  minor_allele[[i]] = Min_allele
  rownames(major_allele_counts[[i]]) = paste(locus, " ", major_allele[[i]], "/", minor_allele[[i]], sep="")
  colnames(major_allele_counts[[i]]) = s
  colnames(major_allele_norm_counts[[i]]) = s
  allele_dist[[i]] = as.matrix(dist(t(Max_norm_counts), method = "euclidean"))
  rownames(allele_dist[[i]]) = colnames(allele_dist[[i]]) = s
  unique_alleles[[i]] = matrix(ncol = length(s), nrow = length(s))
  rownames(unique_alleles[[i]]) = colnames(unique_alleles[[i]]) = s
  tot_alleles[[i]] = unique_alleles[[i]]
  for (j in 1:length(s)) {
    inj = which(allele_counts[,((j*4)-3):(j*4)] > 0)
    if (j < length(s)) {
      for (k in (j+1):length(s)) {
        ink = which(allele_counts[,((k*4)-3):(k*4)] > 0)
        unique_alleles[[i]][j,k] = unique_alleles[[i]][k,j] = 1 - (length(intersect(inj, ink)) / length(union(inj, ink)))
        tot_alleles[[i]][j,k] = tot_alleles[[i]][k,j] = length(union(inj, ink))
      }
    }
  }
}


#############################
### finished reading data ###
#############################

### num of loci ###

par(mfrow = c(3,3), mar = c(2,2,1,1))
for (i in 1:length(genome)) {
  tab = locuspi_matr[[i]]
  tab[!is.na(tab)] = 1
  tab[is.na(tab)] = 0
  hist(apply(tab, 1, sum), breaks = 100, main = genome[i], cex.main = 1.0, xlim = c(0, ncol(tab)))
  #pheatmap(tab, show_rownames = F, main = mag[i])
}

par(mfrow = c(3,3), mar = c(2,2,1,1))
for (i in 1:length(genome)) {
  tab = major_allele_counts[[i]]
  these_sample = colnames(tab)
  these_region = region[match(colnames(tab), sample)]
  these_col = col_region[match(colnames(tab), sample)]
  tab[!is.na(tab)] = 1
  tab[is.na(tab)] = 0
  tab2 = major_allele_counts[[i]]
  ranked_loci = sort(apply(tab, 1, sum), decreasing = T, index.return = T)$ix
  ix = which(apply(tab, 1, sum) > 0)
  pheatmap(tab2[ranked_loci[ix],], show_rownames = F, main = genome[i], cluster_cols = F, cluster_rows = F, labels_col = these_region)
}

par(mfrow = c(3,3), mar = c(2,2,1,1))
for (i in 1:length(genome)) {
  these_sample = colnames(tab)
  tab = major_allele_counts[[i]]
  tab2 = major_allele_counts[[i]]
  tab[!is.na(tab)] = 1
  tab[is.na(tab)] = 0
  ranked_samples = sort(apply(tab, 2, sum), decreasing = T, index.return = T)$ix
  matr = matrix(ncol = 2, nrow = length(ranked_samples))
  ix = which(tab[,ranked_samples[1]] == 1) # loci present in top-rakning sample (sample with most loci)
  matr[1,] = c(1, length(ix))
  these_region = region[match(colnames(tab), sample)]
  these_col = col_region[match(colnames(tab), sample)]
  for (j in 2:length(ranked_samples)) {
    ix = which(apply(tab[,ranked_samples[1:j]], 1, sum) == j)
    matr[j,] = c(j, length(ix))
  }
  ix1 = which(matr[,2] >= 100) # the samples that share >= 100 loci with the sample with most loci 
  if (length(ix1) >= 4) {
    ix2 = which(apply(tab[,ranked_samples[ix1]], 1, sum) == length(ix1))
    #pheatmap(tab2[ix2,ranked_samples[ix1]], show_rownames = F, main = genome[i], cluster_cols = T, labels_col = these_region[ranked_samples[ix1]])
    pcoa <- pcoa(dist(t(tab2[ix2,ranked_samples[ix1]])), correction = "cailliez")
    plot(pcoa$vectors[,1], pcoa$vectors[,2], cex = 2, pch = 21, col = these_col[ranked_samples[ix1]], main = genome[i])
  }
}



i = 3
tab = major_allele_counts[[i]]
tab[!is.na(tab)] = 1
tab[is.na(tab)] = 0
ix = which(apply(tab, 1, sum) >= 15)
tab = major_allele_counts[[i]]
tab[is.na(tab)] = -10
pheatmap(tab[,], show_rownames = F, main = genome[i], cluster_cols = T)
pheatmap(tab[ix,], show_rownames = F, main = genome[i], cluster_cols = T)

i = 2
tab = major_allele_counts[[i]]
tab[!is.na(tab)] = 1
tab[is.na(tab)] = 0
ix = which(apply(tab, 1, sum) >= 6)
pheatmap(tab[,], show_rownames = F, main = mag[i], labels_col = region[])
pheatmap(tab[ix,], show_rownames = F, main = mag[i], labels_col = region[ix])
tab = major_allele_counts[[i]]
tab[is.na(tab)] = -10
ix2 = c(1,4,6,8,9,12)
ix = which(apply(tab[,ix2], 1, min) > -10)
pheatmap(tab[ix,ix2], show_rownames = F, main = paste(mag[i], bacl[i]), cluster_cols = T)

i = 3
tab = major_allele_counts[[i]]
tab[is.na(tab)] = -10
pheatmap(tab[], show_rownames = F, main = genome[i], cluster_cols = T)

#pdf(file = "major_allele_clust.pdf")
for (i in c(1,3:5,7:14)) {
  tab = major_allele_counts[[i]]
  tab[is.na(tab)] = -10
  ix = match(colnames(tab), sample)
  pheatmap(tab, show_rownames = F, main = paste(mag[i], bacl[i], sep = "|"), cluster_cols = T, labels_col = sample[ix])
}
#dev.off()


i = 3
tab = major_allele_counts[[i]]
tab[is.na(tab)] = -10
ix = match(colnames(tab), sample)
#pheatmap(tab, show_rownames = F, main = genome[i], cluster_cols = T, labels_col = region[ix])
ix2 = hm$tree_col$order[2:16]
ix3 = which(apply(tab[,ix2], 1, min) > -10) # i.e. which are present in all ix2 samples
pheatmap(tab[ix3,ix2], show_rownames = F, main = genome[i], cluster_cols = T, labels_col = sample[ix])
pheatmap(tab[ix3,ix2], show_rownames = F, main = genome[i], cluster_cols = T, labels_col = region[ix[ix2]])
pheatmap(tab[ix3,ix2], show_rownames = F, main = genome[i], cluster_cols = T, labels_col = station[ix[ix2]])

ix4 = match(sample[ix[ix2]], colnames(fst_matr[[1]]))
plot(as.dist(fst_matr[[1]][ix4,ix4]), dist(t(tab[ix3,ix2])))  # cor fst vs allele count euclidean dist

pcoa <- pcoa(dist(t(tab[ix3,ix2])), correction = "cailliez")
#pcoa <- pcoa(as.dist(fst_matr[[1]][ix4,ix4]), correction = "cailliez")
pcoa$values[c(1,2),3]
par(mfrow=c(1,1), mar=c(2,2,1,1), xpd = TRUE)
xlab = paste("PC1 (", 100*round(pcoa$values[1,3], 2), "%)", sep = "")
ylab = paste("PC2 (", 100*round(pcoa$values[2,3], 2), "%)", sep = "")
plot(pcoa$vectors[,1], pcoa$vectors[,2], cex = 2, pch = 21, col = col_region[ix[ix2]], bg = , xlab = xlab, ylab = ylab)
text(pcoa$vectors[,1], pcoa$vectors[,2], labels = station[ix[ix2]])

par(mfrow=c(1,1), mar=c(2,2,1,1), xpd = TRUE)
mapdata <- getNOAA.bathy(lat1 = 10, lat2 = 80, lon1 = -100, lon2 = 30, res=20, keep=TRUE)
blues <- colorRampPalette(c("blue","cadetblue1","white"), bias = 0.5)
#plot(mapdata, image = TRUE, bpal = blues(100), asp = 1.48, # to paint depth in blue scale
plot(mapdata, image = TRUE, bpal = NA, asp = 1.48,
     deep = c(-10000, -100, 0), 
     shallow = c(-100, -10, 0),
     step = c(50, 50, 0),
     lwd = c(0.8, 0.8, 1), lty = c(0, 0, 1),
     #col = c("lightgrey", "darkgrey", "black"),
     drawlabel = c(FALSE, FALSE, FALSE),)
text(lon[ix[ix2]], lat[ix[ix2]], cex = 1, col = col_region[ix[ix2]], labels = station[ix[ix2]])

par(mfrow=c(1,1), mar=c(2,2,1,1), xpd = TRUE)
mapdata <- getNOAA.bathy(lat1 = 24, lat2 = 46, lon1 = -10, lon2 = 36, res=5, keep=TRUE)
blues <- colorRampPalette(c("blue","cadetblue1","white"), bias = 0.5)
#plot(mapdata, image = TRUE, bpal = blues(100), asp = 1.48, # to paint depth in blue scale
plot(mapdata, image = TRUE, bpal = NA, asp = 1.48,
     deep = c(-10000, -100, 0), 
     shallow = c(-100, -10, 0),
     step = c(50, 50, 0),
     lwd = c(0.8, 0.8, 1), lty = c(0, 0, 1),
     #col = c("lightgrey", "darkgrey", "black"),
     drawlabel = c(FALSE, FALSE, FALSE),)
text(lon[ix[ix2]], lat[ix[ix2]], cex = 1, col = col_region[ix[ix2]], labels = station[ix[ix2]])
