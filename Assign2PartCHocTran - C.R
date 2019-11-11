#Primary Author: Hoc Tran
#Editor: Braden Judson (B)
#BINF6210 - Assignment 2 Part C

#I chose to ask the question of whether 18S and 16S genes within Daphnia generate the same or different phylogenetic hypotheses? I thought that this would be an interesting question to explore as I had already looked at the phylogenetic hypotheses of the 16S genes and wanted to compare what I saw visually with another gene, such as 18S. This also led me to ask another question, that if they do differ, is there some factor that determines which phylogenetic hypotheses is more correct? If such a factor exists, how must the dendograms differ visually to receive different scores? Would it be apparent visually? I also wanted to know if there was some way to visually compare two dendograms, as looking at them seperately would not reveal much about each other. I decided to answer this question by generating a dataframe that contained 2 sequences per species, one 18S and one 16S nucleotide sequence and creating dendograms for each sequence individually. I made several checks to the data before generating the dendograms and removed any sequences that I thought to be outliers, either by looking at sequence length compared to the mean sequence length or by looking at them after alignment. Any sequences that had missing data in the gene marker column were also removed as it was critical to know where each sequence belonged for generating the dendograms. The rest of the dataset was found to not have any more missing data after this removal, as the remaining number of samples was not very large.  From there, I was able to compare the two dendograms by using the tangogram function, which created a tangogram for a visual comparison.

#B: I would elaborate on the function of both 16 and 18S genes, e.g. both are mitochondrial rRNA coding genes. Is there a reason they should differ?

library(rentrez)
library(seqinr)
library(Biostrings)
library(stringr)
library(tidyverse)
library(stringi)
library(ape)
library(muscle)
library(msa)
library(DECIPHER)
library(dendextend)

#Getting 16S Daphnia DATA----
daphnia16S_search <- entrez_search(db = "nuccore", term = "(Daphnia[ORGN] AND 16S[TITL]) NOT (genome[TITL])", retmax = 1000)#Searching NCBI 
daphnia16S_search #Check search hits
daphnia16S_fetch <- entrez_fetch(db = "nuccore", id = daphnia16S_search$ids, rettype = "fasta") #Gettting it as FASTA format
write(daphnia16S_fetch, "daphnia_fetch16S.fasta", sep = "\n") #Writing to disk

DaphniaStringSet <- readDNAStringSet("daphnia_fetch16S.fasta")#Convert to DNA stringset

Daphnia16S <- data.frame(Sequence_Labels = names(DaphniaStringSet), Daphnia_Nucleotide_Sequence = paste(DaphniaStringSet))#Making a dataframe

Daphnia16S$unique_identifier <- word(Daphnia16S$Sequence_Labels, 1L,) #Getting identifiers from first "word" of label     #B: I believe this is called the accession number. 
Daphnia16S$species_name <- str_extract(Daphnia16S$Sequence_Labels, "Daphnia [a-zA-Z]+\\.* |Daphniopsis [a-zA-Z]+\\.*") #Getting species name from label, Daphniopsis present in this dataset as well, 
Daphnia16S$gene_name <- str_extract(Daphnia16S$Sequence_Labels, "16S [a-zA-Z]+\\.*") #Getting the gene name from the label
Daphnia16S <- Daphnia16S[, c("unique_identifier", "species_name", "gene_name", "Daphnia_Nucleotide_Sequence", "Sequence_Labels")]#Organizing columns

#B: I would recommend some data checking here! You've collected species names from string information, but maybe some aren't what you want?
unique(Daphnia16S$species_name)
#B: All appear valid with the exception of "Daphnia sp.", which is not super informative. I would remove entries of unknown species as your tanglegram may say two Daphnia sp. have dissimilar 16S and/or 18S gene sequence, but that may be expected if they are actually different species. 
Daphnia16S <- Daphnia16S %>% 
  filter(species_name != "Daphnia sp. ")
#B: Again, check to see it works.
unique(Daphnia16S$species_name) #Looks good.


Daphnia.16S.per.species <- Daphnia16S %>% #Getting only one sequence from each species
  group_by(species_name) %>%
  filter(!is.na(species_name))%>% #Remove NAs #B: You can also use the function na.omit
  sample_n(1)

#Getting 18S Daphnia DATA----
Daphnia18S_search <- entrez_search(db = "nuccore", term = "(Daphnia[ORGN] AND 18S[TITL]) NOT (genome[TITL])", retmax = 1000)#Searching NCBI
Daphnia18S_search #Check search hits
Daphnia18S_fetch <- entrez_fetch(db = "nuccore", id = Daphnia18S_search$ids, rettype = "fasta") #Getting it as FASTA format
write(Daphnia18S_fetch, "Daphnia18S_fetch.fasta", sep = "\n") #Writing to disk
Daphnia18SstringSet <- readDNAStringSet("Daphnia18S_fetch.fasta") #Convert to DNA Stringset
Daphnia18Sdf <- data.frame(Sequence_Labels18S = names(Daphnia18SstringSet), Daphnia18S_Nucleotide_Sequence = paste(Daphnia18SstringSet))#Making a dataframe


Daphnia18Sdf$unique_identifier <- word(Daphnia18Sdf$Sequence_Labels18S, 1L,)#Getting unique identifiers from label
Daphnia18Sdf$species_name <- str_extract(Daphnia18Sdf$Sequence_Labels18S, "Daphnia [a-zA-Z]+\\.* |Daphniopsis [a-zA-Z]+\\.*") #Getting species name from label, Daphniopsis present in this dataset as well
Daphnia18Sdf$gene_name <- str_extract(Daphnia18Sdf$Sequence_Labels18S, "18S [a-zA-Z]+\\.*") #Getting gene name from label
Daphnia18Sdf <- Daphnia18Sdf[, c("unique_identifier", "species_name", "gene_name", "Daphnia18S_Nucleotide_Sequence", "Sequence_Labels18S")] #Reorganizing columns

#B: Similar comments as above (line 36-42).
unique(Daphnia18Sdf$species_name)
#B: Here, you have NAs and "Daphnia cf. ", which will only complicate downstream processing and thus should be removed (I see you remove NAs on Line 73). 
Daphnia18Sdf <- Daphnia18Sdf %>% 
  filter(species_name != "Daphnia cf. ")
unique(Daphnia18Sdf$species_name) #Looks good!

Daphnia18Sdf.by.species <- Daphnia18Sdf %>% #Getting only one sequence from each species
  group_by(species_name) %>%
  filter(!is.na(species_name))%>%#Remove NAs
  sample_n(1)

#B: I would also recommend checking data here. For example, are the number of rows (nrow) and number of unique species the same (length(unique(species_name)))? Should be. 

#Combining the two datasets so that each species will have both a 16S sequence and an 18S sequence
#DaphniaOverlap <-merge(Daphnia.16S.per.species, Daphnia18Sdf, by = "unique_identifier", all = F) Couldn't merge by unique_identifiers, so no same individuals with both 16S and 18S data
DaphniaOverlap <-merge(Daphnia.16S.per.species, Daphnia18Sdf.by.species, by = "species_name", all = F)#Used species name instead
#B: Change "Daphnia_Nucleotide_Sequence" to a more informative name using dplyr::rename. 
library(dplyr)
DaphniaOverlap <- DaphniaOverlap %>% 
  rename(Daphnia_16S_Nucleotide_Sequence = Daphnia_Nucleotide_Sequence) #B: Clearer formatting. You could also probably remove the gene_name.x and gene_name.y columns as they are not informative (unless you had a single column called "Gene_Name" that had both genes).

unique(DaphniaOverlap$species_name)#Checks unique species
length(unique(DaphniaOverlap$species_name)) #16
sum(is.na(DaphniaOverlap$Daphnia_Nucleotide_Sequence)) #0
sum(is.na(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #0

#Checks for 18S Nucleotide Sequences----
mean(str_count(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #2024.75, very large
min(str_count(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #435
max(str_count(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #6705, very very big
median(str_count(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #605
DaphniaOverlap$Daphnia18S_Nucleotide_Sequence<- as.character(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)#Change 18S nucleotide sequence to character for use
DaphniaOverlap$Daphnia_Nucleotide_Sequence<- as.character(DaphniaOverlap$Daphnia_16S_Nucleotide_Sequence)#Change 16S nucleotide to character for use
class(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)

#B: For lines 92-95 I get different (but kind of close) values for each.

hist((nchar(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)), xlab = "Nucleotide Sequence Length", main = "Gene Length Frequency for 16S and 18S")
#Histogram for 18S, majority of data is between 0 and 1000, 1 between 2000 and 3000, 3 between 5000 and 6000, 1 between 6000 and 7000.
#B: Are the genes simply different lengths? Needs elaboration. Also, title and x-axis names are unclear.Probably would have made more sense to make a histogram of gene length for EACH gene to check for outliers. 

DaphniaOverlap$Daphnia18S_Nucleotide_Sequence <- DNAStringSet(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence) #Convert to DNA stringset for alignment

Daphnia.18S.alignment <- DNAStringSet(muscle::muscle(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence, maxiters = 2,)) #Run alignment with potential outliers
BrowseSeqs(Daphnia.18S.alignment, htmlFile = paste(tempdir(), "/Daphnia18.S.Alignment.html", sep = "")) #Alignment not very good, very long sequences causes big gaps

DaphniaOverlap$Daphnia18S_Nucleotide_Sequence<- as.character(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)#Convert to character to remove outliers

DaphniaOverlap <- DaphniaOverlap %>% #Remove outliers with lenght greater than 1000
  filter((nchar(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence) < 1000))

hist(nchar(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence)) #Histogram with outliers removed, looks alot better
DaphniaOverlap$Daphnia18S_Nucleotide_Sequence <- DNAStringSet(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence) #Convert to DNA Stringset to run alignment again

DaphniaOverlap.18S.alignment <- DNAStringSet(muscle::muscle(DaphniaOverlap$Daphnia18S_Nucleotide_Sequence, maxiters = 2,)) #Run alignment with outliers removed
BrowseSeqs(DaphniaOverlap.18S.alignment, htmlFile = paste(tempdir(), "/Daphnia18Overlap.S.Alignment.html", sep = "")) #Alignment looks a lot better
#B: Still some pretty gappy regions. I would suggest playing around the gap-penalty to see if that improves it. But yes, removing the >1000bp samples did help a lot! The gappy regions may influence your downstream "distance" measurements between samples.


#Generating a dendogram for 18S----
dnaBin.DaphniaOverlap.18S<- as.DNAbin(DaphniaOverlap.18S.alignment) #Convert to DNABin to use in clustering

distanceMatrix18SOverlap <- dist.dna(dnaBin.DaphniaOverlap.18S, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE) #Generating distance matrix
#Clustering and generating a dendogram
clusters.DaphniaOverlap.18S <- IdClusters(distanceMatrix18SOverlap, method = "single", cutoff= 0.02, showPlot = TRUE, type = "both", verbose = TRUE) 
#B: Keep in mind your ouput here is in a list format, and requires indexing to withdraw the dendrogram itself.
Daphnia18SDendrogram <- clusters.DaphniaOverlap.18S[[2]]
class(Daphnia18SDendrogram) #B: Dendrogram!


#B: You need to comment on your dendrogram! The y-axis goes to 0.3 - is that expected for an rRNA gene? Seems a bit high. This may be due to sequence gaps above. Also the branch labels are unclear (this will be addressed later). 


#Checks for 16S----
mean(str_count(DaphniaOverlap$Daphnia_Nucleotide_Sequence)) #475.94
min(str_count(DaphniaOverlap$Daphnia_Nucleotide_Sequence)) #448
max(str_count(DaphniaOverlap$Daphnia_Nucleotide_Sequence)) #492
median(str_count(DaphniaOverlap$Daphnia_Nucleotide_Sequence)) #482
DaphniaOverlap$Daphnia_Nucleotide_Sequence<- as.character(DaphniaOverlap$Daphnia_Nucleotide_Sequence)#Change 16S nucleotide to character for use
class(DaphniaOverlap$Daphnia_Nucleotide_Sequence)#Check that it is character

#B: For lines 148-151 I get different values.

hist(nchar(DaphniaOverlap$Daphnia_Nucleotide_Sequence))#Histogram, no outliers appear to be present unlike 18S
DaphniaOverlap$Daphnia_Nucleotide_Sequence<- DNAStringSet(DaphniaOverlap$Daphnia_Nucleotide_Sequence)#Convert to DNA string to run alignment

DaphniaOverlap.16S.alignment <- DNAStringSet(muscle::muscle(DaphniaOverlap$Daphnia_Nucleotide_Sequence, maxiters = 2,)) #Run alignment
BrowseSeqs(DaphniaOverlap.16S.alignment, htmlFile = paste(tempdir(), "/Daphnia16Overlap.S.Alignment.html", sep = "")) #Alignment looks good, nothing too strange  #B: Looks decent

#Generating a dendogram for 16S----
dnaBin.DaphniaOverlap.16S<- as.DNAbin(DaphniaOverlap.16S.alignment)#Converting to DNAbin

distanceMatrix16SOverlap <- dist.dna(dnaBin.DaphniaOverlap.16S, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)#Generating a distance matrix

#Clustering and generating a dendogram
clusters.DaphniaOverlap.16S <- IdClusters(distanceMatrix16SOverlap, method = "single", cutoff= 0.02, showPlot = TRUE, type = "both", verbose = TRUE) 

#Seperating cluster data and dendogram data for 18S----
DaphniaOverlap18S_Dendogram <- clusters.DaphniaOverlap.18S[[2]]
class(DaphniaOverlap18S_Dendogram) #Checking that it is a dendogram



#Seperating cluster data and dendogram data for 16S----
DaphniaOverlap16S_Dendogram <- clusters.DaphniaOverlap.16S[[2]]
class(DaphniaOverlap16S_Dendogram)


#B: Also it's important to name the nodes of your dendrogram branches. Currently they are just numbers, but that doesn't mean much. To change them to your species names:
rownames(DaphniaOverlap) <- c(DaphniaOverlap$species_name)
#B: Then, use the package dendextend to change the node labels of your dendrogram.
dendextend::labels(Daphnia18SDendrogram) <- c(rownames(DaphniaOverlap))
class(clusters.DaphniaOverlap.18S[[2]])
plot(Daphnia18SDendrogram) #B: Now species names are correctly assigned to the branches! But they're too big. 
Daphnia18SDendrogram <- set(Daphnia18SDendrogram, "labels_cex", 0.65)
plot(Daphnia18SDendrogram) #B: Readable but not cuttoff! 
#B: If this was going to be the end of your analysis I would recommend naming your plot.

#B: Change your branch tip labels to something meaningful (species names in this case). I already did this above for the 18S dataset. To do this, rename the original file row names.
rownames(DaphniaOverlap) <- c(DaphniaOverlap$species_name)

dendextend::labels(DaphniaOverlap18S_Dendogram) <- c(rownames(DaphniaOverlap)) 
plot(DaphniaOverlap16S_Dendogram) #B: Species as branch tips! Text is too big, but we will address that when constructing the tanglegram. 
dendextend::labels(DaphniaOverlap16S_Dendogram) <- c(rownames(DaphniaOverlap))
plot(DaphniaOverlap16S_Dendogram) #B: Similar comments to above.

#Generating a tanglegram composed of the 18S dendogram and 16S denogram in Daphnia for comparison----
tanglegram (DaphniaOverlap18S_Dendogram, DaphniaOverlap16S_Dendogram)

#Looking at edges that are unique to each tree
dend_diff(DaphniaOverlap18S_Dendogram, DaphniaOverlap16S_Dendogram) #B: I'm not really sure what this function conributes to your hypothesis. Seems like the opposite of a tanglegram almost? 

#B: Tangelgram suggests equally accurate assignment with both 16S and 18S rRNA genes (no crossing lines!). You can also calculate how well the dendrograms agree by using the "entanglement" function.
entanglement(tanglegram(DaphniaOverlap16S_Dendogram, DaphniaOverlap18S_Dendogram)) 
#B: Yields the value 0, which is expected after observing total concordance between dendrograms. 
#B: Tanglegram is difficult to read as the text overlaps the lines in the middle and it isn't immediately clear which side of the tanglegram is which gene. To clean it up, convert it to a dendlist. 
Tanglegram <- dendextend::dendlist(DaphniaOverlap16S_Dendogram, DaphniaOverlap18S_Dendogram) %>% 
  tanglegram(main_left="16S Gene", main_right="18S Gene", bottom_margin = 6, columns_width = c(5,1,5), margin_inner = 11, margin_outer = 2)
#B: Looks way better! Can also change lines to a single colour if you prefer.
#B: I get a different output than what is described below (line 206). 

#The tanglegram produced from the two dendograms of 18S and 16S in Daphnia showed that there were a few unique nodes between the two. The 18S dendogram has less unique nodes than the 16S dendogram. This is most likely a result of individuals randomly being selected for each species, with the two sequences belonging to the same species but from different individuals. There are also two blue lines that connect the two dendograms, which indicates that these two sub trees were present in both dendograms. The evolutionary distance of the two dendograms are slightly different as well, with 18S ranging from 0.0 - 0.4 and 16S ranging from 0.0-0.07. This indicates that 16S may have more genetic variation between species than that of 18S. Further work that I would perform on these dendograms would be to look at the entanglement between them. This would allow me to determine how well these two dendograms are aligned with each other within this tangogram, and if it is not aligned well, I would look to change the configuration to one that does has higher alignment and then analyze the data from there.


#Braden - Random comments:
#Good script commenting! Well done. 
#You don't check the class of your data often - highly recommend doing this after most conversion-type steps.

