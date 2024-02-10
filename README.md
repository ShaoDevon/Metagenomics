# Metagenomics

The goal of this project is to perform metagenomics analysis on a data sample. 

After taking a sample from a live subject, such as a swab from a human, there are many possible genomes present: those from the human host, but also from the various bacterial and even viral organisms present. The sample is sequenced into reads, and while the scientist, researcher or doctor has a set of reference genomes for organisms, they do not know which organisms are present in the sample, let alone which reads map to which organism's genome. 

In order to solve this problem, this project builds a hashtable from each genome using 16-mers of the genome. It then splits each 50-length read into three 16-length sections and checks where in the genome these 16-length sections match. Using Smith-Waterman algorithm, it finds the best genome match for each read. 
