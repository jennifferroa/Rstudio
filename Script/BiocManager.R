if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
#Bioconductor#
#1.upload seqinr
library("seqinr")
help("seqinr-package")
#2. upload the files
choose.files()
dengue<-read.fasta(file = "C:\\Users\\USER\\Documents\\NCBI\\den1.fasta.txt")
#see what type of file is
class(dengue)
#we can access the element of R list object using double square brackets
dengueseq <- dengue[[1]]
#Now we have a vector 
class(dengueseq) #seqfastadna
denq<-dengueseq[1:50]

#Now we can obtain some simple statistics to describe that sequence.
#I can know the sequences' length
length(dengueseq) # print 10735
length(denq) # print 50
# I can count the number of occurrences of the four different nucleotides
#in the sequence.
table(dengueseq)# a(3426), c(2240), g(2770), t(2299)
#I can calculating the GC content of DNA
GC(dengueseq)
#DNA words, is interesting to know the frequency of longer DNA "words"
#calculating 2 words or three or four as a long.
count(dengueseq, 2)
count(dengueseq,3)
count(dengueseq,4)
help(count)# this returns a table whose dimnames are all possible oligomers.
#can use double brackets to extract the values of elements from table.
denguetable<-count(dengueseq,1)# a table of one monomer
print(denguetable) #a(3426), b(2240), g(2770), t(2299)
#to extract the values of elements from tables.
denguetable[[3]] #g(2770)
denguetable[["g"]]#g(2770)
########################################################################
library(seqinr)
dengue<-read.fasta(file = "C:\\Users\\USER\\Documents\\NCBI\\den1.fasta.txt")
dengueseq<-dengue[[1]]
length(dengueseq)
dengueseq[10716:10735]
########################################################################
help.search("complement")#this function complement a nucleic acid sequence.
comp<-comp(dengueseq)
table(comp)#a(2299), c(2770),g(2240), t(3426)
table(dengueseq)#a(3426), c(2240), g(2770), t(2299)
#########################################################################
# how to understand a loop.
for (i in 1:10) {print(i*i)
  
}

for (i in 2:20) {print(i-i)
  
}

Vector<-seq(1:10)
j<-seq(1, 10, by=2)

for (i in 1:4) {print(i*i)
  
}

########################################################################
myfunction <- function(x) {return(20+(x*x))
  
}
myfunction(10)
myfunction(25)

funcion <- function(x) {print(20+(x*x))
  
}
funcion(5)
#########################################################################
dengueseq[452:535]
GC(dengueseq)
#Local fluctuations in GC content within the genome sequence, can reveal
#cases of horizontal transfer or biases in mutation
GC(dengueseq[1:2000])#0.465
GC(dengueseq[2001:4000])#0.4525
GC(dengueseq[4001:6000])#0.4705
GC(dengueseq[6001:8000])#0.479
GC(dengueseq[8001:10000])#0.4545
GC(dengueseq[10001:10735])#0,4993197

#using a for loop I can do the same calculations
#1. I will create a vector with sequence and length
starts<-seq(1, length(dengueseq)-2000, by = 2000)
#2. I will create a new vector to find the length of the vector "starts"
n<-length(starts)
#we set the variable n to be equal to the number of elements in the 
#vector starts,
#now we going to use for loop
#the line for (i in 1:n) means that i will take values of 1-5 in cycles
#chunk vector variable to store the region from nucleotides 1-2000
#chunkGC variable calculate and store  and printed out.


