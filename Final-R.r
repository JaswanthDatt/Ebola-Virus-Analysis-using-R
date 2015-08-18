
#####reading sequences

library(Biostrings)

####Align CDS with Zaire 1976


Seq1=readDNAStringSet("C:/Users/jasda_000/Desktop/Project-1-BIO/KR105238.1/CDS1_NP.txt")

####Align complete sequence with Sudan 1976

Seq2=readDNAStringSet("C:/Users/jasda_000/Desktop/Project-1-BIO/sequence_Sudan.fasta")

####sequence for Zaire 1976

Seq3=readDNAStringSet("C:/Users/jasda_000/Desktop/Project-1-BIO/sequence_Zaire.fasta")


mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)

output_sudan = pairwiseAlignment(pattern=Seq1, subject=Seq2, type = "local", substitutionMatrix = mat,
                      gapOpening = -2, gapExtension = 0)

seq1Aligned_sudan=toString(subject(output_sudan)) 
seq2Aligned_sudan=toString(pattern(output_sudan))

seq1AlignedV_sudan = strsplit(seq1Aligned_sudan, "")[[1]]
seq2AlignedV_sudan = strsplit(seq2Aligned_sudan, "")[[1]] 

mismatch_sudan = 0 
gapcount_sudan = 0
matchcnt_sudan = 0
for (i in 1:length(seq1AlignedV_sudan)) 
{ 

if (seq1AlignedV_sudan[i]!=seq2AlignedV_sudan[i]) 
  { 
if(seq1AlignedV_sudan[i]=="-" || seq2AlignedV_sudan[i]=="-" )
{ gapcount_sudan=gapcount_sudan+1
}

else
{	mismatch_sudan = mismatch_sudan+1 }

} 
}


##################################################################################

output_Zaire = pairwiseAlignment(pattern=Seq1, subject=Seq3, type = "local", substitutionMatrix = mat,
                      gapOpening = -2, gapExtension = 0)

seq1Aligned_Zaire=toString(subject(output_Zaire)) 
seq2Aligned_Zaire=toString(pattern(output_Zaire))

seq1AlignedV_Zaire = strsplit(seq1Aligned_Zaire, "")[[1]]
seq2AlignedV_Zaire = strsplit(seq2Aligned_Zaire, "")[[1]] 

mismatch_Zaire = 0 
gapcount_Zaire = 0
matchcnt_Zaire = 0
for (i in 1:length(seq1AlignedV_Zaire)) 
{ 

if (seq1AlignedV_Zaire[i]!=seq2AlignedV_Zaire[i]) 
  { 
if(seq1AlignedV_Zaire[i]=="-" || seq2AlignedV_Zaire[i]=="-" )
{ gapcount_Zaire=gapcount_Zaire+1
} 
else
{	mismatch_Zaire = mismatch_Zaire+1
}

} 
}

 print(paste("mismatchs with sudan = ", toString(mismatch_sudan))) 
 print(paste("gapcount with sudan = ", toString(gapcount_sudan))) 

 print(paste("mismatchs with Zaire = ", toString(mismatch_Zaire))) 
 print(paste("gapcount with Zaire = ", toString(gapcount_Zaire))) 

nindel(pattern(output_Zaire))
