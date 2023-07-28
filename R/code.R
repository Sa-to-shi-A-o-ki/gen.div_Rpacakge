#DNAbinList can be obtained using read.FASTA function in ape package

#delete columns containing gaps, mixed bases and unknown bases.
deleteGapUnknownColumn<-function(DNAbinList)
{
  deleteColumn<-logical(length(DNAbinList[[1]]))
  for(i in 1:length(DNAbinList))
  {
    #make TRUE for the column containing a gap or a mixed base
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "04")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "02")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "c0")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "a0")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "90")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "60")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "50")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "30")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "e0")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "b0")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "d0")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "70")
    deleteColumn<-deleteColumn | (DNAbinList[[i]] %in% "f0")
  }
  #reverse the bool
  deleteColumn<-!deleteColumn
  #delete the columns
  for(i in 1:length(DNAbinList))
  {
    DNAbinList[[i]]<-DNAbinList[[i]][deleteColumn]
  }
  return(DNAbinList)
}
#count substitution number between two alleles assuming the two have the same length
#assuming allele2 is a single allele, and allele1 is all the alleles
seq.dist<-function(allele1,allele2)
{
  mat<-do.call("rbind",strsplit(allele1,""))
  vec<-strsplit(allele2,"")[[1]]
  return(rowSums(t(t(mat)!=vec)))
}

#calculate nucleotide diversity, expected heterozygosity, Watterson's theta and sigma diversity
calc.diversity<-function(DNAbinList)
{
  #check whether all the sequence lengths are the same
  if(length(unique(summary(DNAbinList)[,1]))!=1)
  {
    return("The sequence length must be the same.")
  }
  alleleSeq<-character(0)
  alleleFreq<-numeric(0)
  charDNA<-as.character(DNAbinList)
  #count and record observed alleles
  for(i in 1:length(DNAbinList))
  {
    currentAlleleSeq<-paste(charDNA[[i]],collapse="")

    matchIndex<-grep(currentAlleleSeq,alleleSeq)
    #when the allele is already recorded
    if(length(matchIndex)>0)
    {
      alleleFreq[matchIndex]<-alleleFreq[matchIndex]+1
    }
    #when the allele is new
    else
    {
      alleleSeq<-c(alleleSeq,currentAlleleSeq)
      alleleFreq<-c(alleleFreq,1)
    }
  }
  seqLength<-nchar(alleleSeq[1])
  alleleFreq<-alleleFreq/sum(alleleFreq)

  #calculate expected heterozygosity
  rawH<-1-sum(alleleFreq^2)
  #calculate nucleotide diversity
  tempPi<-numeric((length(alleleSeq)^2-length(alleleSeq))/2)
  tempSigma<-length(tempPi)
  counter<-1
  isSegregating<-logical(seqLength)
  for(i in 1:length(alleleSeq))
  {
    for(j in 1:length(alleleSeq))
    {
      if(i>=j)
      {
        next
      }
      distance<-seq.dist(alleleSeq[i],alleleSeq[j])
      tempPi[counter]<-alleleFreq[i]*alleleFreq[j]*distance
      tempSigma[counter]<-distance
      counter<-counter+1
      for(k in 1:seqLength)
      {
        if(isSegregating[k])
        {
          next
        }
        if(substring(alleleSeq[i],k,k)!=substring(alleleSeq[j],k,k))
        {
          isSegregating[k]<-TRUE
        }
      }
    }
  }
  rawPi<-sum(tempPi)/seqLength*2
  if(length(alleleSeq)==1)
  {
    rawSigma<-0
  }
  else
  {
    rawSigma<-sum(tempSigma)/seqLength/length(tempSigma)
  }
  H<-rawH*2*length(DNAbinList)/(2*length(DNAbinList)-1)
  if(length(DNAbinList)==1)
  {
    pi<-rawPi
  }
  else
  {
    pi<-rawPi*length(DNAbinList)/(length(DNAbinList)-1)
  }
  theta<-sum(isSegregating)/sum(1/1:(length(DNAbinList)-1))
  psTheta<-theta/seqLength
  result<-matrix(1,nrow=7,ncol=2)
  result[1,1]<-"Expected heterozygosity of population"
  result[1,2]<-rawH
  result[2,1]<-"Expected heterozygosity of sample"
  result[2,2]<-H
  result[3,1]<-"Nucleotide diversity of population"
  result[3,2]<-rawPi
  result[4,1]<-"Nucleotide diversity of sample"
  result[4,2]<-pi
  result[5,1]<-"Watterson's theta"
  result[5,2]<-theta
  result[6,1]<-"per-site Watterson's theta"
  result[6,2]<-psTheta
  result[7,1]<-"sigma (mean substitution rate of all allele types)"
  result[7,2]<-rawSigma
  return(result)
}
