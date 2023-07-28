\name{calc.diversity}
\alias{calc.diversity}
\title{calc.diversity}
\usage{
calc.diversity(DNAbinList)
}
\description{
Calculate genetic diversities from binary DNA data in list format.
}
\arguments{
\item{DNAbinList}{Binary DNA data in List format. This is obtained by read.FASTA function in ape package. All sequences must have the same length.}
}
\value{This function returns seven types of genetic diversities as a matrix.
    1: Expected heterozygosity calculated for a population.
    2: Expected heterozygosity calculated for a sample. (Corrected by 2n/(2n-1).)
    3: Nucleotide diversity calculated for a population.
    4: Nucleotide diversity calculated for a sample. (Corrected by n/(n-1).)
    5: Watterson's theta.
    6: Per-site Watterson's theta. (Divided by the alignment length.)
    7: Final sigma or mean substitution rate between all allele types.
    Wcalhere n is the number of sequences.
}
\examples{
#This example requires ape package.
#library(ape)
#seq<-read.FASTA(system.file("alignment.fasta",package="gen.div"))
#seq<-deleteGapUnknownColumn(seq)
#calc.diversity(seq)
}
\references{
Nei M, Tajima F. 1981. DNA polymorphism detectable by restriction endonucleases. Genetics 97(1):145–163.

Nei M, Li WH. 1979. Mathematical model for studying genetic variation in terms of restriction endonucleases. Proc Natl Acad Sci U.S.A. 76(10):5269–5273.doi.org/10.1073/pnas.76.10.5269

Nei M. 1977. Estimation of average heterozygosity and genetic distance from a small number of individuals. Genetics 89: 583–590. doi.org/10.1093/genetics/89.3.583

Watterson GA. 1975. On the number of segregating sites in genetical models without recombination. Theor Popul Biol. 7(2): 256–276. doi.org/10.1016/0040-5809(75)90020-9
}
