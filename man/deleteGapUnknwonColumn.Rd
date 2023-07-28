\name{deleteGapUnknownColumn}
\alias{deleteGapUnknownColumn}
\title{deleteGapUnknownColumn}
\usage{
deleteGapUnknownColumn(DNAbinList)
}
\description{
Delete columns including a gap, mixed base or unknown base from binary DNA data in list format. The returned data is in the same format as the argument.
}
\arguments{
\item{DNAbinList}{Binary DNA data in List format. This is obtained by read.FASTA function in ape package. All sequences must have the same length.}
}
