# getSeq

Extract DNA sequences from larger contexts like genome sequences in the form of a FASTA file. 
The user can supply a GFF3 file containing the FASTA file's annotation in order to then 
extract DNA sequences by their gene symbols (as found in the GFF file). It is possible to 
specify which annotated element of a given gene shall be extracted (e.g. gene, mRNA, CDS).
Alternatively, the user can manually set start and end positions of a sequence and its 
strandedness in order to retrieve that sequence from the FASTA file.

In both usage scenarios it is possible to define a number of nucleotides to be extracted 
up- and/or downstream of the specified annotated element. Additionally, it is also 
possible to just extract up- or downstream sequences without the related annotated element's 
sequence itself (e.g. in order to create promoteromes).

