# motif-finder
A MarkovCluster algorithm seeking for the motifs in DNA-containing, Roche454-derived FASTA files.

**Keep in mind:**
* You must provide both files: a FASTA file and a corresponding quality scores file.
* Having provided those files, you can choose whether to validate sequences by selecting the confidence level.
* Window size **is not** the length of the motif (sometimes can be, tough). It indirectly indicates a number of substrings obtained from sequence decomposition,
* Altough MarkovCluster is a polynomial-time algorithm, a size of the matrices, many matrix multiplication and squaring operations are factors that can consume time. Be patient!
