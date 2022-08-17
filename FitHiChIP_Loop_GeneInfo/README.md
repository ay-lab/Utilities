
The code is for identifying the genes overlapping FitHiChIP loop endpoints.

In the parameter section, specify the FitHiChIP loops, reference GTF files containing gene annotations.

Specifically, set the columns containing chromosome, TSS and gene names in the GTF file.

Also specify the offset (slack) allowed to compute the overlap between gene TSS and the interacting bins.

For example, offset = 5000 means 5 Kb slack is allowed.
