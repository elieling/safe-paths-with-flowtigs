## Additional software to verify that the "# contig" given by the quast report matches the the number of contigs in the fasta files and to discover the contig length of "# unaligned contig".

Doing these using E.coli bacteria and k = 31.


### Calculating the number of contigs with length greater than x and trying to match that number to the number of # contigs.

Values of #contigs:

    - unitigs: 2,168 
    - trivial omitigs: 1,235 
    - multi safe: 1,240 
    - flowtigs: 1,461 
    - omnitigs: 1,219


For each algorithm, the program gets the value 30 as well as all the lesser values. So, quast reports the number of contigs with lenght greater than 30, which matches with k = 31.


### Calculating the number of contigs with length greater than x and trying to match that number to the number of # unaligned contigs.

Values of #unaligned contigs:

    - unitigs: 1,314 
    - trivial omitigs: 135
    - multi safe: 137
    - flowtigs: 241 
    - omnitigs: 133


#### Results:
    - Unitigs, unaligned contigs: 1314
        > 44 : 1324
        > 45 : 1311
        NO EXACT MATCH!

    - Trivial omnitigs, unaligned contigs: 135
        > 11079 : 135
           ...
        > 11088 : 135

    - multi-safe, unaligned contigs: 137
        > 11000 : 137
           ...
        > 11049 : 137

    - flowtigs, unaligned contigs: 241
        > 4838 : 241
           ...
        > 4859 : 241

    - omnitigs, unaligned contigs: 133
        > 11240 : 133
           ...
        > 11268 : 133








unitigs 1314
66 1310
67 1327
######################################################################
trivial omnitigs 135
65 130
66 136
######################################################################
muti-safe 137
65 132
66 138
######################################################################
flowtigs 241
65 235
66 243
######################################################################
exact 66 133
omnitigs 133
66 133
67 136
