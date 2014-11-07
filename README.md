# pattern_finder - basic usage

pattern_finder identifies patterns a,b,c and d of triplet and quadruplet sites (Padhukasahasram 2004) from simulated data in ms format.

pattern_finder was originally written and conceived by Donna Toleno and Peter Morrell, with later modifications by Jeff Ross-Ibarra. See Morrell et al. 2006 and Toleno et al. 2007 for more information and additional software.

## Usage

Make patterns.pl executable: 

`chmod +x patterns.pl`

Pipe ms simulation into patterns.pl:

`ms 10 1 -t 200 -r 200 10000 -c 1 1000 | ./patterns.pl` 

and it prints out:

triplets acount quads bcount ccount dcount
5984 0.0598262032085562 92900 0.071302475780409 0.0302583423035522 0.00933261571582347

Which correspond to the total number of triplets of parsimony informative SNPs, the proportion of those triplets in pattern a, the number of quadruplets of parsimony informative SNPs, and then the proportion of those quadruplets in b, c, and d, followed by the total number of quadruplets.

## References

Morrell PL, Toleno DM, Lundy KE, Clegg MT (2006) Estimating the contribution of mutation, recombination and gene conversion in the generation of haplotypic diversity. Genetics 173: 1705-1723.

Padhukasahasram B, Marjoram P, Nordborg M (2004) Estimating the rate of gene conversion on human chromosome 21. Am J Hum Genet 75: 386-397.

Toleno DM, Morrell PL, Clegg MT (2007) Error detection in SNP data by considering the likelihood of recombinational history implied by three-site combinations. Bioinformatics 23: 1807-1814.
