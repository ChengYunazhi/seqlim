## SEQLIM
Concatenate and Convert Multiple Sequence Alignments

# Description
SEQLIM is a python script for manipulating biological sequences. It concatenates multiple sequence alignments (MSAs) horizontally or vertically, and converts MSAs into various formats (fasta, phylip, nexus, msf, tsv, and csv). The horizontal concatenation of MSAs is often used for multi-loci/multi-gene phylogenetic analysis and phylogenomics.
 
# Installation

* Install Python 2.7 or higher, Python installers are available at https://www.python.org/.
* Clone or download this repo, make `seqlim.py` executable, add its path to your `PATH`, and test it.
```
$ seqlim.py -h
```
* Or go to a directory with seqlim.py, and use with `python` command.
```
$ python seqlim.py -h
```

# Examples

* Suppose two sequence files in FASTA format in `./test/fasta`.
 
 
`Locus1.fasta`

    >Escheri1
    CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAAC
    >Enteroc1
    UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGC
 
 
`Locus2.fasta`

    >Escheri2
    UAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
    >Enteroc2
    UAGCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
 

* Concatenate these files horizontally.

```
$ seqlim.py cath ./test/fasta
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACUA
GCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGCUA
GCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 

* Concatenate the files vertically.

```
$ seqlim.py catv ./test/fasta
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAAC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGC
>Escheri2
UAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc2
UAGCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 

* Set alternative input sequence format after `-infmt`. SEQLIM accepts `phylip`, `phy`, or `ph` for PHYLIP format; `nexus`, `nex`, or `nxs` for NEXUS format.

```
$ seqlim.py -infmt phylip cath ./test/phylip
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAAC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGC
>Escheri2
UAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc2
UAGCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 
 
* Set output format after `-outfmt`.
``` 
$ seqlim.py -outfmt phylip cath ./test/fasta
 2 114
Escheri1     CCUGGCGGCC GUAGCGCGGU GGUCCCACCU GACCCCAUGC CGAACUCAGA AGUGAAACUA
Enteroc1     UGUGGUGGCG AUAGCGAGAA GGAUACACCU GUUCCCAUGC CGAACACAGA AGUUAAGCUA

             GCGCCGAUGG UAGUGUGGGG UCUCCCCAUG CGAGAGUAGG GAACU--GCC AGGC
             GCGCCGAUUG UAGUGAAGGG UUUCCCUUUG UGAGAGUAGG --ACGUCGCC ACGC
```
 
* The line and block lengths of sequences can be adjusted using `-line\_length` and `-block\_length`, respectively.
```
$ python seqlim.py -outfmt phylip -line_length 50 -block_length 5 cath ./test/fasta
 2 114
Escheri1     CCUGG CGGCC GUAGC GCGGU GGUCC CACCU GACCC CAUGC CGAAC UCAGA
Enteroc1     UGUGG UGGCG AUAGC GAGAA GGAUA CACCU GUUCC CAUGC CGAAC ACAGA

             AGUGA AACUA GCGCC GAUGG UAGUG UGGGG UCUCC CCAUG CGAGA GUAGG
             AGUUA AGCUA GCGCC GAUUG UAGUG AAGGG UUUCC CUUUG UGAGA GUAGG

             GAACU --GCC AGGC
             --ACG UCGCC ACGC
```

* Save an output.
``` 
$ seqlim.py -o ./test/temp/concatenated.fasta cath ./test/fasta
```
 
* Just format conversion.
```
$ seqlim.py -outfmt phylip -o ./test/temp/converted.phylip cnvt ./test/fasta/locus1.fasta
``` 
 
* Convert all sequence files in `./test/fasta` to another format (phylip) and save them in `./test/phylip`.
```
$ seqlim.py -o ./test/phylip -outfmt phylip cnvt ./test/fasta
``` 
 
