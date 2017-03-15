# seqlim
Concatenate and Convert Multiple Sequence Alignments

# Description
`seqlim` includes a python library and an executable for manipulating biological sequences. It concatenates multiple sequence alignments (MSAs) horizontally or vertically, and converts MSAs into various formats (fasta, phylip, nexus, msf, tsv, and csv). The horizontal concatenation of MSAs is often used for multi-loci/multi-gene phylogenetic analysis and phylogenomics.
 
# Installation

* Install Python 2.7 or higher, Python installers are available at https://www.python.org/.
* Clone or download this repo and install using setup.py.

```
$ python setup.py install
```

* Confirm the installation of `seqlim` executable.
```
$ seqlim -h
```

* Confirm the installation of `seqlim` library.

```
$ python
>>> from seqlim import Seq
```


# Executable examples

* Suppose two sequence files in FASTA format in `./test/fasta`.
 
 ```
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
 ```

* Concatenate these files horizontally.

```
$ seqlim cath ./test/fasta
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACUA
GCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGCUA
GCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 

* Concatenate the files vertically.

```
$ seqlim catv ./test/fasta
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAAC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGC
>Escheri2
UAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc2
UAGCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 

* Set an input sequence format after `-infmt`. `seqlim` accepts 'fasta', 'fas', 'mfa', 'fna', 'fsa' or 'fa' for FASTA format, 'phylip' or 'phy' for PHYLIP format and 'msf' for MSF format.

```
$ seqlim -infmt phylip cath ./test/phylip
>Escheri1
CCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAAC
>Enteroc1
UGUGGUGGCGAUAGCGAGAAGGAUACACCUGUUCCCAUGCCGAACACAGAAGUUAAGC
>Escheri2
UAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACU--GCCAGGC
>Enteroc2
UAGCGCCGAUUGUAGUGAAGGGUUUCCCUUUGUGAGAGUAGG--ACGUCGCCACGC
``` 
 
* Set an output sequence format after `-outfmt`. `seqlim` accepts 'fasta', 'fas', 'mfa', 'fna', 'fsa' or 'fa' for FASTA format, 'phylip' or 'phy' for PHYLIP format, 'nexus', 'nex' or 'nxs' for NEXUS format, 'msf' for MSF format, 'csv' for CSV format and 'tsv' for TSV format.
``` 
$ seqlim -outfmt phylip cath ./test/fasta
 2 114
Escheri1     CCUGGCGGCC GUAGCGCGGU GGUCCCACCU GACCCCAUGC CGAACUCAGA AGUGAAACUA
Enteroc1     UGUGGUGGCG AUAGCGAGAA GGAUACACCU GUUCCCAUGC CGAACACAGA AGUUAAGCUA

             GCGCCGAUGG UAGUGUGGGG UCUCCCCAUG CGAGAGUAGG GAACU--GCC AGGC
             GCGCCGAUUG UAGUGAAGGG UUUCCCUUUG UGAGAGUAGG --ACGUCGCC ACGC
```
 
* The line and block lengths of sequences can be adjusted using `-line_length` and `-block_length`, respectively.
```
$ seqlim -outfmt phylip -line_length 50 -block_length 5 cath ./test/fasta
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
$ seqlim -o ./test/temp/concatenated.fasta cath ./test/fasta
```
 
* Just format conversion.
```
$ seqlim -outfmt phylip -o ./test/temp/converted.phylip cnvt ./test/fasta/locus1.fasta
``` 
 
* Convert all sequence files in `./test/fasta` to another format (phylip) and save them in `./test/phylip`.
```
$ seqlim -o ./test/phylip -outfmt phylip cnvt ./test/fasta
``` 
 
