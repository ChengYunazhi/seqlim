seqlim cath ./fasta/
seqlim cath ./fasta/ -outfmt fasta
seqlim cath ./fasta/ -outfmt nex
seqlim cath ./fasta/ -outfmt phylip
seqlim cath ./fasta/ -outfmt tsv
seqlim cath ./fasta/ -outfmt nex -o ./temp/new_locus1.nex 
seqlim cath ./fasta/ -outfmt msf -o ./temp/new_locus1.msf
seqlim cath ./fasta/ -outfmt phylip -o ./temp/new_locus1.phy
seqlim cath ./fasta/ -outfmt phylip -line_length 50 -block_length 10
seqlim catv ./fasta/
seqlim catv ./fasta/ -outfmt fasta
seqlim catv ./fasta/ -outfmt nex
seqlim catv ./fasta/ -outfmt phylip
seqlim catv ./fasta/ -outfmt nex -o ./temp/new_locus1-v.nex 
seqlim catv ./fasta/ -outfmt phylip -o ./temp/new_locus1-v.phy
seqlim cnvt ./fasta/ -outfmt nex -o ./temp
seqlim cnvt ./fasta/ -outfmt phylip -o ./temp
seqlim cnvt ./fasta/ -outfmt csv -o ./temp
seqlim cath ./fasta_gz/ --read_gzipped
seqlim catv ./fasta_gz/ --read_gzipped
seqlim cnvt ./samples/test.phylip -infmt phylip -outfmt fasta
seqlim cnvt ./samples/test.phylip -infmt phylip -outfmt fasta -remove -
