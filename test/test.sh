python ../seqlim.py cath ./fasta/
python ../seqlim.py cath ./fasta/ -outfmt fasta
python ../seqlim.py cath ./fasta/ -outfmt nex
python ../seqlim.py cath ./fasta/ -outfmt phylip
python ../seqlim.py cath ./fasta/ -outfmt tsv
python ../seqlim.py cath ./fasta/ -outfmt nex -o ./temp/new_locus1.nex 
python ../seqlim.py cath ./fasta/ -outfmt phylip -o ./temp/new_locus1.phy
python ../seqlim.py cath ./fasta/ -outfmt phylip -line_length 50 -block_length 10
python ../seqlim.py catv ./fasta/
python ../seqlim.py catv ./fasta/ -outfmt fasta
python ../seqlim.py catv ./fasta/ -outfmt nex
python ../seqlim.py catv ./fasta/ -outfmt phylip
python ../seqlim.py catv ./fasta/ -outfmt nex -o ./temp/new_locus1-v.nex 
python ../seqlim.py catv ./fasta/ -outfmt phylip -o ./temp/new_locus1-v.phy
python ../seqlim.py cnvt ./fasta/ -outfmt nex -o ./temp
python ../seqlim.py cnvt ./fasta/ -outfmt phylip -o ./temp
python ../seqlim.py cnvt ./fasta/ -outfmt csv -o ./temp
python ../seqlim.py cnvt ./samples/test.phylip -infmt phylip -outfmt fasta
python ../seqlim.py cnvt ./samples/test.phylip -infmt phylip -outfmt fasta -remove -
