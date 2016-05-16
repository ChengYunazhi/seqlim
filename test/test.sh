python ../seqlim.py cath ./samples/fasta/
python ../seqlim.py cath ./samples/fasta/ -outfmt phylip
python ../seqlim.py catv ./samples/fasta/
python ../seqlim.py catv ./samples/fasta/ -outfmt nex
python ../seqlim.py cath ./samples/fasta/ -outfmt nex -o ./samples/new_locus1.nex 
python ../seqlim.py cath ./samples/fasta/ -outfmt tsv
python ../seqlim.py cnvt ./samples/fasta/ -outfmt nex -d ./temp
python ../seqlim.py cnvt ./samples/fasta/ -outfmt csv -d ./temp
python ../seqlim.py cnvt ./samples/test.phylip -infmt phylip -outfmt fasta
 
