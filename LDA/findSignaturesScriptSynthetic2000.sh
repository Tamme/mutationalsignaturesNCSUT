#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c

for sigs in {3..7}
	do
	./lda est 0.5 $sigs settings.txt ../data/lda_synthetic_2000_1.txt  random "2000synthetic_rs_"$sigs"sigs_full_moreconv"
	   
	for i in {2..40}
	do
	   echo "Iteration"$i" "$sigs" signatures"
	   ./lda est 0.5 $sigs settings.txt "../data/lda_synthetic_2000_"$i".txt"  "2000synthetic_rs_"$sigs"sigs_full_moreconv/final" "2000synthetic_rs_"$sigs"sigs_full_moreconv"
	done
	python topics.py "2000synthetic_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
