#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c

for sigs in {3..7}
	do
	./lda est 0.5 $sigs settings.txt ../data/lda_synthetic_20_1.txt  random "20synthetic_rs_"$sigs"sigs_full_moreconv"
	
	for i in {1..99}
	do
		echo "Iteration"$i" "$sigs" signatures"
	   ./lda est 0.5 $sigs settings.txt "../data/lda_synthetic_20_"$i".txt"  "20synthetic_rs_"$sigs"sigs_full_moreconv/final" "20synthetic_rs_"$sigs"sigs_full_moreconv"
	done
	python topics.py "20synthetic_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
