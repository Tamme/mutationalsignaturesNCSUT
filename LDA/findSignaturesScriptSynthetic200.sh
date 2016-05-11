#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c

for sigs in {3..7}
	do
	./lda est 0.5 $sigs settings.txt ../data/lda_synthetic_200_1.txt  random "200synthetic_rs_"$sigs"sigs_full_moreconv"
	   
	for i in {2..40}
	do
	   echo "Iteration"$i" "$sigs" signatures"
	   ./lda est 0.5 $sigs settings.txt "../data/lda_synthetic_200_"$i".txt"  "200synthetic_rs_"$sigs"sigs_full_moreconv/final" "200synthetic_rs_"$sigs"sigs_full_moreconv"
	done
	python topics.py "200synthetic_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
