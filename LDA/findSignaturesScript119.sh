#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c

for sigs in {4..12}
	do
	./lda est 0.5 $sigs settings2.txt ../data/lda_original_119_breast1.txt  random "119breast_rs_"$sigs"sigs_full_moreconv"
	   
	for i in {2..20}
	do
	   echo "Welcome"$i" "$sigs"times"

	   ./lda est 0.5 $sigs settings2.txt "../data/lda_original_119_breast"$i".txt"  "119breast_rs_"$sigs"sigs_full_moreconv/final" "119breast_rs_"$sigs"sigs_full_moreconv"
    
	done
	python topics.py "119breast_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
