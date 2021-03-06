#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c-506

for sigs in {4..12}
	do
	./lda2 est 0.5 $sigs settings3.txt ../data/lda_original_506_breast1.txt  random "506breast_rs_"$sigs"sigs_full_moreconv"
	   
	for i in {2..20}
	do
	   echo "Welcome"$i" "$sigs"times"

	   ./lda2 est 0.5 $sigs settings3.txt "../data/lda_original_506_breast"$i".txt"  "506breast_rs_"$sigs"sigs_full_moreconv/final" "506breast_rs_"$sigs"sigs_full_moreconv"
    
	done
	python topics.py "506breast_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
