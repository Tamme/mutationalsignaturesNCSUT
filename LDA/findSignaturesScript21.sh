#!/bin/bash
cd /home/tamme/Desktop/Masters/lda-c

for sigs in {8..12}
	do
	./lda est 0.5 $sigs settings.txt ../data/lda_original_21_breast1.txt  random "21breast_rs_"$sigs"sigs_full_moreconv"
	   
	for i in {2..40}
	do
	   echo "Welcome"$i$sigs"times"
	#   var="21breast_yo_"$i"sigs"
	#   echo $var
	   ./lda est 0.5 $sigs settings.txt "../data/lda_original_21_breast"$i".txt"  "21breast_rs_"$sigs"sigs_full_moreconv/final" "21breast_rs_"$sigs"sigs_full_moreconv"
    
	done
	python topics.py "21breast_rs_"$sigs"sigs_full_moreconv/final.beta" ../data/gen_vocab.txt 96
done
