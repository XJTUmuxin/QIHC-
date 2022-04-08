#!/bin/bash
echo -e "assemble L and get contigs by using canu\n"

echo -e "command:\n"

echo -e "canu -s ./../spec.txt -d ./../run -p test genomeSize=100k -pacbio ./../raw_L/raw_L.fasta\n"

canu -s ./../spec.txt -d ./../run -p test genomeSize=100k -pacbio ./../raw_L/raw_L.fasta

echo -e "get the reference-seq by linking the contigs\n"

./get_ref

echo -e "map L to ref by using blasr\n"

echo -e "command:\n"

echo -e "blasr   ./../raw_L/raw_L.fasta\n
        ./../ref/ref.fasta\n
        --header\n
        --bestn 1\n
        -m 5\n
        --out\n
        ./../map_result/L_ref.out\n"

blasr   ./../raw_L/raw_L.fasta\
        ./../ref/ref.fasta\
        --header\
        --bestn 1\
        -m 5\
        --out\
        ./../map_result/L_ref.out

echo -e "get Lm and Lu according to the map result\n"

./get_Lm_Lu

echo -e "map S to Lu by using blasr\n"

echo -e "command\n"

echo -e "blasr   ./../raw_S/raw_S.fasta\n
        ./../raw_L/Lu.fasta\n
        --header\n
        -m 5\n
        --out\n
        ./../map_result/S_Lu.out\n"

blasr   ./../raw_S/raw_S.fasta\
        ./../raw_L/Lu.fasta\
        --header\
        -m 5\
        --out\
        ./../map_result/S_Lu.out

echo -e "judge the heterozygosity of Lu and correct the Lu\n"

./judge_correct_Lu

echo -e "map S to ref and map Lm to ref by using blasr\n"

echo -e "command:\n"

echo -e "blasr   ./../raw_S/raw_S.fasta\n
        ./../ref/ref.fasta\n
        --header\n
        --bestn 1\n
        -m 5\n
        --out\n
        ./../map_result/S_ref.out\n

blasr   ./../raw_L/Lm.fasta\n
        ./../ref/ref.fasta\n
        --header\n
        --bestn 1\n
        -m 5\n
        --out\n
        ./../map_result/Lm_ref.out\n"

blasr   ./../raw_S/raw_S.fasta\
        ./../ref/ref.fasta\
        --header\
        --bestn 1\
        -m 5\
        --out\
        ./../map_result/S_ref.out

blasr   ./../raw_L/Lm.fasta\
        ./../ref/ref.fasta\
        --header\
        --bestn 1\
        -m 5\
        --out\
        ./../map_result/Lm_ref.out

echo -e "judge the heterozygosity of Lm and correct the Lm\n"

./judge_correct_Lm