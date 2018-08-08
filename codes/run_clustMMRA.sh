#!/bin/bash
echo -e "You are running clustMMRA pipeline. The pipeline is structured in 3 steps: \n "
echo -e "The sthird step due to the network analysis is computationally expensive, we suggest for this reason to run it on a cluster for more details look at http://wiki.c2b2.columbia.edu/califanolab/index.php/Software/ARACNE. \n" 
echo -e "The parameters needed for the analysis will be asked while the program is running. For more informations look at ReadMe.pdf \n"


#######  first step
echo -e "STEP1: microRNA clusters differential expression analysis \n "
if [ ! -f ../results/step1_results.txt ]; then
	mkdir ../results
	echo -n "Enter the number of subtypes > "
	read text2
	echo "You entered: $text2"
	echo -n "Enter the path to the microRNA dataset file > "
	read text3
	echo "You entered: $text3"
	declare -a my_array
	for i in $(seq 1 $text2)
	do
		echo -n "Enter the starting column of class $i > "
		read text
		echo "You entered: $text"
		my_array=("${my_array[@]}" $text) 
		echo  ${my_array[@]}
		echo -n "Enter the ending column of class $i > "
		read text
		echo "You entered: $text" 
		my_array=("${my_array[@]}" $text) 
		echo ${my_array[@]}

	done
	
	echo -n "Enter the Kolmogorov-Smirnov pvalue threshold that you want to use/test (suggested 0.05)> "
	read text4
	echo "You entered: $text4"
	echo -n "Enter the log fold change threshold that you want to use/test ( if you don't want to use a log fold change threshold just set it to 0) > "
	read text5
	echo "You entered: $text5"
	Rscript step1.R $text2 $text3 $text4 $text5 ${my_array[@]}

	if [ ! -f ../results/step1_results.txt ]; then
		exit
	fi
else
	echo -e "step 1 previously executed go directly to step 2. If you haven't executed step 1 delete file  ../results/step1_results.txt \n"
	echo -n "Enter the number of subtypes > "
	read text2
	echo "You entered: $text2"
	echo -n "Enter the path to the microRNA dataset file > "
	read text3
	echo "You entered: $text3"
fi	
################ second step
if [ ! -f ../results/step2_results.txt ]; then
	echo -n "Step 2: Target transcripts enrichment analysis"
		Rscript step2.R 
	if [ ! -f ../results/step2_results.txt ]; then
		exit
	fi
else
	echo -e "Step 2 previously executed go directly to step 3. If you haven't executed step 2 delete file  ../results/step2_results.txt \n"
fi
############### third step
echo -e "Step2: consolidation of candidate microRNAs by network analysis \n"
echo -e "Network reconstruction using ARACNE \n"
mkdir ../results/boot
echo -n "Enter the path to the mRNA dataset file > "
read text9
echo "You entered: $text9"
Rscript preproc_phase2.R  $text3 $text9
FILES=../results/boot/
for NAME in ../results/boot/*.txt
do
	sed -i 's/"//g' $NAME
	perl ./aracne/Scripts/qsubbootstrap.pl BootstrapNetworks 1 100 -i $NAME -p 1e-7 –e 0 -h 
	perl Scripts/getconsensusnet.pl BoostrapNetworks 1e-12
	b=_k*.adj
	d=${NAME%.txt} 
	f=${NAME#../results/boot/} 
	a="$NAME$b"
	c=net_

	perl conversion_network1.pl -g1 $a -g2 "../results/$c$f"
done
echo -e "FET test of the obtained networks \n"
echo -n "Enter the Fisher pvalue threshold that you want to use (we suggest to estimate the threshold through a null model for more details look at ReadMe.pdf)> "
read text10
echo "You entered: $text10"
Rscript FET.R $text2 $text10