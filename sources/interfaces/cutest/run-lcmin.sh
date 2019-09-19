#!/bin/bash

lcmin=$HOME/Dropbox/USP/mestrado/codigo/lcmin
bin=$lcmin/bin/cutest

outputs_lcmin=$HOME/Dropbox/USP/mestrado/codigo/lcmin/outputs_lcmin

if [ ! -d $outputs_lcmin ]
then
    mkdir $outputs_lcmin
fi

results=lcmin-table.out

lcmin_output=lcmin-tabline.out

# Null table line (update limit time)
nulltabline=" 00000    00000     9.9999999999999999D+99     9.9999999999999999D+99    9.9999D+99    9.9999D+99    9.9999D+99         0         0         0         0         0         0         0     7200.00        0.00   &"

bkp_ifs=$IFS
IFS=$'\n'

printf "\n  %10s   %6s   %6s  |  %4s   %6s   %6s   %24s   %24s   %11s   %11s   %11s   %7s   %7s   %7s   %7s   %7s   %7s   %7s   %9s   %9s\n" \
    "prob" "var" "constr" "flag" "var" "const" "sca obj funct" "unsca obj funct" "scafeas" "feas" "nlp err" "#f" "#g" "#h" "#bf" "#bg" "#iit" "#oit" "pre time" "opt time" >> $results

while read line
do

    problem=$(echo $line | awk -F" " '{print $1}' | sed -e 's/^ *//g' -e 's/ *$//g')
    n=$(echo $line | awk -F" " '{print $2}' | sed -e 's/^ *//g' -e 's/ *$//g')
    m=$(echo $line | awk -F" " '{print $3}' | sed -e 's/^ *//g' -e 's/ *$//g')

    if make -C $lcmin cutest PROBNAME=$problem
    then

	mv $bin/OUTSDIF.d .
	printf "%s\n%s\n" $problem $(date +%c) > runtime.txt
	ulimit -St 7200
	$bin/lcmin
	ulimit -St unlimited
	printf "%s\n" $(date +%c) >> runtime.txt

	mv lcmin.out $outputs_lcmin/$problem-lcmin.out

	if [ -e $lcmin_output ]
	then

	    flag=$(sed -n '1p' $lcmin_output | sed -e 's/^ *//g' -e 's/ *$//g')
	    tabline=$(sed -n '2p' $lcmin_output)

	    if [ $flag -ge -6 ]
	    then
		printf "  %10s   %6d   %6d  |  %4d   %s\n" $problem $n $m $flag $tabline >> $results
	    else
		printf "  %10s   %6d   %6d  |  %4d   %s\n" $problem $n $m $flag $nulltabline >> $results
	    fi

	else
	    flag=-50
	    printf "  %10s   %6d   %6d  |  %4d   %s\n" $problem $n $m $flag $nulltabline >> $results
	fi

	rm -f fort.20 ma41.out ma57.out mc58.out mc59.out AUTOMAT.d OUTSDIF.d $lcmin_output initial_point.txt solution.txt runtime.txt

    else
	flag=-100
	printf "  %10s   %6d   %6d  |  %4d   %s\n" $problem $n $m $flag $nulltabline >> $results
    fi

done < $1

IFS=$bkp_ifs