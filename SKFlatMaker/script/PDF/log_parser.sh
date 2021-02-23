#!/bin/bash

if [ -z "$1" ]
then
    echo "Usage: $0 logfile"
    exit 1
fi

filename=$1

readarray -t lines< <(cat $filename|grep "\[SKFlatMaker::endRun,lines\]"|sed 's/^\[SKFlatMaker::endRun,lines\] //'|sed 's/&lt;/</g'|sed 's/&gt;/>/g'|grep "<weight.*</weight>" -A1 -B1)

scaleA=0
scaleB=0
combine=""
pdfA=0
pdfB=0
asA=0
asB=0
asScaleA=0
asScaleB=0
depth=0
candidate=()
for i in $(seq 0 $((${#lines[@]}-1)));
do
    line=${lines[$i]}
    if [[ "$line" == *"<weightgroup"* ]]
    then
	printf "%"$depth"s";echo "$line"
	depth=$(($depth+1))
	printf "%"$depth"s"
    elif echo $line|grep -q "<weight.*</weight>"
    then
	id=$(echo ${lines[$i]}|grep -o " id=['\"][0-9]*['\"]"|grep -o '[0-9]*')
	##echo -n $id" "
	candidate+=("$id")
    elif [[ "$line" == *"</weightgroup>"* ]]
    then
	depth=$((depth-1))
	printf "%"$depth"s";echo "$line"
	if [ "$depth" -eq 0 ]
	then
	    candidate=( $(printf "%s\n" "${candidate[@]}"|sort -n) )
	    echo "${candidate[@]}" "total=${#candidate[@]}" 
	    if [ "${#candidate[@]}" -eq 9 ] && [ "$scaleA" -eq 0 ] && [ "$scaleB" -eq 0 ]
	    then
		scaleA=${candidate[0]}
		scaleB=${candidate[${#candidate[@]}-1]}
	    elif [ "${#candidate[@]}" -eq 100 ] && [ "$pdfA" -eq 0 ] && [ "$pdfB" -eq 0 ]
	    then
		pdfA=${candidate[0]}
		pdfB=${candidate[${#candidate[@]}-1]}
	    elif [ "${#candidate[@]}" -eq 102 ] && [ "$pdfA" -eq 0 ] && [ "$pdfB" -eq 0 ]
	    then
		pdfA=${candidate[0]}
		pdfB=${candidate[${#candidate[@]}-3]}
		asA=${candidate[${#candidate[@]}-2]}
		asB=${candidate[${#candidate[@]}-1]}
	    elif [ "${#candidate[@]}" -eq 103 ] && [ "$pdfA" -eq 0 ] && [ "$pdfB" -eq 0 ]
	    then
		pdfA=${candidate[1]}
		pdfB=${candidate[${#candidate[@]}-3]}
		asA=${candidate[${#candidate[@]}-2]}
		asB=${candidate[${#candidate[@]}-1]}
	    fi
	    candidate=()
	fi
    fi
done
lhaid=$(cat $filename|grep "=[ ]*lhaid"|awk '{print $2}'|sort|uniq)
[ -z "$lhaid" ] && { echo "WARNING cannot find default lhapdf id. use lhaid=306000";lhaid=306000; }
echo lhaid=$lhaid
lhapdf=$(grep $lhaid $LHAPDF_DATA_PATH/pdfsets.index|awk '{print $2}')
echo lhapdf=$lhapdf
pdfinfo=$(head -n1 $LHAPDF_DATA_PATH/$lhapdf/${lhapdf}.info)
echo $pdfinfo|grep --color replica && combine=gaussian
echo $pdfinfo|grep -i --color hessian && combine=hessian
test $((($pdfB-$pdfA)%2)) -eq 0 && pdfA=$(($pdfA+1))
echo $pdfinfo|grep --color=always 0\.117|grep --color 0\.119 && asScaleA=1.5 && asScaleB=1.5
echo $pdfinfo|grep --color=always 0\.116|grep --color 0\.120 && asScaleA=0.75 && asScaleB=0.75

echo "#####DEBUG#####"
echo $pdfinfo
echo $scaleA $scaleB
echo $combine
echo $pdfA $pdfB
echo $asA $asB
echo $asScaleA $asScaleB
echo "####AUTO####"
printf "echo \"%d,%d\t%s\t%d,%d\t%d,%d\t%s,%s\" > ../MCPDFInfo/%s.txt\n" "$scaleA" "$scaleB" "$combine" "$pdfA" "$pdfB" "$asA" "$asB" "$asScaleA" "$asScaleB" $(echo $filename|sed 's/logs_//'|sed 's/.log$//')
read -p "exec?(y/n): " YES
if [ "$YES" = "y" ] || [ "$YES" = "Y" ];then
    mkdir -p ../MCPDFInfo/$(echo $filename|sed 's/logs_//'|sed 's/.log$//'|xargs -i dirname {})
    printf "echo \"%d,%d\t%s\t%d,%d\t%d,%d\t%s,%s\" > ../MCPDFInfo/%s.txt\n" "$scaleA" "$scaleB" "$combine" "$pdfA" "$pdfB" "$asA" "$asB" "$asScaleA" "$asScaleB" $(echo $filename|sed 's/logs_//'|sed 's/.log$//')|xargs -i sh -c '{}'
fi
