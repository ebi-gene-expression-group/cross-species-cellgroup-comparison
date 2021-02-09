all_top=''

while read -r d; do 
    best_result=$(grep "known intersecting cell types were predicted as top match by marker gene composition" $d/*.report.md | sed 's/:/ /' | sort -nrk2 | head -n 1 | grep -oP '\d+ of \d+')    
    this_top=$(grep "$best_result known intersecting cell types were predicted as top match by marker gene composition" $d/*.report.md| awk -F'[/:]' '{print $4}' | sed 's/.report.md//')

    echo -e "Top for $d is $best_result, found in: \n$this_top"
    all_top="$all_top\n$this_top"
done <<< "$(ls -d output/*/ | grep E-)"


echo "Overall top parameter sets:"

echo -e "$all_top" | sort |  uniq -c | sort 
