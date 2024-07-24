#!/bin/bash
##
mkdir submit_split
COUNT=0
echo "#!/bin/bash" > clean_splitted.sh
while read p; do
        COUNT=$(( $COUNT + 1 ))
        echo "#!/bin/bash" > submit_split/submit_clean_$COUNT
        echo "#SBATCH --account=rrg-shaferab" >> submit_split/submit_clean_$COUNT
        echo "#SBATCH --mem-per-cpu 20G" >> submit_split/submit_clean_$COUNT
        echo "#SBATCH --cpus-per-task=10" >> submit_split/submit_clean_$COUNT
        echo "#SBATCH --time=12:00:00" >> submit_split/submit_clean_$COUNT
        echo $p >> submit_split/submit_clean_$COUNT
        ##
        echo "sbatch submit_split/submit_clean_$COUNT" >> clean_splitted.sh;
done < clean.sh
echo "Split is done on $COUNT sequences"
##
