while read INTIME; do

        echo "submit job for time segment $INTIME"

        sbatch --job-name=dissdiag_$INTIME.run --output=io/dissdiag_$INTIME.out --error=io/dissdiag_$INTIME.err --export=code=$INTIME dissipation_diagnostics.sh

done < ../download_time_list.txt