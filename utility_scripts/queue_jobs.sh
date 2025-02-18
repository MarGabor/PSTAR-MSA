#!/bin/bash

dir="../jobs_to_run_no_ref"

counter=0
for job_name in "$dir"/*; do
	counter=$((counter+1))
	job_name=$(basename $job_name)
	echo $job_name
	if [ $counter -lt 6 ]; then
		continue
	fi
	if [ $counter -gt 10 ]; then
		exit 0
	fi
	short_job_name=$job_name
	if [ ${#job_name} -gt 11 ]; then
		short_job_name=${job_name:0:11}
	fi

	input=$(printf "#!/bin/bash\npython3 -m PSTAR-MSA --wrapjob --msadir "$dir/$job_name" --outputdir ../results/ --jobname $short_job_name --diamondfile ../Tools/DIAMOND/diamond --diamonddbfile ../atom_diamond_db/diamond_db.dmnd --verbose --verbose --threads 20 --baseline --dalibindir ../Tools/DaliLite.v5/bin/ --locpdbdb ../pdb_database/ --samplesize 1000")
	echo "$input"
	echo "$input" | sbatch --nodes=1 --tasks-per-node=1 --cpus-per-task=10 --mem=20000
done
