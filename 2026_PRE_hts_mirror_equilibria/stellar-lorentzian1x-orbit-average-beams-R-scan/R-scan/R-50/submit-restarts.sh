#!/bin/bash

# Array to store job IDs
job_ids=()

# Ask the user how many jobs to submit
# read -p "Enter the number of jobs to submit: " num_jobs
num_jobs=5

# Submit the first job and capture its job ID
first_job_id=$(sbatch jobscript-gkyl-stellar-amd | awk '{print $4}')
job_ids+=($first_job_id)

# Submit subsequent jobs with dependency on the previous job
for i in $(seq 2 $num_jobs)
do
  previous_job_id=${job_ids[-1]}
  next_job_id=$(sbatch --dependency=afterany:$previous_job_id jobscript-gkyl-stellar-amd | awk '{print $4}')
  job_ids+=($next_job_id)
done

# Print all job IDs
echo "Submitted jobs with IDs: ${job_ids[@]}"