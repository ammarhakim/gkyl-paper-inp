#!/bin/bash

#==============================================================================
# submit-restarts-perlmutter.sh
# 
# Submit a chain of dependent jobs to the SLURM scheduler. Each job depends
# on the previous job completing, allowing for automatic restart chains.
#==============================================================================

set -o pipefail

# Default values
num_jobs=1
restart_from=""
jobscript="jobscript-gkyl-perlmutter"

#==============================================================================
# Help function
#==============================================================================
show_help() {
  cat <<EOF
Usage: submit-restarts-perlmutter.sh [OPTIONS]

Options:
  --jobs NUM, -j NUM              Number of jobs to submit in the chain
                                  (default: 1)
  
  --restart-from JOB_ID, -r JOB_ID
                                  Start the chain from an existing job ID.
                                  The first new job will depend on JOB_ID.
  
  --jobscript FILE                Path to jobscript to submit
                                  (default: jobscript-gkyl-perlmutter)
  
  --help, -h                      Show this help message and exit
EOF
}

#==============================================================================
# Argument parsing
#==============================================================================
while [[ $# -gt 0 ]]; do
  case $1 in
    --jobs|-j)
      num_jobs="$2"
      shift 2
      ;;
    --restart-from|-r)
      restart_from="$2"
      shift 2
      ;;
    --jobscript)
      jobscript="$2"
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1" >&2
      echo "Use --help for usage information." >&2
      exit 1
      ;;
  esac
done

#==============================================================================
# Validation
#==============================================================================
if ! [[ "$num_jobs" =~ ^[0-9]+$ ]] || [ "$num_jobs" -lt 1 ]; then
  echo "Error: Number of jobs must be a positive integer." >&2
  exit 1
fi

if [ -n "$restart_from" ] && ! [[ "$restart_from" =~ ^[0-9]+$ ]]; then
  echo "Error: Job ID must be a positive integer." >&2
  exit 1
fi

if [ ! -f "$jobscript" ]; then
  echo "Error: Jobscript not found: $jobscript" >&2
  exit 1
fi

#==============================================================================
# Submit jobs
#==============================================================================
job_ids=()

if [ -z "$restart_from" ]; then
  # Submit the first job with no dependency
  first_job_id=$(sbatch "$jobscript" | awk '{print $4}')
  if [ -z "$first_job_id" ]; then
    echo "Error: Failed to submit first job." >&2
    exit 1
  fi
  job_ids+=($first_job_id)
  echo "Submitted initial job: $first_job_id"
else
  # Start from existing job
  job_ids+=($restart_from)
  echo "Starting chain from job: $restart_from"
fi

# Submit subsequent jobs with dependency on the previous job
for i in $(seq 1 $((num_jobs - 1)))
do
  previous_job_id=${job_ids[-1]}
  next_job_id=$(sbatch --dependency=afterany:$previous_job_id "$jobscript" | awk '{print $4}')
  if [ -z "$next_job_id" ]; then
    echo "Error: Failed to submit job $((i + 1))." >&2
    exit 1
  fi
  job_ids+=($next_job_id)
  echo "Submitted job $((i + 1)): $next_job_id (depends on $previous_job_id)"
done

#==============================================================================
# Summary
#==============================================================================
echo ""
echo "Successfully submitted job chain:"
echo "${job_ids[*]}"