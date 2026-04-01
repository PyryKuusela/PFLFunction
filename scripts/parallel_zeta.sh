#!/bin/bash

# Get number of CPU cores
NUM_CORES=$(nproc)
#NUM_CORES=1
echo $NUM_CORES

# Create file:
operator=$1
primes=$2
scaling=$3 # scaling C
label=$4
acc=${5:-0} # p-adic accuracy A; if not set, will use default value
nadd=${6:-0} # number of additional terms; if set, should include acc

python3 create_inputs.py "$primes" "$operator" "$label" "$scaling" "$acc"

# File containing input arguments (space-separated)
INPUT_FILE="inputs/inputs_${4}.txt"

# Function to run each job
run_job() {
    local arg1=$1
    local arg2=$2
    local arg3=$3
    local arg4=$4
    local arg5=$5
    python3 generate_zeta.py "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" "$nadd" "$label"
}

# Job counter
job_count=0

echo "Starting computations."
echo $(date)
#start=`date +%s`

# Read each line from the input file
while IFS= read -r line || [ -n "$line" ]; do
    read -r arg1 arg2 arg3 arg4 arg5 <<< "$line"
    # Start the job in the background
    run_job "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" &
    ((job_count++))

    # If we've reached the max number of parallel jobs, wait for one to finish
    if [[ "$job_count" -ge "$NUM_CORES" ]]; then
        wait -n
        ((job_count--))
    fi
done < "$INPUT_FILE"

# Wait for all remaining background jobs
wait

echo "All jobs finished."
echo $(date)

python3 delete_inputs.py "$label"
#end=`date +%s`

#runtime=$((end-start))
#echo runtime
