#!/bin/bash

submit_jobs() {
    
    GS_DATA=$1
    NARROW_DATA=$2
    WIDE_DATA=$3
    TYPE=$4

    ids=()
    
    # Define your species list here
    SPECIES_LIST=("Abies" "Acer" "Alnus" "Betula" "Ulmus" "Tilia" "rubra" "robur" "Pseudotsuga" "Pinus" "Picea" "Larix" "Fagus" "Carpinus")

    # Submit the generic analysis job for each species
    for species in "${SPECIES_LIST[@]}"; do
        echo "Submitting $species"
        
        # Pass the species name to the job name and export it as an environment variable
        jobid=$(sbatch --partition=snowball --parsable \
               --job-name="${species}_boot" \
               --export=GS_DATA=$GS_DATA,NARROW_DATA=$NARROW_DATA,WIDE_DATA=$WIDE_DATA,TYPE=$TYPE,SPECIES=$species \
               hpc/scripts/generic_boot.sh | cut -d';' -f1)
               
        if [ -z "$jobid" ]; then
            echo "ERROR: submission failed for $species"
            exit 1
        fi
        echo "  Job ID: $jobid"
        ids+=("$jobid")
    done

    # Join job IDs with ':'
    dependency=$(IFS=:; echo "${ids[*]}")
    echo "Submitting compile job with dependency: $dependency"
    sbatch --dependency=afterok:"$dependency" hpc/scripts/compile_boots.sh "$TYPE"

}

submit_jobs \
    ./data/gs_data.csv \
    ./data/narrow_data.csv \
    ./data/wide_data.csv \
    full

submit_jobs \
    ./data/gs_data_restricted.csv \
    ./data/narrow_data_restricted.csv \
    ./data/wide_data_restricted.csv \
    restricted
  