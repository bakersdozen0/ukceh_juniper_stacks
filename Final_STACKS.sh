#!/bin/bash
export PATH="$HOME/stacks-2.68-install/bin:$PATH"

# Define datasets
datasets=("all_UK_samples" "combined_samples")
#"combined_samples"
 

declare -A fastq_dirs
fastq_dirs["all_UK_samples"]="/ssd0/jbaker/stacks/combined_samples/concat_fastq_UK/" 
fastq_dirs["combined_samples"]="/ssd0/jbaker/stacks/combined_samples/concat_fastq"


# Define popmaps and catalog files
declare -A popmaps
popmaps["all_UK_samples"]="/ssd0/jbaker/stacks/UK_popmap.txt"
popmaps["combined_samples"]="/ssd0/jbaker/stacks/all_sample_pop_map"

declare -A catalogs
catalogs["all_UK_samples"]="/ssd0/jbaker/stacks/UK_subset_cat_updated.txt"
catalogs["combined_samples"]="/ssd0/jbaker/stacks/All_Samples_cat_subset.txt"

# Define truncate sizes
truncate_sizes=("100")
# "50" or "100"

# Define parameter sets as: "label M m n"
parameter_sets=(
   #"global_strict 2 5 3"
    "strict2 2 4 2"
    #"strict1 2 4 1"
    #"medium 2 3 1"
    #"lenient 4 2 1"
   #"very_lenient 6 2 2"   
)

# set ulimit to increase memory 

ulimit -n 4096

# Base directory
BASE_DIR="/ssd0/jbaker/stacks/final_param_results"

echo "Verifying and formatting input text files..."

# Create a list containing all the popmap and catalog file paths
files_to_format=("${popmaps[@]}" "${catalogs[@]}")

for file_path in "${files_to_format[@]}"; do
    # First, check if the file actually exists
    if [ -f "$file_path" ]; then
        echo "  - Formatting: $file_path"
        # The sed command finds any carriage return (\r) at the very end ($)
        # of a line and replaces it with nothing, effectively deleting it.
        # The '-i' flag means the file is edited "in-place".
        sed -i 's/\r$//' "$file_path"
    else
        # Print a warning if a defined file is not found
        echo "  - WARNING: File not found, cannot format: $file_path"
    fi
done
echo "File formatting check complete."

# Loop through each dataset
for dataset in "${datasets[@]}"; do
    
    # Get the correct input directory and popmap for the current dataset
    input_dir=${fastq_dirs[$dataset]}
    popmap_file=${popmaps[$dataset]}

    # Check if the input directory exists
    if [ ! -d "$input_dir" ]; then
        echo "Error: Input FASTQ directory not found for dataset '$dataset': $input_dir"
        continue # Skip to the next dataset
    fi

    # Loop through each truncation size
    for truncate_size in "${truncate_sizes[@]}"; do
        
        # Define the output directory for this run's processed files
        proc_output="$BASE_DIR/proc_radtag_${dataset}_t${truncate_size}"
        
        # Run process_radtags on each file individually
        # This step will only run if the output directory doesn't already exist.

        if [ ! -d "$proc_output" ]; then
            mkdir -p "$proc_output"
            echo "Running process_radtags individually for all files in $input_dir (t=$truncate_size)..."
            
            # Loop through every FASTQ file in the input directory
            for fastq_file in "$input_dir"/*.fastq.gz; do
                # Print a message to show progress
                echo "  -> Processing: $(basename "$fastq_file")"
                
                process_radtags \
                    -f "$fastq_file" \
                    -o "$proc_output" \
                    --renz-1 pstI --renz-2 mseI \
                    -c \
                    -q \
                    -r \
                    --truncate "$truncate_size"
            done
        else
            echo "Skipping process_radtags for $dataset (t=$truncate_size), output directory already exists."
        fi


        # for param in "${parameter_sets[@]}"; do
        #     read label M m n <<< "$param"

        #     ustacks_output="$BASE_DIR/ustacks_out_${dataset}_t${truncate_size}_${label}"
        #     mkdir -p "$ustacks_output"

        #     echo "Running denovo_map.pl for $dataset (t=$truncate_size, label=$label, M=$M, m=$m, n=$n)..."

        #     denovo_map.pl --samples "$proc_output"\
        #         --popmap "${popmaps[$dataset]}" \
        #         -o "$ustacks_output" -M "$M" -n "$n" -r 0.8 \
        #         --catalog-popmap "${catalogs[$dataset]}" -T 54 \
        #         -X "ustacks: -m $m" -X "populations: --vcf"

        #     touch "$BASE_DIR/populations_${dataset}_t${truncate_size}_${label}.done"
        # done

        #Take 2 with "populations" filtering steps: 
        for param in "${parameter_sets[@]}"; do
            read label M m n <<< "$param"

            # EDIT 1: Add "_filtered" to the output directory name to save results separately
            ustacks_output="$BASE_DIR/ustacks_out_${dataset}_t${truncate_size}_${label}_filtered"
            mkdir -p "$ustacks_output"

            echo "Running denovo_map.pl for $dataset (t=$truncate_size, label=$label, M=$M, m=$m, n=$n)..."

            denovo_map.pl --samples "$proc_output"\
                --popmap "${popmaps[$dataset]}" \
                -o "$ustacks_output" -M "$M" -n "$n" -r 0.8 \
                --catalog-popmap "${catalogs[$dataset]}" -T 54 \
                -X "ustacks: -m $m" \
                -X "populations: --vcf -r 0.75 -p 3 --min-mac 2 --hwe" # <-- EDIT 2: Add filtering flags

            touch "$BASE_DIR/populations_${dataset}_t${truncate_size}_${label}_filtered.done"
        done


    done
done

echo "Denovo mapping pipeline completed!"
