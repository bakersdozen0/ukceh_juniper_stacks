# Concat FASTQ files from same indivudals:  
# The '-p' flag means it won't complain if the directory already exists.
mkdir -p ./stacks_ready_fastq

echo "Combining FASTQ files based on the master mapping file..."

# 1. Read the master mapping file, skipping the header line (`tail -n +2`).
# 2. For each line, read the two comma-separated values into the variables TARGETID and GENOTYPE.
tail -n +2 master_genotype_to_targetid_mapping.csv | while IFS=, read -r TARGETID GENOTYPE
do
  # Define the name of the source file we are looking for (e.g., "4458805.fastq.gz")
  SOURCE_FILE="${TARGETID}.fastq.gz"
  
  # Define the name of the destination file we want to create/add to (e.g., "NUU-JC-12.fastq.gz")
  DEST_FILE="./stacks_ready_fastq/${GENOTYPE}.fastq.gz"

  # Safety check: make sure the source file actually exists
  if [ -f "$SOURCE_FILE" ]; then
    echo "Appending ${SOURCE_FILE}  >>  ${DEST_FILE}"
    
    # Use `cat` to read the source file and `>>` to APPEND its contents 
    # to the destination file. This is what combines the replicates.
    cat "$SOURCE_FILE" >> "$DEST_FILE"
  else
    echo "Warning: Source file ${SOURCE_FILE} not found. Skipping."
  fi
done

echo "Done! Your combined files are in the stacks_ready_fastq directory."