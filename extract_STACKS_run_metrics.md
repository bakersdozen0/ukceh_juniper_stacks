#!/bin/bash

# --- CONFIG ---
BASE_DIR="/ssd0/jbaker/stacks/final_param_results"
output_dir="./run_summary_output"  # where to store combined results
mkdir -p "$output_dir"

combined_gstacks_file="$output_dir/combined_gstacks.csv"
combined_summary_file="$output_dir/combined_summary.csv"
combined_loci_metrics_file="$output_dir/combined_loci_metrics.csv"
combined_sample_distance_file="$output_dir/combined_sample_distances.tsv" # Renamed for clarity

# Clear old files
> "$combined_gstacks_file"
> "$combined_summary_file"
> "$combined_loci_metrics_file"
> "$combined_sample_distance_file"

# --- Parameter sets from the first script ---
datasets=("all_UK_samples" "combined_samples")
truncate_sizes=("50" "100")
parameter_sets=(
    "medium 2 3 1"
    "lenient 4 2 1"
    "very_lenient 6 2 2"
    "strict1 2 4 1"
    "strict2 2 4 2"
    "global_strict 2 5 3"
)

# --- Loop through runs ---
for dataset in "${datasets[@]}"; do
    for truncate_size in "${truncate_sizes[@]}"; do
        for param in "${parameter_sets[@]}"; do
            read label M m n <<< "$param"

            dir="${BASE_DIR}/ustacks_out_${dataset}_t${truncate_size}_${label}"
            if [ ! -d "$dir" ]; then
                echo "Directory $dir does not exist. Skipping."
                continue
            fi

            echo "Processing: $dataset, t=$truncate_size, label=$label"

            # ---------- GSTACKS ----------
            gstacks_file="${dir}/gstacks.log.distribs"
            if [ -f "$gstacks_file" ]; then
                awk -v t="$truncate_size" -v M="$M" -v m="$m" -v n="$n" -v dataset="$dataset" -v label="$label" '
                  BEGIN {OFS="\t"}
                  /^BEGIN effective_coverages_per_sample/ {section="effective_coverages_per_sample"; next}
                  /^END effective_coverages_per_sample/ {section=""; next}
                  /^BEGIN phasing_rates_per_sample/ {section="phasing_rates_per_sample"; next}
                  /^END phasing_rates_per_sample/ {section=""; next}
                  /^[^#]/ && section != "" {
                    print $0, t, M, m, n, dataset, label, section
                  }
                ' "$gstacks_file" >> "$combined_gstacks_file"
            fi

            # ---------- SUMMARY ----------
            summary_file="${dir}/populations.sumstats_summary.tsv"
            if [ -f "$summary_file" ]; then
                awk -v t="$truncate_size" -v M="$M" -v m="$m" -v n="$n" -v dataset="$dataset" -v label="$label" '
                  BEGIN {OFS="\t"}
                  /^# Variant positions/ {section="variant_positions"; next}
                  /^# All positions/ {section="all_positions"; next}
                  /^[^#]/ && section != "" {
                    print $0, t, M, m, n, dataset, label, section
                  }
                ' "$summary_file" >> "$combined_summary_file"
            fi

            # ---------- LOCI METRICS FROM VCF ----------
            vcf_file="${dir}/populations.snps.vcf"
            if [ -f "$vcf_file" ]; then
                tmp_prefix=$(mktemp) 

                vcftools --vcf "$vcf_file" --site-mean-depth --out "$tmp_prefix"
                vcftools --vcf "$vcf_file" --site-quality --out "$tmp_prefix"
                vcftools --vcf "$vcf_file" --missing-site --out "$tmp_prefix"

                if [ -f "$tmp_prefix.ldepth.mean" ] && [ -f "$tmp_prefix.lqual" ] && [ -f "$tmp_prefix.lmiss" ]; then
                    paste <(awk 'NR>1{print $1, $2, $3}' "$tmp_prefix.ldepth.mean") \
                          <(awk 'NR>1{print $3}' "$tmp_prefix.lqual") \
                          <(awk 'NR>1{print $5}' "$tmp_prefix.lmiss") \
                    | awk -v t="$truncate_size" -v M="$M" -v m="$m" -v n="$n" -v dataset="$dataset" -v label="$label" '
                        BEGIN {OFS="\t"}
                        {
                            chrom=$1; pos=$2; meandepth=$3; phred=$4; missrate=$5;
                            reproducibility=1-missrate;
                            print chrom, pos, meandepth, phred, missrate, reproducibility, t, M, m, n, dataset, label
                        }
                    ' >> "$combined_loci_metrics_file"
                else
                    echo "Warning: One or more vcftools output files not found for $dir"
                fi
                rm -f "$tmp_prefix".*
            fi

            # ---------- SAMPLE DISTANCE ----------
            if [ -f "$vcf_file" ]; then
                # --- FIX: Replaced deprecated '+distance' with the modern 'gtcheck -G 1' command ---
                bcftools gtcheck -G 1 "$vcf_file" \
                | awk -v t="$truncate_size" -v M="$M" -v m="$m" -v n="$n" -v dataset="$dataset" -v label="$label" \
                'NR>1 {print $0, t, M, m, n, dataset, label}' \
                >> "$combined_sample_distance_file"
            fi
        done
    done
done

echo "Data processing complete."