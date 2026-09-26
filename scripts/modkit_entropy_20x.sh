#!/usr/bin/env bash
set -Eeuo pipefail
#
# Rapid methylation-entropy workflow:
#   1. Estimate whole-reference mean depth for each input BAM.
#   2. Randomly subsample high-depth BAMs to approximately 15x each.
#   3. Index the temporary BAMs.
#   4. Run modkit entropy on all nine inputs.
#   5. Delete temporary BAMs automatically, including after failure/interruption.
#
# Requirements: bash, awk, samtools, modkit
#
# Notes:
# - Subsampling is global and random, so 15x is an expected whole-genome mean,
#   not an exact cap of 15 reads at every position.
# - MM/ML/MN modification tags are retained by samtools view.
# - Only mapped primary alignments are written to temporary BAMs, matching the
#   primary-read analysis intended here and reducing temporary I/O.
# - --num-reads is deliberately not set: it only controls Modkit's preliminary
#   filtering-threshold estimation, not entropy-window coverage.

REFERENCE="/home/luke/Documents/data/E.faecium_C68/E.faecium-C68-UoB2025.fasta"

TARGET_COV_PER_BAM=20
SAMTOOLS_THREADS=30
MODKIT_THREADS=30
SEED_BASE=260713

# Change this if /tmp lacks space, for example:
# TMP_ROOT="/home/luke/Documents/data/tmp"
TMP_ROOT="${TMPDIR:-/tmp}"

SAMPLES=(HC1 HC2 HC3 VC1 VC2 VC3 AC1 AC2 AC3)


die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}


have_bam_index() {
    local bam=$1
    [[ -s "${bam}.bai" ||
       -s "${bam%.bam}.bai" ||
       -s "${bam}.csi" ||
       -s "${bam%.bam}.csi" ]]
}


# Estimate mean primary-alignment depth over the complete reference.
# samtools depth is used for compatibility with older samtools releases.
#
# Only mapped primary alignments are counted, matching the temporary-BAM filter.
# Zero-coverage positions are represented through the full-reference denominator,
# so emitting every zero-depth position is unnecessary.
estimate_mean_depth() {
    local bam=$1
    local genome_length=$2

    samtools depth \
        -d 0 \
        -G 0x904 \
        "$bam" |
    awk -v genome_length="$genome_length" '
        {
            depth_bases += $3
        }
        END {
            if (genome_length <= 0) {
                exit 2
            }
            printf "%.10f\n", depth_bases / genome_length
        }
    '
}


# Arguments:
#   1: modification label, e.g. 6mA, 5mC, or 4mC
#   2: input path
#   3: input layout: "nested" or "flat"
#   4: primary base: A or C
run_entropy_15x() (
    set -Eeuo pipefail

    local mod=${1:?Missing modification label}
    local input_path=${2:?Missing input path}
    local layout=${3:?Missing input layout}
    local base=${4:?Missing primary base}

    command -v samtools >/dev/null 2>&1 || die "samtools is not in PATH"
    command -v modkit   >/dev/null 2>&1 || die "modkit is not in PATH"
    command -v awk      >/dev/null 2>&1 || die "awk is not in PATH"

    [[ -r "$REFERENCE" ]] || die "Reference FASTA is not readable: $REFERENCE"
    [[ "$base" == "A" || "$base" == "C" ]] ||
        die "Base must be A or C, received: $base"
    [[ "$layout" == "nested" || "$layout" == "flat" ]] ||
        die "Layout must be nested or flat, received: $layout"
    [[ -d "$TMP_ROOT" && -w "$TMP_ROOT" ]] ||
        die "Temporary root is not a writable directory: $TMP_ROOT"

    local out="/home/luke/Documents/data/E.faecium_C68/methyl/#NEW/UPGRADE/${mod}_ME"
    local bed="${out}/50bp_regions.bed"
    local sampling_log="${out}/subsampling_15x.tsv"
    local tmpdir

    mkdir -p "$out"

    if [[ ! -s "${REFERENCE}.fai" ]]; then
        samtools faidx "$REFERENCE"
    fi

    awk 'BEGIN { OFS="\t" }
        {
            for (start = 0; start < $2; start += 50) {
                end = start + 50
                if (end > $2) {
                    end = $2
                }
                print $1, start, end
            }
        }
    ' "${REFERENCE}.fai" > "$bed"

    local genome_length
    genome_length=$(
        awk '{ total += $2 }
             END {
                 if (total <= 0) exit 1
                 printf "%.0f\n", total
             }' "${REFERENCE}.fai"
    ) || die "Could not calculate reference length from ${REFERENCE}.fai"

    tmpdir=$(mktemp -d -p "$TMP_ROOT" "modkit_${mod}_15x.XXXXXXXX")
    cleanup() {
        rm -rf -- "$tmpdir"
    }
    trap cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP

    printf 'Temporary directory: %s\n' "$tmpdir" >&2
    printf 'sample\tinput_bam\toriginal_mean_depth\ttarget_depth\tsampling_fraction\texpected_depth\tmodkit_input\n' \
        > "$sampling_log"

    local -a in_args=()
    local sample bam mean_depth fraction expected_depth tmpbam seed sample_spec
    local index=0

    for sample in "${SAMPLES[@]}"; do
        case "$layout" in
            nested)
                bam="${input_path}/${sample}/${mod}_ref_align-sort.bam"
                ;;
            flat)
                bam="${input_path}/${sample}.bam"
                ;;
        esac

        [[ -r "$bam" ]] || die "Input BAM is not readable: $bam"
        samtools quickcheck -v "$bam" ||
            die "samtools quickcheck failed for: $bam"

        printf '[%s] Estimating mean depth...\n' "$sample" >&2
        mean_depth=$(estimate_mean_depth "$bam" "$genome_length") ||
            die "Coverage estimation failed for: $bam"

        # Reject malformed output before using it in Bash printf or arithmetic.
        [[ "$mean_depth" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
            die "Coverage estimator returned a non-numeric value for $bam: <$mean_depth>"

        awk -v depth="$mean_depth" 'BEGIN { exit !(depth > 0) }' ||
            die "Estimated depth is zero for: $bam"

        fraction=$(
            awk -v target="$TARGET_COV_PER_BAM" -v depth="$mean_depth" '
                BEGIN {
                    fraction = target / depth
                    if (fraction > 1) fraction = 1
                    if (fraction < 0) fraction = 0
                    printf "%.12f\n", fraction
                }
            '
        )

        expected_depth=$(
            awk -v depth="$mean_depth" -v fraction="$fraction" '
                BEGIN { printf "%.4f\n", depth * fraction }
            '
        )

        seed=$((SEED_BASE + index))
        index=$((index + 1))

        # Avoid copying a BAM that is already at or below the requested depth.
        if awk -v fraction="$fraction" 'BEGIN { exit !(fraction >= 0.999999999999) }'
        then
            printf '[%s] %.2fx is <= target %.2fx; using the original BAM.\n' \
                "$sample" "$mean_depth" "$TARGET_COV_PER_BAM" >&2

            if ! have_bam_index "$bam"; then
                printf '[%s] Creating missing BAM index...\n' "$sample" >&2
                samtools index -@ "$SAMTOOLS_THREADS" "$bam"
            fi

            in_args+=(--in-bam "$bam")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$sample" "$bam" "$mean_depth" "$TARGET_COV_PER_BAM" \
                "$fraction" "$expected_depth" "$bam" >> "$sampling_log"
            continue
        fi

        tmpbam="${tmpdir}/${sample}.${mod}.approximately_15x.bam"

        printf '[%s] %.2fx -> expected %.2fx; retaining fraction %s\n' \
            "$sample" "$mean_depth" "$expected_depth" "$fraction" >&2

        # -1 uses fast BAM compression.
        # -F 0x904 removes unmapped, secondary and supplementary records.
        # Subsampling retains a random fraction of templates using a fixed seed.
        # Older samtools releases use:
        #     -s INT.FRAC
        # where INT is the random seed and 0.FRAC is the retained fraction.
        # Example: seed=260713 and fraction=0.0769 -> -s 260713.0769
        sample_spec="${seed}.${fraction#0.}"

        samtools view \
            -1 \
            -@ "$SAMTOOLS_THREADS" \
            -F 0x904 \
            -s "$sample_spec" \
            -o "$tmpbam" \
            "$bam"

        samtools quickcheck -v "$tmpbam" ||
            die "Temporary BAM failed quickcheck: $tmpbam"

        samtools index -@ "$SAMTOOLS_THREADS" "$tmpbam"
        in_args+=(--in-bam "$tmpbam")

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$sample" "$bam" "$mean_depth" "$TARGET_COV_PER_BAM" \
            "$fraction" "$expected_depth" "$tmpbam" >> "$sampling_log"
    done

    printf 'Running Modkit on %d inputs...\n' "${#SAMPLES[@]}" >&2

    # Remove only prior Modkit result files so a rerun cannot be confused
    # with stale output. The BED and subsampling audit TSV are preserved.
    rm -f -- "${out}/regions.bed" "${out}/windows.bedgraph"

    modkit entropy \
        "${in_args[@]}" \
        --ref "$REFERENCE" \
        --regions "$bed" \
        --out-bed "$out" \
        --threads "$MODKIT_THREADS" \
        --base "$base"

    printf 'Completed: %s\n' "$out" >&2
    printf 'Subsampling record: %s\n' "$sampling_log" >&2
    printf 'Temporary BAMs are now being deleted.\n' >&2
)


# 6mA: BAMs are stored as:
#   <path>/<sample>/6mA_ref_align-sort.bam
run_entropy_15x \
    6mA \
    "/home/luke/Documents/data/E.faecium_C68/methyl" \
    nested \
    A

# 5mC and 4mC: BAMs are stored as:
#   <path>/<sample>.bam
run_entropy_15x \
    5mC \
    "/home/luke/Documents/data/E.faecium_C68/methyl/#NEW/pileup/5mC" \
    flat \
    C

run_entropy_15x \
    4mC \
    "/home/luke/Documents/data/E.faecium_C68/methyl/#NEW/pileup/4mC" \
    flat \
    C
