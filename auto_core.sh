#!/bin/bash
# Author: Luke B Jones
#

#############################################
# PROJECT / INPUT PATHS
#############################################
# working directory
WD=""
SCRIPT_DIR="${WD}/scripts"
THREADS=32

# fasta
GENOME="${WD}/"
fai="${GENOME}.fai"

# Input modification BAM folders.
A_bam="${WD}/demux/6mA/merged_bams"
C_bam="${WD}/demux/4mC_5mC/merged_bams"

# Transcriptomic FASTQs.
trans="${WD}/transcript/raw"
loc="${WD}"


# Genome-feature / annotation files.
PGAP="${WD}/genomic-features/PGAP/annot.gff"
ncRNA="${WD}/genomic-features/cmsearch/ncRNA_all.bed"
gene_an="${WD}/genomic-features/hybrid_genes.gff"
hybrid_bed="${WD}/genomic-features/hybrid_genes.bed"
operons="${WD}/genomic-features/operon-mapper/operons.tsv"
operon_gff="${WD}/genomic-features/operon-mapper/operon-mapper_annot.gff"
operon_list="${WD}/genomic-features/operon-mapper/operon_list.tsv"
dna_annot="${WD}/genomic-features/operon-mapper/DNA_annot.fasta"

# Other study inputs / resources.
PubMETA="${WD}/MTase/MTase-conservation/PubMLST_metadata.tsv"
MTase_final="${WD}/MTase/MTase-search/3.LOG2FC.tsv"
KO="${WD}/genomic-features/KO_KEGG_BlastKOALA.tsv"
RFAM_CM="${WD}/resources/Rfam.cm"
REBASE="${WD}/resources/both_rebase.fasta"
PUBMLST_CONTIGS="${WD}/MTase/MTase-conservation/pubmlst_contigs"

# Tn analysis inputs.
tn_free="${WD}/Tn-free-genome.fasta"
tn_contig="${WD}/Tn6085.fasta"
tn_name="Tn6085"
control_fastq="${WD}/complete_control.fastq"

# Plot colours retained from the later 4mC/5mC analysis.
H_col="#1b9e77,#8dcfbb"
V_col="#7570b3,#bab8d9"
A_col="#d95f02,#ec8780"


#############################################
# SCRIPTS USED BY THIS ANALYSIS
# All local helper scripts are expected in ${WD}/scripts/
#############################################

pgap_runner="${SCRIPT_DIR}/pgap.py"
operon_mapper="${SCRIPT_DIR}/PGAP_GFF_to_operon-mapper.py"
ncRNA_bed_script="${SCRIPT_DIR}/ncRNA_bed.py"
gff_hybrid="${SCRIPT_DIR}/hybridise_gffs.py"
operon_mapper_conv="${SCRIPT_DIR}/operon-mapper_to_bed.py"
upstream_find="${SCRIPT_DIR}/find_promoters.py"

methyl_stats="${SCRIPT_DIR}/methyl_stats.py"
methylome_plot="${SCRIPT_DIR}/methylome_plot.py"   # original pre-split global plot helper
fig3="${SCRIPT_DIR}/fig.3.py"
fig4="${SCRIPT_DIR}/genome_methy_graph.py"
overlap_bases="${SCRIPT_DIR}/sig_mods.py"
meta_py="${SCRIPT_DIR}/maisem_table.py"
conver_OM="${SCRIPT_DIR}/convert_OM_to_index.py"

deseq2="${SCRIPT_DIR}/run_deseq2_featurecounts.R"
volcano="${SCRIPT_DIR}/volcano_plots.py"
single_volcano="${SCRIPT_DIR}/single_volcano.py"  # defined in the original main script
ko_log="${SCRIPT_DIR}/ko_log2fc_kegg_merge.py"
cluster_kegg="${SCRIPT_DIR}/cluster_KEGG.R"
kegg_plot="${SCRIPT_DIR}/KEGG.py"
kegg_distro="${SCRIPT_DIR}/plot_log_ditr-kegg.py"

consen_prot="${SCRIPT_DIR}/blastp_tsv_consensus.py"
pgap_check="${SCRIPT_DIR}/PGAP_methyl_filter.py"
interpro="${SCRIPT_DIR}/intrpo_web-submission.py"
inter_consen="${SCRIPT_DIR}/mtase_consensus.py"
mt_con="${SCRIPT_DIR}/mtase_consensus_rebase_inter_pgap.py"

sum_VCF="${SCRIPT_DIR}/sumarise_clair3.py"
intersect_VCF="${SCRIPT_DIR}/intersect_VCF.py"

gff_bed="${SCRIPT_DIR}/gff-to-bed.py"
calc_2="${SCRIPT_DIR}/entropy_express_calc.py"
calc_3="${SCRIPT_DIR}/entropy_stats_calc_3.py"
base_en_stat="${SCRIPT_DIR}/entro_stat.py"
zero_mod="${SCRIPT_DIR}/count_zero_mod.py"
count_unmethyl="${SCRIPT_DIR}/never_methylated.py"

background="${SCRIPT_DIR}/background_regions.py"
hypergeom="${SCRIPT_DIR}/hypergeom.py"
depth="${SCRIPT_DIR}/DOM_STAN_coverage.py"

pub_extract="${SCRIPT_DIR}/PubMLST_to_database.py"
qurey_maker="${SCRIPT_DIR}/blastn_mtase_qurey.py"
bivariat="${SCRIPT_DIR}/bivariat_plot.py"

tn_extr="${SCRIPT_DIR}/extract_overhangs_TE.py"
tn_count="${SCRIPT_DIR}/tn-search.py"
tree="${SCRIPT_DIR}/tree_phylo.py"

# Referenced by the supplied control-control validation shell after its DMR runs.
control_validation="${SCRIPT_DIR}/section2_replicate_validation.py"

# Later control-state ME / subsequent methylation-response analysis.
me_control_entropy="${SCRIPT_DIR}/modkit_entropy_20x.sh"
me_node_enrichment="${SCRIPT_DIR}/me_node_enrichment_v4.py"


#############################################
# (1) SETUP METHYLATION DATA
#############################################

# Combined 4mC/5mC calls were the original input. The later correction split the
# resulting cytosine BED by modification identity for DMR analysis.

mkdir -p \
  "${loc}/methyl/4mC_5mC" \
  "${loc}/methyl/6mA" \
  "${loc}/methyl/4mC/bin" \
  "${loc}/methyl/5mC/bin" \
  "${loc}/methyl/#NEW/4mC" \
  "${loc}/methyl/#NEW/5mC" \
  "${loc}/methyl/#NEW/SIG/4mC" \
  "${loc}/methyl/#NEW/SIG/5mC"

meth_fun(){
    local name="$1"
    local base_dir="${loc}"

    mkdir -p "${base_dir}/methyl/${name}"

    # Align reads using the modification-aware BAMs used in the study.
    dorado aligner "${GENOME}" "${C_bam}/${name}.bam" \
      -o "${base_dir}/methyl/${name}/4mC_5mC_align" --threads "${THREADS}"

    dorado aligner "${GENOME}" "${A_bam}/${name}.bam" \
      -o "${base_dir}/methyl/${name}/6mA_align" --threads "${THREADS}"

    # Sort and index aligned BAMs.
    samtools sort -@ "${THREADS}" \
      -o "${base_dir}/methyl/${name}/4mC_5mC_ref_align-sort.bam" \
      "${base_dir}/methyl/${name}/4mC_5mC_align/${name}.bam"

    samtools sort -@ "${THREADS}" \
      -o "${base_dir}/methyl/${name}/6mA_ref_align-sort.bam" \
      "${base_dir}/methyl/${name}/6mA_align/${name}.bam"

    samtools index -@ "${THREADS}" -b \
      "${base_dir}/methyl/${name}/4mC_5mC_ref_align-sort.bam"

    samtools index -@ "${THREADS}" -b \
      "${base_dir}/methyl/${name}/6mA_ref_align-sort.bam"

    # Original pileups.
    modkit pileup --ref "${GENOME}" --threads "${THREADS}" \
      "${base_dir}/methyl/${name}/4mC_5mC_ref_align-sort.bam" \
      "${base_dir}/methyl/${name}/4mC_5mC.bed"

    modkit pileup --ref "${GENOME}" --threads "${THREADS}" \
      "${base_dir}/methyl/${name}/6mA_ref_align-sort.bam" \
      "${base_dir}/methyl/${name}/6mA.bed"

    # Preserve the original consolidated outputs.
    cp "${base_dir}/methyl/${name}/4mC_5mC.bed" \
       "${base_dir}/methyl/4mC_5mC/${name}.bed"
    cp "${base_dir}/methyl/${name}/6mA.bed" \
       "${base_dir}/methyl/6mA/${name}.bed"

    bgzip -@ "${THREADS}" "${base_dir}/methyl/4mC_5mC/${name}.bed"
    bgzip -@ "${THREADS}" "${base_dir}/methyl/6mA/${name}.bed"
    tabix -@ "${THREADS}" -p bed "${base_dir}/methyl/4mC_5mC/${name}.bed.gz"
    tabix -@ "${THREADS}" -p bed "${base_dir}/methyl/6mA/${name}.bed.gz"

    # Later 4mC / 5mC correction used for DMR analysis.
    # 5mC
    awk -F'\t' '$4=="m"' \
      "${base_dir}/methyl/${name}/4mC_5mC.bed" \
      > "${base_dir}/methyl/5mC/bin/${name}.bed"

    # 4mC
    awk -F'\t' '$4=="21839"' \
      "${base_dir}/methyl/${name}/4mC_5mC.bed" \
      > "${base_dir}/methyl/4mC/bin/${name}.bed"

    bgzip -@ "${THREADS}" "${base_dir}/methyl/5mC/bin/${name}.bed"
    bgzip -@ "${THREADS}" "${base_dir}/methyl/4mC/bin/${name}.bed"
    tabix -p bed "${base_dir}/methyl/5mC/bin/${name}.bed.gz"
    tabix -p bed "${base_dir}/methyl/4mC/bin/${name}.bed.gz"

    # Coverage review used in the original analysis.
    samtools depth -aa -d 0 \
      "${base_dir}/methyl/${name}/4mC_5mC_ref_align-sort.bam" \
      > "${base_dir}/methyl/${name}/coverage_4mC_5mC.txt"

    samtools depth -aa -d 0 \
      "${base_dir}/methyl/${name}/6mA_ref_align-sort.bam" \
      > "${base_dir}/methyl/${name}/coverage_6mA.txt"
}

# Run manually when required:
#for sample in HT1 HT2 HT3 HC1 HC2 HC3 VT1 VT2 VT3 VC1 VC2 VC3 AT1 AT2 AT3 AC1 AC2 AC3; do
#    meth_fun "$sample"
#done


#############################################
# BASE-LEVEL RESOLUTION BEDS
#############################################

# These A.bed and C.bed files are the same region files used for:
#   1. the original treatment-vs-control DMR analysis,
#   2. the later split 4mC / 5mC DMR analysis,
#   3. the control-vs-control validation.

#seqkit locate -r -p A "${GENOME}" > "${loc}/methyl/A_gff.bed"
#seqkit locate -r -p C "${GENOME}" > "${loc}/methyl/C_gff.bed"

#awk 'BEGIN {OFS="\t"}
#     /^seqID/ || NF < 7 { next }
#     {
#       if ($5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) { next }
#       chrom       = $1
#       strand      = $4
#       chromStart  = $5 - 1
#       chromEnd    = $6
#       name        = $7
#       score       = "."
#       print chrom, chromStart, chromEnd, name, score, strand
#     }' "${loc}/methyl/A_gff.bed" | sort -k1,1 -k2,2n \
#     > "${loc}/methyl/1DMbase/A.bed"

#awk 'BEGIN {OFS="\t"}
#     /^seqID/ || NF < 7 { next }
#     {
#       if ($5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) { next }
#       chrom       = $1
#       strand      = $4
#       chromStart  = $5 - 1
#       chromEnd    = $6
#       name        = $7
#       score       = "."
#       print chrom, chromStart, chromEnd, name, score, strand
#     }' "${loc}/methyl/C_gff.bed" | sort -k1,1 -k2,2n \
#     > "${loc}/methyl/1DMbase/C.bed"

#rm "${loc}/methyl/A_gff.bed" "${loc}/methyl/C_gff.bed"

# Optional check retained from the original script:
#numA=$(wc -l < "${loc}/methyl/1DMbase/A.bed")
#numC=$(wc -l < "${loc}/methyl/1DMbase/C.bed")
#if [ "$numA" -lt 5000 ] || [ "$numC" -lt 5000 ]; then
#    echo -e "\e[31mError: Insufficient base calls (A: $numA, C: $numC). Exiting.\e[0m"
#    exit 1
#fi


#############################################
# BASE-LEVEL DMR ANALYSIS: 6mA
# Original A-modification analysis retained.
#############################################

base_dmr_6mA(){
    local gp="$1"
    local group=""

    if [ "$gp" = "H" ]; then
        group="H202"
    elif [ "$gp" = "V" ]; then
        group="VANC"
    elif [ "$gp" = "A" ]; then
        group="AMPI"
    else
        echo "Error: group parameter must be H, V, or A"
        return 1
    fi

    local dir="${loc}/methyl/6mA"
    local out="${loc}/methyl/1DMbase/${group}_6mA"
    local region="${loc}/methyl/1DMbase/A.bed"

    mkdir -p "${loc}/methyl/1DMbase/20_pct"

    modkit dmr multi \
      -s "${dir}/${gp}C1.bed.gz" control \
      -s "${dir}/${gp}C2.bed.gz" control \
      -s "${dir}/${gp}C3.bed.gz" control \
      -s "${dir}/${gp}T1.bed.gz" treatment \
      -s "${dir}/${gp}T2.bed.gz" treatment \
      -s "${dir}/${gp}T3.bed.gz" treatment \
      --base A -o "${out}" --ref "${GENOME}" \
      --regions-bed "${region}" --header --missing warn -t "${THREADS}"

    awk 'BEGIN {
           FS = "\t"; OFS = "\t"
         }
         NR == 1 {
           print $0, "differential_pct_modification"; next
         }
         {
           control_pct_modified = $13 + 0;
           treatment_pct_modified = $14 + 0;
           differential_pct_modification = treatment_pct_modified - control_pct_modified;
           print $0, differential_pct_modification
         }' "${out}/control_treatment.bed" \
         > "${out}/${group}_6mA_DIFF.bed"

    awk 'BEGIN{FS=OFS="\t"}
         NR==1 {print $0; next}
         {
           ctot = $8 + 0;
           ttot = $10 + 0;
           diff = $NF + 0;
           if (ctot >= 10 && ttot >= 10 && (diff > 0.2 || diff < -0.2)) {
             print $0
           }
         }' "${out}/${group}_6mA_DIFF.bed" \
         > "${loc}/methyl/1DMbase/20_pct/${group}_6mA.bed"
}

# Run manually:
#base_dmr_6mA H
#base_dmr_6mA V
#base_dmr_6mA A


#############################################
# BASE-LEVEL DMR ANALYSIS: 4mC / 5mC
# Later split correction retained as performed.
#############################################

base_dmr_C(){
    local gp="$1"
    local group=""

    if [ "$gp" = "H" ]; then
        group="H202"
    elif [ "$gp" = "V" ]; then
        group="VANC"
    elif [ "$gp" = "A" ]; then
        group="AMPI"
    else
        echo "Error: group parameter must be H, V, or A"
        return 1
    fi

    dmr_fun(){
        local mod="$1"
        local out="${loc}/methyl/#NEW/${mod}/${group}"

        mkdir -p "${out}" "${loc}/methyl/#NEW/SIG/${mod}"

        modkit dmr multi \
          -s "${loc}/methyl/${mod}/bin/${gp}C1.bed.gz" control \
          -s "${loc}/methyl/${mod}/bin/${gp}C2.bed.gz" control \
          -s "${loc}/methyl/${mod}/bin/${gp}C3.bed.gz" control \
          -s "${loc}/methyl/${mod}/bin/${gp}T1.bed.gz" treatment \
          -s "${loc}/methyl/${mod}/bin/${gp}T2.bed.gz" treatment \
          -s "${loc}/methyl/${mod}/bin/${gp}T3.bed.gz" treatment \
          --base C -o "${out}" --ref "${GENOME}" \
          --regions-bed "${loc}/methyl/1DMbase/C.bed" \
          --header --missing warn -t "${THREADS}"

        awk 'BEGIN {
               FS = "\t"; OFS = "\t"
             }
             NR == 1 {
               print $0, "differential_pct_modification"; next
             }
             {
               control_pct_modified = $13 + 0;
               treatment_pct_modified = $14 + 0;
               differential_pct_modification = treatment_pct_modified - control_pct_modified;
               print $0, differential_pct_modification
             }' "${out}/control_treatment.bed" \
             > "${out}/control_treatment_DIFF.bed"

        awk 'BEGIN{FS=OFS="\t"}
             NR==1 {print $0; next}
             {
               ctot = $8 + 0;
               ttot = $10 + 0;
               diff = $NF + 0;
               if (ctot >= 10 && ttot >= 10 && (diff > 0.2 || diff < -0.2)) {
                 print $0
               }
             }' "${out}/control_treatment_DIFF.bed" \
             > "${loc}/methyl/#NEW/SIG/${mod}/${group}.bed"
    }

    dmr_fun 4mC
    dmr_fun 5mC
}

# Run manually:
#base_dmr_C H
#base_dmr_C V
#base_dmr_C A


#############################################
# CONTROL-vs-CONTROL DMR VALIDATION
#############################################
#
# The supplied validation workflow compared all three control groups pairwise.
# For E. faecium these are:
#   H = H2O2 control replicates (HC1-HC3)
#   V = VANC control replicates (VC1-VC3)
#   A = AMPI control replicates (AC1-AC3)
#
# As in the supplied validation script, one control set is passed to modkit as
# "control" and the second control set as "treatment". No new A/C region BEDs are
# generated: the exact A.bed and C.bed made above are reused.

CC_DIR="${loc}/methyl/#NEW/CONTROL_CONTROL"
CC_TMP="${CC_DIR}/test"
mkdir -p "${CC_TMP}"

control_control_dmr(){
    local mod="$1"
    local base="$2"
    local left="$3"
    local right="$4"
    local label="$5"

    local input_dir=""
    local regions=""

    if [ "$mod" = "6mA" ]; then
        input_dir="${loc}/methyl/6mA"
        regions="${loc}/methyl/1DMbase/A.bed"
    else
        input_dir="${loc}/methyl/${mod}/bin"
        regions="${loc}/methyl/1DMbase/C.bed"
    fi

    modkit dmr multi \
      -s "${input_dir}/${left}C1.bed.gz" control \
      -s "${input_dir}/${left}C2.bed.gz" control \
      -s "${input_dir}/${left}C3.bed.gz" control \
      -s "${input_dir}/${right}C1.bed.gz" treatment \
      -s "${input_dir}/${right}C2.bed.gz" treatment \
      -s "${input_dir}/${right}C3.bed.gz" treatment \
      --base "${base}" \
      -o "${CC_TMP}" \
      --ref "${GENOME}" \
      --regions-bed "${regions}" \
      --header \
      --missing warn \
      -t "${THREADS}"

    mv "${CC_TMP}/control_treatment.bed" \
       "${CC_DIR}/${label}_${mod}.bed"
}

# Pair order follows the three pairwise comparisons used in the supplied
# control-control validation: first-vs-second, first-vs-third, third-vs-second.
#
# 4mC
#control_control_dmr 4mC C H V HxV
#control_control_dmr 4mC C H A HxA
#control_control_dmr 4mC C A V AxV
#
# 5mC
#control_control_dmr 5mC C H V HxV
#control_control_dmr 5mC C H A HxA
#control_control_dmr 5mC C A V AxV
#
# 6mA
#control_control_dmr 6mA A H V HxV
#control_control_dmr 6mA A H A HxA
#control_control_dmr 6mA A A V AxV

# The supplied S. aureus control-validation shell then called
# section2_replicate_validation.py. That helper was not present in the supplied
# archive, and the only saved invocation is specific to H202/METH/SNAP, so an
# Enterococcus invocation is not invented here.


#############################################
# UPDATED METHYLATION FIGURES / SUMMARY STATS
#############################################

mkdir -p \
  "${loc}/methyl/#NEW/fig.3" \
  "${loc}/methyl/#NEW/fig.4/5mC" \
  "${loc}/methyl/#NEW/fig.4/4mC" \
  "${loc}/methyl/#NEW/fig.4/6mA" \
  "${loc}/methyl/#NEW/STAT" \
  "${loc}/methyl/#NEW/SIG/6mA"

# Figure 3 rerun after the cytosine split.
#python "$fig3" \
#  -png "${loc}/methyl/#NEW/fig.3/summary.png" \
#  -o "${loc}/methyl/#NEW/fig.3" \
#  -i "${loc}/methyl/#NEW/SIG/4mC/H202.bed" H202 4mC "$H_col" \
#  -i "${loc}/methyl/#NEW/SIG/5mC/H202.bed" H202 5mC "$H_col" \
#  -i "${loc}/methyl/#NEW/SIG/5mC/VANC.bed" VANC 5mC "$A_col" \
#  -i "${loc}/methyl/#NEW/SIG/4mC/VANC.bed" VANC 4mC "$A_col" \
#  -i "${loc}/methyl/#NEW/SIG/5mC/AMPI.bed" AMPI 5mC "$V_col" \
#  -i "${loc}/methyl/#NEW/SIG/4mC/AMPI.bed" AMPI 4mC "$V_col"

# The saved correction notes that the already-filtered 6mA files were moved into
# the new SIG tree for Figure 4. These are the same >20% files from the original DMR.
#cp "${loc}/methyl/1DMbase/20_pct/H202_6mA.bed" "${loc}/methyl/#NEW/SIG/6mA/H202.bed"
#cp "${loc}/methyl/1DMbase/20_pct/VANC_6mA.bed" "${loc}/methyl/#NEW/SIG/6mA/VANC.bed"
#cp "${loc}/methyl/1DMbase/20_pct/AMPI_6mA.bed" "${loc}/methyl/#NEW/SIG/6mA/AMPI.bed"

fig4_fun(){
    local mod="$1"
    python "$fig4" \
      -i "${loc}/methyl/#NEW/SIG/${mod}/H202.bed" H202 "#1b9e77" \
      -i "${loc}/methyl/#NEW/SIG/${mod}/VANC.bed" VANC "#d95f02" \
      -i "${loc}/methyl/#NEW/SIG/${mod}/AMPI.bed" AMPI "#7570b3" \
      -o "${loc}/methyl/#NEW/fig.4/${mod}" \
      --fai "${GENOME}.fai" \
      -b 50
}

#fig4_fun 4mC
#fig4_fun 5mC
#fig4_fun 6mA

# Methylation summary statistics after the split.
# 6mA continues to use its original DMR outputs; 4mC and 5mC use #NEW.
#
# NOTE: argument ordering is retained from the saved later correction script.
#python "$methyl_stats" \
#  -i "${loc}/methyl/#NEW/SIG/4mC/H202.bed" "${loc}/methyl/#NEW/4mC/H202/control_treatment_DIFF.bed" H202_4mC \
#  -i "${loc}/methyl/#NEW/SIG/4mC/VANC.bed" "${loc}/methyl/#NEW/4mC/VANC/control_treatment_DIFF.bed" VANC_4mC \
#  -i "${loc}/methyl/#NEW/SIG/4mC/AMPI.bed" "${loc}/methyl/#NEW/4mC/AMPI/control_treatment_DIFF.bed" AMPI_4mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/H202.bed" "${loc}/methyl/#NEW/5mC/H202/control_treatment_DIFF.bed" H202_5mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/AMPI.bed" "${loc}/methyl/#NEW/5mC/VANC/control_treatment_DIFF.bed" VANC_5mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/VANC.bed" "${loc}/methyl/#NEW/5mC/AMPI/control_treatment_DIFF.bed" AMPI_5mC \
#  -i "${loc}/methyl/1DMbase/20_pct/H202_6mA.bed" "${loc}/methyl/1DMbase/H202_6mA/H202_6mA_DIFF.bed" H202_6mA \
#  -i "${loc}/methyl/1DMbase/20_pct/VANC_6mA.bed" "${loc}/methyl/1DMbase/VANC_6mA/VANC_6mA_DIFF.bed" VANC_6mA \
#  -i "${loc}/methyl/1DMbase/20_pct/AMPI_6mA.bed" "${loc}/methyl/1DMbase/AMPI_6mA/AMPI_6mA_DIFF.bed" AMPI_6mA \
#  -o "${loc}/methyl/#NEW/STAT/METHY_STATS.tsv" \
#  --gene "$hybrid_bed" \
#  --upstream "${loc}/genomic-features/upstream.bed"

# Overlap of significant sites after split.
#python "$overlap_bases" \
#  -i H202 "${loc}/methyl/#NEW/SIG/5mC/H202.bed" \
#  -i VANC "${loc}/methyl/#NEW/SIG/5mC/VANC.bed" \
#  -i AMPI "${loc}/methyl/#NEW/SIG/5mC/AMPI.bed" \
#  -o "${loc}/methyl/#NEW/STAT/5mC_overlap_bases.bed"

#python "$overlap_bases" \
#  -i H202 "${loc}/methyl/#NEW/SIG/4mC/H202.bed" \
#  -i VANC "${loc}/methyl/#NEW/SIG/4mC/VANC.bed" \
#  -i AMPI "${loc}/methyl/#NEW/SIG/4mC/AMPI.bed" \
#  -o "${loc}/methyl/#NEW/STAT/4mC_overlap_bases.bed"


#############################################
# (2) GENOMIC FEATURES
#############################################

# Count contigs if required.
#contigs=$(grep -c '^>' "$GENOME")
#echo "${contigs} contigs found in the genome."

# PGAP.
#python "$pgap_runner" -D podman \
#  -g "$GENOME" \
#  -s "Enterococcus faecium" \
#  -o "${loc}/genomic-features/PGAP" \
#  -c 32 -m 90g -n

# Prepare files for Operon-mapper.
#python "$operon_mapper" \
#  --protien "$PGAP" \
#  --genome "$GENOME" \
#  --gff "$PGAP" \
#  --o-gff "${loc}/genomic-features/operon-mapper/operon-mapper_annot.gff" \
#  --o-nucl "${loc}/genomic-features/operon-mapper/DNA_annot.fasta"

# CMSEARCH.
#cmsearch --cpu 32 --cut_ga \
#  --tblout "${loc}/genomic-features/cmsearch/rfam_hits.tbl" \
#  "$RFAM_CM" "$GENOME" \
#  > "${loc}/genomic-features/cmsearch/rfam_hits.out"

#python "$ncRNA_bed_script" \
#  "${loc}/genomic-features/cmsearch/rfam_hits.tbl" \
#  > "$ncRNA"

# Hybridise PGAP and cmsearch results.
#python "$gff_hybrid" --PGAP "$PGAP" --RNA "$ncRNA" -o "$gene_an"

o_dir="${loc}/genomic-features"

# Standardise Operon-mapper output.
#python "$operon_mapper_conv" \
#  -i "$operons" \
#  -o "${o_dir}/operon.bed" \
#  --gff "$operon_gff"

# Promoters / upstream regions.
#python "$upstream_find" \
#  -i "${o_dir}/operon.bed" \
#  -o "${o_dir}/upstream.bed" \
#  -f "$gene_an" \
#  -fai "$fai"


#############################################
# (3) DESEQ2 TRANSCRIPTOME ANALYSIS
#############################################

index_dir="${loc}/transcript/genome_index"

# Build HISAT2 genome index once.
#mkdir -p "${index_dir}"
#hisat2-build "${GENOME}" "${index_dir}/genome_index"

tran_fun(){
    local prefix="$1"
    local out_dir="${loc}/transcript"
    mkdir -p "${out_dir}/${prefix}"
    local loc_sample="${out_dir}/${prefix}"

    cutadapt -j 32 -q 35 -m 140 \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o "${loc_sample}/${prefix}_trimmed_R1.fastq" \
      -p "${loc_sample}/${prefix}_trimmed_R2.fastq" \
      "${trans}/${prefix}.R1.raw.fastq.gz" \
      "${trans}/${prefix}.R2.raw.fastq.gz"

    hisat2 -x "${index_dir}/genome_index" \
      -1 "${loc_sample}/${prefix}_trimmed_R1.fastq" \
      -2 "${loc_sample}/${prefix}_trimmed_R2.fastq" \
      -S "${loc_sample}/${prefix}_aligned_reads.sam" \
      -p 30 --phred33 --dta --no-mixed --no-discordant --score-min L,0,-0.05

    samtools view -Sb -q 40 "${loc_sample}/${prefix}_aligned_reads.sam" | \
      samtools sort -o "${loc_sample}/${prefix}_aligned_reads_sorted.bam"

    samtools view -c -F 260 "${loc_sample}/${prefix}_aligned_reads_sorted.bam"

    featureCounts \
      -a "${gene_an}" \
      -F GTF \
      -t .,gene,riboswitch,pseudogene,binding_site,sequence_feature \
      -o "${out_dir}/${prefix}_gene_counts.txt" \
      --ignoreDup -g ID -T 32 -p -C --primary \
      "${loc_sample}/${prefix}_aligned_reads_sorted.bam"
}

#for sample in H202_T1 H202_T2 H202_T3 H202_C1 H202_C2 H202_C3 \
#              VANC_T1 VANC_T2 VANC_T3 VANC_C1 VANC_C2 VANC_C3 \
#              AMPI_T1 AMPI_T2 AMPI_T3 AMPI_C1 AMPI_C2 AMPI_C3; do
#    tran_fun "$sample"
#done

out_dir="${loc}/transcript"
deseq2_outdir="${out_dir}/DEseq2"

run_deseq_fun(){
    local treat="$1"

    "$deseq2" \
      "${out_dir}/${treat}_C1_gene_counts.txt" \
      "${out_dir}/${treat}_C2_gene_counts.txt" \
      "${out_dir}/${treat}_C3_gene_counts.txt" \
      "${out_dir}/${treat}_T1_gene_counts.txt" \
      "${out_dir}/${treat}_T2_gene_counts.txt" \
      "${out_dir}/${treat}_T3_gene_counts.txt" \
      "control,control,control,treatment,treatment,treatment" \
      "${deseq2_outdir}/${treat}.csv"
}

#mkdir -p "${deseq2_outdir}"
#for treatment in VANC H202 AMPI; do
#    run_deseq_fun "$treatment"
#done

#python "$volcano" \
#  -i "${deseq2_outdir}/H202.csv" H2O2 \
#     "${deseq2_outdir}/VANC.csv" VANC \
#     "${deseq2_outdir}/AMPI.csv" AMPI \
#  -o "$deseq2_outdir" --y 120 --x 7


#############################################
# KEGG / GSEA
#############################################

organisom="ko"
kegg_o="${loc}/transcript/KEGG"

#mkdir -p "$kegg_o"

# Merge KEGG KO annotations and log2FC.
#python "$ko_log" --ko "$KO" --deseq2 "${deseq2_outdir}/AMPI.csv" -o "${kegg_o}/AMPI_KO.tsv"
#python "$ko_log" --ko "$KO" --deseq2 "${deseq2_outdir}/VANC.csv" -o "${kegg_o}/VANC_KO.tsv"
#python "$ko_log" --ko "$KO" --deseq2 "${deseq2_outdir}/H202.csv" -o "${kegg_o}/H202_KO.tsv"

# Cluster GSEA.
#Rscript "$cluster_kegg" -i "${kegg_o}/AMPI_KO.tsv" -o "${kegg_o}/AMPI_GSAE" -c "$organisom"
#Rscript "$cluster_kegg" -i "${kegg_o}/VANC_KO.tsv" -o "${kegg_o}/VANC_GSAE" -c "$organisom"
#Rscript "$cluster_kegg" -i "${kegg_o}/H202_KO.tsv" -o "${kegg_o}/H202_GSAE" -c "$organisom"

# Plot.
#python "$kegg_plot" \
#  --data "${kegg_o}/H202_GSAE_results_ordered.tsv" H2O2 \
#  --data "${kegg_o}/VANC_GSAE_results_ordered.tsv" VANC \
#  --data "${kegg_o}/AMPI_GSAE_results_ordered.tsv" AMPI \
#  -o "${kegg_o}/KEGG"

#python "$kegg_distro" \
#  -i "$KO" \
#  --de "${deseq2_outdir}/AMPI.csv" AMPI \
#  --de "${deseq2_outdir}/VANC.csv" VANC \
#  --de "${deseq2_outdir}/H202.csv" H2O2 \
#  -o "${kegg_o}/KEGG_distrobution.png"


#############################################
# (4) MTASE SEARCH
#############################################

mtase_loc="${loc}/MTase/MTase-search"
gene_prot="${loc}/genomic-features/PGAP/annot.faa"

# PGAP candidate check.
#python "$pgap_check" \
#  -i "${gene_prot}" \
#  -o "${mtase_loc}/PGAP_filter_PASS.tsv"

# REBASE blastp.
#mkdir -p "${mtase_loc}/REBASE_local"
#cd "${mtase_loc}/REBASE_local"
#makeblastdb -in "$REBASE" -dbtype prot -out REBASE_db
#blastp \
#  -query "${gene_prot}" \
#  -db REBASE_db \
#  -out REBASE_blastp.tsv \
#  -num_threads 32 \
#  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
#  -max_target_seqs 5 \
#  -evalue 1e-10
#python "$consen_prot" \
#  -i REBASE_blastp.tsv \
#  -o "${mtase_loc}/REBASE_hit.tsv"

# InterPro scan.
#mkdir -p "${mtase_loc}/InterPro"
#cd "${mtase_loc}/InterPro"
#python "$interpro"
#cd ../
#python "$inter_consen" \
#  -i intrpro_results.tsv \
#  -o "${mtase_loc}/Inter_out.tsv"
#mv intrpro_results.tsv "${mtase_loc}/InterPro/intrpro_results.tsv"

# Final consensus.
#python "$mt_con" \
#  -p "${mtase_loc}/PGAP_filter_PASS.tsv" \
#  -i "${mtase_loc}/Inter_out.tsv" \
#  -r "${mtase_loc}/REBASE_hit.tsv" \
#  -o "${mtase_loc}/1.MTase.CONSENSUS.tsv"


#############################################
# (5) PSEUDOGENE RE-ARRANGEMENTS / CLAIR3
#############################################

dock_clair3(){
    local name="$1"
    local MODEL_NAME="r1041_e82_400bps_sup_v500"

    # The historical Clair3 analysis used the aligned combined cytosine BAM.
    # It was not a 4mC-vs-5mC DMR calculation, so this input is retained.
    local BAM="${loc}/methyl/${name}/4mC_5mC_ref_align-sort.bam"
    local REF="${GENOME}"
    local OUT="${loc}/MTase/psuedogenes/${name}"

    mkdir -p "$OUT"

    if [ ! -f "${BAM}.bai" ]; then
        samtools index "$BAM" || {
            echo "Failed to index $BAM"
            return 1
        }
    fi

    docker run --rm -u "$(id -u):$(id -g)" \
      -v "$(dirname "$BAM")":"$(dirname "$BAM")" \
      -v "$(dirname "$REF")":"$(dirname "$REF")" \
      -v "$OUT":"$OUT" \
      hkubal/clair3:latest \
      /opt/bin/run_clair3.sh \
        --bam_fn="$BAM" \
        --ref_fn="$REF" \
        --threads=32 \
        --platform=ont \
        --model_path="/opt/models/${MODEL_NAME}" \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --min_mq=20 \
        --min_coverage=5 \
        --enable_long_indel \
        --output="$OUT"
}

#for sample in HT1 HT2 HT3 HC1 HC2 HC3 VT1 VT2 VT3 VC1 VC2 VC3 AT1 AT2 AT3 AC1 AC2 AC3; do
#    dock_clair3 "$sample" || echo "Sample $sample failed"
#done

#python "$sum_VCF" \
#  --in "${loc}/MTase/psuedogenes/" \
#  --fai "$fai" \
#  > "${loc}/MTase/psuedogenes/summary.tsv"

#python "$intersect_VCF" \
#  --in "$MTase_final" \
#  --bed "$hybrid_bed" \
#  --vcf-dir "${loc}/MTase/psuedogenes/" \
#  > "${loc}/MTase/psuedogenes/matrix.tsv"


#############################################
# (6) FINAL METHYLATION METADATA TABLE
#############################################

# Operon index used by the metadata table.
#python "$conver_OM" \
#  -g "$operon_gff" \
#  -i "$operons" \
#  -o "${loc}/genomic-features/operon-mapper/operon_index.tsv"

# Final table combines split 4mC/5mC DMRs with the original 6mA DMRs.
#python "$meta_py" \
#  -i "${loc}/methyl/#NEW/SIG/4mC/H202.bed" H2O2 4mC \
#  -i "${loc}/methyl/#NEW/SIG/4mC/VANC.bed" VANC 4mC \
#  -i "${loc}/methyl/#NEW/SIG/4mC/AMPI.bed" AMPI 4mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/H202.bed" H2O2 5mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/VANC.bed" VANC 5mC \
#  -i "${loc}/methyl/#NEW/SIG/5mC/AMPI.bed" AMPI 5mC \
#  -i "${loc}/methyl/1DMbase/20_pct/H202_6mA.bed" H2O2 6mA \
#  -i "${loc}/methyl/1DMbase/20_pct/VANC_6mA.bed" VANC 6mA \
#  -i "${loc}/methyl/1DMbase/20_pct/AMPI_6mA.bed" AMPI 6mA \
#  -r "$gene_an" gene \
#  -r "${loc}/genomic-features/upstream.bed" upstream \
#  -d "${deseq2_outdir}/AMPI.csv" AMPI \
#  -d "${deseq2_outdir}/VANC.csv" VANC \
#  -d "${deseq2_outdir}/H202.csv" H2O2 \
#  --gff "$PGAP" \
#  --op "${loc}/genomic-features/operon-mapper/operon_index.tsv" \
#  -o "${loc}/NEW_META_report.tsv"


#############################################
# (7A) METHYLOMIC ENTROPY: ORIGINAL 6mA
#############################################

entro_A="${loc}/methyl/1DMbase/2entropy"
promoters="${o_dir}/upstream.bed"
gene_an_bed="${loc}/genomic-features/hybrid_genes.bed"
genes="$gene_an_bed"

# Convert hybrid GFF to BED if required.
#python "$gff_bed" -i "$gene_an" -o "$gene_an_bed"

# Original 6mA entropy analysis retained at its original location.
groups=(HT HC VC VT AC AT)
sites=(promoter coding)

bam_triple_6mA(){
    local gp="$1"
    echo "${loc}/methyl/${gp}1/6mA_ref_align-sort.bam" \
         "${loc}/methyl/${gp}2/6mA_ref_align-sort.bam" \
         "${loc}/methyl/${gp}3/6mA_ref_align-sort.bam"
}

run_entropy_6mA(){
    local gp="$1"
    local site="$2"
    local region=""

    if [[ "$site" == "promoter" ]]; then
        region="$promoters"
    else
        region="$genes"
    fi

    read -r bam1 bam2 bam3 <<<"$(bam_triple_6mA "$gp")"

    for b in "$bam1" "$bam2" "$bam3"; do
        [[ -s "$b" ]] || {
            echo "[warn] Missing BAM: $b — skipping ${gp}/6mA/${site}" >&2
            return
        }
        [[ -s "${b}.bai" ]] || {
            echo "[warn] Missing index: ${b}.bai — run 'samtools index' first; skipping ${gp}/6mA/${site}" >&2
            return
        }
    done

    modkit entropy \
      --in-bam "$bam1" \
      --in-bam "$bam2" \
      --in-bam "$bam3" \
      --ref "$GENOME" \
      --regions "$region" \
      --out-bed "${entro_A}/${site}/6mA_${gp}.bed" \
      --base A \
      --threads 32 \
      --log-filepath "${entro_A}/${site}/6mA_${gp}.log" \
      --force
}

#mkdir -p "${entro_A}/promoter" "${entro_A}/coding" "${entro_A}/stats"
#for gp in "${groups[@]}"; do
#  for site in "${sites[@]}"; do
#    run_entropy_6mA "$gp" "$site"
#  done
#done

calc_fun_6mA(){
    local region="$1"

    python "$calc_2" -p "6mA_${region}" --gff "$gene_an" \
      -i "${entro_A}/${region}/6mA_HT.bed/regions.bed" H2O2 \
      -i "${entro_A}/${region}/6mA_AT.bed/regions.bed" AMPI \
      -i "${entro_A}/${region}/6mA_VT.bed/regions.bed" VANC \
      -c "${entro_A}/${region}/6mA_HC.bed/regions.bed" H2O2 \
      -c "${entro_A}/${region}/6mA_AC.bed/regions.bed" AMPI \
      -c "${entro_A}/${region}/6mA_VC.bed/regions.bed" VANC \
      -o "${entro_A}/stats" \
      -p "6mA_${region}" \
      -d "${deseq2_outdir}/H202.csv" H2O2 \
      -d "${deseq2_outdir}/AMPI.csv" AMPI \
      -d "${deseq2_outdir}/VANC.csv" VANC \
      --operon "${loc}/genomic-features/operon-mapper/operon_index.tsv"

    python "$calc_3" -p "6mA_${region}" --gff "$gene_an" \
      -i "${entro_A}/${region}/6mA_HT.bed/regions.bed" H2O2 \
      -i "${entro_A}/${region}/6mA_AT.bed/regions.bed" AMPI \
      -i "${entro_A}/${region}/6mA_VT.bed/regions.bed" VANC \
      -c "${entro_A}/${region}/6mA_HC.bed/regions.bed" H2O2 \
      -c "${entro_A}/${region}/6mA_AC.bed/regions.bed" AMPI \
      -c "${entro_A}/${region}/6mA_VC.bed/regions.bed" VANC \
      -o "${entro_A}/stats" \
      -p "6mA_${region}" \
      -d "${deseq2_outdir}/H202.csv" H2O2 \
      -d "${deseq2_outdir}/AMPI.csv" AMPI \
      -d "${deseq2_outdir}/VANC.csv" VANC \
      --operon "${loc}/genomic-features/operon-mapper/operon_index.tsv" \
      >> "${entro_A}/stats/DELTA-ENTRO.stat.tsv"
}

#> "${entro_A}/stats/DELTA-ENTRO.stat.tsv"
#for region in "${sites[@]}"; do
#    calc_fun_6mA "$region"
#done

ent_stat_6mA(){
    local region="$1"
    local beds=""

    if [[ "$region" == "promoter" ]]; then
        beds="$promoters"
    elif [[ "$region" == "coding" ]]; then
        beds="$hybrid_bed"
    else
        echo "ERROR: region '$region' does not match expected values {promoter|gene} in ent_stat_6mA()" >&2
        return 1
    fi

    python "$base_en_stat" \
      -i "${entro_A}/${region}/6mA_HC.bed/regions.bed" \
      -i "${entro_A}/${region}/6mA_AC.bed/regions.bed" \
      -i "${entro_A}/${region}/6mA_VC.bed/regions.bed" \
      --bed "$beds" \
      --name "control_6mA_${region}" \
      >> "${entro_A}/entropy_stats.tsv"

    # Retained from the original analysis command record.
    python "$base_en_stat" \
      -i "${entro_A}/${region}/6mA_HC.bed/regions.bed" \
      -i "${entro_A}/${region}/6mA_AC.bed/regions.bed" \
      -i "${entro_A}/${region}/6mA_VC.bed/regions.bed" \
      --bed "$beds" \
      --name "treatment_6mA_${region}" \
      >> "${entro_A}/entropy_stats.tsv"
}

#: > "${entro_A}/entropy_stats.tsv"
#for region in "${sites[@]}"; do
#    ent_stat_6mA "$region"
#done


#############################################
# (7B) METHYLOMIC ENTROPY: LATER 4mC / 5mC SPLIT
#############################################

entro_C="${loc}/methyl/#NEW/entropy"

mkdir -p \
  "${loc}/methyl/#NEW/pileup/4mC" \
  "${loc}/methyl/#NEW/pileup/5mC"

# The saved correction separated the combined cytosine BAM specifically for entropy.
pileup_fun_C(){
    local sample="$1"

    # 4mC: ignore 5mC.
    modkit adjust-mods \
      "${loc}/methyl/${sample}/4mC_5mC_align/${sample}.bam" \
      "${loc}/methyl/#NEW/pileup/4mC/${sample}.bam" \
      --ignore m

    samtools index -@ 32 \
      "${loc}/methyl/#NEW/pileup/4mC/${sample}.bam"

    # 5mC: ignore 4mC.
    modkit adjust-mods \
      "${loc}/methyl/${sample}/4mC_5mC_align/${sample}.bam" \
      "${loc}/methyl/#NEW/pileup/5mC/${sample}.bam" \
      --ignore 21839

    samtools index -@ 32 \
      "${loc}/methyl/#NEW/pileup/5mC/${sample}.bam"
}

#for sample in HT1 HT2 HT3 HC1 HC2 HC3 VT1 VT2 VT3 VC1 VC2 VC3 AT1 AT2 AT3 AC1 AC2 AC3; do
#    pileup_fun_C "$sample"
#done

mods_C=(4mC 5mC)

bam_triple_C(){
    local gp="$1"
    local mod="$2"

    if [[ "$mod" == "4mC" ]]; then
        echo "${loc}/methyl/#NEW/pileup/4mC/${gp}1.bam" \
             "${loc}/methyl/#NEW/pileup/4mC/${gp}2.bam" \
             "${loc}/methyl/#NEW/pileup/4mC/${gp}3.bam"
    elif [[ "$mod" == "5mC" ]]; then
        echo "${loc}/methyl/#NEW/pileup/5mC/${gp}1.bam" \
             "${loc}/methyl/#NEW/pileup/5mC/${gp}2.bam" \
             "${loc}/methyl/#NEW/pileup/5mC/${gp}3.bam"
    fi
}

run_entropy_C(){
    local gp="$1"
    local mod="$2"
    local site="$3"
    local region=""

    if [[ "$site" == "promoter" ]]; then
        region="$promoters"
    else
        region="$genes"
    fi

    read -r bam1 bam2 bam3 <<<"$(bam_triple_C "$gp" "$mod")"

    for b in "$bam1" "$bam2" "$bam3"; do
        [[ -s "$b" ]] || {
            echo "[warn] Missing BAM: $b — skipping ${gp}/${mod}/${site}" >&2
            return
        }
        [[ -s "${b}.bai" ]] || {
            echo "[warn] Missing index: ${b}.bai — run 'samtools index' first; skipping ${gp}/${mod}/${site}" >&2
            return
        }
    done

    modkit entropy \
      --in-bam "$bam1" \
      --in-bam "$bam2" \
      --in-bam "$bam3" \
      --ref "$GENOME" \
      --regions "$region" \
      --out-bed "${entro_C}/${site}/${mod}_${gp}.bed" \
      --base C \
      --threads 32 \
      --log-filepath "${entro_C}/${site}/${mod}_${gp}.log" \
      --force
}

#mkdir -p "${entro_C}/promoter" "${entro_C}/coding" "${entro_C}/stats"
#for gp in "${groups[@]}"; do
#  for mod in "${mods_C[@]}"; do
#    for site in "${sites[@]}"; do
#      run_entropy_C "$gp" "$mod" "$site"
#    done
#  done
#done

calc_fun_C(){
    local mod="$1"
    local region="$2"

    python "$calc_2" -p "${mod}_${region}" --gff "$gene_an" \
      -i "${entro_C}/${region}/${mod}_HT.bed/regions.bed" H2O2 \
      -i "${entro_C}/${region}/${mod}_AT.bed/regions.bed" AMPI \
      -i "${entro_C}/${region}/${mod}_VT.bed/regions.bed" VANC \
      -c "${entro_C}/${region}/${mod}_HC.bed/regions.bed" H2O2 \
      -c "${entro_C}/${region}/${mod}_AC.bed/regions.bed" AMPI \
      -c "${entro_C}/${region}/${mod}_VC.bed/regions.bed" VANC \
      -o "${entro_C}/stats" \
      -p "${mod}_${region}" \
      -d "${deseq2_outdir}/H202.csv" H2O2 \
      -d "${deseq2_outdir}/AMPI.csv" AMPI \
      -d "${deseq2_outdir}/VANC.csv" VANC \
      --operon "${loc}/genomic-features/operon-mapper/operon_index.tsv"

    echo "########## ${mod} ${region}" \
      >> "${entro_C}/stats/DELTA-ENTRO.stat.tsv"

    python "$calc_3" -p "${mod}_${region}" --gff "$gene_an" \
      -i "${entro_C}/${region}/${mod}_HT.bed/regions.bed" H2O2 \
      -i "${entro_C}/${region}/${mod}_AT.bed/regions.bed" AMPI \
      -i "${entro_C}/${region}/${mod}_VT.bed/regions.bed" VANC \
      -c "${entro_C}/${region}/${mod}_HC.bed/regions.bed" H2O2 \
      -c "${entro_C}/${region}/${mod}_AC.bed/regions.bed" AMPI \
      -c "${entro_C}/${region}/${mod}_VC.bed/regions.bed" VANC \
      -o "${entro_C}/stats" \
      -p "${mod}_${region}" \
      -d "${deseq2_outdir}/H202.csv" H2O2 \
      -d "${deseq2_outdir}/AMPI.csv" AMPI \
      -d "${deseq2_outdir}/VANC.csv" VANC \
      --operon "${loc}/genomic-features/operon-mapper/operon_index.tsv" \
      >> "${entro_C}/stats/DELTA-ENTRO.stat.tsv"
}

#> "${entro_C}/stats/DELTA-ENTRO.stat.tsv"
#for mod in "${mods_C[@]}"; do
#  for region in "${sites[@]}"; do
#    calc_fun_C "$mod" "$region"
#  done
#done

ent_stat_C(){
    local mod="$1"
    local region="$2"
    local beds=""

    if [[ "$region" == "promoter" ]]; then
        beds="$promoters"
    elif [[ "$region" == "coding" ]]; then
        beds="$hybrid_bed"
    else
        echo "ERROR: region '$region' does not match expected values {promoter|gene} in ent_stat_C()" >&2
        return 1
    fi

    python "$base_en_stat" \
      -i "${entro_C}/${region}/${mod}_HC.bed/regions.bed" \
      -i "${entro_C}/${region}/${mod}_AC.bed/regions.bed" \
      -i "${entro_C}/${region}/${mod}_VC.bed/regions.bed" \
      --bed "$beds" \
      --name "control_${mod}_${region}" \
      >> "${entro_C}/entropy_stats.tsv"

    # Retained from the later correction command record.
    python "$base_en_stat" \
      -i "${entro_C}/${region}/${mod}_HC.bed/regions.bed" \
      -i "${entro_C}/${region}/${mod}_AC.bed/regions.bed" \
      -i "${entro_C}/${region}/${mod}_VC.bed/regions.bed" \
      --bed "$beds" \
      --name "treatment_${mod}_${region}" \
      >> "${entro_C}/entropy_stats.tsv"
}

#: > "${entro_C}/entropy_stats.tsv"
#for mod in "${mods_C[@]}"; do
#  for region in "${sites[@]}"; do
#    ent_stat_C "$mod" "$region"
#  done
#done

# Zero-modification summaries retained across the old 6mA and new cytosine entropy trees.
#python "$zero_mod" \
#  -i "${entro_C}/" \
#  -o "${entro_C}/entropy_zero_summary_C.tsv"

#python "$zero_mod" \
#  -i "${entro_A}" \
#  -o "${entro_C}/entropy_zero_summary_A.tsv"

# Never-methylated summary from the later corrected analysis tree.
#python "$count_unmethyl" \
#  -i "${loc}/methyl/#NEW" \
#  -o "${loc}/methyl/#NEW/STAT/NEVER_METHYLATED.TSV"


#############################################
# (7C) CONTROL-STATE ME / LATER METHYLATION RESPONSE
#############################################
#
# Later analysis tested whether methylomic-entropy structure in the pooled
# controls was associated with sites that subsequently showed the methylation
# responses represented in NEW_META_report.tsv.
#
# The saved Enterococcus analysis first generated pooled control ME tracks from
# HC1-HC3, VC1-VC3 and AC1-AC3 for 6mA, 5mC and 4mC, then ran
# me_node_enrichment_v4.py against the final methylation metadata and #NEW DMR
# tree. Commands below retain the saved Enterococcus invocation structure.

me_upgrade="${loc}/methyl/#NEW/UPGRADE"

# Generate the pooled control-state ME tracks:
#   ${me_upgrade}/4mC_ME/
#   ${me_upgrade}/5mC_ME/
#   ${me_upgrade}/6mA_ME/
#
# This helper retains the historical entropy/subsampling procedure used for
# this later analysis.
#bash "$me_control_entropy"

# Test whether control-state ME is associated with later methylation response.
#python "$me_node_enrichment" \
#  --metadata "${loc}/NEW_META_report.tsv" \
#  --fasta "$GENOME" \
#  --fai "$fai" \
#  --gff "$PGAP" \
#  --operons "$operon_list" \
#  --dmr-root "${loc}/methyl/#NEW/" \
#  --me 4mC x "${me_upgrade}/4mC_ME/" \
#  --me 5mC x "${me_upgrade}/5mC_ME/" \
#  --me 6mA x "${me_upgrade}/6mA_ME/" \
#  -o "${me_upgrade}/ME_enritchment"


#############################################
# (8) METHYLATION ENRICHMENT
#############################################

# Background regions.
#python "$background" \
#  -p "${o_dir}/operon.bed" \
#  -u "$promoters" \
#  -g "$gene_an_bed" \
#  -f "$fai" \
#  -o "${o_dir}/background.bed"

# 4mC.
#python "$hypergeom" \
#  --sig "${loc}/methyl/#NEW/SIG/4mC/H202.bed" H202 \
#  --full "${loc}/methyl/#NEW/4mC/H202/control_treatment_DIFF.bed" H202 \
#  --sig "${loc}/methyl/#NEW/SIG/4mC/VANC.bed" VANC \
#  --full "${loc}/methyl/#NEW/4mC/VANC/control_treatment_DIFF.bed" VANC \
#  --sig "${loc}/methyl/#NEW/SIG/4mC/AMPI.bed" AMPI \
#  --full "${loc}/methyl/#NEW/4mC/AMPI/control_treatment_DIFF.bed" AMPI \
#  --promoters "$promoters" \
#  --operons "${o_dir}/operon.bed" \
#  --genes "$hybrid_bed" \
#  -o "${loc}/methyl/#NEW/STAT" \
#  -p 4mC \
#  -b "${o_dir}/background.bed"

# 5mC.
#python "$hypergeom" \
#  --sig "${loc}/methyl/#NEW/SIG/5mC/H202.bed" H202 \
#  --full "${loc}/methyl/#NEW/5mC/H202/control_treatment_DIFF.bed" H202 \
#  --sig "${loc}/methyl/#NEW/SIG/5mC/VANC.bed" VANC \
#  --full "${loc}/methyl/#NEW/5mC/VANC/control_treatment_DIFF.bed" VANC \
#  --sig "${loc}/methyl/#NEW/SIG/5mC/AMPI.bed" AMPI \
#  --full "${loc}/methyl/#NEW/5mC/AMPI/control_treatment_DIFF.bed" AMPI \
#  --promoters "$promoters" \
#  --operons "${o_dir}/operon.bed" \
#  --genes "$hybrid_bed" \
#  -o "${loc}/methyl/#NEW/STAT" \
#  -p 5mC \
#  -b "${o_dir}/background.bed"

# 6mA (original output tree retained).
#python "$hypergeom" \
#  --sig "${loc}/methyl/1DMbase/20_pct/H202_6mA.bed" H202 \
#  --full "${loc}/methyl/1DMbase/H202_6mA/H202_6mA_DIFF.bed" H202 \
#  --sig "${loc}/methyl/1DMbase/20_pct/VANC_6mA.bed" VANC \
#  --full "${loc}/methyl/1DMbase/VANC_6mA/VANC_6mA_DIFF.bed" VANC \
#  --sig "${loc}/methyl/1DMbase/20_pct/AMPI_6mA.bed" AMPI \
#  --full "${loc}/methyl/1DMbase/AMPI_6mA/AMPI_6mA_DIFF.bed" AMPI \
#  --promoters "$promoters" \
#  --operons "${o_dir}/operon.bed" \
#  --genes "$hybrid_bed" \
#  -o "${loc}/methyl/2stats" \
#  -p 6mA \
#  -b "${o_dir}/background.bed"


#############################################
# (9) DATA SANITY / COVERAGE
#############################################

#python "$depth" \
#  -i "${loc}/methyl/HT1/coverage_6mA.txt" Treatment_H202 \
#  -i "${loc}/methyl/HT2/coverage_6mA.txt" Treatment_H202 \
#  -i "${loc}/methyl/HT3/coverage_6mA.txt" Treatment_H202 \
#  -i "${loc}/methyl/HC1/coverage_6mA.txt" Control_H202 \
#  -i "${loc}/methyl/HC2/coverage_6mA.txt" Control_H202 \
#  -i "${loc}/methyl/HC3/coverage_6mA.txt" Control_H202 \
#  -i "${loc}/methyl/VT1/coverage_6mA.txt" Treatment_VANC \
#  -i "${loc}/methyl/VT2/coverage_6mA.txt" Treatment_VANC \
#  -i "${loc}/methyl/VT3/coverage_6mA.txt" Treatment_VANC \
#  -i "${loc}/methyl/VC1/coverage_6mA.txt" Control_VANC \
#  -i "${loc}/methyl/VC2/coverage_6mA.txt" Control_VANC \
#  -i "${loc}/methyl/VC3/coverage_6mA.txt" Control_VANC \
#  -i "${loc}/methyl/AT1/coverage_6mA.txt" Treatment_AMPI \
#  -i "${loc}/methyl/AT2/coverage_6mA.txt" Treatment_AMPI \
#  -i "${loc}/methyl/AT3/coverage_6mA.txt" Treatment_AMPI \
#  -i "${loc}/methyl/AC1/coverage_6mA.txt" Control_AMPI \
#  -i "${loc}/methyl/AC2/coverage_6mA.txt" Control_AMPI \
#  -i "${loc}/methyl/AC3/coverage_6mA.txt" Control_AMPI \
#  --outdir "${loc}" \
#  --smooth 200 --one-figure


#############################################
# (10) MTASE CONSERVATION
#############################################

cons_mt="${loc}/MTase/MTase-conservation"
MTase_qurey="${cons_mt}/MTase_query.fasta"

# The PubMLST metadata used here had headers:
# id, isolate, year, age_yr
#
#python "$pub_extract" \
#  -i "$PUBMLST_CONTIGS" \
#  -t "$PubMETA" \
#  -o "${cons_mt}/database.fasta" \
#  --first-only --max 250

#python "$qurey_maker" \
#  -i "$MTase_final" \
#  -g "$dna_annot" \
#  -o "$MTase_qurey"

#makeblastdb \
#  -in "${cons_mt}/database.fasta" \
#  -dbtype nucl \
#  -out "${cons_mt}/PubMLST"

#blastn \
#  -query "$MTase_qurey" \
#  -db "${cons_mt}/PubMLST" \
#  -out "${cons_mt}/blast_PubMLST.out" \
#  -num_threads 32 \
#  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
#  -max_target_seqs 100000 \
#  -evalue 1e-10 \
#  -qcov_hsp_perc 70

#python "$bivariat" \
#  -meta "$PubMETA" \
#  -in "${cons_mt}/blast_PubMLST.out" \
#  -o "${cons_mt}/bivarient_density" \
#  --index "$MTase_qurey"


#############################################
# (11) MOTIF SEARCH
#############################################
#
# This section is retained as recorded in the original main script. The supplied
# later 4mC/5mC correction does not contain a replacement motif-search rerun, so
# no new split-motif procedure is invented here.

mot_dir="${loc}/methyl/3motif"
motif_mods=(4mC_5mC 6mA)
motif_groups=(HT HC AC AT VC VT)

motif_fun(){
    local name="$1"
    local mod="$2"
    local f="${loc}/methyl/${name}/${mod}.bed"

    modkit motif search \
      --in-bedmethyl "$f" \
      --ref "$GENOME" \
      --threads 32 \
      --out-table "${mot_dir}/${name}_${mod}_table.tsv" \
      --log-filepath "${mot_dir}/${name}_${mod}.log" \
      --min-sites 20
}

#mkdir -p "${mot_dir}"
#for g in "${motif_groups[@]}"; do
#  for i in 1 2 3; do
#    name="${g}${i}"
#    for mod in "${motif_mods[@]}"; do
#      motif_fun "$name" "$mod"
#    done
#  done
#done


#############################################
# (12) Tn INSERTION SEARCH
#############################################

fastq="$control_fastq"
tn_search_name="all_control"
BIN="${loc}/Tn-search/BIN"

#mkdir -p "$BIN"
#samtools faidx "$tn_contig" || true
#samtools faidx "$tn_free" || true

#minimap2 -t 32 -ax map-ont --secondary=yes \
#  "$tn_contig" \
#  "$fastq" \
#  | samtools sort -@ 32 -o "${BIN}/${tn_search_name}.tn.bam"

#samtools index -@ 32 "${BIN}/${tn_search_name}.tn.bam"

#python "$tn_extr" \
#  --bam "${BIN}/${tn_search_name}.tn.bam" \
#  --tn-contig "$tn_name" \
#  --out-fasta "${BIN}/${tn_search_name}.overhangs.fa" \
#  --out-meta "${BIN}/${tn_search_name}.overhangs.meta.tsv" \
#  --min-mapq-tn 20 \
#  --min-overhang 800 \
#  --tn-end-window 800

#minimap2 -t 32 -ax map-ont --secondary=no \
#  "$tn_free" \
#  "${BIN}/${tn_search_name}.overhangs.fa" \
#  | samtools sort -@ 32 \
#      -o "${BIN}/${tn_search_name}.overhangs2genome.bam"

#samtools index -@ 32 \
#  "${BIN}/${tn_search_name}.overhangs2genome.bam"

#minimap2 -t 32 -ax map-ont --secondary=yes \
#  "$tn_free" "$fastq" \
#  -o "${BIN}/${tn_search_name}_ALL-READS-ON-GENOME.bam"

#samtools sort -@ 32 \
#  "${BIN}/${tn_search_name}_ALL-READS-ON-GENOME.bam" \
#  -o "${BIN}/${tn_search_name}_ALL-READS-ON-GENOME.sort.bam"

#samtools index -@ 32 \
#  "${BIN}/${tn_search_name}_ALL-READS-ON-GENOME.sort.bam"

#python "$tn_count" \
#  --bam_over "${BIN}/${tn_search_name}.overhangs2genome.bam" \
#  --bam_all "${BIN}/${tn_search_name}_ALL-READS-ON-GENOME.sort.bam" \
#  --meta "${BIN}/${tn_search_name}.overhangs.meta.tsv" \
#  --out "${loc}/Tn-search/${tn_search_name}.tn_sites.tsv" \
#  --min-mapq-genome 20 \
#  --min-anchor-bp 60 \
#  --cluster 1 \
#  --min-support 20 \
#  --fai "${tn_free}.fai" \
#  --tn-contig "Tn6085" \
#  --tn-min-mapq 20 \
#  --tn-min-align-bp 100 \
#  --tn-summary-out "${loc}/Tn-search/${tn_search_name}.tn_global_summary.tsv" \
#  --tn-write-readlists-prefix "${BIN}/${tn_search_name}"


#############################################
# (13) PHYLOGENETIC TREE
#############################################

MTase_core="${loc}/MTase/MTase-search/tmp_MTase.tsv"
MTase_ref="${loc}/MTase/MTase-search/homology/REFERENCES/combined_output.fasta"
tree_out="${loc}/MTase/MTase-search/"
annot_2="${loc}/MTase/MTase-search/homology/REFERENCES/annot.faa"

# The annot.faa used here included genes not present in the original file
# where required (for example pseudogene candidates).
#
#python "$tree" \
#  -i "$MTase_core" \
#  -fa "$annot_2" \
#  -re "$MTase_ref" \
#  -o1 "${tree_out}TREE_1.png" \
#  -o2 "${tree_out}TREE_2.png" \
#  --hm1 "${tree_out}HEAT_1.png" \
#  --hm2 "${tree_out}HEAT_2.png"

