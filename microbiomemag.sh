#!/usr/bin/env bash
###############################################################################
# Metagenomics pipeline: raw reads -> QC -> assembly -> binning -> MAGs ->
#                         taxonomy -> quality -> functional annotation
#
# Assumes paired-end Illumina short reads. Designed to be run sample-by-sample
# (set SAMPLE/R1/R2) or looped over a sample sheet — see bottom of script.
#
# Tool versions are pinned loosely via conda/mamba environments so the
# pipeline is reproducible. Install everything up front (see ENV SETUP),
# then run the numbered STEP functions in order.
###############################################################################

# Gaurav Sablok
# gsablok@proton.me

set -euo pipefail
IFS=$'\n\t'

########################################
# 0. USER CONFIG — edit per run
########################################
SAMPLE="sample1"               # sample name / prefix
R1="raw/${SAMPLE}_R1.fastq.gz" # raw forward reads
R2="raw/${SAMPLE}_R2.fastq.gz" # raw reverse reads
THREADS=16
MEMGB=128
OUTDIR="results/${SAMPLE}"
HOST_BOWTIE2_INDEX="ref/host_genome"    # e.g. human genome index; leave "" to skip host removal
GTDBTK_DATA_PATH="/data/gtdbtk_release" # set GTDBTK_DATA_PATH env var, downloaded separately
CHECKM2_DB="/data/checkm2_db/uniref100.KO.1.dmnd"
EGGNOG_DB_DIR="/data/eggnog_db"
ASSEMBLER="metaspades" # metaspades | megahit

mkdir -p "${OUTDIR}"/{00_qc,01_hostfree,02_assembly,03_mapping,04_binning,05_refine,06_checkm,07_gtdbtk,08_annotation,logs}

########################################
# ENV SETUP (run once, not part of the pipeline loop)
########################################
: <<'ENV_SETUP'
mamba create -n mg_qc      -c bioconda -c conda-forge fastp bowtie2 samtools -y
mamba create -n mg_asm     -c bioconda -c conda-forge spades megahit -y
mamba create -n mg_map     -c bioconda -c conda-forge bwa samtools -y
mamba create -n mg_bin     -c bioconda -c conda-forge metabat2 maxbin2 concoct das_tool -y
mamba create -n mg_qual    -c bioconda -c conda-forge checkm2 -y
mamba create -n mg_tax     -c bioconda -c conda-forge gtdbtk -y
mamba create -n mg_annot   -c bioconda -c conda-forge prokka eggnog-mapper -y
mamba create -n mg_derep   -c bioconda -c conda-forge drep coverm -y
ENV_SETUP

########################################
# STEP 1 — Quality control & trimming
########################################
step1_qc() {
	echo "[STEP1] QC/trimming ${SAMPLE}"
	conda run -n mg_qc fastp \
		-i "${R1}" -I "${R2}" \
		-o "${OUTDIR}/00_qc/${SAMPLE}_R1.trim.fastq.gz" \
		-O "${OUTDIR}/00_qc/${SAMPLE}_R2.trim.fastq.gz" \
		--detect_adapter_for_pe \
		--thread "${THREADS}" \
		-j "${OUTDIR}/00_qc/${SAMPLE}.fastp.json" \
		-h "${OUTDIR}/00_qc/${SAMPLE}.fastp.html" \
		2>"${OUTDIR}/logs/${SAMPLE}.fastp.log"
}

########################################
# STEP 2 — Host / contaminant read removal (optional)
########################################
step2_host_removal() {
	if [[ -z "${HOST_BOWTIE2_INDEX}" ]]; then
		echo "[STEP2] Skipping host removal (no index configured)"
		cp "${OUTDIR}/00_qc/${SAMPLE}_R1.trim.fastq.gz" "${OUTDIR}/01_hostfree/${SAMPLE}_R1.clean.fastq.gz"
		cp "${OUTDIR}/00_qc/${SAMPLE}_R2.trim.fastq.gz" "${OUTDIR}/01_hostfree/${SAMPLE}_R2.clean.fastq.gz"
		return
	fi
	echo "[STEP2] Removing host reads for ${SAMPLE}"
	conda run -n mg_qc bowtie2 -x "${HOST_BOWTIE2_INDEX}" \
		-1 "${OUTDIR}/00_qc/${SAMPLE}_R1.trim.fastq.gz" \
		-2 "${OUTDIR}/00_qc/${SAMPLE}_R2.trim.fastq.gz" \
		-p "${THREADS}" --very-sensitive -S /dev/null \
		--un-conc-gz "${OUTDIR}/01_hostfree/${SAMPLE}_R%.clean.fastq.gz" \
		2>"${OUTDIR}/logs/${SAMPLE}.bowtie2_host.log"
	mv "${OUTDIR}/01_hostfree/${SAMPLE}_R1.clean.fastq.gz" "${OUTDIR}/01_hostfree/${SAMPLE}_R1.clean.fastq.gz" 2>/dev/null || true
}

########################################
# STEP 3 — Assembly (metaSPAdes or MEGAHIT)
########################################
step3_assembly() {
	echo "[STEP3] Assembling ${SAMPLE} with ${ASSEMBLER}"
	local R1C="${OUTDIR}/01_hostfree/${SAMPLE}_R1.clean.fastq.gz"
	local R2C="${OUTDIR}/01_hostfree/${SAMPLE}_R2.clean.fastq.gz"
	if [[ "${ASSEMBLER}" == "metaspades" ]]; then
		conda run -n mg_asm metaspades.py \
			-1 "${R1C}" -2 "${R2C}" \
			-o "${OUTDIR}/02_assembly" \
			-t "${THREADS}" -m "${MEMGB}" \
			2>"${OUTDIR}/logs/${SAMPLE}.metaspades.log"
		ASSEMBLY_FASTA="${OUTDIR}/02_assembly/contigs.fasta"
	else
		rm -rf "${OUTDIR}/02_assembly/megahit_out"
		conda run -n mg_asm megahit \
			-1 "${R1C}" -2 "${R2C}" \
			-o "${OUTDIR}/02_assembly/megahit_out" \
			-t "${THREADS}" --memory 0.9 \
			2>"${OUTDIR}/logs/${SAMPLE}.megahit.log"
		ASSEMBLY_FASTA="${OUTDIR}/02_assembly/megahit_out/final.contigs.fa"
	fi
	# length-filter contigs (>=1500 bp is standard for binning)
	conda run -n mg_qc seqkit seq -m 1500 "${ASSEMBLY_FASTA}" \
		>"${OUTDIR}/02_assembly/${SAMPLE}.contigs.filtered.fasta"
}

########################################
# STEP 4 — Map reads back to contigs (for coverage/binning)
########################################
step4_mapping() {
	echo "[STEP4] Mapping reads to contigs for ${SAMPLE}"
	local ASM="${OUTDIR}/02_assembly/${SAMPLE}.contigs.filtered.fasta"
	conda run -n mg_map bwa index "${ASM}"
	conda run -n mg_map bwa mem -t "${THREADS}" "${ASM}" \
		"${OUTDIR}/01_hostfree/${SAMPLE}_R1.clean.fastq.gz" \
		"${OUTDIR}/01_hostfree/${SAMPLE}_R2.clean.fastq.gz" \
		2>"${OUTDIR}/logs/${SAMPLE}.bwa.log" |
		conda run -n mg_map samtools sort -@ "${THREADS}" -o "${OUTDIR}/03_mapping/${SAMPLE}.sorted.bam" -
	conda run -n mg_map samtools index "${OUTDIR}/03_mapping/${SAMPLE}.sorted.bam"
}

########################################
# STEP 5 — Binning: MetaBAT2 + MaxBin2 + CONCOCT, then DAS_Tool consensus
########################################
step5_binning() {
	echo "[STEP5] Binning contigs for ${SAMPLE}"
	local ASM="${OUTDIR}/02_assembly/${SAMPLE}.contigs.filtered.fasta"
	local BAM="${OUTDIR}/03_mapping/${SAMPLE}.sorted.bam"
	local BINDIR="${OUTDIR}/04_binning"

	# depth file shared by MetaBAT2 / MaxBin2
	conda run -n mg_bin jgi_summarize_bam_contig_depths \
		--outputDepth "${BINDIR}/depth.txt" "${BAM}"

	# MetaBAT2
	mkdir -p "${BINDIR}/metabat2"
	conda run -n mg_bin metabat2 -i "${ASM}" -a "${BINDIR}/depth.txt" \
		-o "${BINDIR}/metabat2/bin" -t "${THREADS}" -m 1500

	# MaxBin2 (wants a simplified abundance file: contig \t depth)
	awk 'NR>1{print $1"\t"$3}' "${BINDIR}/depth.txt" >"${BINDIR}/maxbin_abund.txt"
	mkdir -p "${BINDIR}/maxbin2"
	conda run -n mg_bin run_MaxBin.pl -contig "${ASM}" \
		-abund "${BINDIR}/maxbin_abund.txt" \
		-out "${BINDIR}/maxbin2/bin" -thread "${THREADS}"

	# CONCOCT
	mkdir -p "${BINDIR}/concoct"
	conda run -n mg_bin cut_up_fasta.py "${ASM}" -c 10000 -o 0 --merge_last \
		-b "${BINDIR}/concoct/contigs_10K.bed" >"${BINDIR}/concoct/contigs_10K.fasta"
	conda run -n mg_bin concoct_coverage_table.py \
		"${BINDIR}/concoct/contigs_10K.bed" "${BAM}" >"${BINDIR}/concoct/coverage_table.tsv"
	conda run -n mg_bin concoct \
		--composition_file "${BINDIR}/concoct/contigs_10K.fasta" \
		--coverage_file "${BINDIR}/concoct/coverage_table.tsv" \
		-b "${BINDIR}/concoct/" -t "${THREADS}"
	conda run -n mg_bin merge_cutup_clustering.py \
		"${BINDIR}/concoct/clustering_gt1000.csv" >"${BINDIR}/concoct/clustering_merged.csv"
	mkdir -p "${BINDIR}/concoct/fasta_bins"
	conda run -n mg_bin extract_fasta_bins.py "${ASM}" \
		"${BINDIR}/concoct/clustering_merged.csv" --output_path "${BINDIR}/concoct/fasta_bins"

	# DAS_Tool: reconcile the three bin sets into one consensus set
	conda run -n mg_bin Fasta_to_Scaffolds2Bin.sh -i "${BINDIR}/metabat2" -e fa >"${BINDIR}/metabat2.scaffolds2bin.tsv"
	conda run -n mg_bin Fasta_to_Scaffolds2Bin.sh -i "${BINDIR}/maxbin2" -e fasta >"${BINDIR}/maxbin2.scaffolds2bin.tsv"
	conda run -n mg_bin Fasta_to_Scaffolds2Bin.sh -i "${BINDIR}/concoct/fasta_bins" -e fa >"${BINDIR}/concoct.scaffolds2bin.tsv"

	mkdir -p "${OUTDIR}/05_refine"
	conda run -n mg_bin DAS_Tool \
		-i "${BINDIR}/metabat2.scaffolds2bin.tsv,${BINDIR}/maxbin2.scaffolds2bin.tsv,${BINDIR}/concoct.scaffolds2bin.tsv" \
		-l metabat2,maxbin2,concoct \
		-c "${ASM}" \
		-o "${OUTDIR}/05_refine/${SAMPLE}_DASTool" \
		--write_bins -t "${THREADS}" \
		2>"${OUTDIR}/logs/${SAMPLE}.dastool.log"
	# final MAGs land in: ${OUTDIR}/05_refine/${SAMPLE}_DASTool_DASTool_bins/*.fa
}

########################################
# STEP 6 — Bin quality: completeness / contamination (CheckM2)
########################################
step6_checkm2() {
	echo "[STEP6] CheckM2 quality assessment for ${SAMPLE}"
	conda run -n mg_qual checkm2 predict \
		--input "${OUTDIR}/05_refine/${SAMPLE}_DASTool_DASTool_bins" \
		--output-directory "${OUTDIR}/06_checkm" \
		-x fa --threads "${THREADS}" --database_path "${CHECKM2_DB}"
	# Filter to "good" MAGs: completeness >=50%, contamination <=10% (MIMAG medium quality)
	# Adjust thresholds for high-quality MAGs: completeness >=90%, contamination <=5%
	awk -F'\t' 'NR==1 || ($2>=50 && $3<=10)' \
		"${OUTDIR}/06_checkm/quality_report.tsv" >"${OUTDIR}/06_checkm/passed_bins.tsv"
}

########################################
# STEP 7 — Taxonomic classification (GTDB-Tk)
########################################
step7_gtdbtk() {
	echo "[STEP7] GTDB-Tk taxonomic classification for ${SAMPLE}"
	export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH}"
	conda run -n mg_tax gtdbtk classify_wf \
		--genome_dir "${OUTDIR}/05_refine/${SAMPLE}_DASTool_DASTool_bins" \
		--out_dir "${OUTDIR}/07_gtdbtk" \
		-x fa --cpus "${THREADS}" --mash_db "${GTDBTK_DATA_PATH}/mash_db"
}

########################################
# STEP 8 — Functional annotation of each MAG (Prokka + eggNOG-mapper)
########################################
step8_annotation() {
	echo "[STEP8] Annotating MAGs for ${SAMPLE}"
	local BINS="${OUTDIR}/05_refine/${SAMPLE}_DASTool_DASTool_bins"
	for bin_fa in "${BINS}"/*.fa; do
		local bin_name
		bin_name=$(basename "${bin_fa}" .fa)
		local ANNOT_OUT="${OUTDIR}/08_annotation/${bin_name}"
		mkdir -p "${ANNOT_OUT}"

		# Gene calling + basic annotation (fast, per-genome)
		conda run -n mg_annot prokka --outdir "${ANNOT_OUT}/prokka" \
			--prefix "${bin_name}" --cpus "${THREADS}" --metagenome \
			--force "${bin_fa}" \
			2>"${OUTDIR}/logs/${SAMPLE}.${bin_name}.prokka.log"

		# Deeper functional annotation: KEGG/COG/GO via eggNOG-mapper on predicted proteins
		conda run -n mg_annot emapper.py \
			-i "${ANNOT_OUT}/prokka/${bin_name}.faa" \
			--itype proteins \
			-o "${bin_name}" \
			--output_dir "${ANNOT_OUT}/eggnog" \
			--data_dir "${EGGNOG_DB_DIR}" \
			--cpu "${THREADS}" \
			2>"${OUTDIR}/logs/${SAMPLE}.${bin_name}.eggnog.log"
	done
}

########################################
# STEP 9 (optional) — Cross-sample dereplication + relative abundance
# Run only after processing ALL samples, pooling every sample's passed MAGs
# into one directory (e.g. all_mags/) first.
########################################
step9_derep_and_abundance() {
	echo "[STEP9] Dereplicating MAGs across samples and computing abundance"
	conda run -n mg_derep dRep dereplicate results/dRep_output \
		-g all_mags/*.fa -p "${THREADS}" -comp 50 -con 10
	conda run -n mg_derep coverm genome \
		-1 raw/*_R1.fastq.gz -2 raw/*_R2.fastq.gz \
		--genome-fasta-directory results/dRep_output/dereplicated_genomes \
		-x fa -t "${THREADS}" -m relative_abundance \
		>results/mag_relative_abundance.tsv
}

########################################
# RUN — comment out steps you don't need
########################################
step1_qc
step2_host_removal
step3_assembly
step4_mapping
step5_binning
step6_checkm2
step7_gtdbtk
step8_annotation
# step9_derep_and_abundance   # run once, after looping this script over all samples

echo "[DONE] Pipeline complete for ${SAMPLE}. MAGs: ${OUTDIR}/05_refine/${SAMPLE}_DASTool_DASTool_bins"
