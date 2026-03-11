#!/usr/bin/env bash
# =============================================================================
# Step 9 — netMHCpan Binding Prediction
# Predicts H-2 K^b and H-2 D^b binding for 8/9/10-mer peptides.
# Alleles: H-2-Kb, H-2-Db (C57BL/6 H-2b haplotype)
# Threshold: Rank_EL < 0.5% = strong binder (SB)
#            Rank_EL < 2.0% = weak binder   (WB)
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

PEP_DIR="${RESULTS_DIR}/peptides"
OUTDIR="${RESULTS_DIR}/netmhcpan"
S3_PEP="${S3_RESULTS}/peptides"
S3_OUT="${S3_RESULTS}/netmhcpan"
mkdir -p "${OUTDIR}"

# Locate netMHCpan binary
if command -v netMHCpan &>/dev/null; then
    NMHC="netMHCpan"
elif [[ -f "${NETMHCPAN}" ]]; then
    NMHC="${NETMHCPAN}"
else
    log "[ERROR] netMHCpan not found."
    log "  Install from: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    log "  Or via conda: mamba install -c bioconda netmhcpan"
    exit 1
fi
log "netMHCpan: $(${NMHC} -h 2>&1 | grep 'NetMHCpan' | head -1 || echo 'found')"

for PLEN in 8 9 10; do
    FASTA="${PEP_DIR}/peptides_${PLEN}mer.fasta"

    if [[ ! -f "${FASTA}" ]]; then
        aws s3 cp "${S3_PEP}/peptides_${PLEN}mer.fasta" "${FASTA}" 2>/dev/null || {
            log "[WARN] No ${PLEN}-mer FASTA — skipping"
            continue
        }
    fi

    N_PEPS=$(grep -c "^>" "${FASTA}" || echo 0)
    log "========================================"
    log "netMHCpan: ${PLEN}-mers (${N_PEPS} peptides)"
    log "========================================"

    for ALLELE in H-2-Kb H-2-Db; do
        ALLELE_SAFE="${ALLELE//-/_}"
        OUT_TXT="${OUTDIR}/${PLEN}mer_${ALLELE_SAFE}.txt"
        OUT_TSV="${OUTDIR}/${PLEN}mer_${ALLELE_SAFE}.tsv"

        if aws s3 ls "${S3_OUT}/${PLEN}mer_${ALLELE_SAFE}.tsv" &>/dev/null; then
            skip "${PLEN}mer ${ALLELE} predictions (S3)"
            continue
        fi

        log "  Running: ${ALLELE} ${PLEN}-mers..."
        ${NMHC} \
            -a "${ALLELE}" \
            -l "${PLEN}" \
            -f "${FASTA}" \
            -BA \
            -xls \
            -xlsfile "${OUT_TXT}" \
            -inptype 0 \
            > "${OUT_TXT}.raw" 2>&1

        # Parse output to TSV
        if [[ -f "${OUT_TXT}" ]]; then
            # netMHCpan XLS output: tab-separated, parse header + data rows
            awk 'NR>1 && /^[0-9]/' "${OUT_TXT}" \
            | awk -v allele="${ALLELE}" -v plen="${PLEN}" \
                'BEGIN{OFS="\t"; print "allele","plen","pos","peptide","core","el_rank","ba_rank","binder_class"}
                {
                    pep=$3; el=$12; ba=$14
                    cls = (el+0 < 0.5) ? "SB" : (el+0 < 2.0) ? "WB" : "NB"
                    if (cls != "NB") print allele,plen,$2,pep,$5,el,ba,cls
                }' > "${OUT_TSV}"
            SB=$(awk '$8=="SB"' "${OUT_TSV}" | wc -l)
            WB=$(awk '$8=="WB"' "${OUT_TSV}" | wc -l)
            log "    ${ALLELE}: SB=${SB}  WB=${WB}"
        else
            # Fallback: parse raw text output
            grep -v "^#\|^-\|^$\|^Pos\|^Protein\|^Number" "${OUT_TXT}.raw" \
            | awk -v allele="${ALLELE}" -v plen="${PLEN}" \
                'NF>=12 {
                    pep=$3; el=$12
                    cls = (el+0 < 0.5) ? "SB" : (el+0 < 2.0) ? "WB" : "NB"
                    if (cls != "NB") print allele"\t"plen"\t"$1"\t"pep"\t"$6"\t"el"\t.\t"cls
                }' > "${OUT_TSV}"
            SB=$(awk '$8=="SB"' "${OUT_TSV}" | wc -l)
            WB=$(awk '$8=="WB"' "${OUT_TSV}" | wc -l)
            log "    ${ALLELE}: SB=${SB}  WB=${WB} (raw parse)"
        fi

        s3up "${OUT_TSV}"     "netmhcpan/$(basename ${OUT_TSV})"
        rm -f "${OUT_TXT}.raw"
    done
done

# Merge all predictions into one table
MERGED="${OUTDIR}/all_predictions.tsv"
echo -e "allele\tplen\tpos\tpeptide\tcore\tel_rank\tba_rank\tbinder_class" > "${MERGED}"
cat "${OUTDIR}"/*mer_*.tsv | grep -v "^allele" >> "${MERGED}" || true
s3up "${MERGED}" "netmhcpan/all_predictions.tsv"

TOTAL_SB=$(awk -F'\t' '$8=="SB"' "${MERGED}" | wc -l)
TOTAL_WB=$(awk -F'\t' '$8=="WB"' "${MERGED}" | wc -l)
log "All predictions: Strong Binders=${TOTAL_SB}  Weak Binders=${TOTAL_WB}"
log "Merged: ${MERGED}"
log "Step 9 complete."
