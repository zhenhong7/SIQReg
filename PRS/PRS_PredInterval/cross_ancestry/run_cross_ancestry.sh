#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config_cross_ancestry.sh"

ANCESTRIES=("asn" "afr" "white_euro")
SCALES=("original" "bc")

mkdir -p "${OUTBASE}/logs"

for anc in "${ANCESTRIES[@]}"; do
    for scale in "${SCALES[@]}"; do
        mkdir -p "${OUTBASE}/${anc}/${scale}_scale"
    done
done

wait_for_jobs() {
    local job_ids=("$@")
    for jid in "${job_ids[@]}"; do
        while squeue -j "${jid}" -h 2>/dev/null | grep -q "${jid}"; do
            sleep 30
        done
        STATE=$(sacct -j "${jid}" --format=State --noheader -P | head -1 | tr -d ' ')
        if [[ "${STATE}" != "COMPLETED" ]]; then
            echo "Job ${jid} failed: ${STATE}" >&2
            exit 1
        fi
    done
}

for anc in "${ANCESTRIES[@]}"; do
    for scale in "${SCALES[@]}"; do
        echo "${anc} / ${scale}"

        JOB_ID=$(sbatch --parsable --array=1-25 "${SCRIPT_DIR}/cross_setup.sh" "${anc}" "${scale}")
        wait_for_jobs "${JOB_ID}"

        JOB_ID=$(sbatch --parsable --array=1-25 "${SCRIPT_DIR}/cross_score_pgs.sh" "${anc}" "${scale}")
        wait_for_jobs "${JOB_ID}"

        JOB_ID=$(sbatch --parsable --array=1-25 "${SCRIPT_DIR}/cross_predinterval.sh" "${anc}" "${scale}")
        wait_for_jobs "${JOB_ID}"
    done
done

module load gcc/12.1.0
module load R/4.4.1
Rscript "${SCRIPT_DIR}/cross_results.R" "${OUTBASE}"
