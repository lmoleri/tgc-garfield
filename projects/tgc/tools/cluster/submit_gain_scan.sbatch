#!/bin/bash
# SLURM array template for the microscopic gain-vs-position scan.
#
# Each array task runs one job from a manifest produced by:
#     python3 tools/gain_scan.py emit --config config/scan_gain.json \
#             --out $SCRATCH/gain_scan [--xchunks K]
# which prints the array size N.  Submit with that N:
#     sbatch --array=0-$((N-1)) tools/cluster/submit_gain_scan.sbatch \
#            $SCRATCH/gain_scan/manifest.txt $SCRATCH/gain_scan
# Then merge + plot on a login node:
#     python3 tools/gain_scan.py merge --out $SCRATCH/gain_scan
#     python3 tools/gain_scan.py plot  --csv  $SCRATCH/gain_scan/gain_map.csv
#
# The gain scan has ion drift OFF, so it needs NO ion-mobility table and NO
# GARFIELD_INSTALL at runtime — only the Garfield/ROOT runtime libraries on
# LD_LIBRARY_PATH (set by your site's module/setup below).

#SBATCH --job-name=tgc_gain
#SBATCH --array=0-0            # OVERRIDE on the command line: --array=0-$((N-1))
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00        # per task; size for (x-points-in-job) x n_events x ~10 s
#SBATCH --output=tgc_gain_%A_%a.out
## #SBATCH --account=YOUR_ACCOUNT
## #SBATCH --partition=YOUR_PARTITION

set -euo pipefail

MANIFEST="${1:?usage: sbatch [...] submit_gain_scan.sbatch <manifest.txt> <out_dir>}"
OUT="${2:?missing <out_dir>}"

# ── Configure for your cluster ───────────────────────────────────────────────
# Path to projects/tgc on the cluster (where build/tgc_sim and the .gas table live):
REPO="${TGC_REPO:-$HOME/Garfield++/projects/tgc}"
# Make the Garfield/ROOT runtime libraries findable. Use whatever your site needs,
# e.g. one of:
#   source "$GARFIELD_HOME/share/Garfield/setupGarfield.sh"
#   module load root garfield
# (If you build ion-drift runs instead, also: export GARFIELD_INSTALL=/path/to/garfield)
# ─────────────────────────────────────────────────────────────────────────────

cd "$REPO"   # so tgc_sim finds the cached .gas table in the working directory

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
cfg=$(printf '%s' "$line" | awk '{print $1}')
name=$(printf '%s' "$line" | awk '{print $2}')
echo "[task $SLURM_ARRAY_TASK_ID] config=$cfg run-name=$name out=$OUT"

./build/tgc_sim --config "$cfg" --out "$OUT" --run-name "$name"
