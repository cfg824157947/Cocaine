#! bash
set -euo pipefail

fragment_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/Fixed_Recall_peaks/sample_cluster_fragment/"
cut_dir="/project/zduren/durenlab/palmetto/cham/cocain/Final/Fragment/cut_dir/"

mkdir -p "$cut_dir"
# Function: fragments -> cut sites (+4 / -5); keeps barcode/count if present
mk_cuts() {
    awk 'BEGIN{OFS="\t"}
  {
    chr=$1; s=$2; e=$3;
    bc=(NF>=4 ? $4 : ".");
    cnt=(NF>=5 ? $5 : 1);
    # left cut: start+4 (1 bp)
    print chr, s+4, s+5, bc, cnt;
    # right cut: end-5 (1 bp)
    print chr, e-5, e-4, bc, cnt;
  }' "$1"
}

# Loop over all .bed files
for f in "$fragment_dir"/*.bed; do
    base=$(basename "$f" .bed)
    echo "→ $base"
#    mk_cuts "$f" | bgzip -c > "$cut_dir/${base}.cuts.bed.gz"
    mk_cuts "$f" \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -S 70% --parallel="$(nproc)" \
        | bgzip -@ "$(nproc)" -c > "$cut_dir/${base}.cuts.bed.gz"
    tabix -p bed -f "$cut_dir/${base}.cuts.bed.gz"

done

echo "Done. Cuts written to: $cut_dir"
