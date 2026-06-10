#!/bin/bash
# Script to organize figures into subdirectories
# Run from repository root, e.g.: bash figures/organize_figures.sh
# Usage: bash organize_figures.sh

set -euo pipefail

FDIR="$(pwd)/figures"

echo "=== Creating directory structure ==="
mkdir -p "$FDIR/Main_Figures/Figure1_scRNA_UPR_landscape"
mkdir -p "$FDIR/Main_Figures/Figure2_CellChat_communication"
mkdir -p "$FDIR/Main_Figures/Figure3_consensus_clustering"
mkdir -p "$FDIR/Main_Figures/Figure4_immune_genomic_pathway"
mkdir -p "$FDIR/Main_Figures/Figure5_ML_model_validation"
mkdir -p "$FDIR/Main_Figures/Figure6_nomogram_prognosis"
mkdir -p "$FDIR/Supplementary_Figures"
mkdir -p "$FDIR/QC_Figures"
mkdir -p "$FDIR/Analysis_Figures"

echo "=== Copying Figure 1 files ==="
# Fig1* files
for f in "$FDIR"/Fig1*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure1_scRNA_UPR_landscape/"
done
# Part1 differential analysis related: tam_m1m2, cd8_exhaustion/cd8_exh, dc_antigen/dc_hlaii
for pattern in "*tam_m1m2*" "*cd8_exhaustion*" "*cd8_exh*" "*dc_antigen*" "*dc_hlaii*"; do
    for f in "$FDIR"/${pattern}.pdf; do
        [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure1_scRNA_UPR_landscape/"
    done
done

echo "=== Copying Figure 2 files ==="
# Fig2* (already includes Fig2_chord_*)
for f in "$FDIR"/Fig2*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure2_CellChat_communication/"
done

echo "=== Copying Figure 3 files ==="
for f in "$FDIR"/Fig3*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure3_consensus_clustering/"
done

echo "=== Copying Figure 4 files ==="
for f in "$FDIR"/Fig4*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure4_immune_genomic_pathway/"
done

echo "=== Copying Figure 5 files ==="
for f in "$FDIR"/Fig5*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure5_ML_model_validation/"
done

echo "=== Copying Figure 6 files ==="
for f in "$FDIR"/Fig6*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Main_Figures/Figure6_nomogram_prognosis/"
done

echo "=== Copying Supplementary files ==="
for f in "$FDIR"/FigS*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Supplementary_Figures/"
done
for f in "$FDIR"/SuppFig*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/Supplementary_Figures/"
done

echo "=== Copying QC files ==="
for pattern in "sc_qc*" "sc_elbow*" "sc_clustree*" "preprocessing*"; do
    for f in "$FDIR"/${pattern}.pdf; do
        [ -f "$f" ] && cp "$f" "$FDIR/QC_Figures/"
    done
done
# batch_correction (matches FigS1_batch_correction which is already in Supplementary;
# but standalone batch_correction files go to QC)
for f in "$FDIR"/batch_correction*.pdf; do
    [ -f "$f" ] && cp "$f" "$FDIR/QC_Figures/"
done

echo "=== Copying remaining files to Analysis_Figures ==="
# Collect all PDFs that were already classified
declare -A classified
for d in "$FDIR/Main_Figures"/Figure*/ "$FDIR/Supplementary_Figures/" "$FDIR/QC_Figures/"; do
    for f in "$d"*.pdf; do
        [ -f "$f" ] && classified["$(basename "$f")"]=1
    done
done

# Copy unclassified PDFs to Analysis_Figures
for f in "$FDIR"/*.pdf; do
    [ -f "$f" ] || continue
    bn="$(basename "$f")"
    if [[ -z "${classified[$bn]+x}" ]]; then
        cp "$f" "$FDIR/Analysis_Figures/"
    fi
done

echo ""
echo "=== Summary ==="
echo "Figure1: $(ls "$FDIR/Main_Figures/Figure1_scRNA_UPR_landscape/" 2>/dev/null | wc -l) files"
echo "Figure2: $(ls "$FDIR/Main_Figures/Figure2_CellChat_communication/" 2>/dev/null | wc -l) files"
echo "Figure3: $(ls "$FDIR/Main_Figures/Figure3_consensus_clustering/" 2>/dev/null | wc -l) files"
echo "Figure4: $(ls "$FDIR/Main_Figures/Figure4_immune_genomic_pathway/" 2>/dev/null | wc -l) files"
echo "Figure5: $(ls "$FDIR/Main_Figures/Figure5_ML_model_validation/" 2>/dev/null | wc -l) files"
echo "Figure6: $(ls "$FDIR/Main_Figures/Figure6_nomogram_prognosis/" 2>/dev/null | wc -l) files"
echo "Supplementary: $(ls "$FDIR/Supplementary_Figures/" 2>/dev/null | wc -l) files"
echo "QC: $(ls "$FDIR/QC_Figures/" 2>/dev/null | wc -l) files"
echo "Analysis: $(ls "$FDIR/Analysis_Figures/" 2>/dev/null | wc -l) files"
total=0
for d in "$FDIR/Main_Figures"/Figure*/ "$FDIR/Supplementary_Figures/" "$FDIR/QC_Figures/" "$FDIR/Analysis_Figures/"; do
    count=$(ls "$d"*.pdf 2>/dev/null | wc -l)
    total=$((total + count))
done
echo "Total classified: $total files"
echo ""
echo "=== Done! Original files preserved. ==="
