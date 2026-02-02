#!/usr/bin/env python3
"""
KIR Genotyping Multi-Resolution Evaluation for HPRC 44 Samples (v2)

Evaluates KIR genotyping tools at 3-digit, 5-digit, and 7-digit resolution levels.

New Precision Logic (CDS-aware):
- Pred >= n, GT >= n, match → TP
- Pred >= n, GT >= n, mismatch → FP
- Pred >= n, GT < n, prefix match → Excluded (not in denominator)
- Pred >= n, GT < n, prefix mismatch → FP
- Pred < n → Excluded (not in denominator)
- Unmatched Pred >= n → FP (over-call)
- Unmatched Pred < n → Excluded

Fusion Allele Handling:
- Fusion alleles (e.g., '00101e2DP1*00201') are truncated at 'e' for matching
- If GT is fusion and paired:
  - GT fusion >= n, Pred >= n → FP
  - GT fusion >= n, Pred < n → Excluded
  - GT fusion < n → Excluded

Outputs:
1. Recall table (3/5/7-digit)
2. Precision table (3/5/7-digit)
3. Terminal display, TSV export, and TXT log
"""

import csv
import os
import re
from collections import defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================

SAMPLE_ID_DIR = os.path.join(BASE_DIR, "sampleID")
BASE_DIR = "{your_base_directory}"  # Set your base directory here
DATA_DIR = os.path.join(BASE_DIR, "{your_data_directory}")  # Set your data directory here
RESULT_DIR = os.path.join(BASE_DIR, "{your_result_directory}")  # Set your output directory
TABLE_DIR = os.path.join(RESULT_DIR, "tables")

# Ground truth (skirt annotation, 44 samples)
GROUND_TRUTH_CSV = os.path.join(DATA_DIR, "groundtruth", "hprc_summary_v1_2_e.tsv")

# Tool configurations
TOOL_CONFIGS = {
    "GraphKIR_exonfirst_hg19": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hs37d5-homo-t1-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR-exonfirst-hg38(noAlt)": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hg38noalt-homo-t1-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR-exonfirst-hg38-extract1(mainOnly)": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hg38new-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_exonfirst_hg38_extract2(alt)": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hg38new2-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_exonfirst_hg38_extract3(alt + noInterGenic)": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hg38new3-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_full_hg19": {
        "file": os.path.join(DATA_DIR, "graphkir_hg19_homo_hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_full_hg38": {
        "file": os.path.join(DATA_DIR, "graphkir_hg38_homo_hprc.tsv"),
        "type": "graphkir"
    },
    "PING-wgs": {
        "file": os.path.join(DATA_DIR, "hprc_pingsample.index_hs37d5.bwa.part_strict.result_ping_wgs.merge.tsv"),
        "type": "graphkir"
    },
    "Geny": {
        "file": os.path.join(DATA_DIR, "geny_hprc44.txt"),
        "type": "geny"
    },
}

# HPRC sample IDs
HPRC_44_SAMPLE_IDS = {
    "HG002", "HG00438", "HG005", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741",
    "HG01071", "HG01106", "HG01109", "HG01175", "HG01243", "HG01258", "HG01358", "HG01361",
    "HG01891", "HG01928", "HG01952", "HG01978", "HG02055", "HG02080", "HG02109", "HG02145",
    "HG02148", "HG02257", "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818",
    "HG02886", "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579",
    "NA18906", "NA19240", "NA20129", "NA21309",
}

# Output files
OUTPUT_TSV_RECALL = os.path.join(TABLE_DIR, "evaluation_hprc_alldigit_recall.tsv")
OUTPUT_TSV_PRECISION = os.path.join(TABLE_DIR, "evaluation_hprc_alldigit_precision.tsv")
OUTPUT_TSV_F1 = os.path.join(TABLE_DIR, "evaluation_hprc_alldigit_f1.tsv")
OUTPUT_TXT = os.path.join(TABLE_DIR, "evaluation_hprc_alldigit.txt")

# KIR genes
KIR_GENES = [
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
    "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
    "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1"
]


# ============================================================================
# DATA LOADING
# ============================================================================


def load_ground_truth(csv_file):
    """Load ground truth data from CSV file."""
    data = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row['ID'].replace('.bam', '')
            gene_alleles = {}
            for gene in KIR_GENES:
                alleles_str = row[gene].strip()
                if alleles_str:
                    alleles = [
                        a[1:] if a.startswith('`') else a
                        for a in alleles_str.split(';')
                        if a
                    ]
                    gene_alleles[gene] = alleles
                else:
                    gene_alleles[gene] = []
            data[sample_id] = gene_alleles
    return data


def parse_graphkir_alleles(alleles_string):
    """Parse Graph-KIR alleles string."""
    gene_alleles = defaultdict(list)
    entries = alleles_string.strip().split('_')
    for entry in entries:
        if '*' not in entry:
            continue
        gene, allele = entry.split('*', 1)
        allele = allele.rstrip('e')
        gene_alleles[gene].append(allele)
    return gene_alleles


def load_graphkir_tsv(input_tsv):
    """Load Graph-KIR TSV format."""
    results = {}
    with open(input_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['id']
            alleles_string = row['alleles']
            gene_alleles = parse_graphkir_alleles(alleles_string)
            formatted_alleles = {}
            for gene in KIR_GENES:
                if gene in gene_alleles and gene_alleles[gene]:
                    formatted_alleles[gene] = gene_alleles[gene]
                else:
                    formatted_alleles[gene] = []
            results[sample_id] = formatted_alleles
    return results


def parse_geny_output(input_file):
    """Parse Geny output file."""
    results = {}
    current_sample = None
    current_alleles = defaultdict(list)

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Sample:'):
                if current_sample:
                    results[current_sample] = dict(current_alleles)
                current_sample = line.split('Sample:')[1].strip()
                current_alleles = defaultdict(list)
            elif line.startswith('[ilp] KIR'):
                parts = line.split()
                if len(parts) >= 4:
                    gene = parts[1]
                    allele = parts[2]
                    if gene in KIR_GENES:
                        current_alleles[gene].append(allele)
        if current_sample:
            results[current_sample] = dict(current_alleles)
    return results


def load_geny_data(input_txt):
    """Load Geny output."""
    results = parse_geny_output(input_txt)
    formatted_results = {}
    for sample_id, gene_dict in results.items():
        formatted_alleles = {}
        for gene in KIR_GENES:
            formatted_alleles[gene] = gene_dict.get(gene, [])
        formatted_results[sample_id] = formatted_alleles
    return formatted_results


def load_prediction_data(file_path, tool_type):
    """Load prediction data based on tool type."""
    if tool_type == "graphkir":
        return load_graphkir_tsv(file_path)
    elif tool_type == "geny":
        return load_geny_data(file_path)
    else:
        raise ValueError(f"Unknown tool type: {tool_type}")


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def is_fusion_allele(allele):
    """Check if allele is a fusion gene (contains 'e' followed by digit)."""
    # Match pattern like '00101e3DL1*03501' or '00101e2DP1*00201'
    return bool(re.search(r'e\d', allele))


def get_allele_part_before_fusion(allele):
    """
    Get the allele part before fusion marker 'e'.
    For '00101e2DP1*00201' returns '00101'.
    For non-fusion alleles, returns the original allele.
    """
    if not allele:
        return allele
    # Remove gene prefix if present
    if '*' in allele:
        allele = allele.split('*', 1)[1]
    # Check for fusion marker
    match = re.search(r'e\d', allele)
    if match:
        return allele[:match.start()]
    return allele


def get_digit_length(allele):
    """
    Get the number of digits in an allele string.
    For fusion alleles, only count digits before 'e'.
    """
    if not allele:
        return 0
    # Get the part before fusion marker (if any)
    allele_part = get_allele_part_before_fusion(allele)
    # Extract only digits
    digits_only = re.sub(r'\D', '', allele_part)
    return len(digits_only)


def truncate_allele(allele, n_digits):
    """
    Truncate allele to n digits.
    For fusion alleles, truncate from the part before 'e'.
    """
    if not allele:
        return ''
    # Get the part before fusion marker (if any)
    allele_part = get_allele_part_before_fusion(allele)
    # Extract only digits
    digits_only = re.sub(r'\D', '', allele_part)
    # Truncate
    if len(digits_only) >= n_digits:
        return digits_only[:n_digits]
    else:
        return digits_only


# ============================================================================
# GREEDY MATCHING
# ============================================================================

def greedy_match_for_precision(pred_alleles, gt_alleles, n_digit):
    """
    Perform greedy matching between predictions and ground truth at n-digit level.
    
    Matching priority:
    1. Exact match at n-digit level (both >= n)
    2. Prefix match (one side < n, the other >= n, prefix matches)
    3. Both < n and equal
    
    Returns:
        pairs: list of tuples ((pred_truncated, pred_orig_len), (gt_truncated, gt_orig_len, gt_is_fusion))
        unmatched_pred: list of (pred_truncated, pred_orig_len)
        unmatched_gt: list of (gt_truncated, gt_orig_len, gt_is_fusion)
    """
    # Step 1: Prepare truncated alleles with original length info
    pred_prepared = []
    for p in pred_alleles:
        orig_len = get_digit_length(p)
        truncated = truncate_allele(p, n_digit) if orig_len >= n_digit else truncate_allele(p, orig_len)
        pred_prepared.append((truncated, orig_len))

    gt_prepared = []
    for g in gt_alleles:
        orig_len = get_digit_length(g)
        truncated = truncate_allele(g, n_digit) if orig_len >= n_digit else truncate_allele(g, orig_len)
        fusion = is_fusion_allele(g)
        gt_prepared.append((truncated, orig_len, fusion))

    pairs = []
    used_pred = set()
    used_gt = set()

    # Pass 1: Exact matches at n-digit level (both >= n, truncated values equal)
    for gi, (gt_val, gt_len, gt_fusion) in enumerate(gt_prepared):
        if gi in used_gt:
            continue
        if gt_len < n_digit:
            continue  # Skip GT < n in this pass
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len < n_digit:
                continue  # Skip Pred < n in this pass
            
            # Both >= n-digit, compare truncated values
            if pred_val == gt_val:
                pairs.append(((pred_val, pred_len), (gt_val, gt_len, gt_fusion)))
                used_pred.add(pi)
                used_gt.add(gi)
                break

    # Pass 2: Prefix matches (Pred >= n, GT < n, prefix matches)
    for gi, (gt_val, gt_len, gt_fusion) in enumerate(gt_prepared):
        if gi in used_gt:
            continue
        if gt_len >= n_digit:
            continue  # Skip GT >= n in this pass
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len < n_digit:
                continue  # Skip Pred < n in this pass
            
            # Pred >= n, GT < n: check if GT is prefix of Pred
            if pred_val[:gt_len] == gt_val:
                pairs.append(((pred_val, pred_len), (gt_val, gt_len, gt_fusion)))
                used_pred.add(pi)
                used_gt.add(gi)
                break

    # Pass 3: Prefix matches (Pred < n, GT >= n, prefix matches)
    for gi, (gt_val, gt_len, gt_fusion) in enumerate(gt_prepared):
        if gi in used_gt:
            continue
        if gt_len < n_digit:
            continue  # Skip GT < n in this pass
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len >= n_digit:
                continue  # Skip Pred >= n in this pass
            
            # Pred < n, GT >= n: check if Pred is prefix of GT
            if gt_val[:pred_len] == pred_val:
                pairs.append(((pred_val, pred_len), (gt_val, gt_len, gt_fusion)))
                used_pred.add(pi)
                used_gt.add(gi)
                break

    # Pass 4: Both < n and equal
    for gi, (gt_val, gt_len, gt_fusion) in enumerate(gt_prepared):
        if gi in used_gt:
            continue
        if gt_len >= n_digit:
            continue
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len >= n_digit:
                continue
            
            # Both < n, must be exactly equal
            if pred_val == gt_val:
                pairs.append(((pred_val, pred_len), (gt_val, gt_len, gt_fusion)))
                used_pred.add(pi)
                used_gt.add(gi)
                break

    # Collect unmatched
    unmatched_pred = [pred_prepared[i] for i in range(len(pred_prepared)) if i not in used_pred]
    unmatched_gt = [gt_prepared[i] for i in range(len(gt_prepared)) if i not in used_gt]

    return pairs, unmatched_pred, unmatched_gt


# ============================================================================
# EVALUATION FUNCTIONS
# ============================================================================

def evaluate_recall_at_digit(ground_truth, predicted, sample_ids, n_digit):
    """
    Evaluate recall at n-digit level.

    Recall = TP / (number of GT alleles with >= n digits)

    Note: Fusion alleles are included in denominator if their truncated part
    (before 'e') >= n-digit. However, fusion alleles cannot be TP.
    """
    total_tp = 0
    total_gt_qualified = 0

    gt_samples = set(ground_truth.keys())
    if sample_ids is not None:
        gt_samples = gt_samples & sample_ids

    for sample_id in sorted(gt_samples):
        for gene in KIR_GENES:
            if sample_id in predicted:
                pred_alleles = predicted[sample_id].get(gene, [])
            else:
                pred_alleles = []

            gt_alleles = ground_truth[sample_id].get(gene, [])

            # Count GT alleles with >= n digits (fusion included based on truncated length)
            gt_qualified = [g for g in gt_alleles if get_digit_length(g) >= n_digit]
            total_gt_qualified += len(gt_qualified)

            if not gt_qualified:
                continue

            # Greedy matching
            pairs, _, _ = greedy_match_for_precision(pred_alleles, gt_qualified, n_digit)

            # Count TP: pairs where both >= n digits and match at n-digit level
            for (pred_val, pred_len), (gt_val, gt_len, gt_fusion) in pairs:
                if pred_len >= n_digit and gt_len >= n_digit:
                    if pred_val == gt_val and not gt_fusion:
                        total_tp += 1

    recall = total_tp / total_gt_qualified if total_gt_qualified > 0 else 0.0
    return recall, total_tp, total_gt_qualified


def evaluate_precision_at_digit(ground_truth, predicted, sample_ids, n_digit):
    """
    Evaluate precision at n-digit level with CDS-aware logic.

    Rules for paired alleles:
    - Pred >= n, GT >= n, non-fusion, match → TP
    - Pred >= n, GT >= n, non-fusion, mismatch → FP
    - Pred >= n, GT >= n, fusion → FP
    - Pred >= n, GT < n, non-fusion, prefix match → Excluded
    - Pred >= n, GT < n, non-fusion, prefix mismatch → FP (won't happen due to matching)
    - Pred >= n, GT < n, fusion, prefix match → Excluded (fusion < n)
    - Pred >= n, GT fusion >= n → FP
    - Pred < n → Excluded (regardless of GT)

    Rules for unmatched Pred:
    - Pred >= n → FP (over-call)
    - Pred < n → Excluded

    Precision = TP / (TP + FP)
    """
    total_tp = 0
    total_fp = 0

    gt_samples = set(ground_truth.keys())
    if sample_ids is not None:
        gt_samples = gt_samples & sample_ids

    for sample_id in sorted(gt_samples):
        for gene in KIR_GENES:
            if sample_id in predicted:
                pred_alleles = predicted[sample_id].get(gene, [])
            else:
                pred_alleles = []

            gt_alleles = ground_truth[sample_id].get(gene, [])

            # Greedy matching
            pairs, unmatched_pred, _ = greedy_match_for_precision(pred_alleles, gt_alleles, n_digit)

            # Evaluate pairs
            for (pred_val, pred_len), (gt_val, gt_len, gt_fusion) in pairs:
                # Case: Pred < n → Excluded
                if pred_len < n_digit:
                    continue  # Excluded from this n-digit level
                
                # From here: Pred >= n
                
                if gt_fusion:
                    # GT is fusion allele
                    if gt_len >= n_digit:
                        # GT fusion >= n, Pred >= n → FP
                        total_fp += 1
                    else:
                        # GT fusion < n → Excluded
                        continue
                else:
                    # GT is not fusion
                    if gt_len >= n_digit:
                        # Both >= n: compare at n-digit level
                        if pred_val == gt_val:
                            total_tp += 1
                        else:
                            total_fp += 1
                    else:
                        # GT < n (CDS-only): check prefix match
                        # If paired, it means prefix matched → Excluded
                        # (prefix mismatch wouldn't be paired)
                        continue  # Excluded

            # Evaluate unmatched predictions (over-call)
            for (pred_val, pred_len) in unmatched_pred:
                if pred_len >= n_digit:
                    total_fp += 1
                # else: Pred < n → Excluded

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    return precision, total_tp, total_fp


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def generate_recall_table_text(results, gt_counts):
    """Generate recall table as text string."""
    lines = []
    lines.append("")
    lines.append("=" * 75)
    lines.append("RECALL TABLE (3/5/7-digit) - HPRC 44 Samples [v2 CDS-aware]")
    lines.append("=" * 75)
    lines.append("")

    header1 = f"{'Tool':<30} {'3-digit':<14} {'5-digit':<14} {'7-digit':<14}"
    lines.append(header1)

    header2 = f"{'':<30} {'(GT=' + str(gt_counts[3]) + ')':<14} {'(GT=' + str(gt_counts[5]) + ')':<14} {'(GT=' + str(gt_counts[7]) + ')':<14}"
    lines.append(header2)
    lines.append("-" * 75)

    for row in results:
        line = (
            f"{row['tool']:<30} "
            f"{row['recall_3d']:<14.4f} "
            f"{row['recall_5d']:<14.4f} "
            f"{row['recall_7d']:<14.4f}"
        )
        lines.append(line)

    lines.append("=" * 75)
    lines.append("")
    return "\n".join(lines)


def generate_precision_table_text(results):
    """Generate precision table as text string."""
    lines = []
    lines.append("")
    lines.append("=" * 75)
    lines.append("PRECISION TABLE (3/5/7-digit) - HPRC 44 Samples [v2 CDS-aware]")
    lines.append("=" * 75)
    lines.append("")

    header = f"{'Tool':<30} {'3-digit':<14} {'5-digit':<14} {'7-digit':<14}"
    lines.append(header)
    lines.append("-" * 75)

    for row in results:
        line = (
            f"{row['tool']:<30} "
            f"{row['precision_3d']:<14.4f} "
            f"{row['precision_5d']:<14.4f} "
            f"{row['precision_7d']:<14.4f}"
        )
        lines.append(line)

    lines.append("=" * 75)
    lines.append("")
    return "\n".join(lines)


def generate_precision_table_with_denom_text(results):
    """Generate precision table with denominator counts: 0.9xx(xxx)"""
    lines = []
    lines.append("")
    lines.append("=" * 90)
    lines.append("PRECISION TABLE WITH DENOMINATOR (3/5/7-digit) - HPRC 44 Samples [v2 CDS-aware]")
    lines.append("=" * 90)
    lines.append("")

    header = f"{'Tool':<30} {'3-digit':<18} {'5-digit':<18} {'7-digit':<18}"
    lines.append(header)
    lines.append("-" * 90)

    for row in results:
        val_3d = f"{row['precision_3d']:.4f}({row['denom_3d']})"
        val_5d = f"{row['precision_5d']:.4f}({row['denom_5d']})"
        val_7d = f"{row['precision_7d']:.4f}({row['denom_7d']})"
        line = f"{row['tool']:<30} {val_3d:<18} {val_5d:<18} {val_7d:<18}"
        lines.append(line)

    lines.append("=" * 90)
    lines.append("")
    return "\n".join(lines)


def generate_precision_table_fraction_text(results):
    """Generate precision table as fraction: TP/denominator"""
    lines = []
    lines.append("")
    lines.append("=" * 90)
    lines.append("PRECISION TABLE (TP/Total) (3/5/7-digit) - HPRC 44 Samples [v2 CDS-aware]")
    lines.append("=" * 90)
    lines.append("")

    header = f"{'Tool':<30} {'3-digit':<18} {'5-digit':<18} {'7-digit':<18}"
    lines.append(header)
    lines.append("-" * 90)

    for row in results:
        val_3d = f"{row['tp_3d']}/{row['denom_3d']}"
        val_5d = f"{row['tp_5d']}/{row['denom_5d']}"
        val_7d = f"{row['tp_7d']}/{row['denom_7d']}"
        line = f"{row['tool']:<30} {val_3d:<18} {val_5d:<18} {val_7d:<18}"
        lines.append(line)

    lines.append("=" * 90)
    lines.append("")
    return "\n".join(lines)


def generate_f1_table_text(results):
    """Generate F1 score table as text string."""
    lines = []
    lines.append("")
    lines.append("=" * 75)
    lines.append("F1 SCORE TABLE (3/5/7-digit) - HPRC 44 Samples [v2 CDS-aware]")
    lines.append("=" * 75)
    lines.append("")

    header = f"{'Tool':<30} {'3-digit':<14} {'5-digit':<14} {'7-digit':<14}"
    lines.append(header)
    lines.append("-" * 75)

    for row in results:
        line = (
            f"{row['tool']:<30} "
            f"{row['f1_3d']:<14.4f} "
            f"{row['f1_5d']:<14.4f} "
            f"{row['f1_7d']:<14.4f}"
        )
        lines.append(line)

    lines.append("=" * 75)
    lines.append("")
    return "\n".join(lines)


def save_recall_tsv(results, gt_counts, output_file):
    """Save recall table to TSV."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Tool', f'3-digit (GT={gt_counts[3]})', f'5-digit (GT={gt_counts[5]})', f'7-digit (GT={gt_counts[7]})'])

        for row in results:
            writer.writerow([
                row['tool'],
                f"{row['recall_3d']:.4f}",
                f"{row['recall_5d']:.4f}",
                f"{row['recall_7d']:.4f}"
            ])


def save_precision_tsv(results, output_file):
    """Save precision table to TSV."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Tool', '3-digit', '5-digit', '7-digit'])

        for row in results:
            writer.writerow([
                row['tool'],
                f"{row['precision_3d']:.4f}",
                f"{row['precision_5d']:.4f}",
                f"{row['precision_7d']:.4f}"
            ])


def save_f1_tsv(results, output_file):
    """Save F1 score table to TSV."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Tool', '3-digit', '5-digit', '7-digit'])

        for row in results:
            # Calculate F1 scores: F1 = 2 * (Precision * Recall) / (Precision + Recall)
            f1_3d = 2 * row['precision_3d'] * row['recall_3d'] / (row['precision_3d'] + row['recall_3d']) if (row['precision_3d'] + row['recall_3d']) > 0 else 0.0
            f1_5d = 2 * row['precision_5d'] * row['recall_5d'] / (row['precision_5d'] + row['recall_5d']) if (row['precision_5d'] + row['recall_5d']) > 0 else 0.0
            f1_7d = 2 * row['precision_7d'] * row['recall_7d'] / (row['precision_7d'] + row['recall_7d']) if (row['precision_7d'] + row['recall_7d']) > 0 else 0.0

            writer.writerow([
                row['tool'],
                f"{f1_3d:.4f}",
                f"{f1_5d:.4f}",
                f"{f1_7d:.4f}"
            ])


def save_txt(text_content, output_file):
    """Save text content to file."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(text_content)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("Loading data...")

    # Load sample IDs
    sample_ids = HPRC_44_SAMPLE_IDS
    print(f"  Sample IDs: {len(sample_ids)}")

    # Load ground truth
    ground_truth = load_ground_truth(GROUND_TRUTH_CSV)
    print(f"  Ground truth: {len(ground_truth)} samples")

    # Count GT alleles at each digit level
    # Fusion alleles are included if their truncated part (before 'e') >= n-digit
    gt_counts = {3: 0, 5: 0, 7: 0}
    for sample_id in sample_ids:
        if sample_id not in ground_truth:
            continue
        for gene in KIR_GENES:
            for allele in ground_truth[sample_id].get(gene, []):
                # get_digit_length handles fusion by taking part before 'e'
                digit_len = get_digit_length(allele)
                if digit_len >= 3:
                    gt_counts[3] += 1
                if digit_len >= 5:
                    gt_counts[5] += 1
                if digit_len >= 7:
                    gt_counts[7] += 1

    print(f"  GT allele counts: 3-digit={gt_counts[3]}, 5-digit={gt_counts[5]}, 7-digit={gt_counts[7]}")

    # Evaluate each tool
    results = []

    for tool_name, config in TOOL_CONFIGS.items():
        file_path = config['file']
        tool_type = config['type']

        if not os.path.exists(file_path):
            print(f"  WARNING: File not found: {file_path}")
            continue

        print(f"  Evaluating {tool_name}...")

        # Load predictions
        predicted = load_prediction_data(file_path, tool_type)

        # Evaluate at each digit level
        recall_3d, tp_r3d, total_r3d = evaluate_recall_at_digit(ground_truth, predicted, sample_ids, 3)
        recall_5d, tp_r5d, total_r5d = evaluate_recall_at_digit(ground_truth, predicted, sample_ids, 5)
        recall_7d, tp_r7d, total_r7d = evaluate_recall_at_digit(ground_truth, predicted, sample_ids, 7)

        precision_3d, tp_p3d, fp_p3d = evaluate_precision_at_digit(ground_truth, predicted, sample_ids, 3)
        precision_5d, tp_p5d, fp_p5d = evaluate_precision_at_digit(ground_truth, predicted, sample_ids, 5)
        precision_7d, tp_p7d, fp_p7d = evaluate_precision_at_digit(ground_truth, predicted, sample_ids, 7)

        # Calculate F1 scores
        def calc_f1(p, r):
            return 2 * p * r / (p + r) if (p + r) > 0 else 0.0

        results.append({
            'tool': tool_name,
            'recall_3d': recall_3d,
            'recall_5d': recall_5d,
            'recall_7d': recall_7d,
            'precision_3d': precision_3d,
            'precision_5d': precision_5d,
            'precision_7d': precision_7d,
            'f1_3d': calc_f1(precision_3d, recall_3d),
            'f1_5d': calc_f1(precision_5d, recall_5d),
            'f1_7d': calc_f1(precision_7d, recall_7d),
            # Store TP and denominator (TP+FP) for precision
            'tp_3d': tp_p3d,
            'tp_5d': tp_p5d,
            'tp_7d': tp_p7d,
            'denom_3d': tp_p3d + fp_p3d,
            'denom_5d': tp_p5d + fp_p5d,
            'denom_7d': tp_p7d + fp_p7d,
        })

    # Generate all table texts
    all_tables = []
    all_tables.append(generate_recall_table_text(results, gt_counts))
    all_tables.append(generate_precision_table_text(results))
    all_tables.append(generate_f1_table_text(results))
    all_tables.append(generate_precision_table_with_denom_text(results))
    all_tables.append(generate_precision_table_fraction_text(results))

    # Combine all tables
    all_text = "\n".join(all_tables)

    # Print tables to terminal
    print(all_text)

    # Save to TSV
    save_recall_tsv(results, gt_counts, OUTPUT_TSV_RECALL)
    save_precision_tsv(results, OUTPUT_TSV_PRECISION)
    save_f1_tsv(results, OUTPUT_TSV_F1)

    # Save to TXT
    save_txt(all_text, OUTPUT_TXT)

    print(f"\nTables saved to:")
    print(f"  {OUTPUT_TSV_RECALL}")
    print(f"  {OUTPUT_TSV_PRECISION}")
    print(f"  {OUTPUT_TSV_F1}")
    print(f"  {OUTPUT_TXT}")


if __name__ == "__main__":
    main()
