#!/usr/bin/env python3
"""
KIR Genotyping Evaluation with Novel Allele Filtering for HPRC 44 Samples
V2 - CDS-Aware Evaluation Logic

Evaluates KIR genotyping tools at 3-digit, 5-digit, and 7-digit resolution levels,
with the ability to filter out genes containing novel alleles in the ground truth.

Novel symbols in GT:
    # → CDS nonsynonymous variant (novel allele)
    = → CDS synonymous variant
    $ → non-CDS region variant (intron/UTR)
    + → genomic allele matching CDS-only allele
    e → fusion gene (e.g., KIR2DS4*00101e3DL1*03501)

Filter levels:
    0: Remove fusion genes only
    3: Also remove # (novel CDS)
    5: Also remove + and = (CDS-only, synonymous)
    7: Also remove $ (non-CDS variants)

Evaluation Logic (for non-filtered genes):
    - Precision denominator: only predictions with >= n-digit
    - CDS-only prefix matches (Pred >= n, GT < n) are excluded from Precision
    - 4-pass greedy matching to pair predicted and ground truth alleles
    - Fusion alleles are handled specially because most tools are unable to predict fusion alleles

Outputs:
    1. Recall table
    2. Precision table
    3. F1 score table
"""

import csv
import os
import re
from collections import defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = "{your_base_directory}"  # Set your base directory here
DATA_DIR = os.path.join(BASE_DIR, "{your_data_directory}")  # Set your data directory here
RESULT_DIR = os.path.join(BASE_DIR, "{your_result_directory}")  # Set your output directory

# Ground truth with novel symbols
GROUND_TRUTH_WITH_NOVEL = os.path.join(DATA_DIR, "groundtruth", "hprc_annotation_skirt.tsv")

# Sample ID file (44 HPRC samples)

# Tool configurations
TOOL_CONFIGS = {
    "GraphKIR_exonfirst_hg19": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hs37d5-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_exonfirst_hg38": {
        "file": os.path.join(DATA_DIR, "graph-kir-exonfirst-hg38-hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_full_hg19": {
        "file": os.path.join(DATA_DIR, "graphkir_hg19_hprc44.tsv"),
        "type": "graphkir"
    },
    "GraphKIR_full_hg38": {
        "file": os.path.join(DATA_DIR, "graphkir_hg38_hprc44.tsv"),
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


# Filter levels
FILTER_LEVELS = [
    {"name": "-fusion", "level": 0},
    {"name": "-fusion -#", "level": 3},
    {"name": "-fusion -# -+ -=", "level": 5},
    {"name": "-fusion -# -+ -= -$", "level": 7},
]

# KIR genes
KIR_GENES = [
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
    "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
    "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1"
]


# ============================================================================
# DATA LOADING
# ============================================================================

def clean_novel_symbols(allele):
    """Remove novel symbols from allele for matching: KIR2DL1*00302$ -> 00302"""
    # Remove gene prefix if present
    if '*' in allele:
        allele = allele.split('*', 1)[1]
    # Remove novel symbols
    allele = re.sub(r'[#$+=]', '', allele)
    return allele


def get_gene_name(allele):
    """Extract gene name from allele: KIR2DL1*00302$ -> KIR2DL1"""
    if '*' in allele:
        return allele.split('*')[0]
    return None


def load_ground_truth_with_novel(tsv_file):
    """
    Load ground truth data with novel symbols preserved.

    Returns:
        dict: {sample_id: {gene: [(allele_with_symbols, allele_clean), ...]}}
    """
    data = {}
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['id']
            alleles_str = row['alleles']

            gene_alleles = defaultdict(list)
            for entry in alleles_str.split('_'):
                if '*' not in entry:
                    continue
                gene = get_gene_name(entry)
                allele_part = entry.split('*', 1)[1]
                allele_clean = clean_novel_symbols(entry)
                gene_alleles[gene].append((allele_part, allele_clean))

            result = {}
            for gene in KIR_GENES:
                result[gene] = gene_alleles.get(gene, [])
            data[sample_id] = result

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
# NOVEL FILTERING
# ============================================================================

def has_fusion(allele_with_symbols):
    """Check if allele is a fusion gene (contains 'e' followed by digit)."""
    return bool(re.search(r'e\d', allele_with_symbols))


def get_excluded_genes_for_sample(gt_sample, filter_level):
    """
    Get set of genes to exclude for a single sample based on filter level.

    Args:
        gt_sample: {gene: [(allele_with_symbols, allele_clean), ...]}
        filter_level: 0, 3, 5, or 7

    Returns:
        set of gene names to exclude
    """
    excluded_genes = set()

    for gene, alleles in gt_sample.items():
        for allele_with_symbols, _ in alleles:
            # Level 0: fusion genes
            if filter_level >= 0 and has_fusion(allele_with_symbols):
                excluded_genes.add(gene)
                parts = allele_with_symbols.split('e')
                for part in parts[1:]:
                    if '*' in part:
                        other_gene = 'KIR' + part.split('*')[0]
                        excluded_genes.add(other_gene)

            # Level 3: novel CDS (#)
            if filter_level >= 3 and '#' in allele_with_symbols:
                excluded_genes.add(gene)

            # Level 5: CDS-only (+) and synonymous (=)
            if filter_level >= 5:
                if '+' in allele_with_symbols or '=' in allele_with_symbols:
                    excluded_genes.add(gene)

            # Level 7: non-CDS ($)
            if filter_level >= 7 and '$' in allele_with_symbols:
                excluded_genes.add(gene)

    return excluded_genes


# ============================================================================
# UTILITY FUNCTIONS (V2)
# ============================================================================

def is_fusion_allele(allele):
    """Check if allele is a fusion gene (contains 'e' followed by digit)."""
    return bool(re.search(r'e\d', allele))


def get_allele_part_before_fusion(allele):
    """
    Get the allele part before fusion marker 'e'.
    For '00101e2DP1*00201' returns '00101'.
    For non-fusion alleles, returns the original allele.
    """
    if not allele:
        return allele
    if '*' in allele:
        allele = allele.split('*', 1)[1]
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
    allele_part = get_allele_part_before_fusion(allele)
    digits_only = re.sub(r'\D', '', allele_part)
    return len(digits_only)


def truncate_allele(allele, n_digits):
    """
    Truncate allele to n digits.
    For fusion alleles, truncate from the part before 'e'.
    """
    if not allele:
        return ''
    allele_part = get_allele_part_before_fusion(allele)
    digits_only = re.sub(r'\D', '', allele_part)
    if len(digits_only) >= n_digits:
        return digits_only[:n_digits]
    else:
        return digits_only


# ============================================================================
# GREEDY MATCHING
# ============================================================================

def greedy_match_for_precision_v2(pred_alleles, gt_alleles_clean, n_digit):
    """
    Perform 4-pass greedy matching between predictions and ground truth at n-digit level.
    
    V2 Matching priority:
    1. Exact match at n-digit level (both >= n)
    2. Prefix match (Pred >= n, GT < n, prefix matches)
    3. Prefix match (Pred < n, GT >= n, prefix matches)
    4. Both < n and equal
    
    Returns:
        pairs: list of tuples ((pred_truncated, pred_orig_len), (gt_truncated, gt_orig_len, gt_is_fusion))
        unmatched_pred: list of (pred_truncated, pred_orig_len)
        unmatched_gt: list of (gt_truncated, gt_orig_len, gt_is_fusion)
    """
    # Prepare truncated alleles with original length info
    pred_prepared = []
    for p in pred_alleles:
        orig_len = get_digit_length(p)
        truncated = truncate_allele(p, n_digit) if orig_len >= n_digit else truncate_allele(p, orig_len)
        pred_prepared.append((truncated, orig_len))

    gt_prepared = []
    for g in gt_alleles_clean:
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
            continue
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len < n_digit:
                continue
            
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
            continue
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len < n_digit:
                continue
            
            if gt_len > 0 and pred_val[:gt_len] == gt_val:
                pairs.append(((pred_val, pred_len), (gt_val, gt_len, gt_fusion)))
                used_pred.add(pi)
                used_gt.add(gi)
                break

    # Pass 3: Prefix matches (Pred < n, GT >= n, prefix matches)
    for gi, (gt_val, gt_len, gt_fusion) in enumerate(gt_prepared):
        if gi in used_gt:
            continue
        if gt_len < n_digit:
            continue
        
        for pi, (pred_val, pred_len) in enumerate(pred_prepared):
            if pi in used_pred:
                continue
            if pred_len >= n_digit:
                continue
            
            if pred_len > 0 and gt_val[:pred_len] == pred_val:
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

def evaluate_recall_with_filter_v2(ground_truth, predicted, sample_ids, n_digit, filter_level):
    """
    Evaluate RECALL at n-digit level with novel filtering and V2 logic.

    Recall = TP / (GT alleles with >= n_digit)

    Returns:
        recall, tp, total_gt_qualified
    """
    total_tp = 0
    total_gt_qualified = 0

    for sample_id in sorted(sample_ids):
        if sample_id not in ground_truth:
            continue

        gt_sample = ground_truth[sample_id]
        pred_sample = predicted.get(sample_id, {gene: [] for gene in KIR_GENES})

        # Get excluded genes for this sample
        excluded_genes = get_excluded_genes_for_sample(gt_sample, filter_level)

        for gene in KIR_GENES:
            if gene in excluded_genes:
                continue

            gt_alleles_clean = [clean for _, clean in gt_sample.get(gene, [])]
            pred_alleles = pred_sample.get(gene, [])

            # Count GT alleles with >= n_digit (fusion included based on truncated length)
            gt_qualified = [g for g in gt_alleles_clean if get_digit_length(g) >= n_digit]
            total_gt_qualified += len(gt_qualified)

            if not gt_qualified:
                continue

            # Greedy matching using V2 logic
            pairs, _, _ = greedy_match_for_precision_v2(pred_alleles, gt_qualified, n_digit)

            # Count TP: pairs where both >= n digits and match at n-digit level (non-fusion GT)
            for (pred_val, pred_len), (gt_val, gt_len, gt_fusion) in pairs:
                if pred_len >= n_digit and gt_len >= n_digit:
                    if pred_val == gt_val and not gt_fusion:
                        total_tp += 1

    recall = total_tp / total_gt_qualified if total_gt_qualified > 0 else 0.0
    return recall, total_tp, total_gt_qualified


def evaluate_precision_with_filter_v2(ground_truth, predicted, sample_ids, n_digit, filter_level):
    """
    Evaluate PRECISION at n-digit level with novel filtering and V2 CDS-aware logic.

    V2 Rules for paired alleles:
    - Pred >= n, GT >= n, non-fusion, match → TP
    - Pred >= n, GT >= n, non-fusion, mismatch → FP
    - Pred >= n, GT >= n, fusion → FP
    - Pred >= n, GT < n, non-fusion, prefix match → Excluded (CDS-only)
    - Pred >= n, GT < n, fusion → Excluded (fusion < n)
    - Pred < n → Excluded (regardless of GT)

    Rules for unmatched Pred:
    - Pred >= n → FP (over-call)
    - Pred < n → Excluded

    Precision = TP / (TP + FP)

    Returns:
        precision, tp, denominator (tp + fp)
    """
    total_tp = 0
    total_fp = 0

    for sample_id in sorted(sample_ids):
        if sample_id not in ground_truth:
            continue

        gt_sample = ground_truth[sample_id]
        pred_sample = predicted.get(sample_id, {gene: [] for gene in KIR_GENES})

        # Get excluded genes for this sample
        excluded_genes = get_excluded_genes_for_sample(gt_sample, filter_level)

        for gene in KIR_GENES:
            if gene in excluded_genes:
                continue

            gt_alleles_clean = [clean for _, clean in gt_sample.get(gene, [])]
            pred_alleles = pred_sample.get(gene, [])

            if not pred_alleles:
                continue

            # Greedy matching using V2 logic
            pairs, unmatched_pred, _ = greedy_match_for_precision_v2(pred_alleles, gt_alleles_clean, n_digit)

            # Evaluate pairs
            for (pred_val, pred_len), (gt_val, gt_len, gt_fusion) in pairs:
                # Case: Pred < n → Excluded
                if pred_len < n_digit:
                    continue
                
                # From here: Pred >= n
                if gt_fusion:
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
                        # GT < n (CDS-only): Excluded (prefix match assumed since paired)
                        continue

            # Evaluate unmatched predictions (over-call)
            for (pred_val, pred_len) in unmatched_pred:
                if pred_len >= n_digit:
                    total_fp += 1
                # else: Pred < n → Excluded

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    return precision, total_tp, total_tp + total_fp


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def generate_table_text(results_by_filter, metric_name, metric_key):
    """Generate a metric table as text string."""
    lines = []
    lines.append("")
    lines.append("=" * 100)
    lines.append(f"{metric_name.upper()} TABLE (3/5/7-digit) - HPRC 44 Samples with Novel Filtering [v2 CDS-aware]")
    lines.append("=" * 100)
    lines.append("")

    # Header
    header = f"{'Filter':<25} {'Tool':<30} {'3-digit':<14} {'5-digit':<14} {'7-digit':<14}"
    lines.append(header)
    lines.append("-" * 100)

    for filter_info in FILTER_LEVELS:
        filter_name = filter_info['name']
        filter_results = results_by_filter[filter_name]

        # Only show GT counts for recall table (not for precision/f1)
        if metric_key == 'recall':
            gt_3d = filter_results['gt_counts'][3]
            gt_5d = filter_results['gt_counts'][5]
            gt_7d = filter_results['gt_counts'][7]
            lines.append(f"{filter_name:<25} {'(GT alleles)':<30} {gt_3d:<14} {gt_5d:<14} {gt_7d:<14}")

        # Tool rows
        for i, tool_result in enumerate(filter_results['tools']):
            tool_name = tool_result['tool']
            val_3d = tool_result[f'{metric_key}_3d']
            val_5d = tool_result[f'{metric_key}_5d']
            val_7d = tool_result[f'{metric_key}_7d']
            # For precision/f1, show filter name on first tool row
            filter_col = filter_name if (metric_key != 'recall' and i == 0) else ''
            lines.append(f"{filter_col:<25} {tool_name:<30} {val_3d:<14.4f} {val_5d:<14.4f} {val_7d:<14.4f}")

        lines.append("")

    lines.append("=" * 100)
    return "\n".join(lines)


def generate_precision_with_denom_text(results_by_filter):
    """Generate precision table with denominator as text string."""
    lines = []
    lines.append("")
    lines.append("=" * 115)
    lines.append("PRECISION TABLE WITH DENOMINATOR (3/5/7-digit) - HPRC 44 Samples with Novel Filtering [v2 CDS-aware]")
    lines.append("=" * 115)
    lines.append("")

    header = f"{'Filter':<25} {'Tool':<30} {'3-digit':<18} {'5-digit':<18} {'7-digit':<18}"
    lines.append(header)
    lines.append("-" * 115)

    for filter_info in FILTER_LEVELS:
        filter_name = filter_info['name']
        filter_results = results_by_filter[filter_name]

        for i, tool_result in enumerate(filter_results['tools']):
            tool_name = tool_result['tool']
            val_3d = f"{tool_result['precision_3d']:.4f}({tool_result['denom_3d']})"
            val_5d = f"{tool_result['precision_5d']:.4f}({tool_result['denom_5d']})"
            val_7d = f"{tool_result['precision_7d']:.4f}({tool_result['denom_7d']})"
            # Show filter name on first tool row
            filter_col = filter_name if i == 0 else ''
            lines.append(f"{filter_col:<25} {tool_name:<30} {val_3d:<18} {val_5d:<18} {val_7d:<18}")

        lines.append("")

    lines.append("=" * 115)
    return "\n".join(lines)


def save_table_tsv(results_by_filter, metric_key, output_file):
    """Save metric table to TSV."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Filter', 'Tool', '3-digit', '5-digit', '7-digit'])

        for filter_info in FILTER_LEVELS:
            filter_name = filter_info['name']
            filter_results = results_by_filter[filter_name]

            # Only output GT alleles row for recall (not for precision/f1)
            if metric_key == 'recall':
                gt_3d = filter_results['gt_counts'][3]
                gt_5d = filter_results['gt_counts'][5]
                gt_7d = filter_results['gt_counts'][7]
                writer.writerow([filter_name, '(GT alleles)', gt_3d, gt_5d, gt_7d])

            for i, tool_result in enumerate(filter_results['tools']):
                # For precision/f1, put filter name on first tool row
                filter_col = filter_name if (metric_key != 'recall' and i == 0) else ''
                writer.writerow([
                    filter_col,
                    tool_result['tool'],
                    f"{tool_result[f'{metric_key}_3d']:.4f}",
                    f"{tool_result[f'{metric_key}_5d']:.4f}",
                    f"{tool_result[f'{metric_key}_7d']:.4f}"
                ])

            writer.writerow([])


def save_combined_txt(all_tables_text, output_file):
    """Save all tables to a single TXT file."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(all_tables_text)


def save_precision_with_denom_tsv(results_by_filter, output_file):
    """Save precision table with denominator to TSV."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Filter', 'Tool', '3-digit', '5-digit', '7-digit'])

        for filter_info in FILTER_LEVELS:
            filter_name = filter_info['name']
            filter_results = results_by_filter[filter_name]

            for i, tool_result in enumerate(filter_results['tools']):
                filter_col = filter_name if i == 0 else ''
                val_3d = f"{tool_result['precision_3d']:.4f}({tool_result['denom_3d']})"
                val_5d = f"{tool_result['precision_5d']:.4f}({tool_result['denom_5d']})"
                val_7d = f"{tool_result['precision_7d']:.4f}({tool_result['denom_7d']})"
                writer.writerow([filter_col, tool_result['tool'], val_3d, val_5d, val_7d])

            writer.writerow([])


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("KIR Genotyping Evaluation with Novel Allele Filtering [v2 CDS-aware]")
    print("=" * 70)
    print()

    # Load sample IDs
    print("Loading data...")
    sample_ids = HPRC_44_SAMPLE_IDS
    print(f"  Sample IDs: {len(sample_ids)}")

    # Load ground truth with novel symbols
    ground_truth = load_ground_truth_with_novel(GROUND_TRUTH_WITH_NOVEL)
    print(f"  Ground truth (with novel): {len(ground_truth)} samples")

    # Filter to only samples in sample_ids
    ground_truth = {k: v for k, v in ground_truth.items() if k in sample_ids}
    print(f"  Ground truth (filtered): {len(ground_truth)} samples")

    # Load predictions for each tool
    predictions = {}
    for tool_name, config in TOOL_CONFIGS.items():
        file_path = config['file']
        tool_type = config['type']

        if not os.path.exists(file_path):
            print(f"  WARNING: File not found: {file_path}")
            continue

        predictions[tool_name] = load_prediction_data(file_path, tool_type)
        print(f"  Loaded {tool_name}: {len(predictions[tool_name])} samples")

    print()

    # Evaluate for each filter level
    results_by_filter = {}

    for filter_info in FILTER_LEVELS:
        filter_name = filter_info['name']
        filter_level = filter_info['level']

        print(f"Evaluating filter: {filter_name} (level={filter_level})")

        filter_results = {
            'gt_counts': {3: 0, 5: 0, 7: 0},
            'tools': []
        }

        gt_count_computed = False

        for tool_name, pred_data in predictions.items():
            print(f"  Evaluating {tool_name}...")

            tool_result = {'tool': tool_name}

            for n_digit in [3, 5, 7]:
                # Calculate recall (V2)
                recall, tp_recall, gt_qualified = evaluate_recall_with_filter_v2(
                    ground_truth, pred_data, sample_ids, n_digit, filter_level
                )

                # Calculate precision (V2)
                precision, tp_prec, denom = evaluate_precision_with_filter_v2(
                    ground_truth, pred_data, sample_ids, n_digit, filter_level
                )

                # Calculate F1
                f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

                tool_result[f'recall_{n_digit}d'] = recall
                tool_result[f'precision_{n_digit}d'] = precision
                tool_result[f'f1_{n_digit}d'] = f1
                tool_result[f'tp_{n_digit}d'] = tp_prec
                tool_result[f'gt_{n_digit}d'] = gt_qualified
                tool_result[f'denom_{n_digit}d'] = denom

                if not gt_count_computed:
                    filter_results['gt_counts'][n_digit] = gt_qualified

            gt_count_computed = True
            filter_results['tools'].append(tool_result)

        results_by_filter[filter_name] = filter_results
        print()

    # Generate tables
    print("Generating tables...")

    recall_table = generate_table_text(results_by_filter, "Recall", "recall")
    precision_table = generate_table_text(results_by_filter, "Precision", "precision")
    precision_denom_table = generate_precision_with_denom_text(results_by_filter)
    f1_table = generate_table_text(results_by_filter, "F1 Score", "f1")

    # Print to terminal
    print(recall_table)
    print(precision_table)
    print(precision_denom_table)
    print(f1_table)

    # Save files
    os.makedirs(RESULT_DIR, exist_ok=True)

    save_table_tsv(results_by_filter, 'recall', os.path.join(RESULT_DIR, 'recall_remove_novel.tsv'))
    save_table_tsv(results_by_filter, 'precision', os.path.join(RESULT_DIR, 'precision_remove_novel.tsv'))
    save_precision_with_denom_tsv(results_by_filter, os.path.join(RESULT_DIR, 'precision_with_denom_remove_novel.tsv'))
    save_table_tsv(results_by_filter, 'f1', os.path.join(RESULT_DIR, 'f1_remove_novel.tsv'))

    all_tables = "\n\n".join([recall_table, precision_table, precision_denom_table, f1_table])
    save_combined_txt(all_tables, os.path.join(RESULT_DIR, 'evaluation_remove_novel.txt'))

    print()
    print("Files saved to:")
    print(f"  {os.path.join(RESULT_DIR, 'recall_remove_novel.tsv')}")
    print(f"  {os.path.join(RESULT_DIR, 'precision_remove_novel.tsv')}")
    print(f"  {os.path.join(RESULT_DIR, 'precision_with_denom_remove_novel.tsv')}")
    print(f"  {os.path.join(RESULT_DIR, 'f1_remove_novel.tsv')}")
    print(f"  {os.path.join(RESULT_DIR, 'evaluation_remove_novel.txt')}")


if __name__ == "__main__":
    main()
