import re
import time
import statistics as stats

# Import from your BLAST implementation
from BLAST import (
    load_pickle, load_genome, load_conf, load_queries,
    find_seeds_exact, find_seeds_variants, should_use_variant_fallback,
    ungapped_extend, gapped_extend_flanks,
    compute_max_prob_from_word, WORD_LEN
)

# ----------------------------
# Helpers: parse truth, intervals, metrics
# ----------------------------

def parse_truth_from_qid(qid: str):
    """Parse db_start / db_end from FASTA header."""
    m1 = re.search(r"db_start=(\d+)", qid)
    m2 = re.search(r"db_end=(\d+)", qid)
    if not (m1 and m2):
        return None
    return int(m1.group(1)), int(m2.group(1))  # [S, E)

def pred_interval_from_alignment(aln):
    """Convert alignment traceback to predicted interval [start, end)."""
    gp = aln.get("aligned_gpos", [])
    gpos = [p for p in gp if p is not None]
    if not gpos:
        return None
    return (min(gpos), max(gpos) + 1)

def overlap_len(a, b):
    (a0, a1), (b0, b1) = a, b
    return max(0, min(a1, b1) - max(a0, b0))

def union_len(a, b):
    inter = overlap_len(a, b)
    (a0, a1), (b0, b1) = a, b
    return (a1 - a0) + (b1 - b0) - inter

def iou(a, b):
    u = union_len(a, b)
    return 0.0 if u == 0 else overlap_len(a, b) / u

def hit_with_tolerance(pred, truth, tol=50):
    """Hit if start within tol OR intervals overlap."""
    (ps, pe), (S, E) = pred, truth
    return (abs(ps - S) <= tol) or (overlap_len(pred, truth) > 0)

# ----------------------------
# Runtime / CPU / Memory snapshot
# ----------------------------

def get_resource_snapshot():
    """
    Cross-platform snapshot:
      - cpu_time_s: process CPU time
      - rss_mb: resident set size (memory), if available
    """
    cpu_time_s = time.process_time()

    # Try psutil for memory (works on Windows/macOS/Linux)
    try:
        import psutil
        rss_bytes = psutil.Process().memory_info().rss
        rss_mb = rss_bytes / (1024 * 1024)
    except Exception:
        rss_mb = None

    return cpu_time_s, rss_mb

# ----------------------------
# Your requested pipeline:
# TOP-200 seeds -> ungapped; TOP-100 ungapped -> gapped
# ----------------------------

def run_pipeline_top200_top100(
    query, genome, seed_db, score_db,
    max_prob_threshold=0.0,
    seeds_top_k=200,
    ungapped_top_k=100,
):
    # 1) seeds
    exact = find_seeds_exact(query, seed_db)
    all_seeds = exact[:]

    if should_use_variant_fallback(exact, max_prob_threshold):
        all_seeds.extend(find_seeds_variants(query, seed_db))

    if not all_seeds:
        return None, None, [], []

    # 2) top-200 seeds by probability
    all_seeds.sort(key=lambda x: x[2], reverse=True)
    top_seeds = all_seeds[:seeds_top_k]

    # 3) ungapped extension
    ungapped = []
    for q_pos, g_pos, _, _ in top_seeds:
        if q_pos < 0 or q_pos + WORD_LEN > len(query):
            continue
        if g_pos < 0 or g_pos + WORD_LEN > len(genome):
            continue
        ungapped.append(ungapped_extend(query, genome, score_db, q_pos, g_pos))

    if not ungapped:
        return None, None, [], []

    ungapped.sort(key=lambda x: x["score"], reverse=True)
    best_ungapped = ungapped[0]
    top_ungapped = ungapped[:ungapped_top_k]

    # 4) gapped extension on top-100 ungapped
    finals = []
    for u in top_ungapped:
        core, full = gapped_extend_flanks(query, genome, score_db, u)
        final = full if full["score"] > core["score"] else core
        finals.append(final)

    finals.sort(key=lambda x: x["score"], reverse=True)
    best_final = finals[0] if finals else None
    return best_final, best_ungapped, top_ungapped, finals

# ----------------------------
# Evaluation over all queries
# ----------------------------

def evaluate_queries(
    queries,
    genome,
    seed_db,
    score_db,
    tol_list=(20, 50),
    max_prob_threshold=0.0,
    seeds_top_k=200,
    ungapped_top_k=100,
):
    records = []

    wall0 = time.perf_counter()
    cpu0, mem0 = get_resource_snapshot()

    for qid, query in queries:
        truth = parse_truth_from_qid(qid)
        print("Evaluating query:", qid)
        if truth is None:
            # skip queries without ground truth coordinates
            continue

        t0 = time.perf_counter()
        c0, _ = get_resource_snapshot()

        best_final, best_ungapped, top_ungapped, top_gapped = run_pipeline_top200_top100(
            query=query,
            genome=genome,
            seed_db=seed_db,
            score_db=score_db,
            max_prob_threshold=max_prob_threshold,
            seeds_top_k=seeds_top_k,
            ungapped_top_k=ungapped_top_k,
        )

        c1, _ = get_resource_snapshot()
        t1 = time.perf_counter()

        rec = {
            "qid": qid,
            "truth": truth,
            "found": best_final is not None,
            "wall_s": (t1 - t0),
            "cpu_s": (c1 - c0),
            "n_ungapped": len(top_ungapped),
            "n_gapped": len(top_gapped),
        }

        if best_final is not None:
            pred = pred_interval_from_alignment(best_final)
            rec["pred"] = pred
            if pred is not None:
                rec["iou"] = iou(pred, truth)
                rec["start_err"] = pred[0] - truth[0]
                rec["end_err"] = pred[1] - truth[1]
                for tol in tol_list:
                    rec[f"hit@{tol}"] = hit_with_tolerance(pred, truth, tol=tol)
            else:
                rec["iou"] = 0.0
                rec["start_err"] = None
                rec["end_err"] = None
                for tol in tol_list:
                    rec[f"hit@{tol}"] = False
        else:
            rec["pred"] = None
            rec["iou"] = 0.0
            rec["start_err"] = None
            rec["end_err"] = None
            for tol in tol_list:
                rec[f"hit@{tol}"] = False

        records.append(rec)

    cpu1, mem1 = get_resource_snapshot()
    wall1 = time.perf_counter()

    # ---- Summary ----
    N = len(records)
    summary = {
        "N": N,
        "wall_total_s": wall1 - wall0,
        "cpu_total_s": cpu1 - cpu0,
        "mem_mb_max": mem1 if mem1 is not None else None,
    }

    for tol in tol_list:
        hits = [1.0 if r.get(f"hit@{tol}", False) else 0.0 for r in records]
        summary[f"acc@{tol}"] = (sum(hits) / N) if N else 0.0

    ious = [r["iou"] for r in records]
    summary["iou_mean"] = stats.mean(ious) if N else 0.0
    summary["iou_median"] = stats.median(ious) if N else 0.0

    start_errs = [r["start_err"] for r in records if r["start_err"] is not None]
    end_errs = [r["end_err"] for r in records if r["end_err"] is not None]

    if start_errs:
        start_errs_sorted = sorted(start_errs)
        summary["start_err_mean"] = stats.mean(start_errs)
        summary["start_err_median"] = stats.median(start_errs)
        summary["start_err_p95"] = start_errs_sorted[int(0.95 * (len(start_errs_sorted) - 1))]

    if end_errs:
        end_errs_sorted = sorted(end_errs)
        summary["end_err_mean"] = stats.mean(end_errs)
        summary["end_err_median"] = stats.median(end_errs)
        summary["end_err_p95"] = end_errs_sorted[int(0.95 * (len(end_errs_sorted) - 1))]

    # per-query runtime summary
    wall_times = [r["wall_s"] for r in records]
    cpu_times = [r["cpu_s"] for r in records]
    if wall_times:
        summary["wall_per_query_mean"] = stats.mean(wall_times)
        summary["wall_per_query_median"] = stats.median(wall_times)
    if cpu_times:
        summary["cpu_per_query_mean"] = stats.mean(cpu_times)
        summary["cpu_per_query_median"] = stats.median(cpu_times)

    return records, summary

# ----------------------------
# Main
# ----------------------------

def main():
    seed_db = load_pickle("probabilistic_db.pkl")
    score_db = load_pickle("score_matrix_db.pkl")
    genome = load_genome("chr22_ancestor.fa")
    conf = load_conf("chr22_ancestor.conf")  # not used for metrics; kept for compatibility
    queries = load_queries("test_queries_random_mut_copy.fa")

    # Optional threshold for variant fallback
    max_prob = 0.0
    try:
        with open("max_prob_variant.txt", "r") as f:
            w = f.read().strip()
        max_prob = compute_max_prob_from_word(seed_db, w)
        print("max_prob_threshold:", max_prob, "(from", w, ")")
    except FileNotFoundError:
        print("max_prob_variant.txt not found -> max_prob_threshold = 0.0")

    records, summary = evaluate_queries(
        queries=queries,
        genome=genome,
        seed_db=seed_db,
        score_db=score_db,
        tol_list=(0,0),
        max_prob_threshold=max_prob,
        seeds_top_k=200,
        ungapped_top_k=10,
    )

    print("\n==== OVERALL SUMMARY ====")
    for k, v in summary.items():
        print(f"{k}: {v}")

if __name__ == "__main__":
    main()
