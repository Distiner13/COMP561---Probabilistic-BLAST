# BLAST.py

import pickle
import time
import matplotlib.pyplot as plt

WORD_LEN = 11
GAP_PENALTY = -2.0
UNGAPPED_X_DROP = 10.0
MAX_SEEDS_TO_EXTEND = 200
GAPPED_WINDOW_PAD = 40

BASE_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


# -------------------- I/O --------------------

def load_pickle(path):
    print("Loading", path)
    t0 = time.time()
    with open(path, "rb") as f:
        obj = pickle.load(f)
    print("Loaded in {:.2f} s".format(time.time() - t0))
    return obj


def load_genome(path):
    parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    genome = "".join(parts)
    print("Genome length:", len(genome))
    return genome


def load_conf(path):
    with open(path, "r") as f:
        text = f.read().strip()
    return [float(x) for x in text.split()] if text else []


def load_queries(path):
    queries = []
    qid = None
    buf = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if qid is not None:
                    queries.append((qid, "".join(buf)))
                qid = line[1:]
                buf = []
            else:
                buf.append(line)

    if qid is not None:
        queries.append((qid, "".join(buf)))

    print("Number of queries:", len(queries))
    return queries


# -------------------- Scoring --------------------

def get_score(score_db, g_pos, q_base):
    return score_db[g_pos][BASE_TO_INDEX[q_base]]


# -------------------- Seed handling --------------------

def split_into_words(seq):
    return [(i, seq[i:i + WORD_LEN]) for i in range(len(seq) - WORD_LEN + 1)]


def find_seeds(query, seed_db):
    seeds = []
    for q_pos, word in split_into_words(query):
        hits = seed_db.get(word)
        if hits:
            for g_pos, _ in hits:
                seeds.append((q_pos, g_pos))
    return seeds


# -------------------- Ungapped extension --------------------

def ungapped_extend(query, genome, score_db, q0, g0):
    k = WORD_LEN

    seed_score = sum(
        get_score(score_db, g0 + t, query[q0 + t]) for t in range(k)
    )

    best = seed_score
    best_lq = q0
    best_lg = g0
    best_rq = q0 + k
    best_rg = g0 + k

    run = seed_score
    q = q0 - 1
    g = g0 - 1

    while q >= 0 and g >= 0:
        run += get_score(score_db, g, query[q])
        if run > best:
            best = run
            best_lq = q
            best_lg = g
        if best - run > UNGAPPED_X_DROP:
            break
        q -= 1
        g -= 1

    run = seed_score
    q = q0 + k
    g = g0 + k

    while q < len(query) and g < len(genome):
        run += get_score(score_db, g, query[q])
        if run > best:
            best = run
            best_rq = q + 1
            best_rg = g + 1
        if best - run > UNGAPPED_X_DROP:
            break
        q += 1
        g += 1

    return {
        "q_start": best_lq,
        "q_end": best_rq,
        "g_start": best_lg,
        "g_end": best_rg,
        "score": best,
        "seed_q": q0,
        "seed_g": g0,
        "seed_score": seed_score,
    }


def ungapped_extend_with_trace_coords(query, genome, score_db, q0, g0):
    k = WORD_LEN
    seed_score = sum(
        get_score(score_db, g0 + t, query[q0 + t]) for t in range(k)
    )

    left = []
    run = seed_score
    best = seed_score
    q = q0 - 1
    g = g0 - 1

    while q >= 0 and g >= 0:
        run += get_score(score_db, g, query[q])
        left.append((g, run))
        best = max(best, run)
        if best - run > UNGAPPED_X_DROP:
            break
        q -= 1
        g -= 1

    right = []
    run = seed_score
    best = seed_score
    q = q0 + k
    g = g0 + k

    while q < len(query) and g < len(genome):
        run += get_score(score_db, g, query[q])
        right.append((g, run))
        best = max(best, run)
        if best - run > UNGAPPED_X_DROP:
            break
        q += 1
        g += 1

    return (g0, seed_score), left, right


# -------------------- Gapped alignment --------------------

def needleman_wunsch_traceback(query_seg, genome_seg, genome_start, score_db):
    n = len(query_seg)
    m = len(genome_seg)

    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    bt = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_PENALTY
        bt[i][0] = "U"

    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + GAP_PENALTY
        bt[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + get_score(
                score_db, genome_start + j - 1, query_seg[i - 1]
            )
            up = dp[i - 1][j] + GAP_PENALTY
            left = dp[i][j - 1] + GAP_PENALTY
            dp[i][j] = max(diag, up, left)
            bt[i][j] = "D" if dp[i][j] == diag else ("U" if dp[i][j] == up else "L")

    i, j = n, m
    aq, ag, gp = [], [], []
    gpos = genome_start + m - 1

    while i > 0 or j > 0:
        if bt[i][j] == "D":
            aq.append(query_seg[i - 1])
            ag.append(genome_seg[j - 1])
            gp.append(gpos)
            i -= 1
            j -= 1
            gpos -= 1
        elif bt[i][j] == "U":
            aq.append(query_seg[i - 1])
            ag.append("-")
            gp.append(None)
            i -= 1
        else:
            aq.append("-")
            ag.append(genome_seg[j - 1])
            gp.append(gpos)
            j -= 1
            gpos -= 1

    return {
        "score": dp[n][m],
        "aligned_query": "".join(reversed(aq)),
        "aligned_genome": "".join(reversed(ag)),
        "aligned_gpos": list(reversed(gp)),
    }


# -------------------- Alignment assembly --------------------

def make_ungapped_alignment(query, genome, q0, g0, length):
    return {
        "score": None,
        "aligned_query": query[q0:q0 + length],
        "aligned_genome": genome[g0:g0 + length],
        "aligned_gpos": list(range(g0, g0 + length)),
    }


def score_ungapped(query, score_db, q0, g0, length):
    return sum(get_score(score_db, g0 + i, query[q0 + i]) for i in range(length))


def gapped_extend_flanks(query, genome, score_db, best):
    q0 = best["q_start"]
    g0 = best["g_start"]
    length = min(best["q_end"] - q0, best["g_end"] - g0)

    core_score = score_ungapped(query, score_db, q0, g0, length)
    core = make_ungapped_alignment(query, genome, q0, g0, length)

    left_q = query[:q0]
    left_g = genome[max(0, g0 - len(left_q) - GAPPED_WINDOW_PAD):g0]
    right_q = query[q0 + length:]
    right_g = genome[g0 + length:g0 + length + len(right_q) + GAPPED_WINDOW_PAD]

    left = needleman_wunsch_traceback(
        left_q, left_g, g0 - len(left_g), score_db
    ) if left_q and left_g else {"score": 0, "aligned_query": "", "aligned_genome": "", "aligned_gpos": []}

    right = needleman_wunsch_traceback(
        right_q, right_g, g0 + length, score_db
    ) if right_q and right_g else {"score": 0, "aligned_query": "", "aligned_genome": "", "aligned_gpos": []}

    total = core_score + left["score"] + right["score"]

    full = {
        "score": total,
        "aligned_query": left["aligned_query"] + core["aligned_query"] + right["aligned_query"],
        "aligned_genome": left["aligned_genome"] + core["aligned_genome"] + right["aligned_genome"],
        "aligned_gpos": left["aligned_gpos"] + core["aligned_gpos"] + right["aligned_gpos"],
    }

    core["score"] = core_score
    return core, full


# -------------------- Visualization --------------------

def plot_alignment_score_trace(seed, left, right):
    plt.figure(figsize=(12, 4))

    if left:
        plt.plot([p[0] for p in left], [p[1] for p in left], label="Ungapped left")
    if right:
        plt.plot([p[0] for p in right], [p[1] for p in right], label="Ungapped right")

    plt.scatter([seed[0]], [seed[1]], color="red", label="Seed")
    plt.xlabel("Genome index")
    plt.ylabel("Cumulative alignment score")
    plt.title("Ungapped BLAST extension")
    plt.legend()
    plt.tight_layout()
    plt.show()


def print_alignment_visual(aln, conf, score_db, block=20):
    q = aln["aligned_query"]
    g = aln["aligned_genome"]
    gp = aln["aligned_gpos"]

    for i in range(0, len(q), block):
        j = min(len(q), i + block)
        print("P_max |", " ".join(f"{conf[p]:5.2f}" if p is not None else "     " for p in gp[i:j]))
        print("DB    |", " ".join(f"  {c}  " for c in g[i:j]))
        print("Query |", " ".join(f"  {c}  " for c in q[i:j]))
        print("Score |", " ".join(
            f"{get_score(score_db, p, qc):5.2f}" if p is not None and qc != "-" else f"{GAP_PENALTY:5.2f}"
            for qc, p in zip(q[i:j], gp[i:j])
        ))
        print()


# -------------------- Main --------------------

def main():
    seed_db = load_pickle("probabilistic_db.pkl")
    score_db = load_pickle("score_matrix_db.pkl")
    genome = load_genome("chr22_ancestor.fa")
    conf = load_conf("chr22_ancestor.conf")
    queries = load_queries("test_queries_generated.fa")

    for qid, query in queries:
        print("\nProcessing:", qid)

        seeds = find_seeds(query, seed_db)
        if not seeds:
            print("No seeds found")
            continue

        results = [ungapped_extend(query, genome, score_db, q, g)
                   for q, g in seeds[:MAX_SEEDS_TO_EXTEND]]
        best = max(results, key=lambda x: x["score"])

        print("Best ungapped score:", round(best["score"], 4))

        seed, left, right = ungapped_extend_with_trace_coords(
            query, genome, score_db, best["seed_q"], best["seed_g"]
        )

        plot_alignment_score_trace(seed, left, right)

        core, full = gapped_extend_flanks(query, genome, score_db, best)
        final = full if full["score"] > core["score"] else core

        print("Final alignment score:", round(final["score"], 4))
        print_alignment_visual(final, conf, score_db)


if __name__ == "__main__":
    main()
