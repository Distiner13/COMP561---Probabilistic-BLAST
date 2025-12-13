# BLAST.py

import pickle
import time

WORD_LEN = 11
GAP_PENALTY = -2.0
UNGAPPED_X_DROP = 5.0
MAX_SEEDS_TO_EXTEND = 200
GAPPED_WINDOW_PAD = 40

BASE_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}

def load_pickle(path):
    print("Loading", path)
    t0 = time.time()
    with open(path, "rb") as f:
        obj = pickle.load(f)
    print("Loaded in {:.2f} s".format(time.time() - t0))
    return obj

def load_genome(path):
    print("Loading genome")
    parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    genome = "".join(parts)
    print("Genome length:", len(genome))
    return genome

def load_queries(path):
    print("Loading queries")
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

def get_score(score_db, g_pos, q_base):
    return score_db[g_pos][BASE_TO_INDEX[q_base]]

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

def ungapped_extend(query, genome, score_db, q0, g0):
    k = WORD_LEN

    score = 0.0
    for t in range(k):
        score += get_score(score_db, g0 + t, query[q0 + t])

    best = score
    best_lq = q0
    best_lg = g0
    best_rq = q0 + k
    best_rg = g0 + k

    run = score
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

    run = best
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
        "seed_g": g0,
    }

def needleman_wunsch(query_seg, genome_seg, genome_start, score_db):
    n = len(query_seg)
    m = len(genome_seg)

    dp = [[0.0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP_PENALTY
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + GAP_PENALTY

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + get_score(
                score_db, genome_start + j - 1, query_seg[i - 1]
            )
            up = dp[i - 1][j] + GAP_PENALTY
            left = dp[i][j - 1] + GAP_PENALTY
            dp[i][j] = max(diag, up, left)

    return dp[n][m]

def main():
    seed_db = load_pickle("probabilistic_db.pkl")
    score_db = load_pickle("score_matrix_db.pkl")
    genome = load_genome("chr22_ancestor.fa")
    queries = load_queries("test_queries_random_mut.fa")

    for qi, (qid, query) in enumerate(queries, 1):
        print("\n==============================")
        print("Processing query", qi, "/", len(queries))
        print("Query ID:", qid)
        print("Query length:", len(query))

        t0 = time.time()
        seeds = find_seeds(query, seed_db)
        print("Seeds found:", len(seeds))

        if not seeds:
            print("No perfect seeds found")
            print("Fallback logic placeholder")
            continue

        ungapped_results = []
        for q_pos, g_pos in seeds[:MAX_SEEDS_TO_EXTEND]:
            ungapped_results.append(
                ungapped_extend(query, genome, score_db, q_pos, g_pos)
            )

        ungapped_results.sort(key=lambda x: x["score"], reverse=True)
        best = ungapped_results[0]

        print("Best ungapped score:", round(best["score"], 4))
        print("Ungapped genome range:", best["g_start"], best["g_end"])
        print("Ungapped query range:", best["q_start"], best["q_end"])

        # Incremental gapped extension around ungapped match
        q_start = best["q_start"]
        q_end = best["q_end"]
        g_start = best["g_start"]
        g_end = best["g_end"]

        win_g_start = max(0, g_start - GAPPED_WINDOW_PAD)
        win_g_end = min(len(genome), g_end + GAPPED_WINDOW_PAD)

        genome_seg = genome[win_g_start:win_g_end]
        query_seg = query[q_start:q_end]

        print("Running gapped extension")
        t1 = time.time()
        gapped_score = needleman_wunsch(
            query_seg, genome_seg, win_g_start, score_db
        )
        print("Gapped score:", round(gapped_score, 4))
        print("Gapped extension time: {:.2f} s".format(time.time() - t1))
        print("Total query time: {:.2f} s".format(time.time() - t0))

if __name__ == "__main__":
    main()
