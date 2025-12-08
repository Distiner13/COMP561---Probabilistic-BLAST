from typing import List, Dict, Optional, Tuple
import math

def read_fasta_sequence(path: str) -> str:
    """
    Reads a FASTA file and returns the concatenated sequence (single contiguous string),
    ignoring all header lines and line breaks. Assumes only one sequence.
    """
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Skip FASTA header
                continue
            seq.append(line.upper())
    return ''.join(seq)


def read_confidences(path: str) -> List[float]:
    """
    Reads a .conf file and returns a list of float values.
    Works whether the file has one value per line or several values per line.
    """
    conf = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            parts = line.split()
            for p in parts:
                conf.append(float(p))
    return conf


def build_score_matrix(
    major_base: str,
    p_max: List[float],
    bases: str = "ACGT",
    bg: Optional[Dict[str, float]] = None,
) -> List[Dict[str, float]]:
    """
    Builds a log-odds score matrix for the probabilistic genome.

    major_base : the major base string (e.g., 'ACGT...'), length L
    p_max      : list of length L, where p_max[i] = P(major_base[i] is correct)
    bases      : alphabet (default = "ACGT")
    bg         : background frequencies (default uniform)

    Returns:
        score[i][b] = log( P(b at position i) / bg[b] )
    """
    if bg is None:
        bg = {b: 1.0 / len(bases) for b in bases}

    if len(major_base) != len(p_max):
        raise ValueError(
            f"Length mismatch: major_base={len(major_base)}, p_max={len(p_max)}"
        )

    log_bg = {b: math.log(bg[b]) for b in bases}
    L = len(major_base)

    score: List[Dict[str, float]] = [{b: 0.0 for b in bases} for _ in range(L)]

    for i in range(L):
        b = major_base[i]
        c = p_max[i]
        for x in bases:
            if x == b:
                p = c
            else:
                p = (1.0 - c) / 3.0
            p = max(p, 1e-12)  # Avoid log(0)
            score[i][x] = math.log(p) - log_bg[x]

    return score


def build_kmer_index(major_base: str, k: int) -> Dict[str, List[int]]:
    """
    Builds a k-mer index from the major_base sequence.
    Maps each k-mer string to the list of starting positions where it occurs.
    """
    index: Dict[str, List[int]] = {}
    L = len(major_base)
    for i in range(L - k + 1):
        word = major_base[i:i+k]
        if word not in index:
            index[word] = []
        index[word].append(i)
    return index


def find_probabilistic_seeds(
    query: str,
    k: int,
    index: Dict[str, List[int]],
    score: List[Dict[str, float]],
    seed_threshold: float,
):
    """
    Optimized probabilistic seed finding:

    1. Use k-mer perfect matches to find candidate positions in the database.
    2. For each candidate (db_pos, q_pos), compute the k-mer score using the score matrix.
    3. Keep only seeds whose score >= seed_threshold.

    Returns:
        seeds: list of (db_pos, q_pos, kmer_score)
    """
    seeds = []
    M = len(query)
    L = len(score)

    for q_pos in range(M - k + 1):
        word = query[q_pos:q_pos + k]
        if word not in index:
            continue  # No such k-mer in DB

        # For all perfect matches in DB
        for db_pos in index[word]:
            if db_pos + k > L:
                continue  # Safety check

            # Compute probabilistic k-mer score
            kmer_score = 0.0
            for t in range(k):
                base = query[q_pos + t]
                kmer_score += score[db_pos + t][base]

            if kmer_score >= seed_threshold:
                seeds.append((db_pos, q_pos, kmer_score))

    return seeds


def _extend_one_direction(
    score: List[Dict[str, float]],
    query: str,
    db_pos: int,
    q_pos: int,
    step_db: int,
    step_q: int,
    drop_off: float,
) -> Tuple[float, int]:
    """
    Ungapped extension in one direction (left or right).
    Uses the BLAST-style drop-off rule:
        stop when (best_score - current_score) > drop_off.

    Returns:
        best_increment: best score gained in this direction
        best_len     : number of steps that achieved best_increment
    """
    len_db = len(score)
    len_q = len(query)

    cur = 0.0
    best = 0.0
    best_len = 0
    steps = 0

    i_db = db_pos
    i_q = q_pos

    while 0 <= i_db < len_db and 0 <= i_q < len_q:

        base = query[i_q]
        cur += score[i_db][base]
        steps += 1

        # Update best
        if cur > best:
            best = cur
            best_len = steps

        # Drop-off condition
        if best - cur > drop_off:
            break

        i_db += step_db
        i_q += step_q

    return best, best_len


def ungapped_extension(
    score: List[Dict[str, float]],
    query: str,
    db_start: int,
    q_start: int,
    k: int,
    drop_off: float = 20.0,
):
    """
    Performs ungapped extension around a perfect-match seed.
    Produces an HSP (High-scoring Segment Pair).

    Returns an HSP dictionary containing:
        score        – total HSP score
        db_start/db_end – DB range (inclusive)
        q_start/q_end   – Query range (inclusive)
        seed_db, seed_q – seed starting positions
        seed_len        – k
    """
    len_db = len(score)
    len_q = len(query)

    # Seed score
    seed_score = 0.0
    for t in range(k):
        pos_db = db_start + t
        pos_q = q_start + t
        if pos_db < 0 or pos_db >= len_db or pos_q < 0 or pos_q >= len_q:
            break
        base = query[pos_q]
        seed_score += score[pos_db][base]

    # Extend left
    if db_start > 0 and q_start > 0:
        left_score, left_len = _extend_one_direction(
            score, query,
            db_start - 1, q_start - 1,
            -1, -1,
            drop_off
        )
    else:
        left_score, left_len = 0.0, 0

    # Extend right
    if db_start + k < len_db and q_start + k < len_q:
        right_score, right_len = _extend_one_direction(
            score, query,
            db_start + k, q_start + k,
            +1, +1,
            drop_off
        )
    else:
        right_score, right_len = 0.0, 0

    total_score = seed_score + left_score + right_score

    # Compute HSP coordinates (inclusive)
    hsp_db_start = db_start - left_len
    hsp_db_end   = db_start + k - 1 + right_len
    hsp_q_start  = q_start - left_len
    hsp_q_end    = q_start + k - 1 + right_len

    hsp = {
        "score": total_score,
        "db_start": hsp_db_start,
        "db_end": hsp_db_end,
        "q_start": hsp_q_start,
        "q_end": hsp_q_end,
        "seed_db": db_start,
        "seed_q": q_start,
        "seed_len": k,
    }
    return hsp




def generate_query(
    major_base: str,
    p_max: List[float],
    length: int,
    start: Optional[int] = None,
    mut_rate: float = 0.0,
) -> Tuple[str, int]:
    """
    generate a query from probabilistic genome.
    major_base: ancestor seq
    p_max: read_confidences
    length: query length
    start: query start index
    mut_rate: mutation rate
    ---return 
    (query, start)
    """
    L = len(major_base) 
    if length <= 0:
        raise ValueError("length has to be positive")

    if start is None:
        max_start = L - length
        if max_start < 0:
            raise ValueError()
        start = random.randrange(0, max_start + 1)
    else:
        if start < 0 or start + length > L:
            raise ValueError()

    bases = "ACGT"
    base_to_idx = {b: i for i, b in enumerate(bases)}
    query_chars = []

    # generate the query according to conf prob
    for offset in range(length):
        pos = start + offset
        maj = major_base[pos]
        c = p_max[pos]  # maj prob

        # A/C/G/T probs at the pos
        probs = [(1.0 - c)/3.0]*4
        probs[base_to_idx[maj]] = c

        base = random.choices(bases, weights=probs, k=1)[0]
        query_chars.append(base)

    # maybe we want mutations
    if mut_rate > 0.0:
        for i in range(length):
            if random.random() < mut_rate:
                current = query_chars[i]
                alternatives = [b for b in bases if b != current]
                query_chars[i] = random.choice(alternatives)

    query = "".join(query_chars)
    return query, start

import random
from typing import List



def write_queries_fasta(
    out_path: str,
    major_base: str,
    p_max: List[float],
    num_queries: int,
    length: int,
    mut_rate: float = 0.0,
    line_width: int = 80, # upper limits of length of a row in fasta
) -> None:
    """
    generate a fa file with customized number of queries 
    """

    with open(out_path, "w") as fout:
        for i in range(num_queries):
            query, start = generate_query(
                major_base=major_base,
                p_max=p_max,
                length=length,
                start=None,
                mut_rate=mut_rate,
            )
            db_start = start
            db_end = start + len(query)-1

            header = f">query{i+1}|db_start={db_start}|db_end={db_end}\n"
            fout.write(header)
            for j in range(0, len(query), line_width):
                fout.write(query[j:j+line_width] + "\n")


def write_kmer_index_fasta(
    major_base: str,
    k: int,
    out_path: str,
    line_width: int = 80,
) -> None:
    """
    use build_kmer_index to build k-mer indexed database.
    output in .fa

    i.e.:
        >ACGT|count=12
        10 283 2949 ...
    """

    index = build_kmer_index(major_base, k)
    kmers_sorted = sorted(index.keys())

    with open(out_path, "w") as fout:
        for kmer in kmers_sorted:
            positions = index[kmer]
            count = len(positions)

            fout.write(f">{kmer}|count={count}\n")

            pos_str = " ".join(str(p) for p in positions)
            for i in range(0, len(pos_str), line_width):
                fout.write(pos_str[i:i+line_width] + "\n")






































def main():

    conf_path = '/Users/xiaoyixu/Downloads/COMP561---Probabilistic-BLAST-main/chr22_ancestor.conf'

    fasta_path = "/Users/xiaoyixu/Downloads/COMP561---Probabilistic-BLAST-main/chr22_ancestor.fa"
    p_max = read_confidences(conf_path)
    sequence = read_fasta_sequence(fasta_path)
    score = build_score_matrix(sequence, p_max)

    query, start = generate_query(sequence, p_max, length = 6, mut_rate=0.005)
    out_path_q = "/Users/xiaoyixu/Downloads/COMP561---Probabilistic-BLAST-main/queries.fa"
    out_path_w = "/Users/xiaoyixu/Downloads/COMP561---Probabilistic-BLAST-main/wordsdb.fa"

    # write_queries_fasta(
    #     out_path=out_path_q,
    #     major_base=sequence,
    #     p_max=p_max,
    #     num_queries=10,
    #     length=200,
    #     mut_rate=0.005,
    # )
    write_kmer_index_fasta(sequence, 11, out_path_w)





    # query = "CTTTCTGACTCCTTACGCTGTCCACTCATTCAGTTGATAAAAAGATAGAAACCCAGAATCTGAAGTCTCCTTTTGTCCCCAACATGCCTGATCACGAGGGCAGAACAGGACAACAGAAATGGCTCACACTTGACAGGCCATCTCTCCATGCCTATATTCATAAAGTATGCAGCCAGTTTTTCTTGTGACACACTTGATGTCCCTTGTAACAATGGGAAAATGTGAAGGACCGTGAGAAATCTGGCAGGCTGGCAAATATTTCCACAATTCAGTTCTTTTTTGATCTGTGAGGCTAATCTAATGATGGTGATAGGGCAAGGGAGTTAATGGCCCCTGAAAGGCCTACAGAAAACCCTCAGCCACGGGGTTAACAGACTGCCACAGGAGACTAAAGGAGAGACTGAGGCCCAGTGGGGAAACCAGCCAATGGGGAAGCAACAGGATCACTGGGGAGGACAGATCTGGGAAGGGAGAACACTCCAGAGCAGAGCTCCCTGAAAAGATCATTAATTCACACAATAGCTATTGAGCACCTACTATGTGCCTGACATTGTTCTAGGCACTGGGGGAAAGCATGGTGAGCAGGACGGGTAATTACGGTAGTCAGAGAATAATTATCTCAGAGTTGTAAGAGACTAAAAGGTCACTATTCGACCCCCTATTATTCACTCATATAACAAACATGTATTGAGAACTTACTATGTGCCGGGCCTTATGCTAGGCACTGGGGGTACGATGTGGACAAAACAGATGAGATTCTTGCCCTCATGGAGCTTACAGTCTAGCGGGGAAGGCAGACAAAAAAATCCCATAAACCAATGTATAATTACGAGCTGCGATAAGGGTCTTGCTGGAAAAGTACAGAATGCCATCAGTACATATGACAGGGGACCTTCCCAGTCTGAGGGGTCAGGGCCAAGACTTGCCTGAGAAAGTGACACTGAGCTCAAATCTGAAGGATGGAGAGGCATAACCAGGTGAAGGATCAGGGAGGTGTACTGATGTCAGAGTACAGGTTTAAATGTCCAAGGTGTGTTCAAGGAATGAGAAGAGGGCGGGGGTGGGGAAGGGCACA"
    # # first three bases manually changed from ATT → CTT

    # k = 11
    # index = build_kmer_index(sequence, k)
    # seed_threshold = 5.0

    # seeds = find_probabilistic_seeds(query, k, index, score, seed_threshold)
    # print("found seeds counts:", len(seeds))
    # for db_pos, q_pos, sc in seeds[:10]:
    #     print(f"DB:{db_pos}, Q:{q_pos}, score={sc:.3f}")

    # drop_off = 15.0
    # best_hsp = None
    # best_score = -1e18

    # for db_pos, q_pos, seed_sc in seeds:
    #     hsp = ungapped_extension(score, query, db_pos, q_pos, k, drop_off)
    #     if hsp["score"] > best_score:
    #         best_score = hsp["score"]
    #         best_hsp = hsp

    # if best_hsp is not None:
    #     print("Best HSP:")
    #     print("  score   :", best_hsp["score"])
    #     print("  DB range:", best_hsp["db_start"], "->", best_hsp["db_end"])
    #     print("  Q range :", best_hsp["q_start"], "->", best_hsp["q_end"])
    # else:
    #     print("No HSP found (no seeds passed the threshold).")


if __name__ == "__main__":
    main()
