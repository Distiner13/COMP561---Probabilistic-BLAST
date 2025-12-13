# build_db.py

import pickle
import time
import os

WORD_LEN = 11

def load_fasta(path):
    # Whole file is read and whitespace is removed
    with open(path, "r") as f:
        text = f.read()

    lines = text.splitlines()

    # FASTA headers are ignored if present
    seq_parts = []
    for line in lines:
        if not line:
            continue
        if line.startswith(">"):
            continue
        seq_parts.append(line.strip())

    seq = "".join(seq_parts)
    seq = "".join(seq.split())

    return seq

def load_conf(path):
    # Whole file is read and tokens are parsed as floats
    with open(path, "r") as f:
        text = f.read().strip()

    if not text:
        return []

    tokens = text.split()
    probs = []

    for t in tokens:
        try:
            probs.append(float(t))
        except ValueError as e:
            raise ValueError(f"Probability parse failed near token: {t}") from e

    return probs

def estimate_resources(genome_length):
    windows = genome_length - WORD_LEN + 1
    entries = windows * (WORD_LEN + 1)

    bytes_per_entry = 130
    total_bytes = entries * bytes_per_entry
    total_gb = total_bytes / (1024 ** 3)

    ops_per_window = 20
    est_time_sec = windows * ops_per_window * 2e-7

    print("Genome length:", genome_length)
    print("Word length:", WORD_LEN)
    print("Total windows:", windows)
    print("Total stored entries:", entries)
    print("Estimated memory usage: {:.2f} GB".format(total_gb))
    print("Estimated build time: {:.1f} seconds".format(est_time_sec))
    print()

def build_database(seq, conf):
    db = {}
    n = len(seq)

    for i in range(n - WORD_LEN + 1):
        word = seq[i:i + WORD_LEN]
        base_prob = 1.0
        probs = []

        # Base word probability is computed from dominant probabilities
        for j in range(WORD_LEN):
            p = conf[i + j]
            base_prob *= p
            probs.append(p)

        # Base word entry is stored
        if word not in db:
            db[word] = []
        db[word].append((i, base_prob))

        # One "_" variant per position is stored
        for j in range(WORD_LEN):
            p = probs[j]
            alt_p = (1.0 - p) / 3 
            var_prob = (base_prob/p) * alt_p 
            var_word = word[:j] + "_" + word[j + 1:]

            if var_word not in db:
                db[var_word] = []
            db[var_word].append((i, var_prob))

    return db

def main():
    output_file = "probabilistic_db.pkl"

    if os.path.exists(output_file):
        print("Database already exists:", output_file)
        print("Delete it manually if you want to rebuild.")
        return

    seq = load_fasta("chr22_ancestor.fa")
    conf = load_conf("chr22_ancestor.conf")

    # Input lengths are validated before build
    if len(conf) != len(seq):
        raise ValueError(
            "Length mismatch: conf has {} probabilities but genome has {} bases".format(
                len(conf), len(seq)
            )
        )

    estimate_resources(len(seq))

    start = time.time()
    db = build_database(seq, conf)
    elapsed = time.time() - start

    # Database is written as a real file in the current directory
    with open(output_file, "wb") as f:
        pickle.dump(db, f, protocol=pickle.HIGHEST_PROTOCOL)

    print("Database written to:", output_file)
    print("Actual build time: {:.1f} seconds".format(elapsed))
    print("Build complete. Do not rerun this script.")

if __name__ == "__main__":
    main()
