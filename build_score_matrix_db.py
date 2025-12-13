# build_score_matrix_db.py

import pickle
import time
import os

BASE_ORDER = ("A", "C", "G", "T")

def load_conf(path):
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
    bytes_per_entry = 4 * 8  # four floats
    total_bytes = genome_length * bytes_per_entry
    total_mb = total_bytes / (1024 ** 2)

    est_time_sec = genome_length * 5e-8

    print("Genome length:", genome_length)
    print("Score values per position: 4")
    print("Estimated memory usage: {:.2f} MB".format(total_mb))
    print("Estimated build time: {:.2f} seconds".format(est_time_sec))
    print()

def build_score_database(conf):
    n = len(conf)
    score_db = [None] * n

    for i in range(n):
        p_max = conf[i]
        p_other = (1.0 - p_max) / 3.0

        score_dom = 1.0
        score_other = 1.0 - ((p_max - p_other) * 3.0)

        # Scores are stored in fixed (A, C, G, T) order
        score_db[i] = (
            score_dom,
            score_other,
            score_other,
            score_other,
        )

    return score_db

def main():
    output_file = "score_matrix_db.pkl"

    if os.path.exists(output_file):
        print("Score matrix database already exists:", output_file)
        print("Delete it manually if you want to rebuild.")
        return

    conf = load_conf("chr22_ancestor.conf")

    estimate_resources(len(conf))

    start = time.time()
    score_db = build_score_database(conf)
    elapsed = time.time() - start

    with open(output_file, "wb") as f:
        pickle.dump(score_db, f, protocol=pickle.HIGHEST_PROTOCOL)

    print("Score matrix database written to:", output_file)
    print("Actual build time: {:.2f} seconds".format(elapsed))
    print("Build complete. Do not rerun this script.")

if __name__ == "__main__":
    main()
