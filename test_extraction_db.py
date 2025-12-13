# test_extraction_db.py

import pickle
import time
import random
import string

DB_FILE = "probabilistic_db.pkl"
WORD_LEN = 11
N_TEST_LOOKUPS = 100_000

def load_database(path):
    # Database is loaded from disk into memory
    with open(path, "rb") as f:
        db = pickle.load(f)
    return db

def random_word():
    # Random query word is generated
    alphabet = ["A", "C", "G", "T", "_"]
    return "".join(random.choice(alphabet) for _ in range(WORD_LEN))

def validate_structure(db):
    # Basic structural properties are checked
    if not isinstance(db, dict):
        raise TypeError("Database is not a dictionary")

    sample_key = next(iter(db))
    sample_val = db[sample_key]

    if not isinstance(sample_key, str):
        raise TypeError("Keys are not strings")

    if not isinstance(sample_val, list):
        raise TypeError("Values are not lists")

    if sample_val:
        entry = sample_val[0]
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise TypeError("Entries are not (index, probability) tuples")

    print("Database structure validated")

def single_lookup(db, word):
    # Single dictionary lookup is performed
    return db.get(word, [])

def benchmark_lookups(db):
    words = [random_word() for _ in range(N_TEST_LOOKUPS)]

    start = time.perf_counter()
    for w in words:
        _ = db.get(w, [])
    elapsed = time.perf_counter() - start

    avg_time = elapsed / N_TEST_LOOKUPS

    print("Total lookups:", N_TEST_LOOKUPS)
    print("Total time: {:.4f} seconds".format(elapsed))
    print("Average lookup time: {:.2e} seconds".format(avg_time))

def main():
    print("Loading database...")
    db = load_database(DB_FILE)

    print("Database loaded")
    print("Total unique keys:", len(db))
    print()

    validate_structure(db)
    print()

    # Known-key lookup sanity check
    example_word = '_AAAAAAAAAA'
    hits = single_lookup(db, example_word)

    print("Sanity check lookup")
    print("Word:", example_word)
    print("Number of hits:", len(hits))
    print("First hit:", hits[0] if hits else None)
    print()

    # Performance benchmark
    print("Benchmarking lookup performance...")
    benchmark_lookups(db)

    print()
    print("Expected complexity: O(1) average per lookup")
    print("Test completed successfully")

if __name__ == "__main__":
    main()
