import random

GENOME_FILE = "chr22_ancestor.fa"
OUTPUT_FILE = "test_queries_generated.fa"

QUERY_LENGTH = 50
N_MUTATIONS = 2
N_QUERIES = 10
RAND_ADDITIONAL = 30

BASES = ["A", "C", "G", "T"]

def load_genome(path):
    parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    return "".join(parts)

def mutate_sequence(seq, n_mut):
    seq = list(seq)
    length = len(seq)

    positions = random.sample(range(length), n_mut)

    for pos in positions:
        original = seq[pos]
        choices = [b for b in BASES if b != original]
        seq[pos] = random.choice(choices)

    return "".join(seq)

def random_flank():
    length = random.randint(1, RAND_ADDITIONAL)
    return "".join(random.choice(BASES) for _ in range(length)), length

def main():
    genome = load_genome(GENOME_FILE)
    genome_len = len(genome)

    with open(OUTPUT_FILE, "w") as out:
        for i in range(N_QUERIES):
            start = random.randint(0, genome_len - QUERY_LENGTH)
            end = start + QUERY_LENGTH

            core = genome[start:end]
            mutated = mutate_sequence(core, N_MUTATIONS)

            prefix, pre_len = random_flank()
            suffix, suf_len = random_flank()

            final_seq = prefix + mutated + suffix

            header = (
                f">query{i+1}"
                f"|db_start={start}"
                f"|db_end={end}"
                f"|prefix_len={pre_len}"
                f"|suffix_len={suf_len}"
            )

            out.write(header + "\n")
            out.write(final_seq + "\n")

    print("Generated", N_QUERIES, "query sequences")
    print("Output written to:", OUTPUT_FILE)

if __name__ == "__main__":
    main()
