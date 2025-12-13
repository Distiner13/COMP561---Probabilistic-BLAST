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

seq = load_fasta("chr22_ancestor.fa")
conf = load_conf("chr22_ancestor.conf")

print(seq[425030:(425030+11)])