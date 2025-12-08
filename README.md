# COMP561---Probabilistic-BLAST

Create your own file to do the assigned tasks. We will merge the project by calling on the different files.

Use Python.
Tasks:
TASK1 Anis
- Read files and extract data
- Finding the most likely sequence, a list of the most probable sequences from most likely to least likely
→ O(L!) memory usage, so instead use top 10 or so most likely databases
- Find the maximum likelihood word
→ O(L) memory usage

TASK2 Yilin
- NW algo, make the gapped and ungapped extension algorithms (with the ability to customize all inputs) (score tracking included)
- Further extension will use a score matrix which takes in the probabilities into consideration
→ 4 x 4 x L → 16L → O(L) memory usage

TASK3 Xiaoyi
- Build database
- Make the BLAST algorithm with the ability to customize all inputs. Given database and query sequence, match them using the words
- Score tracking and ungapped extension, and then gapped. Add the ability to modify a single nucleotide at that position and calculate the score again. Set thresholds and optimization 

- Create a list of queries that we can use as examples

Xiaoyi Comments 12.8:
I added 4 method: generate_query, write_queries_fasta, write_kmer_index_fasta, sw_gapped
Basically, for the first 2 methods, they can generate customized query and output them in a .fa fiile correspondingly.
The third one use Yilin's method build_kmer_index, which I found is already very well-written: build_kmer_index can already genrate a dict of indexed k-mer words
write_kmer_index_fasta is based on the method to output a .fa file, which is more like a database.
I provided 2 examples of my .fa outputs for each methods: queries.fa and wordsdb.fa
As to sw_gapped, I am not confident with the codes since it is not tested, since there is no final proper scoring matrix yet. 
This SW does a gapped alignment based on the db range given by the ungapped extension method.
Another thing I have to mention is that I used AI helping in coding, especially in SW.

* I committed another data_process_xiaoyi.py, which keeps everything Yilin had written and my part. Plz download data_process_xiaoyi.py instead of data_process.py
