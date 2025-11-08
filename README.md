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

