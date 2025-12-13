import pickle

with open("score_matrix_db.pkl", "rb") as f:
    score_db = pickle.load(f)

scores_i = score_db[0]          # O(1)
print(score_db[0])
score_A, score_C, score_G, score_T = scores_i

print("Scores at position 0: {'A': %.4f, 'C': %.4f, 'G': %.4f, 'T': %.4f}" % (score_A, score_C, score_G, score_T))