# COMP561---Probabilistic-BLAST

1. Run build_db.py and build_score_matrix_db.py. After running them, you should be able to see probabilistic_db.pkl and score_matrix_db.pkl generated.
2. If you do not have prepared queries, please run generate_samples.py. It should generate a file test_queries_generated.fa. You can modify the parameters of the generated sequences to have a specific length, a number of mutations, a prefix and a suffix of random lengths. 
3. If you are using a query file with a name other than test_queries_generated.fa, please either rename your file to test_queries_generated.fa, or update line 414 in BLAST.py as follows: queries = load_queries("YOUR_QUERY_FILE_NAME").
4. Run BLAST.py, and should can see the results.
