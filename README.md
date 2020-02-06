# Bioinformatics

**_Python programs_** used to solve classical **_Bioinformatics_** problems that involve analyzing DNA sequences given as Fasta Files. Code is programmed using **_Object Oriented_** style for easy use elswhere. **_Dynammic Programming Algorithms_** are extensively used, which are **_algorithms that essentially divide a large problem into a series of smaller problems, and it uses the solutions to the smaller problems to find an optimal solution to the larger problem_**

## _Gene-Finding.py_ : Application of Hidden Markov Models

Programmed a **_Class_** that has functionalities that help find areas of **_Genetic coding_** given a sequence of Nucleotides(**_DNA sequence in Fasta File_**). This Class first needs to be trained by the provided **_class methods_** using **_gff3 files and fasta files_**, which the Class contains methods to process the data. Once data is processed use the **_Hidden Markov Model_** to train the model and call the **_Viterbi Algorithm_** which is a **_Dynamic Programming Algorithm_** for any given DNA Sequence to show **_predicted areas_** of Genetic Coding



## _Needleman-Wunsch.py & NW_AffineGap.py_

Programmed a **_Class_** that has functionalities that help **_find the optimal aligned protein or nucleotide sequences given two sequences_**. The class can take and process **_fasta files_** or sequences can be hand typed. I used **_Dynamic Programming_** to compare biological sequences also known as the **_optimal matching algorithm and the global alignment technique_**, All parts of the algorithms are found in the class as methods including **_filling matrix table, backtracking, etc etc**

Difference between Needleman-Wunsch.py & NW_AffineGap.py include methods of **_penalty_** for gaps that show up in finding optimal sequence alignment

  
