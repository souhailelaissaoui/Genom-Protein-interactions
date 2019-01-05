# genom-2 

This group project was developed by students from Ecole Polytechnique fédérale de Lausanne

RUNNING THE PROGRAM
(from your genom-2 folder)


STEP 1: Make a build folder and build the program
 
 
--> mkdir build 

--> cd build

--> cmake ../

--> make

documentation		--> make doc
tests 			--> make test


STEP 2: Execute the program 


with Qt  			--> ./genom-2

with terminal 		--> ./Main

from the command line --> ./Main --help

tests 				--> ./gtests

==============================================================================
    HOW THE PROGRAM WORKS: 
==============================================================================

 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PROCEDURE: GET MATRIX FROM SEQUENCES
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
If the user wants to calculate a probability weight matrix based on a list of 
sequences, this is done in the following fashion:
If all sequences are of the same length, the amount that each nucleotide is at 
the specific position is counted in a matrix. In order to avoid a probability 
of zero (which will result in problems in later analysis), the user can choose 
to use nucleotide base probabilities to create a background noise. The 
probability of each nucleotide is added to the corrsponding matrix column. 
Finally, each row is divided by its sum, resulting in a probability weight 
matrix. 
The user can choose to weigh every sequence by its score. If he does this, the 
program adds the score of a specific sequence at the position of each 
nucleotide (instead of just 1, as described above). Then, the program iterates 
through the resulting matrix and subtracts the minimum of each line from its 
values if it is below zero. This ensures that the final probabilities will be 
between 0 and 1. As described above, the user then can choose to create some 
background noise using base probabilities. Finally, the program divides every 
line by its sum, yielding the output probability weight matrix.

If the user wants to calculate aprobability weight matrix based on a list of 
sequences with different lengths, the program will use the Expectation 
Maximization Algorithm (short EM-Algorithm) to compute a matrix. This 
algorithm is described in: 
Do, C. B., & Batzoglou, S. (2008). What is the expectation maximization 
algorithm? Nature Biotechnology, 26(8), 897-899. 
The user can determine when to stop the algorithm with the following  
conditions:
He can set the maximum amount of iterations of the algorithm.
He can set the maximal difference between the two iterating matrices, if the 
difference is smaller than the given value, the program stops.
He also chooses the cutoff of the algorithm. 



 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PROCEDURE: GET ENZYME BINDING POSITIONS ON SEQUENCE
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
In order to find all enzyme binding positions on a sequence, the user has to 
provide a .fasta file containing the sequence, and a .mat file containing the 
probability weight matrix. The matrix is converted to a logarithmic matrix 
(if necessary), using base probabilities if wanted by the user. The program 
will then iterate through all positions of the sequence and calculate 
a score corresponding to the sum of the logarithmic probabilities of the 
nucleotides in the given position. A sequence is considered as a binding 
position if its score is bigger or equal as a user defined cutoff. 
In order to be clear and precise the sequence returned is the 
complementary reverse branch when (-) is indicated.



 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PROCEDURE: LOGO
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
This creates a png visualization of a probability weight matrix. It is saved 
as yourlogo.png in the folder above.



 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PROCEDURE: CORRELATE RESULTS TO BINDING POSITIONS
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
It is possible for the user to compare the binding positions found in previous
analysis to a genomic coordinate file. To every binding position, the genomic 
score of its position is added. This additional score is printed as additional
column of the output .csv file. 
