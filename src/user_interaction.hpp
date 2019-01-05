#ifndef user_interaction_hpp
#define user_interaction_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "Matrix.hpp"

#define COORD 0
#define SEQ 1
#define START 2

/*!
 * @class SearchResult
 *
 * @brief Contains a sequence, a position, a score and a direction.
 */
class SearchResult {
public:
    std::string sequence;
    unsigned int position;
    double score;
    char direction;
};
/*!
 * @class SearchResults
 *
 * @brief Vector of SearchResult
 */
class SearchResults {
public:
    std::string description;
    std::vector<SearchResult> searchResults;
};


// Three columned table which associates two objects.
typedef std::vector<std::vector<unsigned int >> Association_Table;

enum PROCEDURE { MatrixFromSequences, SequencesFromMatrix, Logo, CorrelateResults, Exit };
enum BP_FILL { AllEqual, UserDefined, FromSequence, NotUsed };
enum SEQ_SOURCE { CoordAndSeq, OnlySeq, FromSearchResult };
enum LIST_FILE { Fasta, NormalList, SeparatedList };


//-----------------------------------------------------------------------


/*!
 * @brief Asks what the user wants to do
 *
 * @return Answer of type PROCEDURE
 */
PROCEDURE whatToDo();


/*!
 * @brief Ask the user for the cutoff, returns said cutoff
 *
 * @return cutoff value specified by user
 */
double ask_cutoff();

/*!
 * @brief Ask the user for the filename; checks filename; returns filename
 *
 * @return name of the file the user would like to open
 */
std::string ask_name_fasta();


/*!
 * @brief Allows to ask the name of the matrix file we want to open,
 *		  with exception handling
 *
 * @return Name of the file to open
 */
std::string ask_name_matrix();


/*!
 * @brief Checks if the matrix file is valid
 *
 * @param name of the file we want to open
 * @return 0 if format is valid of 1 if not
 */
bool InvalidFormatMat(std::string file_name);


/*!
 * @brief Checks if the file if valid
 *
 * @param name of the file we want to open
 * @return 0 if format is valid of 1 if not
 */
bool InvalidFormat(std::string file_name);

/*!
 * @brief   Prints out current progress
 * @param   Current position, current filesize
 */
void print_progress(int position, int filesize);

/*!
 * @brief   A simple return
 */
void ret();

/*!
 * @brief Prints search results to file
 *
 * @param Search results, filename
 */
void print_results(SearchResults results, std::string filename);

/*!
 * @brief Prints search results to file
 *
 * @param Search results, file
 */
void print_results2(SearchResults results, std::ofstream &outputfile);

/*!
 * @brief Prints search results and corresponding genomic scores to file
 *
 * @param Search results, a vector of doubles with the genomic scores,  filename
 */
void print_results_correlated(SearchResults results, std::vector <double> gen_score, std::string filename);

/*!
 * @brief warns that there is a problematic char
 *
 * @param character causing problem
 */
void nucleotide_warning(char c);


/*!
 * @brief Ask the choice ( 0 or 1 or 2 ) corresponding to the question "what we want to do"
 *
 * @return the Choice made by the user, as enum BP_FILL.
 */
BP_FILL CoutCin_AskBaseProb();

/*!
 * @brief Ask the Prob corresponding to the nucleotide.
 *
 * @param the nucleotide needed (A or C or G or T)
 * @return the Prob enter by the user for this nucleotide
 */
double CoutCin_AskBaseProb0(char letter);

/*!
 * @brief Ask the name for the outputfile
 *
 * @return Name of the outputfile
 */
std::string Ask_Outputfile_Name();

/*!
 * @brief Handles everything concerning base probabilities, returns a vector of doubles
 *          with the base probabilities
 *
 * @return  Vector with base probabilities
 */
std::vector<double> AskBaseProb();

/*!
 * @brief   Lets the user choose the base probabilities, checks for their
 *          validity.
 *
 * @return  Vector with base probabilities
 */
std::vector<double> User_Defined_Base_Prob();

/*!
 * @brief   Asks the user as what kind of matrix he wants to save his matrix as.
 *
 * @return  MATRIX_TYPE with the corresponding type
 */
MATRIX_TYPE Ask_Return_Matrix_Type();

/*!
 * @brief   Asks the user if he wants to weigh each sequence by its score for
 *          the calculation of the matrix
 *
 * @return  1 if yes, else 0
 */
bool ask_matrix_from_sequences_weighed();

/*!
 * @brief   Asks the user if he wants to create a new matrix based on the
 *          matches found in the previous search
 *
 * @return  1 if yes, else 0
 */
bool ask_matrix_from_search_matches();

/*!
 * @brief   Warns the user that there are either no sequences found or that
 *          they are of different length
 *
 */
void error_input_sequence();

/*!
 * @brief   Informs the user that there was an error reading a coordinate file
 *
 * @param   Line on which it failed
 */

void error_reading_coordinates(unsigned int line);

/*!
 * @brief   Ask the name of a coordinate files
 *
 * @return  Name of user input
 */
std::string ask_coordinate_filename();

/*!
 * @brief   Asks the user if a line description is present
 *
 * @return  1 if yes, 0 if no.
 */
bool ask_line_description_present();

/*!
 * @brief   Ask the user to determine which sequence he wants to analyze with which set of
 *          genomic coordinates
 *
 * @param   Vector of strings with all coordinate descriptions, a vector of strings with all
 *          sequence descriptions
 *
 * @return  Association table, with number of coordinates in column 0 (COORD) and number of sequences
 *          in column 1 (SEQ) and the starting position in column 2 (START)
 */
Association_Table associate_genomic_with_sequences(std::vector<std::string> coordinate_description,
                                                   std::vector<std::string> sequence_description);


/*!
 * @brief   Ask the user to determine which searchResult he wants to analyze with which set of
 *          genomic coordinates
 *
 * @param   Vector of strings with all coordinate descriptions, a vector of searchResults
 *
 * @return  Association table, with id of coordinates in column 0 (COORD), id of corresponding
 *          search results in column 1 (SEQ) and the starting position in column 2 (START)
 */
Association_Table associate_genomic_with_result(std::vector<std::string> coordinate_description,
                                                std::vector<SearchResults> result_list);


/*!
 * @brief   Prints warning if sequence to analyze doesnt exist
 */
void error_sequence_doesnt_exist();


/*!
 * @brief   Asks the user, where he wants to get the sequences from that he wants to analyze
 *          (to create a PWM)
 *
 * @return  Answer as type SEQ_SOURCE
 */
SEQ_SOURCE ask_source_sequence();

//-----------------------------------------------------------------------
// 							LOGO USER INTERACTION
//-----------------------------------------------------------------------

/*!
 * @brief   prints the update on creation of the logo
 */
void logo_in_process();

/*!
 * @brief   Asks the user what the name of the file with the probability matrix is
 *
 * @return  file name as a string
 */
std::string ask_logo_matrix();

/*!
 * @brief   prints logo process for user
 *
 * @param   position being processed
 */
void position_in_process(int pos, int size);

//-----------------------------------------------------------------------



/*!
 * @brief   Asks the user what number of iterations he wants to use
 *
 * @return  the number of iterations the user chose
 */
int ask_iterations (int length);

/*!
 * @brief   Print the matrix into an output file
 *
 * @param   an output file and a matrix
 */
void print_into_file(std::ostream & out, Matrix_Neo matrix);

/*!
 * @brief   Asks the user what cutoff he wants to use and takes care of the possible errors
 *
 * @param   a specific number. The user has to give a bigger number than this one
 * @return  the cutoff the user chose
 */
double ask_cutoff2(double max_score);

/*!
 * @brief   tells the user where (in what file) he is
 *
 */
void path ();

/*!
 * @brief   Returns a path
 *
 * @return 	in what directory we are at the moment
 */
std::string get_working_path();


/*!
 * @brief   Asks the user the maximum of EM_algorithm he wants to do
 *
 * @return 	the maximum of times the user wants to do the EM_algorithm
 */
int maximum_EM ();


/*!
 * @brief   Asks the user what is the difference between two matrices to stop the EM_algorithm
 *
 * @return 	the differences between two matrices the user chose
 */
double differences_matrices ();

/*!
 * @brief   If the user has a separate file with the sequence, it asks what this file looks like.
 *
 * @return  Answer as type LIST_FILE.
 */
LIST_FILE ask_list_file_type();


/*!
 * @brief   If the user has a file where the sequences are separated by characters other than \n,
 *          he can specify the separation character with this function.
 *
 * @return  Answer as type char.
 */
char ask_separation_character();


/*!
 * @brief   Asks the user, which file he wants to open (in general, no specific filetype)
 *
 * @return  Answer as string
 */
std::string ask_inputfile_name();

/*!
 * @brief   Prints "DONE" to terminal
 */
void done();

/*!
 * @brief   Prints Error if no search results were read (for procedure 1)
 */
void error_no_search_result_read();

/*!
 * @brief   Asks user if he wants to continue correlating
 * @return  1 if yes, 0 if no.
 */
bool correlate_more();

/*!
 * @brief   Check if the input from the user is in the acceptable range and returns a double
 * @return  the number the user entered as a double
 */
double ask_for_a_number_infinitely();

/*!
 * @brief   Asks user if he wants correlate the found results to the search results
 * @return  1 if yes, 0 if no.
 */
bool ask_correlate_to_search_results();


/*!
 * @brief   Error message if user launches program with invalid flags
 */
void error_invalid_flags();


/*!
 * @brief   Checks if the file can be opened. Prints error if not.
 * @param   Name of file to open as string
 * @return  1 if yes, 0 if no
 */
bool checkfile(std::string filename);

/*!
 * @brief   Checks if the file already exists. If yes, ask if user wants to overwrite it.
 * @param   Name of file to open as string
 * @return  1 to proceed (overwrite or no file of that name), 0 if not
 */
bool overwrite(std::string filename);

/*!
 * @brief   Checks if the user inputs a 1 or a 0 
 * 
 * @return  user input 1 or 0
 */
bool correct_bool();
//-----------------------------------------------------------------------


#endif
