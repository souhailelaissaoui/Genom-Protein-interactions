#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include "Matrix.hpp"
#include "user_interaction.hpp"
#include "genomic_coordinates.hpp"


//==========================================================================================
// SEQUENCE ANALYSIS
/*!
 * @brief Goes through all sequences in the file and returns sequences with a score above 5
 *         Optimized version 2
 *          Estimated time for PPM on human genome: 30min, with flag -O2, MacBook Pro, cutoff 5
 *
 * @param Name of sequence file, Matrix to compare it with, Cutoff
 * @return A vector of all search results (each in one structure SearchResults)
 */
std::vector<SearchResults> analyze_sequence_opt2(std::string filename, Matrix matrix, double cutoff);


/*!
 * @brief   Goes through all sequences in the file and returns sequences with a score above 5
 *          prints directly to output file.
 *          Optimized version 4
 *
 *
 * @param Name of sequence file, Matrix to compare it with, Cutoff
 * @return A vector of all search results (each in one structure SearchResults)
 */
void analyze_sequence_opt4(std::string filename, Matrix matrix, double cutoff, std::string outputfile_name);


/*!
 * @brief Fills up a Search result. This should be made as a constructor of the class
 *
 * @param A nuc-list of the sequence, the position, the score and the direction of the result
 * @return A filled class SearchResult
 */
SearchResult fill_search_result(std::list<nuc> sequence, int position_counter, double score, char direction);


/*!
 * @brief Determines if input character is valid (a nucleotide)
 *
 * @param Input character
 * @return true if valid, false if not valid
 */
bool valid_character(char character);


/*!
 * @brief   Determines size of a file
 * @param   Name of file to analyze
 * @return  Number of bytes of file, only approximatively
 */
int filesize(std::string filename);


/*!
 * @brief Translates a list of the enum "nuc" to a string
 *
 * @param A list of nucleotides to be translated
 * @return Nucleotide sequence as a string
 */
std::string sequence_string_from_nuc_list(std::list<nuc>);

/*!
 * @brief Creates a probability weight matrix from SearchResults with all the same length
 *          Function doesn't weigh the probabilities by the score of each individual sequence
 *
 * @param An instance of SearchResults with all the same length, base probabilities to add as noise
 * @return A regular probability weight matrix (type absoluteMatrix). Matrix of size 0 in case of error.
 */
Matrix_Neo matrix_from_same_length_sequences_not_weighted(std::vector<SearchResults> input, Base_Prob base_prob);


/*!
 * @brief Creates a probability weight matrix from SearchResults with all the same length
 *          Function weighs each probability by the score of each sequence
 *
 * @param An instance of SearchResults with all the same length, base probabilities to add as noise
 * @return A regular probability weight matrix (type absoluteMatrix). Matrix of size 0 in case of error.
 */
Matrix_Neo matrix_from_same_length_sequences_weighted(std::vector<SearchResults> input, Base_Prob base_prob);


/*!
 * @brief Creates a probability weight matrix from a vector of SearchResults with all the same length
 *          Takes bool as argument if it should be weighed by the score or not
 *
 * @param An instance of SearchResults with all the same length
 * @return A regular probability weight matrix (type absoluteMatrix)
 */
Matrix_Neo matrix_from_same_length( std::vector<SearchResults> input, Base_Prob base_prob, bool weighed );


/*!
 * @brief Checks if all search results have the same length, returns this length
 *
 * @param An instance of SearchResults
 * @return The length if yes, false if no or input doesn't contain searchResults
 */
unsigned int searchResults_same_length(std::vector<SearchResults> input);


/*!
 * @brief   Fills a vector of coordinates from a input file
 *
 * @param   The name of the genomic file, a bool saying whether the file contains a column with
 *          information on every sequence.
 * @return  A vector of coordinates.
 */
std::vector <Coordinates> read_coordinates(std::string filename, bool line_description_present);


/*!
 * @brief   Goes through files and extracts all sequence descriptions present
 *
 * @param   The name of the sequence file.
 * @return  A vector of strings with all sequence descriptions.
 */
std::vector<std::string> extract_sequence_descriptions(std::string filename);


/*!
 * @brief   Returns a vector of coordinates
 *
 * @param   The name of the sequence file.
 * @return  A vector of strings with all sequence descriptions.
 */
std::vector<std::string> get_descriptions_from_coord_list(std::vector<Coordinates> coord_list);

/*!
 * @brief   Returns a matrix
 *
 * @param   One sequence and the number of iterations
 * @return  A matrix made from the sequence according to the number of iterations
 */
Matrix_Neo sequence_to_PPM (std::string sequence, int n);

/*!
 * @brief   Returns a matrix
 *
 * @param   A vector of sequences, the number of iterations and the baseprobabilities for each nucleotide
 * @return  A matrix made from the sequences according to the number of iterations and the baseprobabilities
 */
Matrix_Neo sequences_to_PPM (std::vector<std::string> sequences, unsigned int n);

/*!
 * @brief   Returns a score
 *
 * @param   A sequence and a matrix
 * @return  the score of the sequence according to the matrix
 */
double sequence_score (std::string sequence, Matrix_Neo PPM);

/*!
 * @brief   Returns a vector of sequences
 *
 * @param   A vector of sequences, a matrix and a number that represents the cutoff
 * @return  A vector of sequences found from a set of sequences and its matrix, taking into account the cutoff
 */
std::vector<std::string> PPM_to_Sequence (std::vector<std::string> binding_sites, Matrix_Neo PPM,  double cutoff);

/*!
 * @brief   Returns a matrix
 *
 * @param   A vector of sequences, a number of iterations, a number that represents the cutoff, the baseprobabilities, the maximum and the differences
 * @return  A stabilised matrix made from the different sequences
 */
Matrix_Neo EM_algorithm (std::vector<std::string> Sequences, int n, double cutoff, Base_Prob base_probabilities, int max, double diff);

/*!
 * @brief   Returns a vector of sequences
 *
 * @return  a vector of sequences made from a file that the user gives
 */
std::vector <std::string> Initialization_string();

/*!
 * @brief   returns a score
 *
 * @param   a matrix
 * @return  return the maximal score of the matrix
 */
double max_score (Matrix_Neo & PPM);

/*!
 * @brief   Returns a bool
 *
 * @param   two matrices and a number that represents the differences between the 2 matrices
 * @return  if the two matrices have a difference smaller or equal to the number of differences
 */
bool diff_matrices (Matrix_Neo mat1, Matrix_Neo mat2, double difference);

/*!
 * @brief   returns the size of the smallest sequence among a vector of sequences
 *
 * @param   a vector of sequences
 * @return  the size of the smallest sequence among all sequences
 */
int smallest_length (std::vector<std::string> sequences);


/*!
 * @brief   Goes through sequence file, and saves each sequence as an individual search
 *          result. Direction is always '+', score always 1 and startpos always 0.
 *
 * @param   The name of the sequence file.
 * @return  A SearchResults with the result of the reading of the file.
 */
SearchResults read_sequencefile_to_searchresults(std::string filename);


/*!
 * @brief   Goes through sequence file, and saves each sequence as an individual search
 *          result. Direction is always '+', score always 1 and startpos always 0.
 *
 * @param   The name of the sequence file.
 * @return  A SearchResults with the result of the reading of the file.
 */
std::vector <SearchResults> read_searchresult_file(std::string filename);


/*!
 * @brief   Goes through file and saves each line as individual search result.
 *          Direction is always '+', score always 1 and startpos always 0.
 *
 * @param   The name of the sequence file.
 * @return  A SearchResults with the result of the reading of the file.
 */
SearchResults read_sequence_list_to_searchresults(std::string filename);


/*!
 * @brief   Goes through file and saves each sequence as individual search result.
 *          Each sequence is separated by the character defined in the parameters
 *
 * @param   The name of the sequence file, the delimiting character.
 * @return  A SearchResults with the result of the reading of the file.
 */
SearchResults read_char_separated_to_searchresults(std::string filename, char delim);

/*!
 * @brief   Goes through file and savec each sequence as individual string.
 *
 * @param   The name of the sequence file
 * @return  A vector of strings with the result of the reading of the file.
 */
std::vector <std::string> ExtractSequence(std::string const& entry_name);


/*!
 * @brief   Goes through vector of search results and extracts all sequences in it.
 *
 * @param   The vector of search results
 * @return  A vector of strings.
 */
std::vector <std::string> string_list_from_searchResults(std::vector <SearchResults> input);


Matrix_Neo normalized(Matrix_Neo input);

#endif
