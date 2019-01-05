#ifndef procedures_hpp
#define procedures_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

#include "Matrix.hpp"
#include "utility.hpp"
#include "user_interaction.hpp"

//-----------------------------------------------------------------------

/*! 
 * @Brief   Handles everything if the user wants to determine on
 *          what position a given enzyme binds, with given probability
 *          matrix and sequence
 */
void enzyme_on_sequence();


/*!
 * @Brief   Handles everything if the user wants to determine on
 *          what position a given enzyme binds, with given probability
 *          matrix and sequence. DOesnt use as much memory space
 */
void enzyme_on_sequence2();



/*! 
 * @Brief   Handles everything if the user wants to determine
 *          a probability weight matrix from given sequences
 */
void enzyme_from_sequences();


/*!
 * @Brief   Handles everything to load search results into file, when containing score
 *
 * @return  Vector of search results with the input
 */
std::vector<SearchResults> input_search_results();

/*!
 * @Brief   Handles everything if the user wants to compare the genomic coordinates
 *          with
 */
void correlate_coordinates_with_results(std::vector <SearchResults> result_list);

/*!
 * @Brief   Handles everything if the user knows the length of the binding
 *          site
 */
Matrix binding_length_known(std::vector<SearchResults> seq_to_analyze);

/*!
 * @Brief   Handles everything if the user doesn't know the length of the binding
 *          site
 */
void binding_length_unknown(std::vector<std::string> sequence_list, std::string name, Base_Prob base_probabilities  = {1, 1, 1, 1});

/*! 
 * @Brief   Handles everything if the user wants to extract a list of search results
 *          from coordinates and a sequence file
 *
 * @return  Returns a list of SearchResults of the found sequences
 */
std::vector<SearchResults> seq_source_CoordAndSeq();

/*! @Brief   Handles everything if the user wants to extract a list of search results
 *          from a sequence file. All scores will be set to 1
 *
 * @return  Returns a list of SearchResults of the found sequences
 */
std::vector<SearchResults> seq_source_OnlySeq();

/*!
 * @Brief   Handles everything if the user already has a list of search results which
 *          he wants to load into the program
 *
 * @return  Returns a list of SearchResults of the found sequences
 */
std::vector<SearchResults> seq_source_FromSearchResult();


/*!
 * @Brief   Converts char** to vector of strings
 *
 * @param   size of vector, char** as pointer to string
 * @return  vector of strings
 */

std::vector<std::string> flags_to_strings(int argc, char** argv);


//------------------------------------------------------------------------
// NON-VERBOSE (When program is started from command line)

/*!
 * @brief   Handles everything if user starts program with flag -help
 */
void flag_help();

/*!
 * @brief   Handles everything if user starts program with flag -about
 */
void flag_about();

/*!
 * @brief   Handles everything if user starts program with flag -bindingsites
 * @param   Flags as vector of strings
 */
void flag_bindingsites(std::vector<std::string> flags);

/*!
 * @brief   Handles everything if user starts program with flag -getmatrix
 * @param   Flags as vector of strings
 */
void flag_getmatrix(std::vector<std::string> flags);




#endif
