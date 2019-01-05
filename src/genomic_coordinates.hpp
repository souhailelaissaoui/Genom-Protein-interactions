//
//  genomic_coordinates.hpp
//
//
//  Created by Mättu on 29.11.16.
//
//

#ifndef genomic_coordinates_hpp
#define genomic_coordinates_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include "user_interaction.hpp"


class Coordinates {

private:

    std::vector <unsigned int> start_pos;
    std::vector <unsigned int> length;
    std::vector <std::string> description;
    std::vector <double> score;



    std::string file_description;
    std::string location;
    bool line_description_present;

public:

    /*!
     * @brief Extracts sequences from sequence-file with given genomic coordinates.
     *          Each sequence is analyzed with the same genomic coordinates
     *
     * @param The name of the sequence file, the starting position of the sequence
     */
    // DOESNT WORK
    SearchResults get_search_results(std::string seq_filename, unsigned int sequence, unsigned int startpos);
    /*!
     * @brief   Getter funciton for location.
     *
     * @return  Location parameter of class.
     */
    std::string get_location();


    /*!
     * @brief   Creates an empty genomic coordinates, only with description, location and
     *          the information whether there will be a line description.
     *
     * @param   The file description, the description of the location of the genes, and
     *          the information whether there will be a line description
     */
    Coordinates( std::string file_description, std::string location, bool line_description_present);

    /*!
     * @brief   Fills a new line of genomic coordinates
     *
     * @param   Parameters to fill them with, in a string.
     */
    void fillNewLine(std::string line);

    /*!
     * @brief   Goes through SearchResults and Genomic Coordinates, and creates a vector
     *          of doubles with the genomic score of each sequence
     *
     * @param   A SearchResults to analyze, the starting position of the search results
     *          in relation to the genome.
     *
     * @return  Vector of doubles with the genomic scores
     */
    std::vector <double> position_score(SearchResults input, unsigned int pos_zero);
};


#endif /* genomic_coordinates_hpp */
