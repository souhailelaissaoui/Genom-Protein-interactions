#include "utility.hpp"
#include <assert.h>
#include <iostream>
#include <fstream>

#define LINE_SIZE 10000

std::map< char, nuc > charmap{
    {'A', A}, {'a', A},
    {'C', C}, {'c', C},
    {'G', G}, {'g', G},
    {'T', T}, {'t', T},
    {'N', N}, {'n', N},
    {'-', N}, {'.', N}
};

std::map< char, nuc > complementmap{
    {'A', T}, {'a', T},
    {'C', G}, {'c', G},
    {'G', C}, {'g', C},
    {'T', A}, {'t', A},
    {'N', N}, {'n', N},
    {'-', N}, {'.', N}
};

std::map< nuc, char > backwardsmap {
    {A, 'A'},
    {C, 'C'},
    {G, 'G'},
    {T, 'T'},
    {N, 'N'}
};

//==========================================================================================
// SEQUENCE ANALYSIS
//==========================================================================================
// I am very sorry for the uglyness of this function
//
// TODO: Doesn't work if at position 0; reverse strand doesnt work
// TODO: Progress bar

// Souhail : I deleted the not used functions
//==========================================================================================
// Even further optimized version of the function, INPROGRESS
//
// TODO: Progress bar

std::vector<SearchResults> analyze_sequence_opt2(std::string filename, Matrix matrix, double cutoff) {

    // File to read
    std::ifstream entry_file(filename);

    // Search results of all sequences
    std::vector <SearchResults> output;

    // Search results of one sequence
    SearchResults sequence_matches;

    std::list<nuc> forwardSequence;
    std::list<nuc> backwardSequence;

    bool fill(true);
    //The following variable is commented because it's not used and create warnings - SOuhail
   // bool description(true);

    unsigned int char_counter(0);
    unsigned int length(matrix.get_length());
    unsigned int position_counter(1);
    unsigned int line_size;

    unsigned int print_counter(0);

    int size_of_file;
    size_of_file = filesize(filename);

    //bool first_line(true); // not used ? Souhail
    unsigned int index(0);

    double score;

    char line[LINE_SIZE];
    std::string seq_description;

    // Handles first line
    while (entry_file.peek() == '>' or entry_file.peek() == ';') {
        getline(entry_file, seq_description, '\n');
        sequence_matches.description += seq_description;
    }


    while(entry_file.get(line, LINE_SIZE, '>')) {
        // Sets line_size to amount of characters read in get
        line_size = entry_file.gcount();

        //std::cout << line_size << std::endl;

        index = 0;

        if (print_counter > 10) {
            print_progress(entry_file.tellg(), size_of_file);
            print_counter = 0;

        }
        print_counter++;


        while (fill) {
            if (valid_character(line[index])) {
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.push_front(complementmap[line[index]]);
                char_counter++;

            }

            if(char_counter >= length) {
                fill = false;
            }


            index++;


        }

        // For initialized sequence
        // What to do if forward is binding
        score = matrix.sequence_score(forwardSequence);

        if(score >= cutoff) {

            SearchResult sequence_match(fill_search_result(forwardSequence, position_counter, score, '+'));
            sequence_matches.searchResults.push_back(sequence_match);
        }

        // What to do if backward is binding
        score = matrix.sequence_score(backwardSequence);

        if(score >= cutoff) {
            SearchResult sequence_match(fill_search_result(backwardSequence, position_counter, score, '-'));
            sequence_matches.searchResults.push_back(sequence_match);
        }


        // For all following combinations
        while (index < line_size) {

            // Skips if character is not valid
            if(valid_character(line[index]))
            {
                position_counter++;

                // Updates sequence with the new character
                forwardSequence.pop_front();
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.pop_back();
                backwardSequence.push_front(complementmap[line[index]]);

                // What to do if forward is binding
                score = matrix.sequence_score(forwardSequence);

                if(score >= cutoff) {

                    SearchResult sequence_match(fill_search_result(forwardSequence, position_counter, score, '+'));
                    sequence_matches.searchResults.push_back(sequence_match);
                }

                // What to do if backward is binding
                score = matrix.sequence_score(backwardSequence);

                if(score >= cutoff) {
                    SearchResult sequence_match(fill_search_result(backwardSequence, position_counter, score, '-'));
                    sequence_matches.searchResults.push_back(sequence_match);
                }
            }
            index++;
        }


        if (entry_file.peek() == '>' or entry_file.peek() == EOF) {
            output.push_back(sequence_matches);
            sequence_matches = SearchResults();

            if(!forwardSequence.empty())
                forwardSequence.clear();
            if(!backwardSequence.empty())
                backwardSequence.clear();

            position_counter = 1;
            char_counter = 0;
            //first_line = false;
            fill = true;

            while (entry_file.peek() == '>' or entry_file.peek() == ';') {
                getline(entry_file, seq_description);
                sequence_matches.description += seq_description;
            }
        }
    }

    output.push_back(sequence_matches);

	print_progress(size_of_file, size_of_file);
	ret();


    entry_file.close();
    return output;
}

//==========================================================================================

void analyze_sequence_opt3(std::string filename, Matrix matrix, double cutoff, std::string outputfile_name) {

    std::ofstream output_file;

    output_file.open(outputfile_name);

    // File to read
    std::ifstream entry_file(filename);


    std::list<nuc> forwardSequence;
    std::list<nuc> backwardSequence;

    bool fill(true);
    //The following variable is commented because it's not used and create warnings - SOuhail
    // bool description(true);


    unsigned int char_counter(0);
    unsigned int length(matrix.get_length());
    unsigned int position_counter(1);
    unsigned int line_size;

    unsigned int print_counter(0);


    int size_of_file;
    size_of_file = filesize(filename);

    //bool first_line(true); // not used ? Souhail
    unsigned int index(0);

    double score;

    char line[LINE_SIZE];
    std::string seq_description;


    // Handles first line
    while (entry_file.peek() == '>' or entry_file.peek() == ';') {
        getline(entry_file, seq_description);
        output_file << seq_description << std::endl;
    }


    while(entry_file.get(line, LINE_SIZE, '>')) {
        // Sets line_size to amount of characters read in get
        line_size = entry_file.gcount();

        //std::cout << line_size << std::endl;

        index = 0;

        if (print_counter > 10) {
            print_progress(entry_file.tellg(), size_of_file);
            print_counter = 0;

        }
        print_counter++;


        while (fill) {
            if (valid_character(line[index])) {
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.push_front(complementmap[line[index]]);
                char_counter++;

            }

            if(char_counter >= length) {
                fill = false;
            }


            index++;


        }

        // For initialized sequence
        // What to do if forward is binding
        score = matrix.sequence_score(forwardSequence);

        if(score >= cutoff) {
            output_file << sequence_string_from_nuc_list(forwardSequence) << ";\t" << position_counter
            << ";\t" << score << ";\t" << '+' << std::endl;
        }

        // What to do if backward is binding
        score = matrix.sequence_score(backwardSequence);

        if(score >= cutoff) {
            output_file << sequence_string_from_nuc_list(forwardSequence) << ";\t" << position_counter << ";\t" << score << ";\t" << '-' << std::endl;
        }


        // For all following combinations
        while (index < line_size) {

            // Skips if character is not valid
            if(valid_character(line[index]))
            {
                position_counter++;

                // Updates sequence with the new character
                forwardSequence.pop_front();
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.pop_back();
                backwardSequence.push_front(complementmap[line[index]]);

                // What to do if forward is binding
                score = matrix.sequence_score(forwardSequence);

                if(score >= cutoff) {
                    output_file << sequence_string_from_nuc_list(forwardSequence) << ";\t" << position_counter << ";\t" << score << ";\t" << '+' << std::endl;
                }

                // What to do if backward is binding
                score = matrix.sequence_score(backwardSequence);

                if(score >= cutoff) {
                    output_file << sequence_string_from_nuc_list(forwardSequence) << ";\t" << position_counter << ";\t" << score << ";\t" << '-' << std::endl;
                }
            }
            index++;
        }


        if (entry_file.peek() == '>') {
            while (entry_file.peek() == '>' or entry_file.peek() == ';') {
                getline(entry_file, seq_description);
                output_file << seq_description << std::endl;
            }

        }
    }


    print_progress(size_of_file, size_of_file);
    ret();


    entry_file.close();
    output_file.close();
}



void analyze_sequence_opt4(std::string filename, Matrix matrix, double cutoff, std::string outputfile_name) {

    std::ofstream output_file;

    output_file.open(outputfile_name + ".csv");

    // File to read
    std::ifstream entry_file(filename);

    // Search results of one sequence
    SearchResults sequence_matches;

    std::list<nuc> forwardSequence;
    std::list<nuc> backwardSequence;

    bool fill(true);
    //The following variable is commented because it's not used and create warnings - SOuhail
    // bool description(true);


    unsigned int char_counter(0);
    unsigned int length(matrix.get_length());
    unsigned int position_counter(1);
    unsigned int line_size;

    unsigned int print_counter(0);


    int size_of_file;
    size_of_file = filesize(filename);

    //bool first_line(true); // not used ? Souhail
    unsigned int index(0);

    double score;

    char line[LINE_SIZE];
    std::string seq_description;


    // Handles first line
    while (entry_file.peek() == '>' or entry_file.peek() == ';') {
        getline(entry_file, seq_description, '\n');
        output_file << seq_description << ";" << std::endl;
    }


    while(entry_file.get(line, LINE_SIZE, '>')) {

        // Sets line_size to amount of characters read in get
        line_size = entry_file.gcount();

        //std::cout << line_size << std::endl;

        index = 0;

        if (print_counter > 10) {
            print_progress(entry_file.tellg(), size_of_file);
            print_counter = 0;

        }
        print_counter++;


        while (fill) {
            if (valid_character(line[index])) {
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.push_front(complementmap[line[index]]);
                char_counter++;

            }

            if(char_counter >= length) {
                fill = false;
            }


            index++;


        }

        // For initialized sequence
        // What to do if forward is binding
        score = matrix.sequence_score(forwardSequence);

        if(score >= cutoff) {

            SearchResult sequence_match(fill_search_result(forwardSequence, position_counter, score, '+'));
            sequence_matches.searchResults.push_back(sequence_match);
        }

        // What to do if backward is binding
        score = matrix.sequence_score(backwardSequence);

        if(score >= cutoff) {
            SearchResult sequence_match(fill_search_result(backwardSequence, position_counter, score, '-'));
            sequence_matches.searchResults.push_back(sequence_match);
        }


        // For all following combinations
        while (index < line_size) {

            // Skips if character is not valid
            if(valid_character(line[index]))
            {
                position_counter++;

                // Updates sequence with the new character
                forwardSequence.pop_front();
                forwardSequence.push_back(charmap[line[index]]);
                backwardSequence.pop_back();
                backwardSequence.push_front(complementmap[line[index]]);

                // What to do if forward is binding
                score = matrix.sequence_score(forwardSequence);

                if(score >= cutoff) {
                    SearchResult sequence_match(fill_search_result(forwardSequence, position_counter, score, '+'));
                    sequence_matches.searchResults.push_back(sequence_match);
                }

                // What to do if backward is binding
                score = matrix.sequence_score(backwardSequence);

                if(score >= cutoff) {
                    SearchResult sequence_match(fill_search_result(backwardSequence, position_counter, score, '-'));
                    sequence_matches.searchResults.push_back(sequence_match);
                }
            }
            index++;
        }

        print_results2(sequence_matches, output_file);
        sequence_matches = SearchResults();


        if (entry_file.peek() == '>' or entry_file.peek() == EOF) {
            sequence_matches = SearchResults();

            if(!forwardSequence.empty())
                forwardSequence.clear();
            if(!backwardSequence.empty())
                backwardSequence.clear();

            position_counter = 1;
            char_counter = 0;
            fill = true;

            while (entry_file.peek() == '>' or entry_file.peek() == ';') {
                getline(entry_file, seq_description);
                output_file << std::endl << seq_description << ";" << std::endl;
            }
        }


    }


    print_progress(size_of_file, size_of_file);
    ret();


    entry_file.close();
    output_file.close();
}

//==========================================================================================

SearchResult fill_search_result(std::list<nuc> sequence, int position_counter, double score, char direction)
{
    SearchResult sequence_match;

    sequence_match.sequence = sequence_string_from_nuc_list(sequence);
    sequence_match.position = position_counter;
    sequence_match.score = score;
    sequence_match.direction = direction;

    return sequence_match;
}

//==========================================================================================

bool valid_character(char character) {
    auto it = charmap.find(character);
    if (it == charmap.end())
    {
        // std::cout << "Unknown character: " << character << ". Character is skipped and ignored\n";
        return false;
    }
    else
        return true;
}

//==========================================================================================

int filesize(std::string filename) {
    std::ifstream in(filename, std::ios::binary | std::ios::ate);
    return in.tellg();
}

//==========================================================================================

std::string sequence_string_from_nuc_list(std::list<nuc> sequence) {
    std::string output;
    std::list<nuc>::iterator iterator;

    for (iterator = sequence.begin(); iterator != sequence.end(); iterator++)
        output += backwardsmap[*iterator];
    return output;
}

//==========================================================================================

Matrix_Neo matrix_from_same_length_sequences_not_weighted(std::vector<SearchResults>  input, Base_Prob base_prob) {
    Matrix_Neo output_matrix;

    assert(input.size() >= 1);


    if (base_prob[0] == 1 and
        base_prob[1] == 1 and
        base_prob[2] == 1 and
        base_prob[3] == 1)
    {
        base_prob = {0, 0, 0, 0};
    }

    unsigned int nb_search_results;
    unsigned int length_sequence( searchResults_same_length(input) );

    if(length_sequence == 0) {
        error_input_sequence();
        return output_matrix;
    }


    char character;

    // Initialization of matrix with proper size, with corresponding baseprob for
    // every nucleotide position
    std::vector <double> line;

    for (unsigned int index(0); index<NUMBER_NUCLEOTIDES; index++) {
        line.push_back(base_prob[index]);
    }

    for (unsigned int index(0); index<length_sequence; index++) {
        output_matrix.push_back(line);
    }

    // Counter for every search result
    for (unsigned int k(0); k<input.size(); k++) {
        nb_search_results = input[k].searchResults.size();

        for (unsigned int i(0); i<nb_search_results; i++) {

            for (unsigned int j(0); j<length_sequence; j++) {
                character = input[k].searchResults[i].sequence[j];
                output_matrix[j][charmap[character]] += 1;
            }
        }
    }

    // Division by sum of each line
    double line_sum = Matrix::sum_of_line(output_matrix[0]);

    for (unsigned int i(0); i<length_sequence; i++) {
        for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++) {
            output_matrix[i][j] /= line_sum;
        }
    }

    return output_matrix;
}

//==========================================================================================

Matrix_Neo matrix_from_same_length_sequences_weighted(std::vector<SearchResults>  input, Base_Prob base_prob) {
    Matrix_Neo output_matrix;

    assert(input.size() >= 1);


    if (base_prob[0] == 1 and
        base_prob[1] == 1 and
        base_prob[2] == 1 and
        base_prob[3] == 1)
    {
        base_prob = {0, 0, 0, 0};
    }

    unsigned int length_sequence( searchResults_same_length(input) );

    if(length_sequence == 0) {
        error_input_sequence();
        return output_matrix;
    }


    unsigned int nb_search_results;
    double line_min;
    double line_sum;

    char character;

    // Initialization of matrix with proper size, all zero
    for (unsigned int index(0); index<length_sequence; index++) {
        output_matrix.push_back({0,0,0,0});
    }

    // Counter for every search result, score is added for each match
    for (unsigned int k(0); k<input.size(); k++) {

        nb_search_results = input[k].searchResults.size();

        for (unsigned int i(0); i<nb_search_results; i++) {
            for (unsigned int j(0); j<length_sequence; j++) {
                character = input[k].searchResults[i].sequence[j];
                output_matrix[j][charmap[character]] += input[k].searchResults[i].score;
            }
        }
    }


    // If the smallest number is less or equal 0, it is subtracted from all values from the line
    for (unsigned int i(0); i<length_sequence; i++) {
        line_min = Matrix::min_of_line(output_matrix[i]);

        if (line_min<=0) {
            for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++) {
                output_matrix[i][j] -= line_min;
            }
        }
    }


    // To avoid havnig a probability of 0, the base probability for each nucleotide is added to
    // each line
    for (unsigned int i(0); i<length_sequence; i++) {
        for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++) {
            output_matrix[i][j] += base_prob[j];
        }
    }


    // Every value of the line is divided by the sum of each line, yielding a probability
    for (unsigned int i(0); i<length_sequence; i++) {
        line_sum = Matrix::sum_of_line(output_matrix[i]);

        for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++) {
            output_matrix[i][j] /= line_sum;
        }
    }

    return output_matrix;
}

//==========================================================================================

Matrix_Neo matrix_from_same_length( std::vector<SearchResults>  input, Base_Prob base_prob, bool weighed ) {

    if (weighed)
        return matrix_from_same_length_sequences_weighted( input, base_prob);
    else
        return matrix_from_same_length_sequences_not_weighted( input, base_prob);

}

//==========================================================================================

unsigned int searchResults_same_length(std::vector<SearchResults> input) {

    unsigned int nb_sequences(input.size());
    unsigned int nb_results;

    assert(nb_sequences);

    unsigned int length_result_0(0);

    for (unsigned int index(0); index<nb_sequences; index++) {
        if(input[index].searchResults.size()) {
            length_result_0 = input[0].searchResults[0].sequence.size();
            break;
        }
    }

    for (unsigned int k(0); k<nb_sequences; k++) {
        nb_results = input[k].searchResults.size();


        for (unsigned int i(0); i<nb_results; i++) {
            if (input[k].searchResults[i].sequence.size() != length_result_0) {
                return false;
            }
        }
    }

    return length_result_0;
}

//==========================================================================================

std::vector <Coordinates> read_coordinates(std::string filename, bool line_description_present) {
    std::ifstream file;
    file.open(filename);

    std::string file_description;
    std::string seq_description;
    std::string seq_descr_intermediate;
    std::string line;

    std::vector <Coordinates> output;

    file >> seq_description;
    Coordinates intermediate(file_description, seq_description, line_description_present);
    getline(file, line, '\n');
    intermediate.fillNewLine(line);


    while (file >> seq_descr_intermediate) {
        if (seq_description != seq_descr_intermediate) {
            seq_description = seq_descr_intermediate;
            output.push_back(intermediate);
            intermediate = Coordinates(file_description, seq_description, line_description_present);
        }

        getline(file, line, '\n');
        intermediate.fillNewLine(line);
    }

    output.push_back(intermediate);
    return output;
}

//==========================================================================================

std::vector<std::string> extract_sequence_descriptions(std::string filename) {
    std::ifstream file;
    file.open(filename);

    std::vector<std::string> output;
    std::string intermediate;

    //unsigned int streamsize(file.tellg()); unused

    while (getline(file, intermediate)) {
        if (intermediate[0] == '>') {
            output.push_back(intermediate);
        }
    }

    return output;
}

//==========================================================================================

std::vector<std::string> get_descriptions_from_coord_list(std::vector<Coordinates> coord_list) {
    std::vector<std::string> output;
    for (size_t id(0); id<coord_list.size(); id++) {
        output.push_back(coord_list[id].get_location());
    }
    return output;
}

//==========================================================================================

SearchResults read_sequencefile_to_searchresults(std::string filename) {
    SearchResults output;
    SearchResult intermediate;

    std::string line;
    std::string seq_intermediate;

    std::ifstream file;
    file.open(filename);

    // For every line
    while (getline(file, line)) {
        if (line[0] == '>' or line[0] == ';') {
            if (not(seq_intermediate.empty())) {

                intermediate.sequence = seq_intermediate;
                intermediate.score = 1;
                intermediate.position = 0;
                intermediate.direction = '+';

                output.searchResults.push_back(intermediate);
                seq_intermediate.clear();
                intermediate = SearchResult();
            }
        }
        else {
            seq_intermediate += line;
        }
    }

    // For the last line:
    if (not(seq_intermediate.empty())) {

        intermediate.sequence = seq_intermediate;
        intermediate.score = 1;
        intermediate.position = 0;
        intermediate.direction = '+';

        output.searchResults.push_back(intermediate);
        seq_intermediate.clear();
        intermediate = SearchResult();
    }


    return output;
}

//==========================================================================================

std::vector <SearchResults> read_searchresult_file(std::string filename) {
    std::vector <SearchResults> output;
    SearchResults intermed_chromosome;

    std::ifstream file;
    file.open(filename);

    std::string throwaway;
    std::string line;

    while (getline(file, line)) {

        std::replace( line.begin(), line.end(), ';', ' '); // replace all comma separators with space

        if (line[0] == '>') {
            if (not(intermed_chromosome.searchResults.empty())) {
                output.push_back(intermed_chromosome);
                intermed_chromosome = SearchResults();
            }

            intermed_chromosome.description = line;
        }

        else if (line.find("Seq") != std::string::npos);


        else {
            SearchResult intermed_match;

            std::istringstream line_stream(line);
            if(line_stream >> intermed_match.sequence
                           >> intermed_match.position
                           >> intermed_match.direction
                           >> intermed_match.score) {
                intermed_chromosome.searchResults.push_back(intermed_match);
            }

        }

    }

    if (not(intermed_chromosome.searchResults.empty())) {
        output.push_back(intermed_chromosome);
        intermed_chromosome = SearchResults();
    }


    return output;
}

//==========================================================================================

SearchResults read_sequence_list_to_searchresults(std::string filename) {
    SearchResults output;
    SearchResult intermediate;

    std::ifstream file;
    file.open(filename);

    intermediate.score = 1;
    intermediate.position = 0;
    intermediate.direction = '+';

    output.description = filename;

    while (getline(file, intermediate.sequence)) {
        output.searchResults.push_back(intermediate);
    }
    return output;

}

//==========================================================================================

SearchResults read_char_separated_to_searchresults(std::string filename, char delim) {
    SearchResults output;
    SearchResult intermediate;

    std::ifstream file;
    file.open(filename);

    intermediate.score = 1;
    intermediate.position = 0;
    intermediate.direction = '+';
    intermediate.sequence = "";

    output.description = filename;

    while (getline(file, intermediate.sequence, delim)) {

        if (not(intermediate.sequence.empty())) {
            output.searchResults.push_back(intermediate);
            intermediate.sequence = "";
        }


        // Skips to beginning of next text, by ignoring  all whitespace characters (ASCII value
        // less than 15 or 32)
        while ((file.peek() <= 15 or file.peek() == ' ') and file.peek()!=EOF) {
            file.ignore();
        }

    }
    return output;

}

//==========================================================================================

Matrix_Neo sequence_to_PPM (std::string sequence, int n)
{
	Matrix_Neo results;
	std::vector<double> nucleot (4, 0.0);
	int r=0;
	int m=0;

	for (int i=0;i<n;i++)
	{
		results.push_back(nucleot);
	}

	while (r <= n-1)
	{
		for (size_t i(r); i<sequence.size()-n+r+1; ++i)
		{
			switch (charmap[sequence[i]]) {
            case nuc::A:
                results[m][0] += 1.0;
                break;
            case nuc::C:
                results[m][1] += 1.0;
                break;
            case nuc::G:
                results[m][2] += 1.0;
                break;
            case nuc::T:
                results[m][3] += 1.0;
                break;
            case nuc::N:
				break;
			}

		}
		r+=1;
		m+=1;
	}

	return results;
}

//==========================================================================================

Matrix_Neo sequences_to_PPM (std::vector<std::string> sequences, unsigned int n)
{
	Matrix_Neo PPM;
	std::vector<double> nucleot (4, 0);
	Matrix_Neo result;
	std::string sequence;

	for (size_t i=0;i<n;i++)
	{
		result.push_back(nucleot);
		PPM.push_back(nucleot);
	}

	for(size_t i(0); i<sequences.size(); ++i)
	{
		sequence = sequences[i];
		result = sequence_to_PPM(sequence, n);

		for(size_t i(0); i < PPM.size(); ++i)
		{
			for (size_t j(0);j<4;++j)
			{
				PPM[i][j] += result[i][j];
			}
		}
	}


	return PPM;
}

//==========================================================================================

double sequence_score (std::string sequence, Matrix_Neo PPM)
{
	double score=1.0;

	for(size_t i= 0; i < sequence.size(); ++i)
	{
		switch (charmap[sequence[i]]) {
            case nuc::A:
                score = score * (PPM[i][0]);
                break;
            case nuc::C:
                score = score * (PPM[i][1]);
                break;
            case nuc::G:
                score = score * (PPM[i][2]);
                break;
            case nuc::T:
                score = score * (PPM[i][3]);
                break;
            case nuc::N:
				break;
			}
	}
    return score;
}

//==========================================================================================

std::vector<std::string> PPM_to_Sequence (std::vector<std::string> binding_sites, Matrix_Neo PPM,  double cutoff)
{
	std::vector<std::string> output;
	std::string sequence;
	std::string short_sequence;
	double score=0.0;
	int PPMsize = PPM.size();

	for(size_t i(0); i<binding_sites.size(); ++i)
	{
		sequence = binding_sites[i];

		for(size_t j(0); j<=sequence.size()-PPMsize; ++j)
		{
			short_sequence = sequence.substr(j,PPMsize);
			score = sequence_score(short_sequence, PPM);

			if (score >= cutoff)
			{
				output.push_back(short_sequence);
			}
		}
	}
	return output;
}

Matrix_Neo normalized(Matrix_Neo input){

	Matrix_Neo result;

	for (unsigned int i(0); i<input.size(); i++) {

		double line_sum = Matrix::sum_of_line(input[i]);
		std::vector<double> new_line;

		for (unsigned int j(0); j<4; j++) {
			new_line.push_back(input[i][j]/line_sum);
		}
		result.push_back(new_line);
	}

	return result;
}

//==========================================================================================

double max_score (Matrix_Neo & PPM)
{
	double result=1;

	for(size_t i(0); i < PPM.size(); ++i)
	{
		double max_line=0;

		for (size_t j(0);j<4;++j)
		{
			if (PPM[i][j] > max_line)
			{
				max_line = PPM[i][j];
			}
		}
		result = result * max_line;
	}
	return result;
}

//==========================================================================================

Matrix_Neo EM_algorithm (std::vector<std::string> Sequences, int n, double cutoff, Base_Prob base_probabilities, int max, double diff)
{
	Matrix_Neo results;
	Matrix_Neo first;
	std::vector<std::string> second;
	Matrix_Neo third;
	bool result;
	double line_sum = 0.0;
	int counter = 0;

	first = sequences_to_PPM(Sequences, n);
	result = true;
	counter = 1;

	if (cutoff<=max_score(first) and counter<max)
	{
		for(size_t i(0); i<PPM_to_Sequence(Sequences, first, cutoff).size(); ++i)
		{
			second.push_back(PPM_to_Sequence(Sequences, first, cutoff)[i]);
			third = sequences_to_PPM(second, n);
			result = false;
		}
		counter =2;
	}

	while ((cutoff<=max_score(third)) and (diff_matrices(first,third,diff)==false) and (counter<max))
	{
		second = PPM_to_Sequence(second, third, cutoff);

		for(size_t i(0); i<first.size(); ++i)
		{
			for (size_t j(0);j<4;++j)
			{
				first[i][j] = sequences_to_PPM(second, n)[i][j];
				result = true;
			}
		}

		counter += 1;

		if (cutoff<=max_score(first) and diff_matrices(first,third,diff)==false and counter<max)
		{
			second = PPM_to_Sequence(second, first, cutoff);

			for(size_t i(0); i<third.size(); ++i)
			{
				for (size_t j(0);j<4;++j)
				{
					third[i][j] = sequences_to_PPM(second, n)[i][j];
					result = false;
				}
			}
			counter += 1;
		}
	}

	if (base_probabilities[0] == 1 and
		base_probabilities[1] == 1 and
		base_probabilities[2] == 1 and
		base_probabilities[3] == 1) {
		base_probabilities = {0,0,0,0};
	}

	if (result == true)
	{
		for(size_t i(0); i < first.size(); ++i)
		{
			for (size_t j(0);j<4;++j)
			{
				first[i][j]+=base_probabilities[j];
			}
		}

		for (unsigned int i(0); i<first.size(); i++)
		{
			line_sum = Matrix::sum_of_line(first[i]);

			for (unsigned int j(0); j<4; j++)
			{
				first[i][j] = first[i][j]/line_sum;
			}
		}
		return first;

	} else {

		for(size_t i(0); i < third.size(); ++i)
		{
			for (size_t j(0);j<4;++j)
			{
				third[i][j]+=base_probabilities[j];
			}
		}

		for (unsigned int i(0); i<third.size(); i++)
		{
			line_sum = Matrix::sum_of_line(third[i]);

			for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++)
			{
				third[i][j] = third[i][j]/line_sum;
			}
		}
		return third;
	}
}

//==========================================================================================

std::vector <std::string> Initialization_string()
{
    std::string entry_name = ask_name_fasta();

    std::vector <std::string> sequence_list;

    std::vector <std::string> sequences(ExtractSequence(entry_name));
    for(size_t i(0); i<sequences.size(); ++i){
        std::string seq(sequences[i]);
        sequence_list.push_back(seq);
    }

    return sequence_list;
}

//==========================================================================================

bool diff_matrices (Matrix_Neo mat1, Matrix_Neo mat2, double difference)
{
	bool diff = true;

	for(size_t i(0); i<mat1.size() and diff; ++i)
	{
		for (size_t j(0);j<4;++j)
		{
			if ((std::fabs(mat1[i][j]-mat2[i][j])<=difference) and diff)
			{
				diff = true;

			} else if (diff) {

				diff = false;
			}
		}
	}

	return diff;
}

//==========================================================================================

int smallest_length (std::vector<std::string> sequences)
{
	unsigned int result = sequences[0].size();
	unsigned int size = sequences.size();

	for(size_t i(0); i<size; ++i)
	{
		if (sequences[i].size() < result)
		{
			result = sequences[i].size();
		}
	}

	return result;
}

//==========================================================================================

std::vector <std::string> ExtractSequence(std::string const& entry_name)
{
    std::ifstream entry(entry_name.c_str());
    std::string line;

    std::vector <std::string> sequences; // That is the object that will be returned.

    std::string intermediate_value;

    // Reads lines
    while(getline(entry, line)){
        if(line.front()=='>'){
            sequences.push_back(intermediate_value);
            intermediate_value.clear(); }
        //else if(line.front()=='\n'){}
        else {

            // Makes DNA Sequence on line all uppercase, allowing for easier search afterwards
            for (auto & c: line) c = toupper(c);
            intermediate_value+=line;
        }
    }



    // Allow to register the last value even though there's no >...
    sequences.push_back(intermediate_value);

    // Deletes the first vector created, which is a ghost one. This is just a hack cobbled together, it'd probably be best to correct that in a better way at some point.
    sequences.erase(sequences.begin());


    // Testing the values by showing them
    /*for(size_t i(0); i<sequences.size(); ++i){
     std::cout <<sequences[i] <<'\n';
     }*/

    entry.close(); // Don't you have to close it afterwards?
    return sequences;
}

//==========================================================================================

std::vector <std::string> string_list_from_searchResults(std::vector <SearchResults> input) {
    std::vector <std::string> output;

    for (unsigned int i(0); i<input.size(); i++) {
        for (unsigned int j(0); j<input[i].searchResults.size(); j++) {
            output.push_back(input[i].searchResults[j].sequence);
        }
    }

    return output;
}
