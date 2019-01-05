#include "procedures.hpp"


void enzyme_on_sequence() {

    //system("Clear");
    std::vector <double> base_prob(AskBaseProb());



    Matrix enzyme(ask_name_matrix(), base_prob);

	std::string file_name(Ask_Outputfile_Name());

    std::vector<SearchResults> enzyme_matches;
    enzyme_matches = analyze_sequence_opt2(ask_name_fasta(), enzyme, ask_cutoff());

    for (unsigned index(0); index < enzyme_matches.size(); index++) {
        print_results(enzyme_matches[index],file_name);
    }

    if(ask_matrix_from_search_matches()) {
        Matrix enzyme_renewed(matrix_from_same_length(enzyme_matches, base_prob, ask_matrix_from_sequences_weighed()),MATRIX_TYPE::absoluteMatrix);
        enzyme_renewed.save(file_name + ".mat",  Ask_Return_Matrix_Type());
    }

    if (ask_correlate_to_search_results()) {
        correlate_coordinates_with_results(enzyme_matches);
    }

}


void enzyme_on_sequence2() {

    //system("Clear");
    std::vector <double> base_prob(AskBaseProb());


    Matrix enzyme(ask_name_matrix(), base_prob);

    std::string fasta_name = ask_name_fasta();
    double cutoff = ask_cutoff();
    std::string output_name = Ask_Outputfile_Name();

    analyze_sequence_opt4(fasta_name, enzyme, cutoff, output_name);

    if (ask_matrix_from_search_matches()) {
        std::vector<SearchResults> enzyme_matches (read_searchresult_file(output_name + ".csv"));
        Matrix enzyme_renewed(matrix_from_same_length(enzyme_matches, base_prob, ask_matrix_from_sequences_weighed()),MATRIX_TYPE::absoluteMatrix);
        enzyme_renewed.save(output_name + ".mat",  Ask_Return_Matrix_Type());
    }

    if (ask_correlate_to_search_results()) {
        correlate_coordinates_with_results(read_searchresult_file(output_name  + ".csv"));
    }

}


void enzyme_from_sequences() {

    SEQ_SOURCE seq_source(ask_source_sequence());
    std::vector<SearchResults> seq_to_analyze;

    switch (seq_source) {
        case SEQ_SOURCE::OnlySeq:
            seq_to_analyze = seq_source_OnlySeq();
            break;

        case SEQ_SOURCE::FromSearchResult:
            seq_to_analyze = seq_source_FromSearchResult();
            break;

        default:
            break;
    }

    if (seq_to_analyze.size() == 0) {
        error_no_search_result_read();
        return;
    }



    if(searchResults_same_length(seq_to_analyze)) {
        Matrix result(binding_length_known(seq_to_analyze));
        result.save(Ask_Outputfile_Name(), Ask_Return_Matrix_Type());
    }
    else {
        binding_length_unknown(string_list_from_searchResults(seq_to_analyze), Ask_Outputfile_Name(), AskBaseProb());

    }


}


//=======================================================================

std::vector<SearchResults> input_search_results() {
    std::vector <SearchResults> result_list;
    result_list = read_searchresult_file(ask_inputfile_name());
    return result_list;
}


//=======================================================================

void correlate_coordinates_with_results(std::vector <SearchResults> result_list) {

    Association_Table associate;
    unsigned int startpos;
    unsigned int coord_id;
    unsigned int seq_id;

    std::vector <double> corr_values;

    // Initializes vector of coordinates from a file.
    std::vector <Coordinates> coord_list(read_coordinates(ask_coordinate_filename(),
                                                           ask_line_description_present()));



    do {
        associate = associate_genomic_with_result(get_descriptions_from_coord_list(coord_list),
                                                  result_list);

        coord_id = associate[0][COORD];
        seq_id = associate[0][SEQ];
        startpos = associate[0][START] - 1;


        corr_values = coord_list[coord_id].position_score(result_list[seq_id], startpos);

        print_results_correlated(result_list[seq_id], corr_values, Ask_Outputfile_Name());


    } while (correlate_more());

}




//=======================================================================

Matrix binding_length_known(std::vector<SearchResults> seq_to_analyze) {
    Base_Prob base_probabilities = AskBaseProb();

    return Matrix(matrix_from_same_length(
                                          seq_to_analyze,
                                          base_probabilities,
                                          ask_matrix_from_sequences_weighed()
                                          ),
                  MATRIX_TYPE::absoluteMatrix);
}


//=======================================================================

void binding_length_unknown(std::vector<std::string> sequence_list, std::string name, Base_Prob base_probabilities) {
	unsigned int n=0;
	double t=0.0;
	Matrix_Neo results;
	std::ofstream outputfile;
	double max = 0.0;
	int i = 0;
	int f=0;
	double g = 0.0;

    i = smallest_length(sequence_list);
    n = ask_iterations(i);
    Matrix_Neo normd = normalized(sequences_to_PPM(sequence_list,n));
    max = max_score(normd);

	t = ask_cutoff2(max);
	f = maximum_EM();
	g = differences_matrices ();

	results = EM_algorithm(sequence_list, n, t, base_probabilities, f, g);
	path();
	outputfile.open(name);
	print_into_file(outputfile, results);
	outputfile.close();
}




//=======================================================================
std::vector<SearchResults> seq_source_OnlySeq() {

    std::vector<SearchResults> output;
    std::string filename;

    switch (ask_list_file_type()) {
        case LIST_FILE::Fasta:
            filename = ask_name_fasta();
            output.push_back(read_sequencefile_to_searchresults(filename));
            break;

        case LIST_FILE::NormalList:
            filename = ask_inputfile_name();
            output.push_back(read_sequence_list_to_searchresults(filename));
            break;

        case LIST_FILE::SeparatedList:
            char delim(ask_separation_character());
            filename = ask_inputfile_name();
            output.push_back(read_char_separated_to_searchresults(filename, delim));
            break;
    }

    return output;
}

//=======================================================================

std::vector<SearchResults> seq_source_FromSearchResult() {
    return read_searchresult_file(ask_inputfile_name());
}


// NON VERBOSE
//=======================================================================

std::vector<std::string> flags_to_strings(int argc, char** argv) {
    std::vector<std::string> output;
    for ( int i(0); i<argc; i++) {
        output.push_back(std::string(argv[i]));
    }
    return output;
}


//=======================================================================

void flag_help() {
    std::cout <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "The following flags are allowed to launch the program. " <<
    std::endl << "==============================================================================" <<
    std::endl << "For help:" <<
    std::endl <<
    std::endl << "--help (or -h)" <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "For readme (program and algorithm description):" <<
    std::endl <<
    std::endl << "--about (or -a)" <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "To obtain the binding sites of an enzyme on a sequence: " <<
    std::endl <<
    std::endl << "-bindingsites matrix sequence output [cutoff]" <<
    std::endl << "(or short: -s)" <<
    std::endl <<
    std::endl << "Where: "<<
    std::endl << " matrix is the name of the file containing the matrix, " <<
    std::endl << " sequence is the name a .fasta file that contains the sequences to analyze," <<
    std::endl << " output is the name of the file where the results are saved to," <<
    std::endl << " cutoff is the minimal score a sequence has to have to be considered as match." <<
    std::endl << "Note that if no cutoff is specified, the score is set to be just above the " <<
    std::endl << "score of an all-N sequence." <<
    std::endl << "Also note that when launching from terminal, base probabilities are not " <<
    std::endl << "used. Therefore, the score of all sequences will be below 0." <<
    std::endl << "The score without base probabilities can be obtained from the score with " <<
    std::endl << "base probabilities with the following formula:" <<
    std::endl << "score' [without bp] = score [with bp] - 2 * n [length of binding site]" <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "To obtain a probability weight matrix from a list of sequences, " <<
    std::endl << "all separated by a return: " <<
    std::endl <<
    std::endl << "--getmatrix --list input output" <<
    std::endl << "(or short: -m -list)" <<
    std::endl <<
    std::endl << "Where: " <<
    std::endl << " input is the name of the file containing sequences. " <<
    std::endl << " output is the name of the file that the resulting matrix is saved to." <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "To obtain a probability weight matrix from a list of sequences, " <<
    std::endl << "present in a .fasta file: " <<
    std::endl <<
    std::endl << "--getmatrix --fasta input output" <<
    std::endl << "(or short: -m -fas)" <<
    std::endl <<
    std::endl << "Where:" <<
    std::endl << " input is the name of the file containing sequences. " <<
    std::endl << " output is the name of the file that the resulting matrix is saved to." <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "To obtain a probability weight matrix from the result of a previous " <<
    std::endl << "analyisis (in a .csv file): " <<
    std::endl <<
    std::endl << "--getmatrix --result input output " <<
    std::endl << "(or short: -m -res)" <<
    std::endl <<
    std::endl << "Where:" <<
    std::endl << " input is the name of the file containing sequences. " <<
    std::endl << " output is the name of the file that the resulting matrix is saved to." <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "The following functionalities only exist in the GUI or Q&A edition: " <<
    std::endl << " - Correlation of search results of a sequence analysis to genomic coordinates" <<
    std::endl << " - Creating a probability weight matrix based on sequences in a file that are " <<
    std::endl << "   separated by a character other than a return " <<
    std::endl << " - Creating a probability weight matrix where each sequence is weighted by its " <<
    std::endl << "   score." <<
    std::endl << " - Saving matrices in a form other than as probability weight matrix. " <<
    std::endl <<
    std::endl <<
    std::endl;


}


void flag_about() {
    std::cout <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << "==============================================================================" <<
    std::endl << "    HOW THE PROGRAM WORKS: " <<
    std::endl << "==============================================================================" <<
    std::endl <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "    PROCEDURE: GET MATRIX FROM SEQUENCES" <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "If the user wants to calculate a probability weight matrix based on a list of " <<
    std::endl << "sequences, this is done in the following fashion:" <<
    std::endl << "If all sequences are of the same length, the amount that each nucleotide is at " <<
    std::endl << "the specific position is counted in a matrix. In order to avoid a probability " <<
    std::endl << "of zero (which will result in problems in later analysis), the user can choose " <<
    std::endl << "to use nucleotide base probabilities to create a background noise. The " <<
    std::endl << "probability of each nucleotide is added to the corrsponding matrix column. " <<
    std::endl << "Finally, each row is divided by its sum, resulting in a probability weight " <<
    std::endl << "matrix. " <<
    std::endl << "The user can choose to weigh every sequence by its score. If he does this, the " <<
    std::endl << "program adds the score of a specific sequence at the position of each " <<
    std::endl << "nucleotide (instead of just 1, as described above). Then, the program iterates " <<
    std::endl << "through the resulting matrix and subtracts the minimum of each line from its " <<
    std::endl << "values if it is below zero. This ensures that the final probabilities will be " <<
    std::endl << "between 0 and 1. As described above, the user then can choose to create some " <<
    std::endl << "background noise using base probabilities. Finally, the program divides every " <<
    std::endl << "line by its sum, yielding the output probability weight matrix." <<
    std::endl <<
    std::endl << "If the user wants to calculate aprobability weight matrix based on a list of " <<
    std::endl << "sequences with different lengths, the program will use the Expectation " <<
    std::endl << "Maximization Algorithm (short EM-Algorithm) to compute a matrix. This " <<
    std::endl << "algorithm is described in: " <<
    std::endl << "Do, C. B., & Batzoglou, S. (2008). What is the expectation maximization " <<
    std::endl << "algorithm? Nature Biotechnology, 26(8), 897-899. " <<
    std::endl << "The user can determine when to stop the algorithm with the following  " <<
    std::endl << "conditions:" <<
    std::endl << "He can set the maximum amount of iterations of the algorithm." <<
    std::endl << "He can set the maximal difference between the two iterating matrices, if the " <<
    std::endl << "difference is smaller than the given value, the program stops." <<
    std::endl << "He also chooses the cutoff of the algorithm. " <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "    PROCEDURE: GET ENZYME BINDING POSITIONS ON SEQUENCE" <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "In order to find all enzyme binding positions on a sequence, the user has to " <<
    std::endl << "provide a .fasta file containing the sequence, and a .mat file containing the " <<
    std::endl << "probability weight matrix. The matrix is converted to a logarithmic matrix " <<
    std::endl << "(if necessary), using base probabilities if wanted by the user. The program " <<
    std::endl << "will then iterate through all positions of the sequence and calculate " <<
    std::endl << "a score corresponding to the sum of the logarithmic probabilities of the " <<
    std::endl << "nucleotides in the given position. A sequence is considered as a binding " <<
    std::endl << "position if its score is bigger or equal as a user defined cutoff. " <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "    PROCEDURE: LOGO" <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "This creates a png visualization of a probability weight matrix. It is saved " <<
    std::endl << "as yourlogo.png in the folder above." <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "    PROCEDURE: CORRELATE RESULTS TO BINDING POSITIONS" <<
    std::endl << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<
    std::endl << "It is possible for the user to compare the binding positions found in previous" <<
    std::endl << "analysis to a genomic coordinate file. To every binding position, the genomic " <<
    std::endl << "score of its position is added. This additional score is printed as additional" <<
    std::endl << "column of the output .csv file. " <<
    std::endl <<
    std::endl <<
    std::endl <<
    std::endl;
}

//=======================================================================

void flag_bindingsites(std::vector<std::string> flags) {
    if (flags.size() < 5) {
        error_invalid_flags();
        return;
    }

    if (not(checkfile(flags[2]))) {
        flags[2] += ".mat";
        if (not(checkfile(flags[2]))) {
            return;
        }
    }

    if (not(checkfile(flags[3]))) {
        flags[3] += ".fasta";
        if (not(checkfile(flags[3]))) {
            return;
        }
    }

    if (not(overwrite(flags[4] + ".csv"))) {
        return;
    }

    std::vector <SearchResults> output;
    // Syntax: first matrix, then sequence file, then outputfile, then cutoff,
    Matrix enzyme(flags[2]);

    // If no cutoff is given, take score of all-N sequence as cutoff
    if (flags.size() < 6) {
        double score(enzyme.get_length() * enzyme.get_N_score() + 0.001);
        flags.push_back(std::to_string(score));
    }

    analyze_sequence_opt4(flags[3], enzyme, std::stod(flags[5]), flags[4]);
}


//=======================================================================

void flag_getmatrix(std::vector<std::string> flags) {
    if (flags.size() < 5) {
        error_invalid_flags();
        return;
    }

    if (not(checkfile(flags[3]))) {
        return;
    }

    // Checks if output file is already present
    if (not(overwrite(flags[4] + ".mat"))) {
        return;
    }
    if (not(overwrite(flags[4]))) {
        return;
    }



    std::vector <SearchResults> input;
    if (flags[2] == "--fasta" or flags[2] == "-fas") {
        input.push_back(read_sequencefile_to_searchresults(flags[3]));
    }
    else if (flags[2] == "--list" or flags[2] == "-list") {
        input.push_back(read_sequence_list_to_searchresults(flags[3]));
    }
    else if (flags[2] == "--result" or flags[2] == "-res") {
        input = read_searchresult_file(flags[3]);
    }


    if(searchResults_same_length(input)) {
        Matrix result(binding_length_known(input));
        result.save(flags[4], MATRIX_TYPE::absoluteMatrix);
    }
    else {
        binding_length_unknown(string_list_from_searchResults(input), flags[4]);
    }
}
