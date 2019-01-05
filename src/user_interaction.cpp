#include "user_interaction.hpp"
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <chrono>

//-----------------------------------------------------------------------



PROCEDURE whatToDo() {

    //system("Clear");
    PROCEDURE whatToDo;
    unsigned int answer;

    std::cout << "Welcome to our program. What would you like to do? \n";
    std::cout << "Enter 1 to obtain a probability weight matrix from multiple" << std::endl
    << "sequences." << std::endl;
    std::cout << "Enter 2 to obtain all binding positions of a protein from a " << std::endl
    << "probability matrix on a given nucleotide sequence." << std::endl;
    std::cout << "Enter 3 to obtain a logo." << std::endl;
    std::cout << "Enter 4 to compare the score of sequences to the score " <<
    std::endl << "in genomic coordinates." << std::endl;
    std::cout << "Enter 0 to exit the program." << std::endl;
    std::cout << "This program can be launched directly from terminal, for more information " <<
    std::endl << "run ./Main --help" << std::endl;

    while (true) {
		answer = ask_for_a_number_infinitely();

        if (answer == 1) {
            whatToDo = MatrixFromSequences;
            break;
        }
        else if (answer == 2) {
            whatToDo = SequencesFromMatrix;
            break;
        }
        else if (answer == 3) {
            whatToDo = Logo;
            break;
        }
        else if (answer == 4) {
            whatToDo = CorrelateResults;
            break;
        }
        else if (answer == 0) {
            whatToDo = Exit;
            break;
        }  
           std::cout<<"Unrecognized input. Please try again." << std::endl;
    }

    return whatToDo;
}

//-----------------------------------------------------------------------

double ask_cutoff2(double max_score) {


	double cutoff = 0.0;

	do{
		std::cout << "What cutoff would you like to use?" << std::endl;
		std::cout << "Choose a number bigger than 0 and smaller than " << max_score << std::endl;

		cutoff = ask_for_a_number_infinitely();

		if((cutoff<max_score) and (cutoff>0)){
			break; //success
		}
		std::cout << "Invalid input, please try again." << std::endl;

	}while(true);

	return cutoff;
}

//-----------------------------------------------------------------------

double ask_cutoff() {
	
    std::cout << "What cutoff would you like to use?" << std::endl;

    return ask_for_a_number_infinitely();
}

//-----------------------------------------------------------------------

std::string ask_name_fasta()
{
    std::string entry_name;

    while (true) {
        std::cout <<"Please give the name of your sequence file: ";
        std::cin >> entry_name;

        std::ifstream entry(entry_name.c_str());

        if (entry.fail()) {
            std::cout << "Impossible to read the file, please try again. " << std::endl;
            continue;
        }

        if(InvalidFormat(entry_name)) {
            std::cout << "Unknown file format, please try again. " << std::endl;
            continue;
        }

        entry.close();
        break;
    }

    return entry_name;
}

//-------------------------------------------------------------------------

bool InvalidFormatMat(std::string file_name)
{
	if (file_name.find(".mat") != std::string::npos)
	{
		return 0;
	}
	else return 1;
}

//-------------------------------------------------------------------------
// There's a better function for this

std::string ask_name_matrix()
{
    std::string entry_name;

    while (true) {
	    std::cout <<"Please give the name of your matrix file: ";
	    std::cin >> entry_name;

	    std::ifstream entry(entry_name.c_str());


	    if (entry.fail()) {
	        std::cout << "Impossible to read the file, please try again." << std::endl;
	        continue;
	    }

	    if (InvalidFormatMat(entry_name)) {
			std::cout << "Unknown file format, please try again." << std::endl;
			continue;
		}

	    entry.close();
	    break;
	}

    return entry_name;
}

//-----------------------------------------------------------------------

bool InvalidFormat(std::string file_name)
{

    // Defines list with known file formats
    static const std::vector<std::string> validValues {".fasta", ".fas", ".fna", ".ffn", ".fa"};

    for(unsigned int i = 0; i < validValues.size(); i++) {
        if(file_name.find(validValues[i])!= std::string::npos)
            return 0;
        // Returns 0 if the file extension can be found
    }

    return 1;
}

//-----------------------------------------------------------------------

void nucleotide_warning(char c)
{
	std::cout << "WARNING, nucleotide " << c
                          << " not recognized" << std::endl;
}

//-----------------------------------------------------------------------

void print_progress(int position, int filesize) {


    static auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start); // Elapsed is in seconds

    static int nextPrint(0);
    static int increment(filesize/500);
	int barWidth = 50;


    if(position >= nextPrint) {

        double progress (std::abs((double)position / (double)filesize));
        nextPrint += increment;

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
        std::cout << "] " << std::fixed << std::setprecision(2) << (progress * 100.0) << "% "<< std::setprecision(0) << elapsed.count() / progress * (1-progress) << "s left\r" << std::flush;

	}

}

void ret() {
	std::cout << std::endl;
}

//==========================================================================================

void print_results(SearchResults results, std::string filename) {
    filename += ".csv";
    std::ofstream outputfile;
    outputfile.open(filename, std::ios_base::app);
    unsigned int size = results.searchResults.size();

    outputfile << std::endl;
    outputfile << ">" << results.description << std::endl;
    outputfile << "Seq; Pos; Dir; Seq_Score" << std::endl;

    for (unsigned int i(0); i < size; i++) {
        outputfile << results.searchResults[i].sequence << "; "
                   << results.searchResults[i].position << "; "
                   << results.searchResults[i].direction << "; "
                   << results.searchResults[i].score << std::endl;
    }
    outputfile.close();
}

//==========================================================================================

void print_results2(SearchResults results, std::ofstream &outputfile) {
    unsigned int size = results.searchResults.size();

    for (unsigned int i(0); i < size; i++) {
        outputfile << results.searchResults[i].sequence << "; "
        << results.searchResults[i].position << "; "
        << results.searchResults[i].direction << "; "
        << results.searchResults[i].score << std::endl;
    }
}

//==========================================================================================

void print_results_correlated(SearchResults results, std::vector <double> gen_score, std::string filename) {
    filename += ".csv";
    std::ofstream outputfile;
    outputfile.open(filename, std::ios_base::app);
    unsigned int size = results.searchResults.size();

    outputfile << std::endl;
    outputfile << ">" << results.description << std::endl;
    outputfile << "Seq; Pos; Dir; Seq_Score; Gen_Score" << std::endl;


    for (unsigned int i(0); i < size; i++) {
        outputfile  << results.searchResults[i].sequence << "; "
                    << results.searchResults[i].position << "; "
                    << results.searchResults[i].direction << "; "
                    << results.searchResults[i].score << "; "
                    << gen_score[i] << std::endl;
    }
    outputfile.close();
}

//==========================================================================================

std::string Ask_Outputfile_Name() {
	std::string filename;

	std::cout <<"What would you like to call your outputfile? ";
	std::cin>>filename;

	return filename;
}

//----------------------------------------------------------------------

BP_FILL CoutCin_AskBaseProb()
{
	int choice;

    while (true) {
        std::cout << "What would you like to use as the base probabilities for each type of nucleotide in your sequence? " <<
        std::endl << "The base probabilities are used to normalize the logarithmic probability weight matrix. " <<
        std::endl << "Enter 0 to not use any base probabilities"<<
        std::endl << "Enter 1 to use a base probability of 0.25 for all nucleotides"<<
        std::endl << "Enter 2 to use custom base probabilities" <<

        std::endl;

        choice = (int)ask_for_a_number_infinitely();

        switch (choice) {
            case 2:
                return BP_FILL::UserDefined;

            case 1:
                return BP_FILL::AllEqual;

            case 0:
                return BP_FILL::NotUsed;


            case 3:
                return BP_FILL::FromSequence;


            default:
                std::cout << "Invalid input, try again. " << std::endl;
        }

    }
     return BP_FILL::NotUsed;
}

//----------------------------------------------------------------------

double CoutCin_AskBaseProb0(char C)
{
	double baseProb;

    while (true) {
        std::cout << "Enter the base probability for "<< C << " " ;
        baseProb = ask_for_a_number_infinitely();

        if(baseProb < 0 || baseProb > 1) {
            std::cout << "Invalid value for base probability. All values "
            << std::endl << "for base probabilities must be between 0 and 1 "
            << std::endl << "and must add up to one."
            << std::endl << "Please try again." << std::endl;
        }
        else
            return baseProb;
    }

	return baseProb;
}

//----------------------------------------------------------------------

std::vector<double> AskBaseProb()
{
    switch (CoutCin_AskBaseProb()) {
        case UserDefined:
            return User_Defined_Base_Prob();

        case AllEqual:
            return {0.25, 0.25, 0.25, 0.25};

        default:
            return {1, 1, 1, 1};
    }
}

//----------------------------------------------------------------------

std::vector<double> User_Defined_Base_Prob() {
    double A, C, G, T;
    while(true) {

        A = CoutCin_AskBaseProb0('A');
        C = CoutCin_AskBaseProb0('C');
        G = CoutCin_AskBaseProb0('G');

        T = 1 - A - C - G;

        std::cout << "The value for T is automatically set to " << T << std::endl;
        if(T < 0 || T > 1) {
            std::cout << "Invalid inputs for the base probabilities. " << std::endl
            << "The base probabilities must add up to 1. The sum of the inserted " << std::endl
            << "base probabilities must therefore be smaller than 1." << std::endl;
        }
        else
            return {A, C, G, T};

    }

}

//----------------------------------------------------------------------

MATRIX_TYPE Ask_Return_Matrix_Type() {
    unsigned int answer;

    while (true) {
        std::cout << "As what kind of matrix do you want to save your matrix as?" <<
        std::endl << "Enter 1 to save as a probability matrix." <<
        std::endl << "Enter 2 to save as a weighted probability matrix." <<
        std::endl << "Enter 3 to save as a logarithmic matrix." <<
        std::endl << "Enter 4 to save as a weighted logarithmic matrix." <<
        std::endl;

        answer = (unsigned int)ask_for_a_number_infinitely();

        switch (answer) {
            case 1:
                return MATRIX_TYPE::absoluteMatrix;
            case 2:
                return MATRIX_TYPE::relativeMatrix;
            case 3:
                return MATRIX_TYPE::logMatrix;
            case 4:
                return MATRIX_TYPE::logConstMatrix;
        }

        std::cout << "Invalid input. Please try again. " << std::endl;
    }
    return MATRIX_TYPE::logMatrix;
}

//----------------------------------------------------------------------

bool ask_matrix_from_sequences_weighed() {
    bool answer;
    std::cout << "Would you like to weigh each probability by its score " <<
    std::endl << "in order to calculate the probability weight matrix?" <<
    std::endl;

	answer = correct_bool();
	
    std::cout << "You answered ";
    if (answer) {
		std::cout <<"YES";
	}
	else {
		std::cout <<"NO";
    }

    std::cout << std::endl;
    return answer;
}

//----------------------------------------------------------------------

bool ask_matrix_from_search_matches() {
    bool answer;
	
    std::cout << "Would you like to create a new matrix based on the search matches?" <<
    std::endl;

	answer = correct_bool();

    std::cout << "You answered ";
    if (answer) {
		std::cout <<"YES";
    }
	else {
		std::cout <<"NO";
    }

    std::cout << std::endl;
    return answer;
}

//----------------------------------------------------------------------

void error_input_sequence() {
    std::cout << "Error: Either there was no sequence found to analyze, or the " <<
    std::endl << "sequences to analyze are of different length." << std::endl;
}

//----------------------------------------------------------------------

void error_reading_coordiates(unsigned int line) {
    std::cout << "Error: Reading genomic coordinate files failed on line " << line <<
    std::endl;
}

//----------------------------------------------------------------------

std::string ask_coordinate_filename() {
    std::string answer;
    std::ifstream file; 
    
    do {
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "What file would you like to open for the genomic coordinates." << std::endl;
		std::cin >> answer;
		file.open(answer);
		if (file.fail())
		{
			std::cerr << "this file doesn't exist, please try again." << std::endl;
		}
		
	} while (file.fail());
	file.close();
    return answer;
}

//----------------------------------------------------------------------

bool ask_line_description_present() {
    bool answer;
    std::cout << "Is a descritpion of the coordinates present in the third question?" <<
    std::endl;

    answer = correct_bool();

    return answer;
}

//----------------------------------------------------------------------

Association_Table associate_genomic_with_sequences(std::vector<std::string> coordinate_description,
                                                   std::vector<std::string> sequence_description) {
    Association_Table output;
    unsigned int coord_intermed;
    unsigned int seq_intermed;
    unsigned int startpos_intermed;

    if (coordinate_description.size() == 1 and sequence_description.size() == 1) {
        while (true) {
            std::cout << "At what position of the genomic coordinates does the sequence " <<
            std::endl << "start?" << std::endl;

            if (std::cin >> startpos_intermed) {
                break;
            }
            std::cout << "Error, please try again." << std::endl;

        }
        return {{0, 0, startpos_intermed}};
    }

    std::cout << "The following coordinates could be found in the given files: " << std::endl;
    for (unsigned int id(0); id<coordinate_description.size(); id++) {
        std::cout << "No. "<< id+1 << "\t" << coordinate_description[id] << std::endl;
    }

    std::cout << "The following sequences could be found in the given files: " << std::endl;
    for (unsigned int id(0); id<sequence_description.size(); id++) {
        std::cout << "No. "<< id+1 << "\t" << sequence_description[id] << std::endl;
    }

    if (coordinate_description.size() == sequence_description.size()) {
        bool answer;
        std::cout << "You have as many sequences as genomic coordinates. Would you like to " <<
        std::endl << "analyze them in order?" <<
        std::endl;
        
        answer = correct_bool();
	
        if (answer) {
            std::cout << "For each sequence, enter its starting position in the genomic " <<
            std::endl << "coordinates. " <<
            std::endl;
            for (unsigned int id(0); id<coordinate_description.size(); id++) {
                while (not(std::cin >> startpos_intermed)) {
                    std::cout << "Error for position " << id << ". Please try again." << std::endl;
                }
                output.push_back({id, id, startpos_intermed});
            }
            return output;
        }
    }

    std::cout << "Which sequences would you like to analyze with which sequences?" <<
    std::endl << "Associate sequences with genomic coordinates by entering first the " <<
    std::endl << "number of the sequence and then the number of the genomic coordinates." <<
    std::endl << "You can combine them any way you'd like. " <<
    std::endl << "As a third value, give the starting position of the sequence in the " <<
    std::endl << "genomic coordinates." <<
    std::endl << "As soon as you are done, enter 0." <<
    std::endl;

    while (true) {
        std::cout << "Analysis " << output.size() + 1 << ":" << std::endl;
        std::cout << "Sequence ";
        while(not(std::cin >> seq_intermed))
            std::cout << "Invalid input, please try again." << std::endl;

        if(seq_intermed == 0)
            break;

        std::cout << "Coordinates ";
        while(not(std::cin >> coord_intermed))
            std::cout << "Invalid input, please try again." << std::endl;

        std::cout << "Start position of coordinates ";
        while(not(std::cin >> startpos_intermed))
            std::cout << "Invalid input, please try again." << std::endl;

        output.push_back({coord_intermed-1, seq_intermed-1, startpos_intermed});
    }

    return output;


}

//----------------------------------------------------------------------
// DOESNT WORK

Association_Table associate_genomic_with_result(std::vector<std::string> coordinate_description,
                                                std::vector<SearchResults> result_list) {
    Association_Table output;
    unsigned int coord_intermed;
    unsigned int seq_intermed;
    unsigned int startpos_intermed;

    if (coordinate_description.size() == 1 and result_list.size() == 1) {
        while (true) {
            std::cout << "At what position of the genomic coordinates is the position " <<
            std::endl << "0 of the input sequences?" <<
            std::endl << "If there is no shift, this value is 1." << std::endl;

            if (std::cin >> startpos_intermed) {
                break;
            }
            std::cout << "Error, please try again." << std::endl;

        }
        return {{0, 0, startpos_intermed}};
    }

    std::cout << "The following coordinates could be found in the given files: " << std::endl;
    for (unsigned int id(0); id<coordinate_description.size(); id++) {
        std::cout << "No. "<< id+1 << "\t" << coordinate_description[id] << std::endl;
    }

    std::cout << "The following sequence results could be found in the given files: " << std::endl;
    for (unsigned int id(0); id<result_list.size(); id++) {
        std::cout << "No. "<< id+1 << "\t" << result_list[id].description << std::endl;
    }


    std::cout << "Which sequence result would you like to analyze with which genomic coordinate?" <<
    std::endl << "Enter first the identification number of the sequence results followed by the " <<
    std::endl << "number of the genomic coordinates." <<
    std::endl << "As a third value, give the position on the genomic coordinates corresponding " <<
    std::endl << "to position 0 of the sequence results." <<
    std::endl << "If there is no shift, this value is 1." << std::endl;


    std::cout << "Sequence result: ";
    while(not(std::cin >> seq_intermed)
          and seq_intermed<=0
          and seq_intermed>result_list.size())
        std::cout << "Invalid input, please try again." << std::endl;

    std::cout << "Coordinates: ";
    while(not(std::cin >> coord_intermed)
          and coord_intermed<=0
          and coord_intermed>coordinate_description.size())
        std::cout << "Invalid input, please try again." << std::endl;

    std::cout << "Start position of coordinates: ";
    while(not(std::cin >> startpos_intermed))
        std::cout << "Invalid input, please try again." << std::endl;

    output.push_back({coord_intermed-1, seq_intermed-1, startpos_intermed});

    return output;


}

//----------------------------------------------------------------------

void error_sequence_doesnt_exist() {
    std::cout << "Error: The desired sequence doesn't exist" << std::endl;
}

//----------------------------------------------------------------------

SEQ_SOURCE ask_source_sequence() {
	
	unsigned int answer(0);
	
	do{
		std::cout << "How would you like to obtain the binding sequences to analyze?" <<
		std::endl << "Enter 1 if you have the binding sequences as a list in a file." <<
		std::endl << "Enter 2 if you want to analyze the result of a previous analysis." <<
		std::endl;
		
		answer = (unsigned int)ask_for_a_number_infinitely();
		
		switch (answer) {
            case 1:
                return SEQ_SOURCE::OnlySeq;

            case 2:
                return SEQ_SOURCE::FromSearchResult;
        }
        if ((answer!=1) or (answer!=2))
        {
			std::cout << "Unrecognized input, Please try again." << std::endl;	
		}
	
	}while(true);
}


//-----------------------------------------------------------------------
// 							LOGO USER INTERACTION
//-----------------------------------------------------------------------

void logo_in_process()
{
	std::cout << "Your logo is being generated, this should take a couple of seconds" << std::endl;
	std::cout << "Your logo will be saved in genom-2 with the name yourlogo.png" << std::endl;
}

std::string ask_logo_matrix()
{
    std::string entry_name;

    while (true) {
	    std::cout <<"Please give the name of your matrix file: ";
	    std::cin >> entry_name;

	    std::ifstream entry(entry_name.c_str());


	    if (entry.fail()) {
	        std::cout << "Impossible to read the file, please try again." << std::endl;
	        continue;
	    }

	    entry.close();
	    break;
	}

    return entry_name;
}

void position_in_process(int pos, int size)
{
	std::cout << "Drawing position: " << pos << "/" << size << std::endl;
}

//-----------------------------------------------------------------------

double ask_for_a_number_infinitely(){
	double n=0;
	std::string raw_input;

	while(true){

		std::cin >> raw_input;
		try{
			n = std::stod(raw_input);
			return n;

		}catch(...){
			std::cout << "Unrecognized input, Please try again." << std::endl;
		}
	}
}

//----------------------------------------------------------------------

int ask_iterations (int length)
{
	int n;

	do{
		std::cout << "How long do you want the enzyme binding sites to be?" << std::endl;
		std::cout << "Choose a whole number between 1 and the length of the smallest sequence, which is " << length << std::endl;

		n = (int)ask_for_a_number_infinitely();
		if((n<length) and (n>0)){
			break; //success
		}
		std::cout << "Invalid input, please try again." << std::endl;

	}while(true);

	return n;
}

//----------------------------------------------------------------------

void print_into_file(std::ostream & out, Matrix_Neo matrix){

	for(unsigned int i = 0 ; i < matrix.size(); i++) {
		for(unsigned int j = 0 ; j < matrix[i].size(); j++)
			out << std::left << std::setw(16) << matrix[i][j];
		out << std::endl;
	}

	out.flush();

}

//----------------------------------------------------------------------

std::string get_working_path()
{
   char temp[1024];
   return ( getcwd(temp, 1024) ? std::string( temp ) : std::string("") );
}

//----------------------------------------------------------------------

void path ()
{
	std::cout << "We will create a file in: " << get_working_path() << " with the final absolute position-probability matrix" << std::endl;
}

//----------------------------------------------------------------------

int maximum_EM ()
{

	int n = 0;

	do{
		std::cout << "What's the maximum of times do you want to do the EM_algorithm?" << std::endl;
		std::cout << "Choose a whole number bigger than 0" << std::endl;

		n = (int)ask_for_a_number_infinitely();

		if(n>0){
			break; //success
		}
		std::cout << "Invalid input, please try again." << std::endl;

	}while(true);

	return n;
}

//----------------------------------------------------------------------

double differences_matrices ()
{
	double n = 0.0;

	do{
		std::cout << "What is the difference between the two matrices when doing the EM_algorithm you chose to stop the EM_algorithm ?" << std::endl;
		std::cout << "Choose a number bigger than 0" << std::endl;

		n = ask_for_a_number_infinitely();

		if(n>0){
			break; //success
		}
		std::cout << "Invalid input, please try again." << std::endl;

	}while(true);

	return n;

}

//----------------------------------------------------------------------

LIST_FILE ask_list_file_type() {
    std::cout << "What does your list file look like?" <<
    std::endl << "Enter 1 if it is in a .fasta file, with each list separated " <<
    std::endl << "by descriptions beginning with the caracter >" <<
    std::endl << "Enter 2 if it is a simple list, with each sequence separated by returns." <<
    std::endl << "Enter 3 if the sequences are separated by a character other than a return." <<
    std::endl;

    unsigned int answer(0);
    while (true) {
        answer = (unsigned int)ask_for_a_number_infinitely();
        
        switch (answer) {
            case 1:
                return LIST_FILE::Fasta;

            case 2:
                return LIST_FILE::NormalList;

            case 3:
                return LIST_FILE::SeparatedList;

            default:
                std::cout << "Invalid input, please try again." << std::endl;
                break;

        }
    }

}



//----------------------------------------------------------------------

char ask_separation_character() {
    char answer;
    std::cout << "Please make sure that the delimiting character present right " <<
    std::endl << "after each sequence, and that no character is present after the last "<<
    std::endl << "sequence." <<
    std::endl << "By which character are the sequences separated?" <<
    std::endl;
    bool rightanswer(0);

    while (not(rightanswer)) {
        std::cin.get();
        std::cin.get(answer);
        std::cout <<
        std::endl << "You entered the character '" << answer << "' with ASCII number " << (int) answer <<
        std::endl << " (decimal). Is this the right character? " << std::endl;
  
        rightanswer = correct_bool();	
	
    }
    return answer;
}


//----------------------------------------------------------------------

std::string ask_inputfile_name() {
    std::string entry_name;

    while (true) {
        std::cout <<"Please give the name of your sequence file: ";
        std::cin >> entry_name;

        std::ifstream entry(entry_name.c_str());

        if (entry.fail()) {
            std::cout << "Impossible to read the file, please try again. " << std::endl;
            continue;
        }

        entry.close();
        break;
    }

    return entry_name;
}

//----------------------------------------------------------------------

void done() {
    std::cout << "DONE" << std::endl;
}

//----------------------------------------------------------------------

void error_no_search_result_read() {
    std::cout << "Error: No sequences were read into the program." << std::endl;
}

//----------------------------------------------------------------------

bool correlate_more() {
    std::cout << "Would you like to correlate more sequence results to genomic coordinates?" <<std::endl;
    bool answer;
    
    answer = correct_bool();
	
    return answer;
}

//----------------------------------------------------------------------

bool ask_correlate_to_search_results() {
    std::cout << "Would you like to correlate the found scores to a genomic coordinate file?" <<
	std::endl;
    bool answer;
    
    answer = correct_bool();

    return answer;
}

//----------------------------------------------------------------------

void error_invalid_flags() {
    std::cout << "The program was launched with invalid flags. Aborting." <<
    std::endl << "Run the program with -help flag for more information." <<
    std::endl;
}

//----------------------------------------------------------------------

bool checkfile(std::string filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }

    std::cout << "The file " << filename << " could not be found. Aborting." <<
    std::endl;
    return false;
}

//----------------------------------------------------------------------

bool overwrite(std::string filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        bool answer;
        std::cout << "The file " << filename << " exists already. Would you like " <<
        std::endl << "to overwrite it? " << std::endl;
     
         answer = correct_bool();
        
        std::cout << "Your answer is " << answer << "." << std::endl;
        return answer;
    }
    return true;
}

bool correct_bool() {
	bool answer; 
	do
		{	std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			std::cout << "Enter 1 for yes, 0 for no." <<
			std::endl;
			std::cin >> answer;
			if (std::cin.fail())
			{
				std::cout << "1 or 0 are the only valid inputs." << std::endl;
			}
			
		} while (std::cin.fail());
	return answer;
}
