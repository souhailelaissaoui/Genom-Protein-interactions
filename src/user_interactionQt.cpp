#include <chrono>
#include <QMessageBox>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "askBaseProb.hpp"
#include "matrixFromSequence.hpp"
#include "resultsWindow.hpp"
#include "sequenceFromMatrix.hpp"
#include "user_interaction.hpp"
#include "window.hpp"

//-----------------------------------------------------------------------

double ask_cutoff2(double max_score) {
    QWidget newCut;
    bool ok;
    QString sentence = "Please enter a value between 0 and ";
    sentence +=max_score;
    sentence +=".";
    return QInputDialog::getDouble(&newCut, "Invalid Error for your cutoff", sentence,0,0,max_score,1, &ok);
}

//-----------------------------------------------------------------------

double ask_cutoff() {
  return sequenceFromMatrix::getCutoff();
}

//-----------------------------------------------------------------------

bool ask_binding_length_known() {
  return matrixFromSequence::getBool();
}

//-----------------------------------------------------------------------

unsigned int ask_position() {
  return matrixFromSequence::getPosition();
}

//-----------------------------------------------------------------------

unsigned int ask_length() {
    return matrixFromSequence::getLength();
}

//-----------------------------------------------------------------------

std::string ask_name_fasta() {
    return Window::getFastaLocation();
}

//-------------------------------------------------------------------------

bool InvalidFormatMat(std::string file_name) {
	if (file_name.find(".mat") != std::string::npos)
	{
		return 0;
	}
	else return 1;
}

//-------------------------------------------------------------------------

std::string ask_name_matrix() {
    return Window::getMatrixLocation();
}

//-----------------------------------------------------------------------

bool InvalidFormat(std::string file_name) {
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

void nucleotide_warning(char c) {
    QWidget alarm;
    QString error = "Nucleotide ";
    error += c;
    error += QString::fromStdString(" not recognized");
    QMessageBox::critical(&alarm, "Warning", error);
}

//-----------------------------------------------------------------------

void print_progress(int position, int filesize) {
    //QT Implement
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
  return Window::getOutputName();
}

//----------------------------------------------------------------------

BP_FILL CoutCin_AskBaseProb() {
    int choice = askBaseProb::getBaseChoiceFinal();

        switch (choice) {
            case 2:
                return BP_FILL::UserDefined;

            case 1:
                return BP_FILL::AllEqual;

            case 0:
                return BP_FILL::NotUsed;

            case 3:
                return BP_FILL::FromSequence;

        }
        return BP_FILL::AllEqual;
}

//----------------------------------------------------------------------

double CoutCin_AskBaseProb0(char letter) {
    if (letter == A) return askBaseProb::getProbA();
    else if (letter == C) return askBaseProb::getProbC();
    else if (letter == G) return askBaseProb::getProbG();
    else return askBaseProb::getProbT();
}

//----------------------------------------------------------------------

std::vector<double> AskBaseProb() {
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
        T = CoutCin_AskBaseProb0('T');
        }
            return {A, C, G, T};
}

//----------------------------------------------------------------------

MATRIX_TYPE Ask_Return_Matrix_Type() {
    unsigned int answer=askBaseProb::getMatrixChoice();

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
        return MATRIX_TYPE::logMatrix;
}

//----------------------------------------------------------------------

bool ask_matrix_from_sequences_weighed() {
  return sequenceFromMatrix::isEM() or matrixFromSequence::isEM();
}

//----------------------------------------------------------------------

bool ask_matrix_from_search_matches() {
  return  sequenceFromMatrix::getBool();
}

//----------------------------------------------------------------------

void error_input_sequence() {
    QWidget alarm;
    QString NoSeq = "Error: Either there was no sequence found to analyze, or the sequences to analyze are of different length.";
    QMessageBox::critical(&alarm, "Warning", NoSeq);
}

//----------------------------------------------------------------------

void error_reading_coordiates(unsigned int line) {
    QWidget alarm;
    QString Misread = "Error: Reading genomic coordinate files failed on line ";
    Misread += line;
    QMessageBox::critical(&alarm, "Warning", Misread);
}

//----------------------------------------------------------------------

std::string ask_coordinate_filename() {
    std::string answer;
    std::cout << "What file would you like to open for the genomic coordinates." << std::endl;
    std::cin >> answer;
    return answer;
}

//----------------------------------------------------------------------

bool ask_line_description_present() {
    bool answer;
    std::cout << "Is a descritpion of the coordinates present in the third question?" <<
    std::endl << "Enter 1 for yes, 0 for no." <<
    std::endl;

    while (not(std::cin >> answer)) {
        std::cout << "Error: Invalid input. Please try again." << std::endl;
    }
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
        std::endl << "Enter 1 for yes, 0 for no." <<
        std::endl;
        std::cin >> answer;
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
    unsigned int answer=matrixFromSequence::getSequenceSource();
    while (true) {
        switch (answer) {
            case 1:
                return SEQ_SOURCE::OnlySeq;

            case 2:
                return SEQ_SOURCE::FromSearchResult;
        }
    }
}

//-----------------------------------------------------------------------
// 							LOGO USER INTERACTION
//-----------------------------------------------------------------------

void logo_in_process() {
	std::cout << "Your logo is being generated, this should take a couple of seconds" << std::endl;
	std::cout << "Your logo will be saved in genom-2 with the name yourlogo.png" << std::endl;
}

std::string ask_logo_matrix() {
  return Window::getMatrixLocation();
}

void position_in_process(int pos, int size) {
    //QT TO DO
	std::cout << "Drawing position: " << pos << "/" << size << std::endl;
}

//-----------------------------------------------------------------------

double ask_for_a_number_infinitely(){
    //QT
	double n=0;
	std::string raw_input;

	while(true){

		std::cin >> raw_input;
		try{
			n = std::stod(raw_input);
			return n;

		}catch(...){
			std::cout << "Invalid input, Try again." << std::endl;
		}
	}
}

//----------------------------------------------------------------------

int ask_iterations (int length) {
    QWidget newCut;
    bool ok;
    QString sentence = "Choose a whole number between 1 and the length of the smallest sequence, which is";
    sentence +=length;
    sentence +=".";
    return QInputDialog::getInt(&newCut, "How long do you want the enzyme binding sites to be?", sentence,0,0,length,1, &ok);
}

//----------------------------------------------------------------------

void print_into_file(std::ostream & out, Matrix_Neo matrix) {
	for(unsigned int i = 0 ; i < matrix.size(); i++) {
		for(unsigned int j = 0 ; j < matrix[i].size(); j++)
			out << std::left << std::setw(16) << matrix[i][j];
		out << std::endl;
	}

	out.flush();

}

//----------------------------------------------------------------------

std::string get_working_path() {
   char temp[1024];
   return ( getcwd(temp, 1024) ? std::string( temp ) : std::string("") );
}

//----------------------------------------------------------------------

void path () {
    QWidget notice;
    QString pigeon = "We will create a file in: ";
    pigeon += QString::fromStdString(get_working_path());
    pigeon += QString::fromStdString(" with the final absolute position-probability matrix");
    QMessageBox::information(&notice, "Warning", pigeon);
    }

//----------------------------------------------------------------------

int maximum_EM () {
    //QT TO DO

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

double differences_matrices () {
    //QT TO DO
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
    unsigned int answer = matrixFromSequence::getSequenceOrigin();

        switch (answer) {
            case 1:
                return LIST_FILE::Fasta;

            case 2:
                return LIST_FILE::NormalList;

            case 3:
                return LIST_FILE::SeparatedList;
    }
    return LIST_FILE::Fasta;
}

//----------------------------------------------------------------------

char ask_separation_character() {
    //QT TO DO
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
        std::endl << " (decimal). Is this the right character? " <<
        std::endl << "Enter 1 for yes, 0 for no. ";
        std::cin >> rightanswer;
    }
    return answer;
}

//----------------------------------------------------------------------

std::string ask_inputfile_name() {
    std::cout <<std::endl <<"easd" <<std::endl;
    QWidget newFile;
    QString fileLocation = QFileDialog::getOpenFileName(&newFile, "Open a File", fileLocation, "Whichever file you want");
    return fileLocation.toStdString();
}

//----------------------------------------------------------------------

void done() {
    //QT ??
    std::cout << "DONE" << std::endl;
}

//----------------------------------------------------------------------

void error_no_search_result_read() {
    QWidget notice;
    QMessageBox::information(&notice, "Warning", "Error: No sequences were read into the program.");
}

//----------------------------------------------------------------------

bool correlate_more() {
    //QT TO DO
    std::cout << "Would you like to correlate more sequence results to genomic coordinates?" <<
    std::endl << "Enter 1 for yes, 0 for no." << std::endl;
    bool answer;
    while (not(std::cin >> answer)) {
        std::cout << "Invalid input, please try again." << std::endl;
    }
    return answer;
}

//----------------------------------------------------------------------

bool ask_correlate_to_search_results() {
  return false;
}

//----------------------------------------------------------------------

void error_invalid_flags() {
    //QT ??
    std::cout << "The program was launched with invalid flags. Aborting." <<
    std::endl << "Run the program with -help flag for more information." <<
    std::endl;
}

//----------------------------------------------------------------------

bool checkfile(std::string filename) {
    //QT ?? (FLAGS are with -a? )
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
    //QT ?? (FLAGS are with -a? )
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        bool answer;
        std::cout << "The file " << filename << " exists already. Would you like " <<
        std::endl << "to overwrite it? " <<
        std::endl << "Enter 1 for yes, 0 for no. ";
        std::cin >> answer;
        std::cout << "Your answer is " << answer << "." << std::endl;
        return answer;
    }
    return true;
}
