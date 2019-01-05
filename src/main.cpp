#include <stdio.h>
#include <string>

#include "procedures.hpp"
#include "user_interaction.hpp"
#include "logo.hpp"



/*!
* @brief   Handles everything if the user gives the necessary information
*          while calling the program
*/
void non_verbose(int argc, char **argv);


/*!
* @brief   Handles everything if the user doesn't give any information
*          when running the program.
*/
void verbose();


//-------------------------------------------------------------------------

int main(int argc, char **argv) {

  if (argc == 1) {
    verbose();
  }
  else {
    non_verbose(argc, argv);
  }


  return 0;
}


//-------------------------------------------------------------------------

void non_verbose(int argc, char **argv) {

  std::vector<std::string> flags(flags_to_strings(argc, argv));

  if (flags[1] == "--help" or flags[1] == "-h") {
    flag_help();
  }

  else if (flags[1] == "--logo" or flags[1] == "-l") {
    logo(flags[2]);
  }

  else if (flags[1] == "--bindingsites" or flags[1] == "-s") {
    flag_bindingsites(flags);
  }

  else if (flags[1] == "--getmatrix" or flags[1] == "-m") {
    flag_getmatrix(flags);
  }

  else if (flags[1] == "--about" or flags[1] == "-a") {
    flag_about();
  }

  else {
    void error_invalid_flags();
  }

}



//-------------------------------------------------------------------------

void verbose() {
  PROCEDURE procedure;



  do {
    procedure = whatToDo();

    // What to do if user wants to extract Sequences from a given matrix
    if (procedure == SequencesFromMatrix) {
      enzyme_on_sequence2();
    }

    if (procedure == MatrixFromSequences) {
      enzyme_from_sequences();
    }

    if (procedure == Logo) {
      logo();
    }

    if (procedure == CorrelateResults) {
      correlate_coordinates_with_results(input_search_results());
    }

    done();

    // Waits with program until next input of user.
    getchar();
    getchar();
    //system("CLEAR");
  } while (procedure != Exit);

}



/*! \mainpage Welcome!
*
* \section intro_sec Presenting DBSAP
*
* Welcome to the Doxygen documentation for DBSAP!
*
* This program will help you by providing the tools necessary to analyse
* a DNA sequence or a DNA-binding protein of interest.
*
* In our program you will find three useful functionalities:
*
* (1) You can obtain a probability matrix (or a weighted probability,
* logarithmic or weighted logarithmic matrix) from a file containing the
* sequences of binding sites you have found. These sequences can be of
* a fixed length or of variable length as we have implemented an
* algorithm that can identify the consensus sequence length as well as
* the sequence.
*
* needed inputs: file with sequences (.fasta, .fas, .fna, .ffn)
*
* (2) You can obtain all the binding positions of a protein to a DNA sequence
* of your choosing. You may even specify the cutoff value you prefer: 0.25,
* calculated based on the input sequence or customised by you.
* In order to be clear and precise the sequence returned is the complementary reverse branch when (-) is indicated.
*
* The needed inputs are a matrix file (.mat) and
* a sequence file (.fasta, .fas, .fna, .ffn)
*
* (3) Generate a consensus sequence logo for a DNA protein
* of your choosing: just give us the file name and we will take care of
* the rest.
*
* The needed inputs are a probability weight matrix file (.mat)
*
* (4) Our last tool correlates binding affinities between a genome file
* and your experimental values
*
* The needed inputs are your affinities, genomic file.
*
* We hope you enjoy!
*
*
*
*
* \section install_sec Running the program
*
* Our program is easy to use you'll be analysing sequences in no time! Just follow these 3 steps:
*
* STEP 1: Clone our github project
*
* --> git clone https://github.com/EPFL-SV-cpp-projects/genom-2.git
*
* STEP 2: Make a build folder and build the program
*
* --> mkdir build
*
* --> cd build
*
* --> cmake ../
*
* --> make
*
* STEP 3: Execute the program (with Qt)
*
* --> ./genom-2
*
* without Qt (--> ./Main)
*
* from the command line (--> ./Main -help)
*
* and voilà!
*
* TESTS (in build):
*
* --> ./gtests
*
* \section authors Brought to you by
*
* Angela Saudan, Erica Geneletti, Jérémy Alexandre, Katia Schalk, Lucas
* Zweili, Matthias Minder, Marion Perier, Souhail Elaissaoui,
* Tristan Vyvyan-Robinson aaaaand Gokcen Nurlu!
*
*
*/
