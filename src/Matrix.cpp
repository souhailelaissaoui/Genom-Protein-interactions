#include "Matrix.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <assert.h>


#define TOLERANCE 0.02
#define APPROX_ZERO 1.0E-300

//-----------------------------------------------------------------------
// CONSTRUCTORS
//-----------------------------------------------------------------------

Matrix::Matrix(const std::string& fileName,
               const std::vector<double> base_prob)
:base_prob(base_prob)
{
    Matrix_Neo input_matrix ( read_matrix(fileName) );
    MATRIX_TYPE type;

    type = determine_matrix_type(input_matrix);

    logMatrix = matrix_to_log( input_matrix, type );
    length = logMatrix.size();
    
    N_score = log2(1/(double)sum_of_line(base_prob));
}

//-----------------------------------------------------------------------

Matrix::Matrix(Matrix_Neo input_matrix,
               MATRIX_TYPE type,
               const std::vector<double> base_prob)
:base_prob(base_prob)
{
    logMatrix = matrix_to_log( input_matrix, type );
    length = logMatrix.size();

    N_score = log2(1/(double)sum_of_line(base_prob));
}

//-----------------------------------------------------------------------
// GETTER FUNCTIONS
//-----------------------------------------------------------------------

unsigned int Matrix::get_length() {
    return length;
}

//-----------------------------------------------------------------------

std::vector<double> Matrix::get_base_probabilities() {
    return base_prob;
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::get_logMatrix() {
    return logMatrix;
}

//-----------------------------------------------------------------------

double Matrix::get_N_score() {
    return N_score;
}

//-----------------------------------------------------------------------
// COMPUTING FUNCTIONS
//-----------------------------------------------------------------------

double Matrix::sequence_score(std::list<nuc> sequence) {
    double score(0);
    std::list<nuc>::iterator iterator;
    unsigned int index(0);
    for (iterator = sequence.begin(); iterator != sequence.end(); iterator++) {
        if (*iterator < N) {
            score += logMatrix[index][*iterator];
        }
        if (*iterator == N) {
            score += N_score;
        }
        index++;
    }

    return score;
}


//-----------------------------------------------------------------------
// SAVING MATRIX
//-----------------------------------------------------------------------


void Matrix::save(std::string fileName, MATRIX_TYPE type) {
    Matrix_Neo outputmatrix(log_to_matrix(logMatrix, type));

    std::ofstream outputfile;
    outputfile.open(fileName);

    unsigned int size(outputmatrix.size());

    for (unsigned int i(0); i<size; i++) {
        for (unsigned int j(0); j<NUMBER_NUCLEOTIDES; j++) {
            outputfile << outputmatrix[i][j] << "\t";
        }
        outputfile << "\n";
    }
    outputfile.close();
}


//-----------------------------------------------------------------------
// READING MATRIX
//-----------------------------------------------------------------------

Matrix_Neo Matrix::read_matrix(std::string const& fileName){

    Matrix_Neo input_matrix;

    std::ifstream PWM;
    PWM.open(fileName);

    //send an error if there is a problem opening file
    if (PWM.fail()) {

        throw std::string("Error: Cannot read PWM file");

    } else {

        //create all variables to be use later
        std::string line;
        double A, T, C, G;
        std::vector< double > rowi(4);

        //push back of rows:
        while (true) {

            //(1) make a sting containing ith line
            getline(PWM, line);

            //check if at end of file now
            if (PWM.eof()) break;

            //(2) copy line into new stream
            std::istringstream values(line);

            //(3) read values and copy into variables
            values >> A >> C >> G >> T;

            //(4) make a matrix to pushback for line 1
            rowi.clear();
            rowi.push_back(A);
            rowi.push_back(C);
            rowi.push_back(G);
            rowi.push_back(T);


            //(5) Pushback the new row
            input_matrix.push_back(rowi);


        }
        PWM.close();

        return input_matrix;
    }
}

//-----------------------------------------------------------------------
// DETERMINE TYPE
//-----------------------------------------------------------------------

MATRIX_TYPE Matrix::determine_matrix_type(Matrix_Neo input) {

    unsigned int size(input.size());

    unsigned int id(0);

    double line_sum;
    double minimum_of_line;
    double maximum_of_line;

    // Checks if it is a probability weight matrix
    while(id<size) {
        line_sum = sum_of_line(input[id]);
        minimum_of_line = min_of_line(input[id]);
        maximum_of_line = max_of_line(input[id]);

        if(line_is_reg_ppm(minimum_of_line, maximum_of_line, line_sum))
            return MATRIX_TYPE::absoluteMatrix;

        if(line_is_normed_ppm(minimum_of_line, maximum_of_line))//, line_sum
            return MATRIX_TYPE::relativeMatrix;

        id++;
    }

    Matrix_Neo converted(probMatrix_from_logMatrix( input ));


    id = 0;

    // Checks if it is a probability weight matrix
    while(id<size) {
        line_sum = sum_of_line(converted[id]);
        minimum_of_line = min_of_line(converted[id]);
        maximum_of_line = max_of_line(converted[id]);

        if(line_is_reg_ppm(minimum_of_line, maximum_of_line, line_sum))
            return MATRIX_TYPE::logMatrix;

        if(line_is_normed_ppm(minimum_of_line, maximum_of_line))//, line_sum
            return MATRIX_TYPE::logConstMatrix;

        id++;
    }

    return MATRIX_TYPE::ERROR;
}

//-----------------------------------------------------------------------

double Matrix::max_of_line(std::vector<double> line) {
    std::vector<double>::iterator maximum_of_line;
    maximum_of_line = std::max_element(line.begin(), line.end());
    return *maximum_of_line;
}

//-----------------------------------------------------------------------

double Matrix::min_of_line(std::vector<double> line) {
    std::vector<double>::iterator minimum_of_line;
    minimum_of_line = std::min_element(line.begin(), line.end());
    return *minimum_of_line;
}

//-----------------------------------------------------------------------

double Matrix::sum_of_line(std::vector<double> line) {
    double sum(0.0);
    unsigned int size(line.size());

    for (unsigned int index(0); index < size; index++) {
        sum += line[index];
    }
    return sum;
}

//-----------------------------------------------------------------------
// A line is considered as part of a weighed position probability
// matrix if the sum is 1 and all values are between 0 and 1

bool Matrix::line_is_reg_ppm(double min, double max, double sum) {
    return ((sum < 1 + TOLERANCE)
            and (sum > 1 - TOLERANCE)
            and (max <= 1)
            and (min >= 0));
}


//-----------------------------------------------------------------------
// A line is considered as part of a weighed position probability
// matrix if the maximum is 1 and all values are between 0 and 1

bool Matrix::line_is_normed_ppm(double min, double max) { // the following argument : "double sum" was not used and created warnings - Souhail
    return ((min > 0)
            and (max >= 1 - TOLERANCE)
            and (min <= 1 + TOLERANCE));
}

//-----------------------------------------------------------------------
// MATRIX CONVERSIONS
//-----------------------------------------------------------------------

// TODO: Maybe change the default
Matrix_Neo Matrix::matrix_to_log(Matrix_Neo input_matrix, MATRIX_TYPE type) {

    Matrix_Neo outputmatrix(logMatrix_from_probMatrix(input_matrix));

    switch (type) {
        case MATRIX_TYPE::logMatrix:
            return input_matrix;

        case MATRIX_TYPE::absoluteMatrix:
            return logMatrix_from_probMatrix( input_matrix );

        case MATRIX_TYPE::relativeMatrix:
            return logMatrix_from_normed_PPM( input_matrix );

        case MATRIX_TYPE::logConstMatrix:
            return Matrix::logMatrix_from_logConstMatrix( input_matrix );
            
        default:
            return input_matrix;
    }

}

//-----------------------------------------------------------------------

// TODO: Maybe change the default
Matrix_Neo Matrix::log_to_matrix(Matrix_Neo input_matrix, MATRIX_TYPE type) {
    switch (type) {
        case MATRIX_TYPE::logMatrix:
            return input_matrix;

        case MATRIX_TYPE::absoluteMatrix:
            return probMatrix_from_logMatrix( input_matrix );

        case MATRIX_TYPE::relativeMatrix:
            return Matrix::normed_PPM_from_logMatrix( input_matrix );
            
        case MATRIX_TYPE::logConstMatrix:
            return Matrix::logConstMatrix_from_logMatrix( input_matrix );
            
        default:
            return input_matrix;
    }

}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::logMatrix_from_probMatrix( Matrix_Neo input_matrix )
{

    Matrix_Neo logMatrix;

    std::vector<double> new_line;

    double intermediate;
    size_t size(input_matrix.size());

    for(size_t i(0);i < size; ++i)
    {
        new_line.clear();

        for (size_t j(0);j< NUMBER_NUCLEOTIDES; ++j)/*Read the absoluteMatrix*/
        {

            intermediate = input_matrix[i][j];

            if(intermediate < APPROX_ZERO)/* Solve the problem of -infini case*/
            {
                new_line.push_back(MINUSINFINI);
            } else {
                new_line.push_back(log2(intermediate / base_prob[j]));
            }
        }

        logMatrix.push_back(new_line); /*Stock the line of 4 new values in the logMatrix*/

    }
    return logMatrix;
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::absolute_from_normed_PPM( Matrix_Neo input_matrix  )
{

    Matrix_Neo absoluteMatrix;

    std::vector<double> new_line;

    double line_sum;
    size_t size(input_matrix.size());

    for(size_t i(0);i < size; ++i)
    {
        new_line.clear();
        line_sum = sum_of_line(input_matrix[i]);
        
      
        //assert(line_sum >= 1);
        
        for (size_t j(0);j < NUMBER_NUCLEOTIDES; ++j)/*Read the absoluteMatrix*/
        {
            new_line.push_back(input_matrix[i][j]/line_sum);
        }

        absoluteMatrix.push_back(new_line); /*Stock the line of 4 new values in the logMatrix*/
    }
    return absoluteMatrix;
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::logMatrix_from_normed_PPM( Matrix_Neo input_matrix ) {
    Matrix_Neo absoluteMatrix;
    absoluteMatrix = ( absolute_from_normed_PPM( input_matrix ));
    return logMatrix_from_probMatrix( absoluteMatrix );
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::probMatrix_from_logMatrix( Matrix_Neo input_matrix)
{

    Matrix_Neo probability_matrix;

    std::vector<double> new_line;
    double intermediate;

    size_t size(input_matrix.size());

    for(size_t i(0);i < size; ++i)
    {
        new_line.clear();

        for (size_t j(0);j<NUMBER_NUCLEOTIDES;++j)/*Read the logMatrix*/
        {
            intermediate=input_matrix[i][j];
            
            if(intermediate<=MINUSINFINI)// Solve the problem of -infini cases
            {
                new_line.push_back(0.0);

            } else {
                new_line.push_back(base_prob[j]*std::pow(2, intermediate));/*Calcul the new values we need*/
            }
        }

        probability_matrix.push_back(new_line);
    }
    return probability_matrix;
}

Matrix_Neo Matrix::probMatrix_from_logMatrix()
{
	return this->probMatrix_from_logMatrix(this->logMatrix);
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::logMatrix_from_logConstMatrix( Matrix_Neo input_matrix )
{
	
    Matrix_Neo normed_prob_matrix( Matrix::probMatrix_from_logMatrix(input_matrix));
    return logMatrix_from_normed_PPM(normed_prob_matrix);
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::normed_from_absolute_PPM( Matrix_Neo input_matrix  )
{
    Matrix_Neo normedMatrix;

    std::vector<double> new_line;

    double line_sum;
    size_t size(input_matrix.size());

    for(size_t i(0);i < size; ++i)
    {
        new_line.clear();
        line_sum = max_of_line(input_matrix[i]);

        for (size_t j(0);j < NUMBER_NUCLEOTIDES; ++j)/*Read the absoluteMatrix*/
        {
            new_line.push_back(input_matrix[i][j]/line_sum);
        }

        normedMatrix.push_back(new_line); /*Stock the line of 4 new values in the logMatrix*/
    }
    return normedMatrix;
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::logConstMatrix_from_logMatrix( Matrix_Neo input_matrix )
{
    Matrix_Neo absoluteMatrix(probMatrix_from_logMatrix(input_matrix));
    Matrix_Neo normedMatrix(Matrix::normed_from_absolute_PPM(absoluteMatrix));
    return logMatrix_from_probMatrix(normedMatrix);
}

//-----------------------------------------------------------------------

Matrix_Neo Matrix::normed_PPM_from_logMatrix( Matrix_Neo input_matrix )
{
    Matrix_Neo absoluteMatrix(probMatrix_from_logMatrix(input_matrix));
    return normed_from_absolute_PPM(absoluteMatrix);

}
