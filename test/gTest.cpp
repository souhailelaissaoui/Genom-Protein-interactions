#include "gtest/gtest.h"
#include "../src/Matrix.hpp"
#include "../src/utility.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

//=========== Elements we use for Matrix gTests ==============

Matrix_Neo lm_1({{-0.440569,MINUSINFINI,-1.247930,1.506960},
			     {-0.440569,MINUSINFINI,-1.247930,1.506960}});

Matrix_Neo lm_2({{-2.440569,MINUSINFINI,-3.247930,1.506960 -2},
			   {-2.440569,MINUSINFINI,-3.247930,1.506960 -2}});

Matrix_Neo a_1({{0.184211,0.000,0.105263, 0.710526},
			  {0.184211,0.000,0.105263, 0.710526}});

Matrix_Neo r_1({{0.259260,0.000000,0.148148,1.000},
				{0.259260,0.000000,0.148148,1.000}});


Matrix_Neo lcm_1({{-1.947529,MINUSINFINI,-2.754889,0.000},
				{-1.947529,MINUSINFINI,-2.754889,0.000}});

const std::vector<double> base_probabilities ({0.25,0.25,0.25,0.25});

Matrix_Neo r({{0.259260,0.000000,0.148148,1.000}});
Matrix ma_matrice_1(r, relativeMatrix, base_probabilities);

Matrix_Neo a({{0.184211,0.000,0.105263, 0.710526}});
Matrix ma_matrice_2(a, absoluteMatrix, base_probabilities);

Matrix_Neo lm({{-0.440569,MINUSINFINI,-1.247930,1.506960}});
Matrix ma_matrice_3(lm,logMatrix, base_probabilities);

Matrix_Neo lcm({{-1.947529,MINUSINFINI,-2.754889,0.000}});
Matrix ma_matrice_4(lcm,logConstMatrix, base_probabilities);

Matrix_Neo wrong({{1.506960,0.105263,0.148148,-2.754889}});

Matrix_Neo DBP_PPM{ {0.991265586410457  , 0.00188241672188422, 0.00438979579543401, 0.00246220107222457},
					 {0.00454038913318475, 0.00308716342389013, 0.961116800192759  , 0.0312556472501656},
					 {0.00319260279955123, 0.991529060968172  ,0.00188243089596181 , 0.0033959053363151},
					 {0.891738449490994  ,0.00188241672188422 , 0.00883982892596831, 0.0975393048611529},
					 {0.624171736642371  ,0.251091801698693   ,0.00595596650804169 ,0.118780495150895},
					 {0.728886814047346  ,0.10454942473345	  ,0.0115881573399193  ,0.154975603879284},
					 {0.304793385940606  ,0.491122522739594	  ,0.00446509246430938 ,0.199618998855491}   };

Matrix_Neo DBP_PSSM {{1.98734355068112	,-7.05319824318173	,-5.83163045464332	,-6.66583570234573},
					{-5.78296833594373	,-6.33950242833837	,1.94278367099148	,-2.99973931124021},
					{-6.29105121022664	,1.98772696248885	,-7.05318738012775	,-6.20198804153446},
					{1.83469252898669	,-7.05319824318173	,-4.82176583473376	,-1.35787249910015},
					{1.32001493662631	,0.00628682961322058,-5.39144864337102	,-1.07363014308654},
					{1.54376670664635	,-1.25774297113760	,-4.4312050115156	,-0.689886969274964},
					{0.285903501598009	,0.97415489074836	,-5.80709423329781	,-0.324679058329434} };

Matrix ma_matrice_5 (DBP_PPM, absoluteMatrix, base_probabilities);

Matrix ma_matrice_A	(a_1,absoluteMatrix, base_probabilities);

Matrix ma_matrice_R	(r_1,relativeMatrix, base_probabilities);

Matrix ma_matrice_LOG (lm_1,logMatrix,base_probabilities);

Matrix test_matrix(DBP_PSSM, MATRIX_TYPE::logMatrix);

//========================	Matrix_gTests		========================

/*!
 *@brief Function testing if max_of_line(std::vector<double> line) returns as expected the max value of the line.
 */
TEST (max_of_line_Test, Good_MaxValue)
{
	std::vector<double> line_to_be_tested_1 ({0.184211,0.000,0.105263, 0.710526});
	double answer_1(0.710526);

	ASSERT_EQ(answer_1,Matrix::max_of_line(line_to_be_tested_1));
}


/*!
 *@brief Function testing if min_of_line(std::vector<double> line) returns as expected the min value of the line.
 */
TEST (min_of_line_Test, Good_MinValue)
{
	std::vector<double> line_to_be_tested_2 ({0.259260,0.000000,0.148148,1.000});
	double answer_2(0.000000);

	ASSERT_EQ(answer_2,Matrix::min_of_line(line_to_be_tested_2));
}


/*!
 *@brief Function testing if sum_of_line(std::vector<double> line) returns the expected value.
 */

TEST (sum_of_line_Test, Good_Sum)
{
	std::vector<double> line_to_be_tested_3 ({0.259260,0.000000,0.148148,1.000});
	double answer_3(1.407408);

	ASSERT_EQ(answer_3, Matrix::sum_of_line(line_to_be_tested_3));
}


/*!
 *@brief Function testing if line_is_reg_ppm(double min, double max, double sum) returns the great bool
 */
TEST (line_is_reg_ppm_Test, Good_Bool)
{
	ASSERT_TRUE(ma_matrice_A.line_is_reg_ppm(1,0.5,1));
}

/*!
 *@brief Function testing if line_is_normed_ppm(double min, double max, double sum)  returns the good bool
 */
TEST (line_is_normed_ppm_Test, Good_Bool)
{

	ASSERT_TRUE(ma_matrice_A.line_is_normed_ppm(1,1));

}


/*!
 *@brief Function testing if logMatrix_from_logConstMatrix( Matrix_Neo input_matrix )returns the good logMatrix
 */
TEST (logMatrix_from_logConstMatrix_Test, Good_logMatrix)
{

	Matrix ma_matrice_2(lcm_1, logConstMatrix);

	Matrix_Neo result2 = ma_matrice_2.logMatrix_from_logConstMatrix(lcm_1);

	for (int i=0 ; i<2 ; ++i)
		{
			for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(lm_2[i][j] - result2[i][j]) <0.001);
		}

		}
}

/*!
 *@brief Function testing if logMatrix_from_normed_PPM( Matrix_Neo input_matrix ) returns the good logMatrix
 */

TEST (logMatrix_from_normed_PPM_Test, Good_logMatrix)
{

	Matrix_Neo result_3 = ma_matrice_R. logMatrix_from_normed_PPM(r_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(lm_1[i][j] - result_3[i][j]) <0.001);
		}

	}
}


/*!
 *@brief Function testing if logMatrix_from_probMatrix( Matrix_Neo input_matrix ) returns the good logMatrix
 */

TEST (logMatrix_from_probMatrix_Test, Good_logMatrix)
{

	Matrix_Neo result_4 = ma_matrice_A. logMatrix_from_probMatrix(a_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(lm_1[i][j] - result_4[i][j]) <0.001);
		}

	}
}

/*!
 *@brief Function testing if probMatrix_from_logMatrix( Matrix_Neo input_matrix )returns the good absoluteMatrix
 */

TEST (probMatrix_from_logMatrix_Test, Good_absoluteMatrix)
{

	Matrix_Neo result_5 = ma_matrice_A.probMatrix_from_logMatrix(lm_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(a_1[i][j] - result_5[i][j]) <0.001);
		}

	}
}



/*!
 *@brief Function testing if absolute_from_normed_PPM( Matrix_Neo input_matrix  )returns the good absoluteMatrix
 */

TEST (absolute_from_normed_PPM_Test, Good_absoluteMatrix)
{

	Matrix_Neo result_6 = ma_matrice_R.absolute_from_normed_PPM(r_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(a_1[i][j] - result_6[i][j]) <0.001);
		}

	}
}

/*!
 *@brief Function testing if determine_matrix_type(Matrix_Neo input) returns the good MatrixType
 */
TEST (determine_matrix_type_Test, Good_MatrixType)
{

	ASSERT_EQ(MATRIX_TYPE::absoluteMatrix, ma_matrice_A.determine_matrix_type(a_1) );

}



/*!
 *@brief Function testing if log_to_matrix returns the good Matrix
 */
TEST (log_to_matrix_Test, Good_Matrix)
{

	Matrix_Neo result_15 = ma_matrice_LOG.log_to_matrix(lm_1,MATRIX_TYPE::absoluteMatrix);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(a_1[i][j] - result_15[i][j]) <0.001);
		}

	}
}

/*!
 *@brief Function testing if matrix_to_log(Matrix_Neo input_matrix, MATRIX_TYPE type)returns the good logMatrix
 */
TEST (matrix_to_log_Test, Good_Matrix)
{


	Matrix_Neo result_15 = ma_matrice_A.matrix_to_log(a_1, MATRIX_TYPE::absoluteMatrix);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(lm_1[i][j] - result_15[i][j]) <0.001);
		}

	}
}


/*!
 *@brief Function testing if logConstMatrix_from_logMatrix( Matrix_Neo input_matrix ) returns the good logConstMatrix
 */
TEST (logConstMatrix_from_logMatrix_Test, GoodlogConstMatrix)
{
	Matrix ma_matrice_1(lm_1, logMatrix);
	Matrix_Neo result_1 = ma_matrice_1.logConstMatrix_from_logMatrix(lm_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(lcm_1[i][j] - result_1[i][j]) <0.05);
		}
	}
}



/*!
 *@brief Function testing if normed_from_absolute_PPM( Matrix_Neo input_matrix  )returns the good reltiveMatrix
 */
TEST (normed_from_absolute_PPM_Test, Good_relativeMatrix)
{

	Matrix_Neo result_9 = ma_matrice_A.normed_from_absolute_PPM(a_1);


	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(r_1[i][j] - result_9[i][j]) <0.001);
		}

	}
}

/*!
 *@brief Function testing if normed_PPM_from_logMatrix( Matrix_Neo input_matrix ) returns the good reltiveMatrix
 */

TEST (normed_PPM_from_logMatrix_Test, Good_relativeMatrix)
{

	Matrix_Neo result_13 = ma_matrice_LOG.normed_PPM_from_logMatrix( lm_1);

	for (int i=0 ; i<2 ; ++i)
	{
		for(int j = 0 ; j< 4 ;++j)
		{
			ASSERT_TRUE(std::abs(r_1[i][j] - result_13[i][j]) <0.001);
		}

	}
}


/*!
 *@brief Testing the getters funtions get_length(), get_base_probabilities() and get_logMatrix()
 */
TEST (getters_functions, good_initialisation)
{
	for(size_t i(0) ; i<4 ; ++i)
	{
		ASSERT_EQ(base_probabilities[i], ma_matrice_1.get_base_probabilities()[i]);
		ASSERT_EQ(base_probabilities[i], ma_matrice_2.get_base_probabilities()[i]);
		ASSERT_EQ(base_probabilities[i], ma_matrice_3.get_base_probabilities()[i]);
		ASSERT_EQ(base_probabilities[i], ma_matrice_4.get_base_probabilities()[i]);

		ASSERT_EQ(lm.size(), ma_matrice_1.get_length());
		ASSERT_EQ(lm.size(), ma_matrice_2.get_length());
		ASSERT_EQ(lm.size(), ma_matrice_3.get_length());
		ASSERT_EQ(lm.size(), ma_matrice_4.get_length());

		ASSERT_EQ(lm[0][i], ma_matrice_3.get_logMatrix()[0][i]);
	};

}

/*!
 *@brief Testing if the matrix reading is done correctly
 */
TEST (read_matrix_test, good_reading_file)
{
	std::string outputfile("../test/DBP_PPM.mat");

	Matrix_Neo tmp = Matrix::read_matrix(outputfile);

	for (size_t i(0) ; i<DBP_PPM.size() ; ++i)
	{
		for ( size_t j(0) ; j<DBP_PPM[i].size() ; ++j)

		ASSERT_EQ(DBP_PPM[i][j], tmp[i][j]); /// really not sure

	};
}

/*!
 *@brief Testing the the "score" of a possible nucleotide binding site sequence_score(std::list<nuc> sequence);
 */
TEST (sequence_score_test, good_score)
{
	double our_score(0.0);
	std::list<nuc> sequence {nuc::A,nuc::A,nuc::A,nuc::A,nuc::A,nuc::A,nuc::A};
	double calc_score(0.0);



	for(size_t i(0) ; i < DBP_PSSM.size() ; ++i)
	{
		our_score += DBP_PSSM[i][0];

	}
	calc_score =  test_matrix.sequence_score(sequence);
	/*std::cout << our_score << std::endl;
	std::cout << calc_score << std::endl;*/
	ASSERT_TRUE(std::abs(-5.10229832163 - calc_score) < 0.001);
}

//=========== Elements we use for utility gTests ===========================

SearchResult tmp;
std::list<nuc> sequence = {A,G,A,A,A,T,C};
bool y = false;
std::vector<SearchResult> tmp_1;
std::vector<SearchResults> what_we_use;
SearchResult t_1;
SearchResult t_2;
SearchResult t_3;
SearchResult t_4;
SearchResults t_5;


//======================	 Utility_gTests		==========================

/*!
 *@brief Function testing if the function called sequence_to_PPM returns the expected value
 */
TEST (sequence_to_PPM_Test, PositivePPM)
{
	std::string sequence = "ATCCAG";
	int n = 3;
	Matrix_Neo result ({{1.0,2.0,0.0,1.0},{1.0,2.0,0.0,1.0},{1.0,2.0,1.0,0.0}});
	ASSERT_EQ(result,sequence_to_PPM(sequence, n));
}

/*!
 *@brief Function testing if the function called sequences_to_PPM returns the expected value
 */
TEST (sequences_to_PPM_Test, PositivePPM)
{
	std::vector<std::string> sequences ({{"ATCCAG", "AATGTCG", "TCCGTAAG"}});
	int n = 3;
	Matrix_Neo result ({{4.0,4.0,2.0,5.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	ASSERT_EQ(result,sequences_to_PPM(sequences, n));
}

/*!
 *@brief Function testing if the function called sequence_score returns the expected value
 */
TEST (sequence_score_Test, Positivescore)
{
	std::string sequence ("TCG");
	Matrix_Neo PPM ({{4.0,4.0,2.0,5.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	double result = 125.0;
	ASSERT_EQ(result,sequence_score(sequence, PPM));
}

/*!
 *@brief Function testing if the function called PPM_to_Sequence returns the expected value
 */
TEST (PPM_to_Sequence_Test, PositiveSequences)
{
	std::vector<std::string> binding_sites ({{"ATCCAG", "AATGTCG", "TCCGTAAG"}});
	Matrix_Neo PPM ({{4.0,4.0,2.0,5.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	double cutoff = 100;
	std::vector<std::string> result ({{"TCC", "TCG", "TCC", "CCG"}});
	ASSERT_EQ(result,PPM_to_Sequence(binding_sites, PPM, cutoff));
}

/*!
 *@brief Function testing if the function called max_score returns the expected value
 */
TEST (max_score_Test, positivemaxscore)
{
	Matrix_Neo PPM ({{4.0,4.0,2.0,5.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	double result = 125;
	ASSERT_EQ(result,max_score(PPM));
}

/*!
 *@brief Function testing if the function called EM_algorithm returns the expected value
 */
TEST (EM_algorithm_Test, PositiveEM_algorithm)
{
	std::vector<std::string> sequences ({{"ATCCAG", "AATGTCG", "TCCGTAAG"}});
	int n =3;
	double cutoff = 100.0;
	Base_Prob baseprob {0.25,0.25,0.25,0.25};
	int max = 2;
	double diff = 1.0;
	Matrix_Neo result ({ { 0.05, 0.25, 0.05, 0.65 }, { 0.05, 0.85, 0.05, 0.05 }, { 0.05, 0.45, 0.45, 0.05 } });
	ASSERT_EQ(result,EM_algorithm(sequences, n, cutoff, baseprob, max, diff));
}

/*!
 *@brief Function testing if "diff_matrices" returns the expected value
 */
TEST (diff_matrices_Test, Good_bool)
{
	double difference = 1.0;
	Matrix_Neo mat1 ({{2.9,4.0,2.0,3.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	Matrix_Neo mat2 ({{4.0,4.0,2.0,3.0},{4.0,5.0,2.0,4.0},{3.0,4.0,5.0,3.0}});
	bool result = false;
	ASSERT_EQ(result,diff_matrices(mat1, mat2, difference));
}

/*!
 *@brief Function testing if "smallest_length" returns the expected value
 */
TEST (smallest_length_Test, smallest_length)
{
	std::vector<std::string> sequences ({{"AACGTGACG"},{"ACTG"},{"GGCAGTTTGA"}});
	int result = 4;
	ASSERT_EQ(result,smallest_length(sequences));
}

/*
 *@brief Testing if analyze_sequence_opt2 returns the right output
 */
TEST ( analyze_sequence_opt2_test,valid_output)
{
	std::vector <SearchResult> tmp1({tmp});

	SearchResults tmp3;

	tmp3.description = "doesn't mater";
	tmp3.searchResults = tmp1;

	std::vector <SearchResults> results({tmp3});
	std::vector <SearchResults> expected_value = analyze_sequence_opt2("../test/promoters.fasta", ma_matrice_5, 0.1);


	if (
	results[0].searchResults[0].sequence == expected_value[0].searchResults[0].sequence  and
	results[0].searchResults[0].position == expected_value[0].searchResults[0].position  and
	std::abs(results[0].searchResults[0].score - expected_value[0].searchResults[0].score) <0.01 and
	results[0].searchResults[0].direction== expected_value[0].searchResults[0].direction )
	{
		y = true;
	}
	ASSERT_TRUE(y);
}

/*!
 *@brief Testing if file_fill_search_result creates a good SearchResult
 */
TEST ( file_fill_search_result_test,valid_SearchResult)
{
	SearchResult test = fill_search_result(sequence, 64, 1.07805,'-');
	if (
		tmp.direction == test.direction and
		tmp.position == test.position   and
		tmp.sequence == test.sequence   and
		std::abs(tmp.score-test.score) < 0.01)
		{
			y = true;
			}
	ASSERT_TRUE(y);
}

/*!
 *@brief Testing if sequence_string_from_nuc_list converts a list of nuc in a string
 */
TEST (sequence_string_from_nuc_list_test, valid_nuc_conversion)
 {
	 std::list<nuc> sequence = {A,G,A,A,A,T,C};
	 ASSERT_TRUE(tmp.sequence==sequence_string_from_nuc_list(sequence));

 }

 /*!
 *@brief Testing if the character is valid
 */
TEST (valid_character_test, valid_char)
{
	char test = 'f';
	ASSERT_FALSE(valid_character(test));


}

/*!
 *@brief Testing if filesize returns the good number of characters in a file
 */
TEST (filesize_test, good_nb)
{
	ASSERT_EQ(864, filesize("../test/promoters.fasta"));
}

/*!
 *@brief Testing if extract_sequence_descriptions returns the good description
 */
TEST (extract_sequence_descriptions_test, good_description)
{
	std::vector<std::string> descrip({">chr7|chr7:113842207-113842607", ">chr11|chr11:16380767-16381167"});
	std::vector<std::string> test = extract_sequence_descriptions("../test/promoters.fasta");

	for (size_t i = 0; i < test.size() ; ++i)
	{
		ASSERT_EQ(descrip[i], test[i]);
	}

}

/*!
 *@brief Testing if searchResults_same_length returns the number of nucleotides for the enzyme
 */
TEST (searchResults_same_length_test, good_number)
{
	std::vector <SearchResult> tmp1({tmp});

	SearchResults tmp3;

	tmp3.description = "doesn't mater";
	tmp3.searchResults = tmp1;

	std::vector <SearchResults> results({tmp3});

	ASSERT_EQ(7, searchResults_same_length(results));
}

/*!
 *@brief Testing if matrix_from_same_length_sequences_not_weighted returns the right matrix
 */
TEST (matrix_from_same_length_sequences_not_weighted_test, good_matrix)
{
	Matrix_Neo a = matrix_from_same_length_sequences_not_weighted(what_we_use,  base_probabilities);
	Matrix_Neo b { {0.65, 0.05, 0.25, 0.05},
				   {0.45, 0.05, 0.05, 0.45}	};

	for ( size_t i=0 ; i < a.size() ; ++i)
	{
		for ( size_t j=0 ; j < a[i].size() ; ++j)
		{
			ASSERT_EQ(a[i][j], b[i][j]);
		}
	}
}

/*!
 *@brief Testing if matrix_from_same_length_sequences_weighted returns the right matrix
 */
TEST (matrix_from_same_length_sequences_weighted_test, good_matrix)
{
	Matrix_Neo a = matrix_from_same_length_sequences_weighted(what_we_use,  base_probabilities);
	Matrix_Neo b { {0.55, 0.05, 0.35, 0.05},
				   {0.55, 0.05, 0.05, 0.35}	};

	for ( size_t i=0 ; i < a.size() ; ++i)
	{
		for ( size_t j=0 ; j < a[i].size() ; ++j)
		{
			ASSERT_EQ(a[i][j], b[i][j]);
		}
	}
}


int main(int argc, char **argv) {

	tmp.sequence = "AGAAATC";
	tmp.position = 64;
	tmp.score = 1.07805;
	tmp.direction = '-';
	t_1.sequence = "AT";
	t_1.score = 0.5;

	tmp_1.push_back(t_1);
	t_2.sequence = "AA";
	t_2.score = 1;
	tmp_1.push_back(t_2);
	t_3.sequence = "AT";
	t_3.score = 1;
	tmp_1.push_back(t_3);
	t_4.sequence = "GA";
	t_4.score = 1.5;
	tmp_1.push_back(t_4);
	t_5.searchResults = tmp_1;
	what_we_use.push_back(t_5);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
