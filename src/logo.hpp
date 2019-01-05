#ifndef LOGO_HPP
#define LOGO_HPP

#define cimg_use_png
#define cimg_display 0
#include "../logo/CImg-1.7.9/CImg.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "user_interaction.hpp"
#include <algorithm>
#include <cmath>

using namespace cimg_library;

//==================== structure used by logo() =======================
/*!
 * @class nuc_prob_pair
 *
 * @brief A nucleotide and its probability at a corresponding position.
 */
typedef struct nuc_prob_pair{
    double prob;
    CImg<unsigned char> icon;
    int index;
    nuc_prob_pair(CImg<unsigned char> icon, double prob) : prob(prob), icon(icon){

	}
}nuc_prob_pair;

//==================== functions used by logo() =======================

/*!
     * @brief   reads a matrix from file and returns the probability matrix associated
     *
     * @param   file name
     *
     * @return  the probability matrix
     */
std::vector<std::vector<double>> read_logo_matrix(std::string const& fileName)
{
	Matrix mat1(fileName);
	std::vector<std::vector<double>> result(mat1.probMatrix_from_logMatrix());

    return result;
}

/*!
     * @brief   what size to give a logo letter
     *
     * @param   probability matrix, nucleotide and position
     *
     * @return  size of the letter
     */
double size(std::vector<std::vector<double>> const& PWM, unsigned N, unsigned pos)
{
	double info_content(2);

	for (unsigned i(0); i<4; i++)
	{
		info_content += ((PWM[pos][i])*(std::log2(PWM[pos][i])));
	}

	return (PWM[pos][N])*info_content*500;
}

/*!
     * @brief   compares the type nuc_prob_matrix based on probability
     *
     * @param   the two nuc_prob_pairs to compare
     *
     * @return  true if first param bigger than second
     */
bool compareByProb(const nuc_prob_pair &a, const nuc_prob_pair &b)
{
    return a.prob > b.prob;
}

//======================= logo main functions ========================

/*!
     * @brief   creates the logo, calls on all functions needed
     */
void logo() {

	std::string fileName(ask_logo_matrix());
	std::vector<std::vector<double>> PWM(read_logo_matrix(fileName));

	unsigned size_motif(PWM.size());

	CImg<unsigned char> background(1, 1, 1, 3, 255, 255, 255);
	background.resize(size_motif*500+250, 1400);

	const unsigned char black[] = { 0,0,0 };
	background.draw_line(170,1205,size_motif*500+203,1205,black);
	background.draw_line(170,1206,size_motif*500+203,1206,black);
	background.draw_line(170,1207,size_motif*500+203,1207,black);
	background.draw_line(170,1208,size_motif*500+203,1208,black);
	background.draw_line(170,1209,size_motif*500+203,1209,black);
	background.draw_line(170,1210,size_motif*500+203,1210,black);
	background.draw_line(170,1211,size_motif*500+203,1211,black);

	background.draw_line(195,1230,195,200,black);
	background.draw_line(194,1230,194,200,black);
	background.draw_line(193,1230,193,200,black);
	background.draw_line(192,1230,192,200,black);
	background.draw_line(191,1230,191,200,black);
	background.draw_line(190,1230,190,200,black);
	background.draw_line(189,1230,189,200,black);

	background.draw_line(189,700,170,700,black);
	background.draw_line(189,701,170,701,black);
	background.draw_line(189,702,170,702,black);
	background.draw_line(189,703,170,703,black);
	background.draw_line(189,699,170,699,black);
	background.draw_line(189,698,170,698,black);
	background.draw_line(189,697,170,697,black);

	background.draw_line(189,200,170,200,black);
	background.draw_line(189,201,170,201,black);
	background.draw_line(189,202,170,202,black);
	background.draw_line(189,203,170,203,black);
	background.draw_line(189,204,170,204,black);
	background.draw_line(189,205,170,205,black);
	background.draw_line(189,206,170,206,black);

	background.draw_text(60,50,"bits", black, 0, 1, 103);
	background.draw_text(size_motif*250+80,1290,"position", black, 0, 1, 103);
	background.draw_text(100,210,"2", black, 0, 1, 53);
	background.draw_text(100,680,"1", black, 0, 1, 53);
	background.draw_text(100,1185,"0", black, 0, 1, 53);
	background.draw_text(189,1250,"5'", black, 0, 1, 53);
	background.draw_text(size_motif*500+190,1250,"3'", black, 0, 1, 53);

	logo_in_process();


	for (unsigned pos(0); pos<size_motif; ++pos)
	{
		position_in_process(pos+1,PWM.size());

		int n = 1 + pos;
		background.draw_text(pos*500+450,1230, std::to_string(n).c_str(), black, 0, 1, 53);

		CImg<unsigned char> iconA("../logo/icons/A.png");
		CImg<unsigned char> iconC("../logo/icons/C.png");
		CImg<unsigned char> iconG("../logo/icons/G.png");
		CImg<unsigned char> iconT("../logo/icons/T.png");

		CImg<unsigned char> column(1, 1, 1, 3, 255, 255, 255);
		column.resize(500, 1000);

		double height(1000-size(PWM,0,pos)-size(PWM,1,pos)-size(PWM,2,pos)-size(PWM,3,pos));


		std::vector<nuc_prob_pair> nuc_pairs;

		nuc_pairs.push_back(nuc_prob_pair(iconA, size(PWM,0,pos)));
		nuc_pairs.push_back(nuc_prob_pair(iconC, size(PWM,1,pos)));
		nuc_pairs.push_back(nuc_prob_pair(iconG, size(PWM,2,pos)));
		nuc_pairs.push_back(nuc_prob_pair(iconT, size(PWM,3,pos)));

		std::sort(nuc_pairs.begin(), nuc_pairs.end(), compareByProb);

		for(unsigned i = 0; i < nuc_pairs.size(); i++){
			nuc_pairs[i].icon.resize(500, nuc_pairs[i].prob, -100, -100, 2);
			column.draw_image(0, height, 0, 0, nuc_pairs[i].icon, nuc_pairs[i].icon.get_channel(3),1,255);
			height += nuc_pairs[i].prob;
		}

		background.draw_image(pos*500+200, 200, 0, 0, column);

		background.draw_line(700+pos*500,1211,700+pos*500,1230,black);
		background.draw_line(701+pos*500,1211,701+pos*500,1230,black);
		background.draw_line(702+pos*500,1211,702+pos*500,1230,black);
		background.draw_line(703+pos*500,1211,703+pos*500,1230,black);
		background.draw_line(699+pos*500,1211,699+pos*500,1230,black);
		background.draw_line(698+pos*500,1211,698+pos*500,1230,black);
		background.draw_line(697+pos*500,1211,697+pos*500,1230,black);

	}


	background.save_png("../yourlogo.png");
}

//===================================================================

/*!
     * @brief   creates the logo, calls all needed function
     *
     * @param   the filename with matrix in it
     */
void logo(std::string fileName) {

    std::vector<std::vector<double>> PWM(read_logo_matrix(fileName));

    unsigned size_motif(PWM.size());

    CImg<unsigned char> background(1, 1, 1, 3, 255, 255, 255);
    background.resize(size_motif*500+250, 1400);

    const unsigned char black[] = { 0,0,0 };
    background.draw_line(170,1205,size_motif*500+203,1205,black);
    background.draw_line(170,1206,size_motif*500+203,1206,black);
    background.draw_line(170,1207,size_motif*500+203,1207,black);
    background.draw_line(170,1208,size_motif*500+203,1208,black);
    background.draw_line(170,1209,size_motif*500+203,1209,black);
    background.draw_line(170,1210,size_motif*500+203,1210,black);
    background.draw_line(170,1211,size_motif*500+203,1211,black);

    background.draw_line(195,1230,195,200,black);
    background.draw_line(194,1230,194,200,black);
    background.draw_line(193,1230,193,200,black);
    background.draw_line(192,1230,192,200,black);
    background.draw_line(191,1230,191,200,black);
    background.draw_line(190,1230,190,200,black);
    background.draw_line(189,1230,189,200,black);

    background.draw_line(189,700,170,700,black);
    background.draw_line(189,701,170,701,black);
    background.draw_line(189,702,170,702,black);
    background.draw_line(189,703,170,703,black);
    background.draw_line(189,699,170,699,black);
    background.draw_line(189,698,170,698,black);
    background.draw_line(189,697,170,697,black);

    background.draw_line(189,200,170,200,black);
    background.draw_line(189,201,170,201,black);
    background.draw_line(189,202,170,202,black);
    background.draw_line(189,203,170,203,black);
    background.draw_line(189,204,170,204,black);
    background.draw_line(189,205,170,205,black);
    background.draw_line(189,206,170,206,black);

    background.draw_text(60,50,"bits", black, 0, 1, 103);
    background.draw_text(size_motif*250+80,1290,"position", black, 0, 1, 103);
    background.draw_text(100,210,"2", black, 0, 1, 53);
    background.draw_text(100,680,"1", black, 0, 1, 53);
    background.draw_text(100,1185,"0", black, 0, 1, 53);
    background.draw_text(189,1250,"5'", black, 0, 1, 53);
    background.draw_text(size_motif*500+190,1250,"3'", black, 0, 1, 53);

    logo_in_process();


    for (unsigned pos(0); pos<size_motif; ++pos)
    {
        position_in_process(pos+1,PWM.size());

        char n = '1'+ pos;
        background.draw_text(pos*500+450,1230, &n, black, 0, 1, 53);

        CImg<unsigned char> iconA("../logo/icons/A.png");
        CImg<unsigned char> iconC("../logo/icons/C.png");
        CImg<unsigned char> iconG("../logo/icons/G.png");
        CImg<unsigned char> iconT("../logo/icons/T.png");

        CImg<unsigned char> column(1, 1, 1, 3, 255, 255, 255);
        column.resize(500, 1000);

        double height(1000-size(PWM,0,pos)-size(PWM,1,pos)-size(PWM,2,pos)-size(PWM,3,pos));


        std::vector<nuc_prob_pair> nuc_pairs;

        nuc_pairs.push_back(nuc_prob_pair(iconA, size(PWM,0,pos)));
        nuc_pairs.push_back(nuc_prob_pair(iconC, size(PWM,1,pos)));
        nuc_pairs.push_back(nuc_prob_pair(iconG, size(PWM,2,pos)));
        nuc_pairs.push_back(nuc_prob_pair(iconT, size(PWM,3,pos)));

        std::sort(nuc_pairs.begin(), nuc_pairs.end(), compareByProb);

        for(unsigned i = 0; i < nuc_pairs.size(); i++){
            nuc_pairs[i].icon.resize(500, nuc_pairs[i].prob, -100, -100, 2);
            column.draw_image(0, height, 0, 0, nuc_pairs[i].icon, nuc_pairs[i].icon.get_channel(3),1,255);
            height += nuc_pairs[i].prob;
        }

        background.draw_image(pos*500+200, 200, 0, 0, column);

        background.draw_line(700+pos*500,1211,700+pos*500,1230,black);
        background.draw_line(701+pos*500,1211,701+pos*500,1230,black);
        background.draw_line(702+pos*500,1211,702+pos*500,1230,black);
        background.draw_line(703+pos*500,1211,703+pos*500,1230,black);
        background.draw_line(699+pos*500,1211,699+pos*500,1230,black);
        background.draw_line(698+pos*500,1211,698+pos*500,1230,black);
        background.draw_line(697+pos*500,1211,697+pos*500,1230,black);

    }


    background.save_png("../yourlogo.png");
}

#endif
