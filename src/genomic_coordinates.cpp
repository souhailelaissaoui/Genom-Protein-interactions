#include "genomic_coordinates.hpp"
#define MAXLINESIZE 1000

Coordinates::Coordinates( std::string file_description, std::string location, bool line_description_present)
:file_description(file_description),
location(location),
line_description_present(line_description_present)
{

}


//-------------------------------------------------------------------------

void Coordinates::fillNewLine(std::string line) {
    std::stringstream line_stream(line);
    unsigned int start_pos_line;
    unsigned int end_pos_line;
    std::string line_description;
    double line_score;

    if(not(line_stream >> start_pos_line >> end_pos_line)) {
        return;
    };

    if(line_description_present) {
        if(not(getline(line_stream, line_description, ' '))) {
            return;
        }
    }

    if(not(line_stream >> line_score)) {
        return;
    }

    start_pos.push_back(start_pos_line);
    length.push_back(end_pos_line - start_pos_line);

    if(line_description_present)
        description.push_back(line_description);

    score.push_back(line_score);
}

//-------------------------------------------------------------------------

std::string Coordinates::get_location() {
    return location;
}

//-------------------------------------------------------------------------

std::vector <double> Coordinates::position_score(SearchResults input, unsigned int pos_zero)
{
    std::vector <double> output;
    unsigned int gen_id(0);
    unsigned int seq_id(0);

    unsigned int seq_pos;
    unsigned int seq_length;
    //unsigned int gen_pos;
    //unsigned int gen_length;

    double score_intermed(0);


    while (seq_id < input.searchResults.size())
    {
        seq_pos = input.searchResults[seq_id].position + pos_zero;
        seq_length = input.searchResults[seq_id].sequence.size();
        gen_id = 0;
        score_intermed = 0;

        for (unsigned int i(0); i<seq_length; i++, seq_pos++) {
            while (start_pos[gen_id + 1] < seq_pos and gen_id < start_pos.size() - 1) {
                gen_id++;
            }
            if (start_pos[gen_id] + length[gen_id] >= seq_pos) {

                score_intermed += score[gen_id];
            }
        }
        output.push_back(score_intermed / seq_length);
        seq_id++;
    }


    return output;
}
