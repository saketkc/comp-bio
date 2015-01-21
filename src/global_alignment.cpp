#include "fasta_reader.hpp"
#include "global_alignment.hpp"
#include "scoring_matrix.hpp"

using std::cout;



int main(int argc, char **argv){
    int MATCH = 2;
    int MISMATCH = -1;
    int INDEL = -2;


    std::vector<Fasta> fasta_sequences;

    if(argc < 2){
        std::cerr << std::endl << "Usage: global_alignment <fasta-file> [conig_file]\
            \nNote: Only first two sequences are read for alignment \n\n";
        exit(EXIT_FAILURE);
    }
    if(argc>=3){
        std::cout<< "loading file";
        char *config_file = argv[1];
        readConfigFile(config_file, MATCH, MISMATCH, INDEL);
        std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;

    }
    fasta_sequences = FastaReader(argv[1]);
    std::string seq1 = fasta_sequences[0].get_seqString();
    std::string seq2 = fasta_sequences[1].get_seqString();
    ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);

    performGlobalAlignment(SM, MATCH, MISMATCH, INDEL, seq1, seq2);
    vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
    int score = getOptimalScore(SM);

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << "Seq1: " << seq1 << std::endl;
    std::cout << "Seq2: " << seq2 << std::endl;

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << "Seq1 Length: " << seq1.length() << std::endl;
    std::cout <<  "Seq2 Length: " << seq2.length() << std::endl;

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << seqOutput[0] << std::endl;
    std::cout << seqOutput[1] << std::endl;

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;

    std::cout << "Score: " << score << std::endl;

    //std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;
}
