#include "fasta_reader.hpp"
#include "scoring_matrix.hpp"
#include "timer.hpp"
#include <stdlib.h>     /* exit, EXIT_FAILURE */

using std::cout;

int main(int argc, char **argv){
    int MATCH = 0;
    int MISMATCH = 1;
    int INDEL = 1;


    std::vector<Fasta> fasta_sequences;

    if(argc < 2){
        std::cerr << std::endl << "Usage: global_alignment <fasta-file> [conig_file]\
            \nNote: Only first two sequences are read for alignment \n\n";
        exit(EXIT_FAILURE);
    }
    if(argc>=3){
        std::cout<< "loading file";
        char *config_file = argv[2];
        readConfigFile(config_file, MATCH, MISMATCH, INDEL);
        std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;

    }
    //    std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;
    timestamp_t startTime=0, endTime=0, difference=0;
    startTime = getTimeinMilliSeconds();
    fasta_sequences = FastaReader(argv[1]);
    endTime = getTimeinMilliSeconds();
    difference = endTime-startTime;
    std::cout.precision(15);
    std::string seq1 = fasta_sequences[0].get_seqString();
    std::string seq2 = fasta_sequences[1].get_seqString();
    if(seq1.length()<seq2.length()){
        std::string seqTemp = seq1;
        seq1 = seq2;
        seq2 = seqTemp;
    }
    std::cout << std::endl << "Sequence1 Length: " << seq1.length() << std::endl;
    std::cout <<  "Sequence2 Length: " << seq2.length() << std::endl;
    std::cout << std::endl << "[LOG] Reading complete in: " << difference << " ms" << std::endl;



    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << "Sequence 1: " << seq1 << std::endl;
    std::cout << "Sequence 2: " << seq2 << std::endl;

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;

    startTime = getTimeinMilliSeconds();
    ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);
    performGlobalAlignment(SM, MATCH, MISMATCH, INDEL, seq1, seq2, true);
    endTime = getTimeinMilliSeconds();
    difference = endTime-startTime;

    std::cout.precision(15);

    std::cout << std::endl << "[LOG] Aligning complete in: " << difference << " ms" << std::endl;

    vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
    int score = getOptimalScore(SM);
    //std::cout <<  seq1.length() << "," << seq2.length() << std::endl;

    std::cout << std::endl <<  "----------------------Optimal Alignment Start--------------------------" << std::endl;
    std::cout << std::endl << seqOutput[0] << std::endl;
    std::cout << std::endl << seqOutput[1] << std::endl;
    std::cout << std::endl <<  "----------------------Optimal Alignment End--------------------------" << std::endl;
    std::cout << "Score: " << score << std::endl;
   //
    //printScoringMatrix(SM);
    //std::cout << std::endl <<  std::endl;
    //printScoringMatrixType(SM);
   // std::cout <<  seq1.length() << "," << seq2.length()  << "," << difference << "," << score << std::endl;
}
