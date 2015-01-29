#include "fasta_reader.hpp"
#include "scoring_matrix.hpp"
#include "timer.hpp"
#include <cstdlib>

using std::min;

int main(int argc, char **argv){
    int dMATCH = 0;
    int dMISMATCH = 1;
    int dINDEL = 1;

    std::cout << "match: " << dMATCH << " MISMATCH: " << dMISMATCH << " INDEL: " << dINDEL << std::endl;
    std::vector<Fasta> fasta_sequences;

    if(argc < 2){
        std::cerr << std::endl << "Usage: global_alignment <fasta-file> [conig_file]\
            \nNote: Only first two sequences are read for alignment \n\n";
        exit(EXIT_FAILURE);
    }
    if(argc>=3){
        std::cout<< "loading file";
        char *config_file = argv[2];
        readConfigFile(config_file, dMATCH, dMISMATCH, dINDEL);
        std::cout << "match: " << dMATCH << " MISMATCH: " << dMISMATCH << " INDEL: " << dINDEL << std::endl;

    }
    timestamp_t startTime=0, endTime=0, difference=0;
    startTime = getTimeinMilliSeconds();
    fasta_sequences = FastaReader(argv[1]);
    endTime = getTimeinMilliSeconds();
    difference = endTime-startTime;


    std::string seq1 = fasta_sequences[0].get_seqString();
    std::string seq2 = fasta_sequences[1].get_seqString();

    if(seq1.length()<seq2.length()){
        std::string seqTemp = seq1;
        seq1 = seq2;
        seq2 = seqTemp;
    }
    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << "Sequence 1: " << seq1 << std::endl;
    std::cout << "Sequence 2: " << seq2 << std::endl;

    std::cout << std::endl <<  "--------------------------------------------" << std::endl;
    std::cout << std::endl << "Sequence1 Length: " << seq1.length() << std::endl;
    std::cout <<  "Sequence2 Length: " << seq2.length() << std::endl;

    startTime = getTimeinMilliSeconds();
    ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);


    //Best score possible = minimum distance
    //For sequence A,B we would have |i-j| > k
    //such that A would have exactly k+! characters aligned
    //with a gap and so does B. The remaining (n-k-1) characters
    //are aligned to each other
    //Assume m>n |A|=m, |B|=n
    //Hence bestDistance = (n-k-1)dMATCH_or_dMISMATCH + (2(k+1)+m-n)dINDEL

    int k=1;
    int alphaK = 32000;
    int minDistance = (seq2.length()-k-1)*dMATCH + (2*(k+1)+abs(seq1.length()-seq2.length()))*dINDEL;
    std::cout << "alpha: "<< alphaK << "best distance: " << minDistance << std::endl;
    while(alphaK > minDistance ){
        performKBandAlignment(SM, k, dMATCH, dMISMATCH, dINDEL, seq1, seq2);
        k*=2;
        minDistance = (seq2.length()-k-1)*dMATCH + (2*(k+1)+abs(seq1.length()-seq2.length()))*dINDEL;
        alphaK = getOptimalScore(SM);
        //std::cout << "alpha: "<< alphaK << "best distance: " << minDistance << std::endl;
    }
    endTime = getTimeinMilliSeconds();
    std::cout.precision(15);
    std::cout << std::endl << "[LOG] Reading complete in: " << difference << " ms" << std::endl;
    difference = endTime-startTime;

    std::cout.precision(15);

    std::cout << std::endl << "[LOG] Aligning complete in: " << difference << " ms" << std::endl;
    std::cout << std::endl << "K = " << k << std::endl;
    std::cout << std::endl <<  std::endl;
    //printScoringMatrix(SM);
    std::cout << std::endl <<  std::endl;
    //printScoringMatrixType(SM);
    vector<string> seqOutput = getOptimalAlignmentFromKBand(SM, seq1, seq2, k);
    //vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
    int score = getOptimalScore(SM);
    std::cout << std::endl <<  "----------------------Optimal Alignment Start--------------------------" << std::endl;
    std::cout << std::endl << seqOutput[0] << std::endl;
    std::cout << std::endl << seqOutput[1] << std::endl;
    std::cout << std::endl <<  "----------------------Optimal Alignment End--------------------------" << std::endl;
    std::cout << "Score: " << score << std::endl;
}

