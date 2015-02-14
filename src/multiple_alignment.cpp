#include "fasta_reader.hpp"
#include "sequence_preprocessor.hpp"
#include "timer.hpp"
#include <cstdlib>     /* exit, EXIT_FAILURE */

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
        char *config_file = argv[2];
        readConfigFile(config_file, MATCH, MISMATCH, INDEL);
        std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;

    }
    timestamp_t startTime=0, endTime=0, difference=0;
    startTime = getTimeinMilliSeconds();
    fasta_sequences = FastaReader(argv[1]);
    endTime = getTimeinMilliSeconds();
    difference = endTime-startTime;
    std::cout << difference;
    //Determine type of Sequence is DNA or AA.
    //preprocessSequences(fasta_sequences);
    createProfileFromSequences(fasta_sequences);
    int numberOfSequences = fasta_sequences.size();
    return numberOfSequences;
}
