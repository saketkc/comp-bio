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
    //std::cout <<  "[LOG] Loading sequences complete in: " << difference << " ms"<<std::endl;
    //vector<DistanceMatrix> PM = createProfileFromSequences(fasta_sequences);
    vector<Profile> PM = createProfileFromSequences(fasta_sequences);
    //printProfile(PM);
    DistanceMatrix DM = calculateDistanceMatrix(PM);
    //cout.precision (15);
    for (unsigned int i=0; i<DM.getRows(); i++){
        std::cout<<std::endl;
        for (unsigned int j=0; j<DM.getColumns(); j++){
            std::cout<<DM.getValue(i,j)<<" ";
        }
    }
    //DM.print();
    std::cout<<std::endl;
    //ProfileAlignment P = calculatePairwiseAlignment(PM[1], PM[0]);
    string seq1 = fasta_sequences[0].get_seqString();
    string seq2 = fasta_sequences[1].get_seqString();
    //std::cout <<  "[LOG] Loading sequences complete in: " << difference << " ms"<<std::endl;
   // vector<string> x =  getOptimalProfileAlignment(P, seq1, seq2);
    //PM.shiftCountsUp();
    // PM.print();
    //int numberOfSequences = fasta_sequences.size();
    //PM.deleteColumn(1);
    //PM.print();
    return 1;// numberOfSequences;
}
