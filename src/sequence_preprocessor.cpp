#include <vector>
#include "sequence_profile.hpp"
#include "fasta_reader.hpp"
#include <stdio.h>
#include <ctype.h>

bool isDNASequence(Fasta &firstSequence){
//Read first 10 characters to determine if sequence is DNA or AA
int count = 0;
std::string seq = firstSequence.get_seqString();
for(int i=0; i<10;i++){
    if(std::find (DNA, DNA+4, toupper(seq[i]) )){
        count +=1;
    }
}
if (count>=5){
    return true;
}
return false;
}

ProfileMatrix createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
    bool isDNA = isDNASequence(fasta_sequences[0]);
    //FIXME
    unsigned int rows = sizeof(DNA);
    unsigned int columns = fasta_sequences.size();
    const char *seqMap = DNA;
    std::string seq;
    if (!isDNA){
        rows = sizeof(AA);//AA.size();
        seqMap = AA;
    }

    ProfileMatrix SP = ProfileMatrix(rows, columns);
    std::cout << "Rows: " << rows << " Columns: " << columns << std::endl;

    //Loop through the sequences(columns)
    //Increase count vectors for DNA//AA Sequence
    const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
    for (unsigned int j=0; j<columns; j++){
        seq = fasta_sequences[j].get_seqString();
        for (unsigned int x=0; x<seq.length(); x++){
            const char *p = std::find(seqMap, end, seq[x]);
            if (p!=end){
                int dist = std::distance(seqMap, p);
                SP.incrementCount(dist,j);
            }
            else{
                std::cout << "[Error]:: Found" << seq[x]<<std::endl;
            }

        }
    }
    return SP;
}
