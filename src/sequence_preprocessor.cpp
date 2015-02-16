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

//vector<DistanceMatrix> createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
vector<Profile> createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
    bool isDNA = isDNASequence(fasta_sequences[0]);
    //FIXME
    unsigned int rows = 5;//izeof(DNA);
    unsigned int columns = fasta_sequences.size();
   std::cout <<"END: "<<columns<<std::endl;
    const char *seqMap = DNA;
    std::string seq;
    if (!isDNA){
        rows = 21;//AA.size();
        seqMap = AA;
    }
    //std::cout<<"SEQMAP: "<<std::begin(seqMap) << std::endl;
   //const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
   //std::cout <<"END: "<<end<<std::endl;
    std::vector<DistanceMatrix> v;
    for (unsigned int i=0; i<columns; i++){
        seq = fasta_sequences[i].get_seqString();
        std::cout<<"LENGTH" << seq.length() << std::endl;
        DistanceMatrix SP = DistanceMatrix(rows, seq.length());
        v.push_back(SP);
        for (unsigned int x=0; x<seq.length(); x++){

            const char *p = std::find(seqMap, seqMap+sizeof(seqMap), seq[x]);
            std::cout<<" "<<seq[x]<<std::endl;
            if (p!=seqMap+sizeof(seqMap)){
                int dist = std::distance(seqMap, p);
                SP.incrementCount(dist, x);
            }
            else{
                std::cout << "[Error]:: Found" << seq[x]<<std::endl;
            }

        }
       SP.print();
       std::cout<<std::endl;
       std::cout<<"FINE HERE"<<std::endl;
   }
   return v;
}

int createDistanceMatrix(DistanceMatrix &SP){
    /*const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
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
    return SP;*/
    for(int i=0; i<SP.getRows(); i++){

    }
    return 1;

}
