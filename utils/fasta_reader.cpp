#include "fasta_reader.hpp"

int main(int argc, char** argv){
    std::vector<Fasta> fasta_sequences;

    fasta_sequences = FastaReader(argv[1]);

    for (unsigned int i=0; i<fasta_sequences.size(); i++){
        std::cout << fasta_sequences[i].get_seqName() << " : " << fasta_sequences[i].get_seqString() << std::endl;
    }
}
