#include "utils/fasta_reader.hpp"
#include "global_alignment.hpp"

using std::cout;

void readConfigFile(const char* const config_file, int &pMATCH, int &pMISMATCH, int &pINDEL){
    CSimpleIni ini(true, true, true);
    ini.SetUnicode();
    std::cout << "match: " << pMATCH << " MISMATCH: " << pMISMATCH << " INDEL: " << pINDEL << std::endl;
    if (ini.LoadFile(config_file) < 0){
        std::cerr << "Error loading file. Setting default values\n";
    }
    else {
        pMATCH = atoi(ini.GetValue("GlobalAlignment", "match", "2"));
        pMISMATCH = atoi(ini.GetValue("GlobalAlignment", "mismatch", "-1"));
        pINDEL = atoi(ini.GetValue("GlobalAlignment", "indel", "-2"));
    }
}



int main(int argc, char **argv){
    int MATCH = 2;
    int MISMATCH = -1;
    int INDEL = -2;

    int seq1_length, seq2_length;

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
        std::cout << "match: " << MATCH << " MISMATCH: " << MISMATCH << " INDEL: " << INDEL << std::endl;
    fasta_sequences = FastaReader(argv[1]);

    //We assume only two sequences for now
    fasta_sequences[0].set_seqString("_"+fasta_sequences[0].get_seqString());
    fasta_sequences[1].set_seqString("_"+fasta_sequences[1].get_seqString());

    std::string seq1 = fasta_sequences[0].get_seqString();
    std::string seq2 = fasta_sequences[1].get_seqString();
    seq1_length = seq1.length();
    seq2_length = seq2.length();
    ScoringMatrix SM(seq1_length, seq2_length, INDEL);

    ScoringInfo SI;
    int match, del_seq2, del_seq1;
    for (int i=1; i<seq1_length; i++){
       for(int j=1; j<seq2_length; j++){

            match = SM.getMatrixEntry(i-1, j-1).score;
            if (seq1[i]==seq2[j])
                match += MATCH;
            else
                match += MISMATCH;
            del_seq2 = SM.getMatrixEntry(i-1, j).score+INDEL;
            del_seq1 = SM.getMatrixEntry(i, j-1).score+INDEL;
            SM.optimize(i, j, match, del_seq2, del_seq1);
            //std::cout << SM.getMatrixEntry(i, j).type << " ";
        }
    }
    int score = SM.getMatrixEntry(seq1_length-1, seq2_length-1).score;
    std::string seq1Output = "";
    std::string seq2Output = "";
    seq1_length = seq1_length - 1;
    seq2_length = seq2_length - 1;
    while (seq1_length > 0  || seq2_length > 0 ){
        SI = SM.getMatrixEntry(seq1_length, seq2_length);
        if(SI.type=='M'){
            seq1Output = seq1[seq1_length] + seq1Output;
            seq2Output = seq2[seq2_length] + seq2Output;
            seq1_length = seq1_length - 1;
            seq2_length = seq2_length - 1;
        }
        else if (SI.type=='2'){
            seq1Output = seq1[seq1_length] + seq1Output;
            seq2Output = "-" + seq2Output;
            seq1_length = seq1_length - 1;
        }
        else if (SI.type=='1'){
            seq1Output = "-" + seq1Output;
            seq2Output = seq2[seq2_length] + seq2Output;
            seq2_length = seq2_length - 1;
        }
        else{
            std::cerr << "Unknown score. Exiting since this is surely a bug!" << std::endl;
            exit(EXIT_FAILURE);
        }

    }

    std::cout << seq1Output << std::endl;
    std::cout << seq2Output << std::endl;
    std::cout << "Score: " << score << std::endl;
}
