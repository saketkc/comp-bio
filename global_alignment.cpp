#include "utils/fasta_reader.hpp"
using std::cout;


class ScoringInfo{
    public:
        int score;
        //MATCH = A
        //MISMATCH = M
        //INDEL = I
        char type;
};



int main(int argc, char **argv){
    std::vector<Fasta> fasta_sequences;
    fasta_sequences = FastaReader(argv[1]);

    //We assume only two sequences for now
    fasta_sequences[0].set_seqString("_"+fasta_sequences[0].get_seqString());
    fasta_sequences[1].set_seqString("_"+fasta_sequences[1].get_seqString());

    int n = fasta_sequences[0].seqString.length();
    int m = fasta_sequences[1].seqString.length();

    ScoringInfo R[n][m];
    int MATCH = 1;
    int MISMATCH = -22;
    int INDEL = 0;
    std::cout << "n: " << n << std::endl;
    for(int i=0; i<n; i++){
        R[i][0].score = INDEL*i;
        R[i][0].type = 'I';
    }

    for(int i=0; i<m; i++){
        R[0][i].score = INDEL*i;
        R[0][i].type = 'I';

    }

    for (int i=1; i<n; i++){
        for(int j=1; j<m; j++){
            int match = R[i-1][j-1].score;
            if (fasta_sequences[0].seqString[i]==fasta_sequences[1].seqString[j])
                match += MATCH;
            else
                match += MISMATCH;
            int indel_u = R[i-1][j].score+INDEL;
            int indel_v = R[i][j-1].score+INDEL;
            cout << "m: " << match << " u: " << indel_u << " v: " << indel_v << std::endl;
            if (match>indel_u && match > indel_v){
                R[i][j].score = match;
                R[i][j].type = 'A';
            }
            else if(indel_u > indel_v){
                R[i][j].type = 'I';
                R[i][j].score = indel_u;
            }
            else{
                R[i][j].type = 'I';
                R[i][j].score = indel_v;

            }
        }

    }

    for (int i =0; i<n; i++){
        for (int j=0; j<m; j++){
            std::cout << R[i][j].score << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Score: " << R[n-1][m-1].score << std::endl;



}
