#include "utils/fasta_reader.hpp"
using std::cout;

int maxScore(int x, int y, int z){
    cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
    return std::max(std::max(x, y), z);
}

int score(char x, char y){
    if (x==y){
        return 1;
    }
    return -1;
}

int main(int argc, char **argv){
    std::vector<Fasta> fasta_sequences;
    fasta_sequences = FastaReader(argv[1]);

    //We assume only two sequences for now
    fasta_sequences[0].set_seqString("_"+fasta_sequences[0].get_seqString());
    fasta_sequences[1].set_seqString("_"+fasta_sequences[1].get_seqString());

    int n = fasta_sequences[0].seqString.length();
    int m = fasta_sequences[1].seqString.length();

    int R[n][m];

    int INDEL = -2;

    for(int i=0; i<=n; i++){
        R[i][0] = INDEL*i;
    }

    for(int i=0; i<=m; i++){
        R[0][i] = INDEL*i;
    }

    for (int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            R[i][j] = maxScore(R[i-1][j-1]+1, R[i-1][j]+INDEL, R[i][j-1]+INDEL);
        }
    }

    std::cout << "Score: " << R[n][m] << std::endl;



}
