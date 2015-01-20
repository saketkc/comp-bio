#include "SimpleIni.h"
#include <ctime>

class ScoringInfo{
    public:
        //Let the score and type be public to avoid unecessary get, set methods
        int score;
        //MATCH = A
        //MISMATCH = M
        //INDEL = I
        char type;
};

class ScoringMatrix{
    int rows, columns;
    //Pointer to a pointer, since an array <=> pointer
    ScoringInfo **R;
    void initMatrix(){
        // Pointer to pointer p is assigned to an array of integers
        R = new ScoringInfo*[rows];
        for (int i=0; i<rows; i++){
            R[i] = new ScoringInfo [columns];
        }
    }
    public:
        ScoringMatrix(int rows, int columns, int INDEL);
        void optimize(int i, int j, int match_mismatch_score, int indel_seq1, int indel_seq2);
        ScoringInfo getMatrixEntry(int i, int j);
        ~ScoringMatrix();

};


ScoringMatrix::ScoringMatrix(int rows, int columns, int INDEL) :rows(rows), columns(columns){
    initMatrix();
    for (int i=0; i<rows; i++){
        for (int j=0; j<columns; j++){
            R[i][j].score = 0;
            if (i==0){
                R[i][j].score = j*INDEL;
                R[i][j].type = '1';
            }
            if (j==0){
                R[i][j].score = i*INDEL;
                R[i][j].type = '2';
            }
        }
    }
}
ScoringInfo ScoringMatrix::getMatrixEntry(int i, int j){
    return R[i][j];
}
void ScoringMatrix::optimize(int i, int j, int match_mismatch_score, int del_seq2, int del_seq1){
    // Is the match mismatch score greater than the other two scores
    // If the scores are equal, we still prefer a Match/Mismatch over insertion/deletion

    // An optimum way to find the max of three numbers
    // TODO : This optimize function should be generic and not
    // retriscted to just three scores
    int max = match_mismatch_score;
    char type = 'M';

    //deletion on seq2
    (max < del_seq2) && (max = del_seq2) && (type = '2');
    //deletion on seq1
    (max < del_seq1) && (max = del_seq1) && (type = '1');

    R[i][j].score = max;
    R[i][j].type = type;
}

ScoringMatrix::~ScoringMatrix(){
    //Free the matrix
    for (int i=0; i<rows; i++){
        //Delete columns
        delete [] R[i];
    }
    //Delete rows
    delete [] R;
}
