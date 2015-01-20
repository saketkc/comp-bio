#include "SimpleIni.h"
#include <ctime>

class ScoringInfo{
    public:
        //Let the score and type be public to avoid unecessary get, set methods
        int score;
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
        ScoringMatrix(int rows, int columns);
        void optimize(int i, int j, int match_mismatch_score, int indel_seq1, int indel_seq2);
        void initializeIndelPenalties(int INDEL);
        const int  getRowSize() const;
        const int  getColumnSize() const;
        const ScoringInfo getMatrixEntry(int i, int j) const;
        ~ScoringMatrix();

};


ScoringMatrix::ScoringMatrix(int rows, int columns) :rows(rows), columns(columns){
    initMatrix();
    for (int i=0; i<rows; i++){
        for (int j=0; j<columns; j++){
            R[i][j].score = 0;
            if (i==0){
                //Deletion in Seq1
                R[i][j].type = '1';
            }
            if (j==0){
                //Deletion in Seq2
                R[i][j].type = '2';
            }
        }
    }
}

void ScoringMatrix::initializeIndelPenalties(int INDEL){
    for(int i=0; i<rows; i++){
        R[i][0].score = i*INDEL;
    }

    for(int i=0; i<columns; i++){
        R[0][i].score = i*INDEL;
    }
}

const int ScoringMatrix::getColumnSize() const{
    return columns;
}

const int ScoringMatrix::getRowSize() const{
    return rows;
}

const ScoringInfo ScoringMatrix::getMatrixEntry(int i, int j) const{
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
