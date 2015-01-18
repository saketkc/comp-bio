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
    ScoringInfo **p;
    void initMatrix(){
        // Pointer to pointer p is assigned to an array of integers
        p = new ScoringInfo*[rows];
        for (int i=0; i<rows; i++){
            p[i] = new ScoringInfo [columns];
        }
    }
    public:
        ScoringMatrix(int rows, int columns);
        ~ScoringMatrix();

};
ScoringMatrix::ScoringMatrix(int rows, int columns){
    initMatrix();
    for (int i=0; i<rows; i++){
        for (int j=0; i<columns; i++){
            p[i][j].score = 0;
            p[i][j].type = 'I';
        }
    }
}
ScoringMatrix::~ScoringMatrix(){
    for (int i=0; i<rows; i++){
        //Delete columns
        delete [] p[i];
    }
    //Delete rows
    delete [] p;
}
