#include "sequence_profile.hpp"
using std::min;
using std::max;

ProfileMatrix::ProfileMatrix(int rows, int columns) :rows(rows), columns(columns){
    initProfileMatrix();
}

ProfileMatrix::~ProfileMatrix(){
    //Free the matrix
    for (int i=0; i<rows; i++){
        //Delete columns
        delete [] P[i];
    }
    //Delete rows
    delete [] P;
}
