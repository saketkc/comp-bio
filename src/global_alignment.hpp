#include "SimpleIni.h"
#include <ctime>
#include <iostream>
using std::string;
using std::vector;
int Factorial(int n){
    return n;
}
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
void readConfigFile(const char* const config_file, int &pMATCH, int &pMISMATCH, int &pINDEL){
    CSimpleIni ini(true, true, true);
    ini.SetUnicode();
    //std::cout << "match: " << pMATCH << " MISMATCH: " << pMISMATCH << " INDEL: " << pINDEL << std::endl;
    if (ini.LoadFile(config_file) < 0){
        std::cerr << "Error loading file. Setting default values\n";
    }
    else {
        pMATCH = atoi(ini.GetValue("GlobalAlignment", "match", "2"));
        pMISMATCH = atoi(ini.GetValue("GlobalAlignment", "mismatch", "-1"));
        pINDEL = atoi(ini.GetValue("GlobalAlignment", "indel", "-2"));
    }
}

vector<string> getOptimalAlignment(const ScoringMatrix &SM, string &seq1, string &seq2){
    int seq1_length = SM.getRowSize()-1;
    int seq2_length = SM.getColumnSize()-1;

    ScoringInfo SI;

    std::string seq1Output = "";
    std::string seq2Output = "";
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

    vector<string> output;
    output.push_back(seq1Output);
    output.push_back(seq2Output);
    return output;
}

int  getOptimalScore(const ScoringMatrix &SM){
    int seq1_length = SM.getRowSize()-1;
    int seq2_length = SM.getColumnSize()-1;
    int score = SM.getMatrixEntry(seq1_length, seq2_length).score;
    return score;
}

void performGlobalAlignment(ScoringMatrix &SM, const int &MATCH, const int &MISMATCH, const int &INDEL, std::string &seq1, std::string &seq2){

    int match=0, del_seq2=0, del_seq1=0;
    int seq1_length = SM.getRowSize();
    int seq2_length = SM.getColumnSize();

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

}

//TODO This should be generic too for multiple sequences
ScoringMatrix createScoringMatrixFromSequences(string &seq1, string &seq2) {




    //We assume only two sequences for now

    seq1 = "_" + seq1;
    seq2 = "_" + seq2;
    int seq1_length = seq1.length();
    int seq2_length = seq2.length();
    ScoringMatrix SM(seq1_length, seq2_length);
    return SM;


}
