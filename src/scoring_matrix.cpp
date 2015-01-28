#include "scoring_matrix.hpp"
using std::min;
using std::max;
//Put constructors in .cpp unless they are inline implementations
// http://stackoverflow.com/a/4761871/756986

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

//Better to inline this
inline void ScoringMatrix::initializeIndelPenalties(int INDEL){
    for(int i=0; i<rows; i++){
        R[i][0].score = i*INDEL;
    }

    for(int i=0; i<columns; i++){
        R[0][i].score = i*INDEL;
    }
}

//Function where the decision takes place
void ScoringMatrix::optimize(int i, int j, int match_mismatch_score, int del_seq2, int del_seq1, bool isDistance){
    // Is the match mismatch score greater than the other two scores
    // If the scores are equal, we still prefer a Match/Mismatch over insertion/deletion

    // An optimum way to find the max of three numbers
    // TODO : This optimize function should be generic and not
    // retriscted to just three scores
    if(!isDistance){
    int max = match_mismatch_score;
    char type = 'M';

    //deletion on seq2
    (max < del_seq2) && (max = del_seq2) && (type = '2');
    //deletion on seq1
    (max < del_seq1) && (max = del_seq1) && (type = '1');

    R[i][j].score = max;
    R[i][j].type = type;}
    else{
    int min = match_mismatch_score;
    char type = 'M';
    (min > del_seq2) && (min = del_seq2) && (type='2');
    (min > del_seq1) && (min = del_seq1) && (type='1');
    R[i][j].type = type;
    R[i][j].score = min;

    }
}

void ScoringMatrix::reset(int i, int j, int value){
    R[i][j].score = value;
}
void ScoringMatrix::minimumDistance(int i, int j, int match_distance, int del_seq2, int del_seq1){
    int min = match_distance;
    char type = 'M';
    (min > del_seq2) && (min = del_seq2) && (type='2');
    (min > del_seq1) && (min = del_seq1) && (type='1');
    R[i][j].type = type;
    R[i][j].score = min;
}

//Destructor
ScoringMatrix::~ScoringMatrix(){
    //Free the matrix
    for (int i=0; i<rows; i++){
        //Delete columns
        delete [] R[i];
    }
    //Delete rows
    delete [] R;
}


void printScoringMatrix(const ScoringMatrix &SM){

    int seq1_length = SM.getRowSize();
    int seq2_length = SM.getColumnSize();
    std::cout<< "Seq1 length: " << seq1_length << " Seq2 length: " << seq2_length << std::endl;
    for(int i=0; i<seq1_length; i++){
        for(int j=0; j<seq2_length; j++){
            std::cout<< SM.getMatrixEntry(i,j).score << " ";
        }
        std::cout << std::endl;
    }
}

void printScoringMatrixType(const ScoringMatrix &SM){

    int seq1_length = SM.getRowSize();
    int seq2_length = SM.getColumnSize();
    for(int i=0; i<seq1_length; i++){
        for(int j=0; j<seq2_length; j++){
            std::cout<< SM.getMatrixEntry(i,j).type << " ";
        }
        std::cout << std::endl;
    }
}

vector<string> getOptimalAlignment(const ScoringMatrix &SM, string &seq1, string &seq2){

    int seq1_length = SM.getRowSize()-1;
    int seq2_length = SM.getColumnSize()-1;
    ScoringInfo SI;
    vector<string> output;
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

           //std::cerr << "Unknown score. Exiting since this is surely a bug! "  << SI.type << " score: " << SI.score <<  std::endl;
            std::cout << "Unknwon score. I: " << seq1_length << " j " << seq2_length << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    output.push_back(seq1Output);
    output.push_back(seq2Output);
    return output;
}
vector<string> getOptimalAlignmentFromKBand(const ScoringMatrix &SM, string &seq1, string &seq2, int k){

    int seq1_length = SM.getRowSize()-1;
    int seq2_length = SM.getColumnSize()-1;
    ScoringInfo SI;
    vector<string> output;
    std::string seq1Output = "";
    std::string seq2Output = "";

    while (seq1_length > 0  || seq2_length > 0 ){
        SI = SM.getMatrixEntry(seq1_length, seq2_length);
        if (!insideBand(seq1_length, seq2_length, k)){
            std::cout << "Unknwon score so increasing I: " << seq1_length << " j " << seq2_length << std::endl;
            seq2_length+=1;
        }
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
            seq2_length+=1;
           //std::cerr << "Unknown score. Exiting since this is surely a bug! "  << SI.type << " score: " << SI.score <<  std::endl;
            std::cout << "Unknwon score. I: " << seq1_length << " j " << seq2_length << std::endl;
           //exit(EXIT_FAILURE);
        }
    }

    output.push_back(seq1Output);
    output.push_back(seq2Output);
    return output;
}


int getOptimalScore(const ScoringMatrix &SM){
    int seq1_length = SM.getRowSize()-1;
    int seq2_length = SM.getColumnSize()-1;
    int score = SM.getMatrixEntry(seq1_length, seq2_length).score;
    return score;
}

void performGlobalAlignment(ScoringMatrix &SM, const int &MATCH, const int &MISMATCH, const int &INDEL, std::string &seq1, std::string &seq2, bool isDistance){

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
            SM.optimize(i, j, match, del_seq2, del_seq1, isDistance);
        }
    }

}


void performKBandAlignment(ScoringMatrix &SM, int k, const int &dMATCH, const int &dMISMATCH, const int &dINDEL, std::string &seq1, std::string &seq2){
    int match=0, del_seq2=0, del_seq1=0;
    int seq1_length = SM.getRowSize();
    int seq2_length = SM.getColumnSize();
    int left, right;
    for (int j=1; j<seq2_length; j++){
        SM.reset(0,j,j*dINDEL);
    }
    for (int i=1; i<seq1_length; i++){
        left = max(0, i-k);
        right = min(seq2_length, i+k);

        for(int j=left; j<right; j++){
            SM.reset(i,j,0);
        }
        SM.reset(i,0,i*dINDEL);
    }
    for (int i=1; i<seq1_length; i++){
        //left = max(0, i-(k + abs(seq1_length-seq2_length))/2);
        //right = min(seq2_length, i+(k + abs(seq1_length-seq2_length))/2);
        left = max(0, i-k);
        left = min(left, seq2_length-1);
        right = min(seq2_length-1, i+k);
        std::cout <<"k: " << k << "i: " << i << " left: " << left<< " right: "<< right << std::endl;
        for(int j=left; j<=right; j++){

            if(insideBand(i-1,j-1,k)){
                match = SM.getMatrixEntry(i-1, j-1).score;
            }
            else{
                match=0;
            }
            if (seq1[i]==seq2[j])
                match += dMATCH;
            else
                match += dMISMATCH;
            if (insideBand(i-1,j,k)){
                del_seq2 = SM.getMatrixEntry(i-1,j).score;
            }
            else {
                del_seq2 = 0;
            }
            del_seq2 += dINDEL;
            if(insideBand(i, j-1, k)){
            del_seq1 = SM.getMatrixEntry(i, j-1).score;
            }
            else{
                del_seq1 = 0;
            }

            del_seq1 += dINDEL;
            SM.minimumDistance(i, j, match, del_seq2, del_seq1);
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

bool insideBand(int i, int j, int k){
    if(std::abs(i-j)<=k)
        return true;
    return false;
}
