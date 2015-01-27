#include <ctime>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

using std::string;
using std::vector;

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
    public:
        void initMatrix(){
            // Pointer to pointer p is assigned to an array of integers
            R = new ScoringInfo*[rows];
            for (int i=0; i<rows; i++){
                R[i] = new ScoringInfo [columns];
            }
        }
        ScoringMatrix(int rows, int columns);
        void initializeIndelPenalties(int INDEL);
        void reset(int i, int j, int value);
        void optimize(int i, int j, int match_mismatch_score, int indel_seq1, int indel_seq2);
        void minimumDistance(int i, int j, int match_distance, int del_seq2, int del_seq1);
        /**
        * NOTE: inlining is implicit for function implmented
        * inside a class: http://stackoverflow.com/a/86576/756986
        * We are still going to implement the function
        * here too as opposed to http://www.parashift.com/c++-faq-lite/inline-member-fns-more.html
        * since they are mere get functions
        **/


        inline const int  getRowSize() const {return rows;};
        inline const int  getColumnSize() const {return columns;};
        inline const ScoringInfo getMatrixEntry(int i, int j) const {return  R[i][j];};
        ~ScoringMatrix();

};


// Define other methods implemented in scoring_matrix.cpp


vector<string> getOptimalAlignment(const ScoringMatrix &SM, string &seq1, string &seq2);

int getOptimalScore(const ScoringMatrix &SM);
void performGlobalAlignment(ScoringMatrix &SM, const int &MATCH, const int &MISMATCH, const int &INDEL, std::string &seq1, std::string &seq2);
void performKBandAlignment(ScoringMatrix &SM, int k, const int &MATCH, const int &MISMATCH, const int &INDEL, std::string &seq1, std::string &seq2);
ScoringMatrix createScoringMatrixFromSequences(string &seq1, string &seq2);

