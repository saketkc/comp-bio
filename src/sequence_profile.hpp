#include <ctime>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

using std::string;
using std::vector;
const char DNA[] = {'A', 'C', 'G', 'T', '_'};
const char AA[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y', 'V', '_'};

class ProfileAlignment{
    public:
    int seq1Length;
    int seq2Length;
    vector< vector<float> > scores;//(10000);//, vector<float> (seq2Length));
    vector< vector<char> > type;//(1000);//seq1Length, vector<float> (seq2Length));
    ProfileAlignment(int r, int c): seq1Length(r+1), seq2Length(c+1){
                //std::cout<<"I: "<<seq1Length<<std::endl;
                //std::cout<<std::endl;
                //std::cout<<"J: "<<seq2Length<<std::endl;
                //std::cout<<std::endl;
                //scores.resize(seq1Length);
                //type.resize(seq1Length);
                //vector< vector<float> > scores(seq1Length, vector<float> (seq2Length));
                //vector< vector<char> > type(seq1Length, vector<float> (seq2Length));

        for (int i=0; i<seq1Length ;i++){
                //scores[i].resize(seq2Length);
                //type[i].resize(seq2Length);
                vector<float> s;
                vector<char> cc;

            for(int j=0; j<seq2Length; j++){
                s.push_back(0);
                if(i!=0)
                    cc.push_back('X');
                else {
                    cc.push_back('1');
                }
                //std::cout<<"I: "<<i<< " J: "<<j<<std::endl;
            }
            scores.push_back(s);
            type.push_back(cc);
            type[i][0]='2';
        }
        scores[0][0]=0;
        for(int i=0;i<seq1Length;i++){
            scores[i][0]=2*i;
        }
        for(int j=0;j<seq2Length; j++){
            scores[0][j]=2*j;
        }
        for (int i=0; i<seq1Length ;i++){
            for(int j=0; j<seq2Length; j++){
                std::cout<<scores[i][j]<< " ";
            }
            std::cout<<std::endl;
        }
    }
};

class Profile{
    public:
    int rows;
    int columns;
    int seqNumber;
    vector< vector<float> > profile;
    string sequence;
    vector<string> alignment;
    Profile(int seqNumber, int rows, int cols, vector< vector<float> > profile, string sequence, vector<string> alignment): rows(rows),
    columns(cols), seqNumber(seqNumber), profile(profile), sequence(sequence), alignment(alignment) {
    }
};
class DistanceMatrix{
    //Rows correspnoding to 4+1 DNA bases
    //or 20+1 AA.
    //Columns correspond to individual sequence frequency
    int rows, columns;
    //Pointer to a pointer, since an array <=> pointer
    float **P;
    public:
        void initDistanceMatrix(){
            // Pointer to pointer p is assigned to an array of integers
            P = new float*[rows];
            for (int i=0; i<rows; i++){
                P[i] = new float [columns];
                for(int j=0; j<columns; j++){
                    P[i][j] = 0.0;
                }
            }
        }
        void edit(int i, int j, float value){
            P[i][j] = value;
        }
        void print(){
            for (int i=0; i<rows; i++){
                for (int j=0;j<columns;j++){
                    std::cout << P[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }
        DistanceMatrix(int rows, int columns);
        ~DistanceMatrix();
        void deleteColumn(int i){
            columns = columns-1;
            delete [] P[i];

        }
        void incrementCount(int row, int column){
            P[row][column]+=1;

        }
       float getValue(int i,int j){
           return P[i][j];
       }
        void shiftCountsUp(){
            //If rows are 5 shift by 1/4
            //Else shift by 1/20
            float shift = 1.0/4.0;
            if(rows>=5){
                shift = 1.0/20.0;
            }
            for (int i=0;i<rows;i++){
                for(int j=0;j<columns;j++){
                    P[i][j]+=shift;
                }
            }
        }
        int getRows(){
            return rows;
        }
        int getColumns(){
            return columns;
        }

};

