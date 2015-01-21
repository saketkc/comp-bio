#include <ctime>
#include <iostream>
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

