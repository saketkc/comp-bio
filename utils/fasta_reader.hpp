#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::string;
class Fasta{
    public:
        std::string seqName;
        std::string seqString;
        inline void set_seqName(string name) {seqName = name;};
        inline void set_seqString(string str) {seqString = str;};
        inline string get_seqName() {return seqName;};
        inline string get_seqString() {return seqString;};
        inline void clear_seqName() {seqName.clear();};
        inline void clear_seqString() {seqString.clear();};
        inline void clearAll() {clear_seqName(); clear_seqString();};
        inline bool seqName_empty() {return seqName.empty()?true:false;};
        inline bool seqString_empty() {return seqString.empty()?true:false;};
};

std::vector<Fasta> FastaReader(char* path);
void readConfigFile(const char* const config_file, int &pMATCH, int &pMISMATCH, int &pINDEL);
