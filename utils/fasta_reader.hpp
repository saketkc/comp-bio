#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::string;
class Fasta{
    public:
        std::string seqName;
        std::string seqString;
        void set_seqName(string);
        void set_seqString(string);
        string get_seqName();
        string get_seqString();
        void clear_seqName();
        void clear_seqString();
        void clearAll();
        bool seqName_empty();
        bool seqString_empty();
};

void Fasta::set_seqName(string name){
    seqName = name;
}

void Fasta::set_seqString(string str){
    seqString = str;
}

string Fasta::get_seqName(){
    return seqName;
}

string Fasta::get_seqString(){
    return seqString;
}

void Fasta::clear_seqName(){
    seqName.clear();
}

void Fasta::clear_seqString(){
    seqString.clear();
}

void Fasta::clearAll(){
    clear_seqName();
    clear_seqString();
}

bool Fasta::seqName_empty(){
    if (seqName.empty())
        return true;
    return false;
}
bool Fasta::seqString_empty(){
    if (seqString.empty())
        return true;
    return false;
}

std::vector<Fasta> FastaReader(char* path){
    std::ifstream input(path);
    std::string line;
    std::vector<Fasta> fasta_sequences;
    Fasta fasta;
    while (std::getline(input, line).good()){
        //Sequence name starts here
        if(line[0] == '>'){
            if(!fasta.seqName_empty()){
                fasta_sequences.push_back(fasta);
                fasta.clearAll();
            }
            if(fasta.seqName_empty()){
                //Set new fasta seq name and clear earlier fasta string
                fasta.set_seqName(line.substr(1));
            }
                //Push the read seq to stack
        }
        else{
            //Ensure we have a fasta seq name
            if(!fasta.seqName_empty()){
                if (fasta.seqString_empty()){
                    fasta.set_seqString("");
                }
                fasta.set_seqString(fasta.get_seqString()+line);
            }
        }
    }
    if(!fasta.seqName_empty()){
        fasta_sequences.push_back(fasta);
        fasta.clearAll();
    }
    return fasta_sequences;
}
