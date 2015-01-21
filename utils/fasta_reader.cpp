#include "fasta_reader.hpp"

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
