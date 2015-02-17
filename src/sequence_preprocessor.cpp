#include <vector>
#include "sequence_profile.hpp"
#include "fasta_reader.hpp"
#include <stdio.h>
#include <ctype.h>
using std::vector;
bool isDNASequence(Fasta &firstSequence){
//Read first 10 characters to determine if sequence is DNA or AA
int count = 0;
std::string seq = firstSequence.get_seqString();
for(int i=0; i<10;i++){
    if(std::find (DNA, DNA+4, toupper(seq[i]) )){
        count +=1;
    }
}
if (count>=5){
    return true;
}
return false;
}

//vector<DistanceMatrix> createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
vector<Profile> createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
    bool isDNA = isDNASequence(fasta_sequences[0]);
    //FIXME
    unsigned int rows = 5;//izeof(DNA);
    unsigned int count = fasta_sequences.size();
    const char *seqMap = DNA;
    std::string seq;
    if (!isDNA){
        rows = 21;//AA.size();
        seqMap = AA;
    }
    //std::cout<<"SEQMAP: "<<std::begin(seqMap) << std::endl;
   //const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
   //std::cout <<"END: "<<end<<std::endl;
    //std::vector<DistanceMatrix> v;
    std::vector<Profile> v;

    for (unsigned int i=0; i<count; i++){
        seq = fasta_sequences[i].get_seqString();
        //std::cout<<"LENGTH" << seq.length() << std::endl;
        //DistanceMatrix SP = DistanceMatrix(rows, seq.length());
        vector< vector<float> > score(rows, vector<float> (seq.length()));
        for (unsigned int x=0; x<seq.length(); x++){

            const char *p = std::find(seqMap, seqMap+sizeof(seqMap), seq[x]);
            //std::cout<<" "<<seq[x]<<std::endl;
            if (p!=seqMap+sizeof(seqMap)){
                int dist = std::distance(seqMap, p);
                //SP.incrementCount(dist, x);
                score[dist][x]=1;
            }
            else{
                std::cout << "[Error]:: Found" << seq[x]<<std::endl;
            }

        }
        //std::cout<<"rows "<<rows<<std::endl;
        //std::cout<<"i "<<i<<std::endl;

        Profile profile(i, rows, seq.length(), score);
        profile.columns = seq.length();
        //std::cout<<"name: "<<profile.seqNumber<<std::endl;
        v.push_back(profile);
   }
   return v;
}

void printProfile(vector<Profile> &profiles){
    //std::cout<< "Profiles: "<<profiles.size();
    for (unsigned int i=0;i<profiles.size();i++){
        std::cout<<"profile: "<<profiles[i].seqNumber<<std::endl;
        for(int x=0;x<profiles[i].rows ; x++){
            for(int y=0; y<profiles[i].columns; y++){
                std::cout<<profiles[i].profile[x][y]<<" ";
            }
            std::cout<<std::endl;
        }
    }
}
int createDistanceMatrix(DistanceMatrix &SP){
    /*const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
    for (unsigned int j=0; j<columns; j++){
        seq = fasta_sequences[j].get_seqString();
        for (unsigned int x=0; x<seq.length(); x++){
            const char *p = std::find(seqMap, end, seq[x]);
            if (p!=end){
                int dist = std::distance(seqMap, p);
                SP.incrementCount(dist,j);
            }
            else{
                std::cout << "[Error]:: Found" << seq[x]<<std::endl;
            }

        }
    }
    return SP;*/
    for(int i=0; i<SP.getRows(); i++){

    }
    return 1;

}

ProfileAlignment calculatePairwiseAlignment(Profile &seq1Profile, Profile &seq2Profile){
    std::cout<<"seq1 P: "<<seq1Profile.columns<<std::endl;
    std::cout<<"seq2 P: "<<seq2Profile.columns<<std::endl;
    int cols1 = seq1Profile.columns;
    int cols2 = seq2Profile.columns;

    ProfileAlignment P(cols1, cols2);

    float del1 = 0;
    float del2 = 0;
    float min =0;
    char type='M';
    for(int i=1; i<seq1Profile.columns; i++){
        for(int j=1; j<seq2Profile.columns; j++){
            int s = 0;
            for (int x=0; x<seq1Profile.rows; x++){
                for(int y=0; y<seq2Profile.rows; y++){
                    if(x==y){
                        //Match
                        s=s+0*(seq1Profile.profile[x][i]*seq2Profile.profile[y][j]);
                    }
                    else if(x==4 && y==4){
                        //Both '_'
                        s=s+0*(seq1Profile.profile[x][i]*seq2Profile.profile[y][j]);
                    }
                    else if(x==4 || y==4){
                        //Indel
                        s=s+2*(seq1Profile.profile[x][i]*seq2Profile.profile[y][j]);
                    }
                    else{
                        //Mismatch
                        s=s+1*(seq1Profile.profile[x][i]*seq2Profile.profile[y][j]);
                    }
                }
            }
            min = P.scores[i-1][j-1] + s;
            del1 = P.scores[i-1][j] + 2;
            del2 = P.scores[i][j-1] + 2;
            (min > del2) && (min = del2) && (type='2');
            (min > del1) && (min = del1) && (type='1');
            P.type[i][j] = type;
            P.scores[i][j] = min;
        }
    }
    return P;
}

vector<string> getOptimalProfileAlignment(ProfileAlignment P, string &seq1, string &seq2){

    int seq1_length = P.seq1Length;
    int seq2_length = P.seq2Length;
    vector<string> output;
    std::string seq1Output = "";
    std::string seq2Output = "";
    int score;
    char type='X';
    while (seq1_length > 0  || seq2_length > 0 ){
        score  = P.scores[seq1_length][seq2_length];
        if(type=='M'){
            seq1Output = seq1[seq1_length] + seq1Output;
            seq2Output = seq2[seq2_length] + seq2Output;
            seq1_length = seq1_length - 1;
            seq2_length = seq2_length - 1;
        }
        else if (type=='2'){
            seq1Output = seq1[seq1_length] + seq1Output;
            seq2Output = "-" + seq2Output;
            seq1_length = seq1_length - 1;
        }
        else if (type=='1'){
            seq1Output = "-" + seq1Output;
            seq2Output = seq2[seq2_length] + seq2Output;
            seq2_length = seq2_length - 1;
        }
        else{

            std::cerr << "Unknown score. Exiting since this is surely a bug! "  << type << " score: " << score <<  std::endl;
            exit(EXIT_FAILURE);
        }
    }

    output.push_back(seq1Output);
    output.push_back(seq2Output);
    return output;
}
