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

vector<Profile> createProfileFromSequences(std::vector<Fasta> &fasta_sequences){
    bool isDNA = isDNASequence(fasta_sequences[0]);
    //FIXME
    unsigned int count = 5;//izeof(DNA);
    unsigned int rows = fasta_sequences.size();
    const char *seqMap = DNA;
    std::string seq;
    if (!isDNA){
        count = 21;//AA.size();
        seqMap = AA;
    }
    std::vector<Profile> v;
    for(unsigned int i=0;i<rows;i++){
        seq = fasta_sequences[i].get_seqString();
        vector< vector<float> > score(seq.length(), vector<float> (rows));
        for (unsigned int j=0; j<seq.length(); j++){
            for (unsigned int x=0; x<count; x++){
                const char *p = std::find(seqMap, seqMap+sizeof(seqMap), seq[j]);
                if (p!=seqMap+sizeof(seqMap)){
                    int dist = std::distance(seqMap, p);
                    //SP.incrementCount(dist, x);
                    score[j][dist]=1;
                }
                else{
                    //std::cout << "[Error]:: Found" << seq[i]<<std::endl;
                }
            }
        }
        vector<string> alignment;
        alignment.push_back(seq); //Alignment is same for sequences as
        vector<string> name;
        name.push_back(fasta_sequences[i].seqName);
        Profile profile(name, seq.length(), count, score, seq, alignment);
        v.push_back(profile);
    }
   return v;
}

void printProfile(vector<Profile> &profiles){
    //std::cout<< "Profiles: "<<profiles.size();
    for (unsigned int i=0;i<profiles.size();i++){
        //std::cout<<"profile: "<<profiles[i].seqName<<std::endl;
        for(int x=0;x<profiles[i].rows ; x++){
            for(int y=0; y<profiles[i].columns; y++){
                std::cout<<profiles[i].profile[x][y]<<" ";
            }
            std::cout<<std::endl;
        }

    }
}

float calculateHammingDistance(Profile &p1, Profile &p2){
    float comparisons = p1.columns*p2.columns;
    unsigned int maxColumns = p1.columns;
    unsigned int minColumns = p2.columns;
    if(p1.columns<p2.columns){
        maxColumns=p2.columns;
        minColumns=p1.columns;
    }
    float distance = 0.0;

    for(unsigned int i=0;i<minColumns; i++){
        if(p1.sequence[i]==p2.sequence[i]){
            //Fixme doing this to later allow match penalties?
            distance = distance+0;
        }
        else{
            distance = distance+1;
        }
    }
    distance = distance+2*(maxColumns-minColumns);
    return distance*1.0/comparisons;
}

DistanceMatrix calculateDistanceMatrix(vector<Profile> &P){
    DistanceMatrix DM(P.size(), P.size());
    float dist;
    for (unsigned int i=0; i<P.size(); i++){
        for (unsigned int j=0; j<P.size(); j++){
            dist = calculateHammingDistance(P[i], P[j]);
            DM.edit(i,j, dist);
        }
    }
    /*std::cout<<"DistanceMatrix "<<std::endl;
    for (unsigned int i=0; i<P.size(); i++){
        std::cout<<std::endl;
        for (unsigned int j=0; j<P.size(); j++){
            std::cout<<DM.getValue(i,j)<<" ";
        }
    }*/
    std::cout<<std::endl;
    return DM;
}

long optimizer(vector<float> seq1score, vector<float> seq2score, int max){
    float s=0;
    int aln1Size = seq1score.size();
    int aln2Size = seq2score.size();
    for(int i=0;i<aln1Size; i++){
        for(int j=0;j<aln2Size;j++){
            if(i==max || j==max){
                s=s+2*(seq1score[i]*seq2score[j]);
            }
            else if(i==j){
                s=s+0*seq1score[i]*seq2score[j];
            }
            else{
                s=s+1*(seq1score[i]*seq2score[j]);
            }
        }
    }
    return s;
}

vector<string> mergeProfiles(vector<string> &seq1Alignment, vector<string> &seq2Alignment, ProfileAlignment M){
    int seq1Size = seq1Alignment[0].size();
    int seq2Size = seq2Alignment[0].size();

    vector<string> aln1(seq1Alignment.size());
    vector<string> aln2(seq2Alignment.size());
    std::cout<<std::endl;
    //M.print();
    std::cout<<std::endl;
    while (seq1Size >0 || seq2Size > 0){
        if(M.type[seq1Size][seq2Size]=='M'){
            for(unsigned int i=0; i<seq1Alignment.size(); i++){
                aln1[i].append(seq1Alignment[i].substr(seq1Size-1,1));
            }
            for(unsigned int j=0; j<seq2Alignment.size(); j++){
                aln2[j].append(seq2Alignment[j].substr(seq2Size-1,1));
            }
            --seq1Size;
            --seq2Size;
        }
        else if(M.type[seq1Size][seq2Size]=='2'){
            for(unsigned int i=0; i<seq1Alignment.size(); i++){
                aln1[i].append(seq1Alignment[i].substr(seq1Size-1,1));
            }
            for(unsigned int j=0; j<seq2Alignment.size(); j++){
                aln2[j].append("-");
            }
            --seq1Size;
        }
        else if(M.type[seq1Size][seq2Size]=='1'){

            for(unsigned int i=0; i<seq1Alignment.size(); i++){
                aln1[i].append("-");
            }
            for(unsigned int j=0; j<seq2Alignment.size(); j++){
                aln2[j].append(seq2Alignment[j].substr(seq2Size-1,1));
            }
            --seq2Size;
        }
        else {
            std::cout<<"BUG: "<<M.type[seq1Size][seq2Size]<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for(unsigned int i=0;i<aln1.size();i++){
        reverse(aln1[i].begin(), aln1[i].end());
    }
    for(unsigned int j=0;j<aln2.size();j++){
        reverse(aln2[j].begin(), aln2[j].end());
    }
    vector<string> merged;
    for (unsigned int i=0;i<aln1.size();i++){
        merged.push_back(aln1[i]);
    }
    for (unsigned int i=0;i<aln2.size();i++){
        merged.push_back(aln2[i]);
    }
    return merged;
}


Profile aligner(vector<Profile> &Profiles, int minX, int minY){
    int aln1Size = Profiles[minX].alignment[0].size()+1;
    int aln2Size = Profiles[minY].alignment[0].size()+1;
    int columns = Profiles[0].columns;
    int del1;
    int del2;
    vector<float> indelP;
    for(int k=0;k<columns;k++){
        indelP.push_back(0);
    }
    indelP.push_back(1);
    ProfileAlignment M(aln1Size, aln2Size);
    M.scores[0][0]=0;
    M.type[0][0]='X';
    for (int j=1;j<aln2Size;j++){
        M.scores[0][j] = M.scores[0][j-1] + optimizer(Profiles[minY].profile[j-1], indelP, columns);
        M.type[0][j] = '1';
    }
    for (int i=1;i<aln1Size;i++){
        M.scores[i][0] = M.scores[i-1][0] + optimizer(Profiles[minX].profile[i-1], indelP, columns);
        M.type[i][0] = '2';
        char type='M';
        int min=32000;

        for (int j=1;j<aln2Size;j++){
            min = M.scores[i-1][j-1] + optimizer(Profiles[minX].profile[i-1], Profiles[minY].profile[j-1], columns);
            del1 = M.scores[i-1][j] + optimizer(Profiles[minX].profile[i-1], indelP, columns);
            del2 = M.scores[i][j-1] + optimizer(indelP, Profiles[minY].profile[j-1], columns);
            (min > del2) && (min = del2) && (type='2');
            (min > del1) && (min = del1) && (type='1');
            M.type[i][j] = type;
            M.scores[i][j] = min;
        }

    }
    vector<string> mergedA = mergeProfiles(Profiles[minX].alignment, Profiles[minY].alignment, M);
    vector<string> mergedNames;
    for(unsigned int k=0; k<Profiles[minX].seqName.size();k++){
        mergedNames.push_back(Profiles[minX].seqName[k]);
    }
    for(unsigned int k=0; k<Profiles[minY].seqName.size();k++){
        mergedNames.push_back(Profiles[minY].seqName[k]);
    }
    vector< vector<float> > mergedScore(mergedA[0].length(), vector<float> (columns));
    const char *seqMap = DNA;
    if (columns==21){
        seqMap = AA;
    }
    for (unsigned int i=0; i<mergedA[0].size(); i++){
        for (unsigned int j=0; j<mergedA.size(); j++){
            //string seq = mergedA[j][i];
            for (int x=0; x<columns; x++){
                const char *p = std::find(seqMap, seqMap+sizeof(seqMap), mergedA[j][i]);
                if (p!=seqMap+sizeof(seqMap)){
                    int dist = std::distance(seqMap, p);
                    //SP.incrementCount(dist, x);
                    mergedScore[i][dist]=1;
                }
                else{
                    //std::cout << "[Error]:: Found" << seq[i]<<std::endl;
                    //std::cout<<"ERROR"<<std::endl;
                }
            }
        }
    }
    int total = mergedA.size();
    for (unsigned int i=0;i<mergedScore.size();i++){
        for(unsigned int j=0; j<mergedScore[i].size();j++){
            mergedScore[i][j]/=total;
        }
    }
    Profile MP(mergedNames, total, columns, mergedScore, "ddsd", mergedA);
    return MP;

}

void ProfileAligner(vector<Profile> &P){
    int size = P.size();
    while(P.size()>1){
        DistanceMatrix DM = calculateDistanceMatrix(P);
        int minX = -1;
        int minY = -1;
        float minDist = 332000;

        for(int i=0; i<DM.getRows(); i++){
            for(int j=0; j<DM.getColumns(); j++){
                if(i!=j && DM.getValue(i,j)<minDist){
                    minDist = DM.getValue(i, j);
                    minX = i;
                    minY = j;
                }
            }
        }

        Profile new_profile = aligner(P, minX, minY);
        //Delete row column X
        P.erase(P.begin()+minX);
        if(minX<minY)
            minY=minY-1;
        //Delete row,column Y
        P.erase(P.begin()+minY);
        P.push_back(new_profile);
    }

    for(int i=0;i<size; i++){
        std::cout<<"Sequence "<<i+1<<" \t ";
        std::cout<<P[0].alignment[i]<<std::endl;
    }
}

