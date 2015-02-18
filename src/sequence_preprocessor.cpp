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
    unsigned int count = 5;//izeof(DNA);
    unsigned int rows = fasta_sequences.size();
    const char *seqMap = DNA;
    std::string seq;
    if (!isDNA){
        count = 21;//AA.size();
        seqMap = AA;
    }
    //std::cout<<"SEQMAP: "<<std::begin(seqMap) << std::endl;
   //const char* end = seqMap + sizeof(seqMap) / sizeof(seqMap[0]);
   //std::cout <<"END: "<<end<<std::endl;
    //std::vector<DistanceMatrix> v;
    std::vector<Profile> v;
    for(int i=0;i<rows;i++){
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
                    std::cout << "[Error]:: Found" << seq[i]<<std::endl;
                }
            }
        }
        vector<string> alignment;
        alignment.push_back(seq); //Alignment is same for sequences as
        Profile profile(fasta_sequences[i].seqName, seq.length(), count, score, seq, alignment);
        v.push_back(profile);
    }
   return v;
}

void printProfile(vector<Profile> &profiles){
    //std::cout<< "Profiles: "<<profiles.size();
    for (unsigned int i=0;i<profiles.size();i++){
        std::cout<<"profile: "<<profiles[i].seqName<<std::endl;
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
//    std::cout<<"Seq1 "<<p1.sequence<<std::endl;
//    std::cout<<"Seq2 "<<p2.sequence<<std::endl;
//    std::cout<<"Seq1 L "<<p1.columns<<std::endl;
//    std::cout<<"Seq2 L "<<p2.columns<<std::endl;
//    std::cout<<"Distance "<<distance<<std::endl;
//    std::cout<<"XXX: "<<(distance*1.0/comparisons)<<std::endl;
    return distance*1.0/comparisons;
}

Profile performAlignment(vector<Profile> &P, int &minX, int &minY){

    int cols1 = P[minX].columns;
    int cols2 = P[minY].columns;
    Profile seq1Profile = P[minX];
    Profile seq2Profile = P[minY];
    ProfileAlignment PA(cols1, cols2);
    int s=0;
    int del1=0;
    int del2=0;
    char type='M';
    float min=32000;
    for(int i=1; i<cols1; i++){
        for(int j=1; j<cols2; j++){
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
            min = PA.scores[i-1][j-1] + s;
            del1 = PA.scores[i-1][j] + 2;
            del2 = PA.scores[i][j-1] + 2;
            (min > del2) && (min = del2) && (type='2');
            (min > del1) && (min = del1) && (type='1');
            PA.type[i][j] = type;
            PA.scores[i][j] = min;
        }
    }
    //Bakctracing
    //

}
DistanceMatrix calculateDistanceMatrix(vector<Profile> &P){
    DistanceMatrix DM(P.size(), P.size());
    float dist;
    for (unsigned int i=0; i<P.size(); i++){
        for (unsigned int j=0; j<P.size(); j++){
            dist = calculateHammingDistance(P[i], P[j]);
            //std::cout<<"YOOOOO: "<<dist<<std::endl;
            DM.edit(i,j, dist);
        }
    }
    std::cout<<"DistanceMatrix "<<std::endl;
    for (unsigned int i=0; i<P.size(); i++){
        std::cout<<std::endl;
        for (unsigned int j=0; j<P.size(); j++){
            std::cout<<DM.getValue(i,j)<<" ";
        }
    }
    std::cout<<std::endl;
    return DM;
}


ProfileAlignment calculatePairwiseAlignment(Profile &seq1Profile, Profile &seq2Profile){
    int cols1 = seq1Profile.columns;
    int cols2 = seq2Profile.columns;

    ProfileAlignment P(cols1, cols2);

    float del1 = 0;
    float del2 = 0;
    float min =0;
    char type='M';
    for(int i=1; i<cols1; i++){
        for(int j=1; j<cols2; j++){
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
    std::cout<<"Post adjustment"<<std::endl;
    for(int i=0; i<cols1; i++){
        std::cout<<std::endl;
        for(int j=0; j<cols2; j++){
            std::cout<<P.scores[i][j]<< " ";
        }
    }
    for(int i=0; i<cols1; i++){
        std::cout<<std::endl;
        for(int j=0; j<cols2; j++){
            std::cout<<P.type[i][j]<< " ";
        }
    }
    std::cout<<std::endl<<"DDDDD"<<std::endl;
    return P;
}


long optimizer(vector<float> seq1score, vector<float> seq2score, int max){
    float s=0;
    int aln1Size = seq1score.size();
    int aln2Size = seq2score.size();
    for(unsigned int i=0;i<aln1Size; i++){
        for(unsigned int j=0;j<aln2Size;j++){
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
    std::cout<<"Loimit x:"<<M.seq1Length;
   std::cout<<std::endl;
    std::cout<<"Loimit y:"<<M.seq2Length;
 //   std::cout<<std::endl;
 //  std::cout<<"Loimit x:"<<seq1Alignment.size();
   //std::cout<<std::endl;
   //std::cout<<"Loimit y:"<<seq2Alignment.size();
   std::cout<<std::endl;
    M.print();
    while (seq1Size >0 || seq2Size > 0){
        std::cout<<"Se1 Aln size: "<<seq1Size<<std::endl;
        std::cout<<"Se2 Aln size: "<<seq2Size<<std::endl;
        if(M.type[seq1Size][seq2Size]=='M'){
            std::cout<<"CHECK: M"<<std::endl;
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
            std::cout<<"CHECK: 2"<<std::endl;
            for(unsigned int i=0; i<seq1Alignment.size(); i++){
                aln1[i].append(seq1Alignment[i].substr(seq1Size-1,1));
            }
            for(unsigned int j=0; j<seq2Alignment.size(); j++){
                aln2[j].append("-");
            }
            --seq2Size;
        }
        else if(M.type[seq1Size][seq2Size]=='1'){

            std::cout<<"CHECK: 1"<<std::endl;
            for(unsigned int i=0; i<seq1Alignment.size(); i++){
                aln1[i].append("-");
            }
            for(unsigned int j=0; j<seq2Alignment.size(); j++){
                aln2[j].append(seq2Alignment[j].substr(seq2Size-1,1));
            }
            --seq1Size;
        }
        else {
            std::cout<<"ERRRRRRRRRRRRRRRRRRRR: "<<M.type[seq1Size][seq2Size]<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    for(int i=0;i<aln1.size();i++){
        reverse(aln1[i].begin(), aln1[i].end());
    }
    for(int j=0;j<aln2.size();j++){
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
        int min=0;

        for (int j=1;j<aln2Size;j++){
            min = M.scores[i-1][j-1] + optimizer(Profiles[minX].profile[i-1], Profiles[minY].profile[j-1], columns);
            del1 = M.scores[i-1][j] + optimizer(Profiles[minX].profile[i-1], indelP, columns);
            del2 = M.scores[i][j-1] + optimizer(indelP, Profiles[minY].profile[j-1], columns);
            (min > del2) && (min = del2) && (type='1');
            (min > del1) && (min = del1) && (type='2');
            M.type[i][j] = type;
            M.scores[i][j] = min;
        }

    }
    std::cout<<"Starting merging profiles"<<std::endl;
    vector<string> mergedA = mergeProfiles(Profiles[minX].alignment, Profiles[minY].alignment, M);
    std::cout<<"Completed merging profiles"<<std::endl;
    vector< vector<float> > mergedScore(mergedA[0].length(), vector<float> (columns));
    const char *seqMap = DNA;
    if (columns==21){
        seqMap = AA;
    }
    for (unsigned int i=0; i<mergedA[0].size(); i++){
        for (unsigned int j=0; j<mergedA.size(); j++){
            //string seq = mergedA[j][i];
            for (unsigned int x=0; x<columns; x++){
                const char *p = std::find(seqMap, seqMap+sizeof(seqMap), mergedA[j][i]);
                if (p!=seqMap+sizeof(seqMap)){
                    int dist = std::distance(seqMap, p);
                    //SP.incrementCount(dist, x);
                    mergedScore[i][dist]=1;
                }
                else{
                    //std::cout << "[Error]:: Found" << seq[i]<<std::endl;
                    std::cout<<"ERROR"<<std::endl;
                }
            }
        }
    }
    int total = mergedA.size();
    for (int i=0;i<mergedScore.size();i++){
        for(int j=0; j<mergedScore[i].size();j++){
            mergedScore[i][j]/=total;
        }
    }
    Profile MP("test", total, columns, mergedScore, "ddsd", mergedA);
    return MP;

}


vector<string> getOptimalProfileAlignment(ProfileAlignment P, string &seq1, string &seq2){

    int seq1_length = P.seq1Length-2;
    int seq2_length = P.seq2Length-2;
    //std::cout<<"seq1 PP: "<<seq1_length<<std::endl;
    //std::cout<<"seq2 PP: "<<seq2_length<<std::endl;
    vector<string> output;
    std::string seq1Output = "";
    std::string seq2Output = "";
    int score;
    char type='X';
    std::cout<<"seq1 L: "<<seq1_length<<std::endl;
    std::cout<<"seq2 L: "<<seq2_length<<std::endl;

    std::cout<<"seq1 : "<<seq1<<std::endl;
    std::cout<<"seq2 : "<<seq2<<std::endl;
    while (seq1_length > 0  || seq2_length > 0 ){
        score  = P.scores[seq1_length][seq2_length];
        type = P.type[seq1_length][seq2_length];
    std::cout<<"seq1 L: "<<seq1_length<<std::endl;
    std::cout<<"seq2 L: "<<seq2_length<<std::endl;

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
            std::cerr<<" Seq1: " << seq1_length << " seq2_length: " << seq2_length<<std::endl;
            //exit(EXIT_FAILURE);
        }
    std::cout<<"seq1 : "<<seq1Output<<std::endl;
    std::cout<<"seq2 : "<<seq2Output<<std::endl;

    }

    output.push_back(seq1Output);
    output.push_back(seq2Output);
    return output;
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

    for(unsigned int i=0;i<size; i++){
        std::cout<<P[0].alignment[i]<<std::endl;
    }
}

