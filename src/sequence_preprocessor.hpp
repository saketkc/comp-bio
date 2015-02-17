#include "sequence_profile.hpp"
using std::vector;

bool isDNASequence(Fasta firstSequence);
//vector<DistanceMatrix> createProfileFromSequences(std::vector<Fasta> &fasta_sequences);
vector<Profile> createProfileFromSequences(std::vector<Fasta> &fasta_sequences);
void printProfile(std::vector<Profile> &PM);
ProfileAlignment calculatePairwiseAlignment(Profile &seq1Profile, Profile &seq2Profile);
vector<string> getOptimalProfileAlignment(ProfileAlignment P, string &seq1, string &seq2);

DistanceMatrix calculateDistanceMatrix(vector<Profile> &P);
float calculateHammingDistance(Profile &p1, Profile &p2);
void ProfileAligner(vector<Profile> &P);
Profile aligner(vector<Profile> &Profiles, int minX, int minY);
