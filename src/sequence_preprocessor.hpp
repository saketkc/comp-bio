#include "sequence_profile.hpp"
using std::vector;

bool isDNASequence(Fasta firstSequence);
vector<DistanceMatrix> createProfileFromSequences(std::vector<Fasta> &fasta_sequences);
vector<Profile> createProfileFromSequences(std::vector<Fasta> &fasta_sequences);

