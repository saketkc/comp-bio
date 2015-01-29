#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "fasta_reader.cpp"
#include "scoring_matrix.cpp"

class GlobalAlignmentTestFixture {
    public:
        static char input_file[];
        int distanceBasedProcessor(){
            int MATCH = 0;
            int MISMATCH = 1;
            int INDEL = 1;


            std::vector<Fasta> fasta_sequences;

            fasta_sequences = FastaReader(input_file);
            std::string seq1 = fasta_sequences[0].get_seqString();
            std::string seq2 = fasta_sequences[1].get_seqString();
            ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);

            performGlobalAlignment(SM, MATCH, MISMATCH, INDEL, seq1, seq2, true);
            vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
            int score = getOptimalScore(SM);
            return score;
    }
        int scoreBasedProcessor(){
            int MATCH = 2;
            int MISMATCH = -1;
            int INDEL = -2;


            std::vector<Fasta> fasta_sequences;

            fasta_sequences = FastaReader(input_file);
            std::string seq1 = fasta_sequences[0].get_seqString();
            std::string seq2 = fasta_sequences[1].get_seqString();
            ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);

            performGlobalAlignment(SM, MATCH, MISMATCH, INDEL, seq1, seq2, false);
            vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
            int score = getOptimalScore(SM);
            return score;
    }
        int distanceBasedKBandProcessor(){
            int dMATCH = 0;
            int dMISMATCH = 1;
            int dINDEL = 1;


            std::vector<Fasta> fasta_sequences;

            fasta_sequences = FastaReader(input_file);
            std::string seq1 = fasta_sequences[0].get_seqString();
            std::string seq2 = fasta_sequences[1].get_seqString();
            ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);
            int k=1;
            int alphaK = 32000;
            int minDistance = (seq2.length()-k-1)*dMATCH + (2*(k+1)+abs(seq1.length()-seq2.length()))*dINDEL;
            std::cout << "alpha: "<< alphaK << "best distance: " << minDistance << std::endl;
            while(alphaK > minDistance ){
                performKBandAlignment(SM, k, dMATCH, dMISMATCH, dINDEL, seq1, seq2);
                k*=2;
                minDistance = (seq2.length()-k-1)*dMATCH + (2*(k+1)+abs(seq1.length()-seq2.length()))*dINDEL;
                alphaK = getOptimalScore(SM);
            }

            vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
            int score = getOptimalScore(SM);
            return score;
    }
};

//This fasta was created by deleting certain obvious lines from two copied sequences
//Seq1 length = 2210
//Seq2 length = 2175
//Since all are a perfect match: Expected Score: (2210)*2 +(2210-2175)*-2*2=4280
//indel ditance=1 so (2210-2175)*1=35
//
char GlobalAlignmentTestFixture::input_file[] = "./data/align-1000bp-deletions.fasta";
TEST_CASE_METHOD(GlobalAlignmentTestFixture, "Test with score deletions", "[create]" ) {
        REQUIRE( distanceBasedProcessor() == 35 );
        REQUIRE( scoreBasedProcessor() == 4280 );
        REQUIRE( distanceBasedKBandProcessor() == 35);
}


