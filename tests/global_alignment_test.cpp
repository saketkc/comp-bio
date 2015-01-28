#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "fasta_reader.cpp"
#include "scoring_matrix.cpp"

class GlobalAlignmentTestFixture {
    public:
        static char input_file[];
        int sequenceProcessor(){
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
};

//This fasta was created by deleting certain obvious lines from two copied sequences
//Seq1 length = 2210
//Seq2 length = 2175
//Since all are a perfect match: Expected Score: (2210)*2 +(2210-2175)*-2*2=4280
//indel ditance=1 so (2210-2175)*1=35
char GlobalAlignmentTestFixture::input_file[] = "./data/align-1000bp-deletions.fasta";
TEST_CASE_METHOD(GlobalAlignmentTestFixture, "Test with score deletions", "[create]" ) {
        REQUIRE( sequenceProcessor() == 35 );
}


