#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "fasta_reader.hpp"
#include "global_alignment.hpp"

class GlobalAlignmentTestFixture {
    private:
        static char input_file[];
    public:
        int sequenceProcessor(){
            int MATCH = 2;
            int MISMATCH = -1;
            int INDEL = -2;


            std::vector<Fasta> fasta_sequences;

            fasta_sequences = FastaReader(input_file);
            std::string seq1 = fasta_sequences[0].get_seqString();
            std::string seq2 = fasta_sequences[1].get_seqString();
            ScoringMatrix SM = createScoringMatrixFromSequences(seq1, seq2);

            performGlobalAlignment(SM, MATCH, MISMATCH, INDEL, seq1, seq2);
            vector<string> seqOutput = getOptimalAlignment(SM, seq1, seq2);
            int score = getOptimalScore(SM);
            return score;
    }
};
char GlobalAlignmentTestFixture::input_file[] = "../tests/data/test.fasta";
TEST_CASE_METHOD(GlobalAlignmentTestFixture, "Test with score 3", "[create]" ) {
        REQUIRE( sequenceProcessor() == 3 );
}
