#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "fasta_reader.cpp"

class FastaReaderTestFixture{
    public:
        static char input_file[];
        int numberOfSequences(){
            std::vector<Fasta> fasta_sequences;
            fasta_sequences = FastaReader(input_file);
            return fasta_sequences.size();
        }
};

char FastaReaderTestFixture::input_file[] = "./data/multiple_seq2.fasta";
TEST_CASE_METHOD(FastaReaderTestFixture, "Test Fasta Reader", "[create]"){
    REQUIRE(numberOfSequences()==5);
}
