#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <cstdlib>
#include <fstream>


using namespace seqan;


typedef FastFMIndexConfig<void, size_t, 2, 0> TFastConfig;

using read_pos_t = uint16_t;
using read_id_t = unsigned;
using pos_vector_t = std::vector<read_pos_t>;
using read2pos_map_t = std::map<int,std::vector<read_pos_t> >;
using index_t = Index<StringSet<DnaString>, BidirectionalIndex<FMIndex<void,TFastConfig> > >;


const auto boot_time = std::chrono::steady_clock::now();

// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}


void createIndex(std::string & filename, std::string & out){

    print("PARSING FASTA FILE");  

    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(filename));
    readRecords(ids, seqs, seqFileIn);

    print("DONE");

    // INDEX CREATION ----------    
    print("CREATING INDEX");
    index_t index(seqs);

    indexCreate(index);

    print("DONE");
    print("SAVING USING OUTPUT PREFIX");
    save(index,toCString(out));

    // print("CLEARING INDEX");
    // clear(index);
    print("DONE");

    // print("RELOADING INDEX");
    // if(open(index,out.c_str())){
    //     print("DONE");  
    // }
    // else{
    //     print("ERROR");
    // }
 
    // auto delegate= [&](auto & iter, const DnaString & needle, int errors)
    // {
        
    //     for (auto occ : getOccurrences(iter)){
    //         std::cout << needle << " " << getValueI1(occ) << " " << getValueI2(occ)  << std::endl;
    //     }
    // };

    // StringSet<DnaString> kmSet;
    // appendValue(kmSet, DnaString("ATCG"));
    // appendValue(kmSet, DnaString("GTTA"));

    // find<0,0>(delegate, index, kmSet, HammingDistance());

}




int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("create_seq_index");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "output filename"));

    // Parse command line.    &
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::string inFile;
    std::string outFile;
    getArgumentValue(inFile,  parser, 0);
    getArgumentValue(outFile, parser, 1);
    std::cout << outFile << std::endl;
    
    // Create index
    createIndex(inFile, outFile);

    }
