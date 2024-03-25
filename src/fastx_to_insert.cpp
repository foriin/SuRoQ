#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <unordered_map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

using namespace std;

#define GZIP_MAGIC "\037\213"
#define BZIP2_MAGIC "BZ"

int main(int argc, char** argv) {
    string usage = R"(
usage:
	fastxToInsert input.[fq|fa][.gz|.bz2] output.insert
	
	)";
    if (argc != 3) {
        cerr << usage;
        return EXIT_FAILURE;
    }
    string input_fastx_file{argv[1]};
    std::ios::sync_with_stdio(false);
    istream* p_ist_in{&std::cin};
    if (input_fastx_file != "stdin" && input_fastx_file != "-") {
        p_ist_in = new std::ifstream{input_fastx_file};
    }
    if (!(*p_ist_in)) {
        cerr << "error: cannot open file " << argv[1] << "for reading" << endl;
        return EXIT_FAILURE;
    }
    boost::iostreams::filtering_istream in;
    char magic_number[3];
    p_ist_in->get(magic_number, 3);
    if (memcmp(magic_number, GZIP_MAGIC, 2) == 0) {
        in.push(boost::iostreams::gzip_decompressor());
    } else if (memcmp(magic_number, BZIP2_MAGIC, 2) == 0) {
        in.push(boost::iostreams::bzip2_decompressor());
    }
    p_ist_in->seekg(0, ios::beg);
    in.push(*p_ist_in);
    ofstream out{argv[2]};
    if (!out) {
        cerr << "error: cannot open file " << argv[2] << " for writing" << endl;
        return EXIT_FAILURE;
    }
    string line, sequence;
    unordered_map<string, int> counter;
    bool isFasta = false;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '@' || line[0] == '>') {
            isFasta = (line[0] == '>');
            if (!getline(in, sequence)) break;
            ++counter[sequence];
        }
        if (!isFasta) {
            // Skip next two lines for FASTQ
            in.ignore(numeric_limits<streamsize>::max(), '\n');
            in.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    }
    for (const auto& seq : counter) {
        out << seq.first << '\t' << seq.second << '\n';
    }
    out.close();
}

