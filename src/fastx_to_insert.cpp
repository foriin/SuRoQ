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

bool isValidExtension(const std::string& filename) {
	size_t dotPos = filename.find_last_of(".");
	if (dotPos == std::string::npos) {
		return false;
	}
	
	std::string ext = filename.substr(dotPos + 1);
	std::string baseExt; // Base extension without possible compression extension

	// Check for compression extension
	if (ext == "gz" || ext == "bz2") {
		size_t secondDotPos = filename.find_last_of(".", dotPos - 1);
		if (secondDotPos == std::string::npos) {
			return false;
		}
		baseExt = filename.substr(secondDotPos + 1, dotPos - secondDotPos - 1);
	} else {
		baseExt = ext;
	}

	// Convert baseExt to lowercase just in case
	std::transform(baseExt.begin(), baseExt.end(), baseExt.begin(),
			[](unsigned char c) { return std::tolower(c); });
	// Check if baseExt if one of the supported extensions
	return baseExt == "fq" || baseExt == "fastq" || baseExt == "fa" || baseExt == "fasta";
}

bool hasRepeatingNucleotides(const std::string& sequence, int threshold) {
    if (sequence.length() < threshold) return false;
    
    char lastChar = sequence[0];
    int count = 1;
    
    for (size_t i = 1; i < sequence.length(); ++i) {
        if (sequence[i] == lastChar) {
            ++count;
            if (count >= threshold) return true;
        } else {
            lastChar = sequence[i];
            count = 1;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    string usage = R"(
usage:
	fastxToInsert input.[fastq|fq|fasta|fa][.gz|.bz2] output.insert
	
	)";
    if (argc != 3) {
        cerr << usage;
        return EXIT_FAILURE;
    }
    string input_fastx_file{argv[1]};
    if (!isValidExtension(input_fastx_file)) {
	    cerr << "Input file has an invalid extension!" << endl;
	    return EXIT_FAILURE;
    }
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
            // Skip sequences with 10 or more identical consecutive nucleotides
            if (hasRepeatingNucleotides(sequence, 10)) continue;
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

