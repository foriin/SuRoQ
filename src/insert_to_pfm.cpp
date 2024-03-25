#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

const int NUCLEOTIDE_LENGTH = 18;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " file.insert\n";
        return EXIT_FAILURE;
    }

    std::ifstream file(argv[1]);
    if (!file) {
        std::cerr << "Error opening file.\n";
        return EXIT_FAILURE;
    }

    // Matrix to store nucleotide counts: [nucleotide][position]
    std::map<char, std::vector<int>> freq_matrix = {
        {'A', std::vector<int>(NUCLEOTIDE_LENGTH, 0)},
        {'C', std::vector<int>(NUCLEOTIDE_LENGTH, 0)},
        {'G', std::vector<int>(NUCLEOTIDE_LENGTH, 0)},
        {'T', std::vector<int>(NUCLEOTIDE_LENGTH, 0)}
    };

    std::string sequence;
    int count;
    while (file >> sequence >> count) {
        for (int i = 0; i < NUCLEOTIDE_LENGTH && i < sequence.length(); ++i) {
            char nucleotide = sequence[i];
            if (freq_matrix.find(nucleotide) != freq_matrix.end()) {
                freq_matrix[nucleotide][i] += count;
            }
        }
    }

    file.close();

    // Output the frequency matrix
    for (const auto& pair : freq_matrix) {
        std::cout << pair.first << "\t";
        for (int val : pair.second) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}

