//
// Created by Dennis Hecker on 21.05.21.
//

# include <iostream>
# include <string>
# include <fstream>
# include <sstream>
# include <unordered_set>
# include <getopt.h>
# include <vector>

/* Goes through the extracted sequences in a file and replaces the invalid characters with an N. It additionally writes
 * a txt-file with the length of the longest sequence in the file. Required for TRAP to read the sequence file as
 * stream.
 *
 * g++ ReplaceInvalidChars.cpp -std=c++11 -O3 -o ReplaceInvalidChars
 *  ./ReplaceInvalidChars -i file -o output -d path_for_maxRox_file
 *
 *  Part of STARE: https://github.com/SchulzLab/STARE
 */

int main(int argc, char **argv) {
    using namespace std;

    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-i input sequence file\n-o output file\n-d file to write the maximal row_length to";

    vector<string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
    for (int i = 1; i < argc; ++i) {
        string current_arg;
        current_arg = argv[i];
        for (string h_flag : h_flags) {
            if (!current_arg.compare(h_flag)) {
                cout << parameter_help;
                return 1;
            }
        }
    }
    if (argc < 5){
        cout << parameter_help;
        return 1;
    }

    // Reading the input args:
    string i_input_file, o_output_file, d_row_file;
    static struct option long_options[] =
            {
                    {"i", required_argument, NULL, 'i'},
                    {"o",   required_argument, NULL, 'o'},
                    {"d",   required_argument, NULL, 'd'},
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "i:o:d:", long_options, NULL)) != -1) {
        switch (arg) {
            case 'i':
                i_input_file = optarg;
                break;
            case 'o':
                o_output_file = optarg;
                break;
            case 'd':
                d_row_file = optarg;
                break;
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for the required parameters.
    for (string p: {i_input_file, o_output_file, d_row_file}){
        if (p.size() < 1){
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    unordered_set<char> allowed_chars = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
    string row;
    int max_len = 0;
    ifstream Read_Sequence(i_input_file);
    ofstream write_output(o_output_file);
    if (Read_Sequence) {
        while (!Read_Sequence.eof()) {
            getline(Read_Sequence, row);
            if (row.substr(0, 1) == ">") {
                write_output << row << "\n";
            } else {
                if (row.size() > max_len) {
                    max_len = row.size();
                }
                for (char &c : row) {
                    if (allowed_chars.find(c) != allowed_chars.end()) {
                        write_output << c;
                    } else {
                        write_output << "N";
                    }
                }
                write_output << "\n";
            }
        }
    }
    else {
        cout << "ERROR could not open sequence file to filter invalid characters:\n" << i_input_file << endl;
        return 1;
    }

    // Write file with the maximal row length, needed by TRAP in the subsequent step.
    ofstream write_row_len(d_row_file);
    write_row_len << to_string(max_len + 10);  // Plus a small buffer, just in case.
}
