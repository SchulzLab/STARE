//
// Created by Dennis Hecker on 09.07.21.
//
# include <iostream>
# include <string>
# include <vector>
# include <getopt.h>
# include <ctype.h>
# include <sstream>
# include <fstream>
# include <chrono>
# include "STARE_MiscFunctions.h"

/*
 * Used if one wants to reshape the Gene-TF matrix that the TF_Gene_Scorer produces per cell into a Cell-TF matrix per
 * gene (horizontally concatenated), for example as input for GAZE. Will only be executed if the respective flag is set
 * when calling STARE (-z). Outputs a metadata file with the order of genes and TFs to index the columns.
 *
 * g++ Reshape_toCellTF_perGene.cpp STARE_MiscFunctions.cpp -std=c++11 -o Reshape_toCellTF_perGene
 * ./Reshape_toCellTF_perGene -i input_folder -o output_folder
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 */


int main(int argc, char **argv) {
    using namespace std;
    auto start0 = chrono::high_resolution_clock::now();

    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-i folder with all files to reshape\n-o output-folder where to write the new files to\n";

    vector <string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
    for (int i = 1; i < argc; ++i) {
        string current_arg;
        current_arg = argv[i];
        for (string h_flag : h_flags) {
            if (!current_arg.compare(h_flag)) {
                cout << parameter_help << endl;
                return 1;
            }
        }
    }

    // Reading the input args:
    string i_input_folder, o_output_folder, c_cores;
    static struct option long_options[] =
            {
                    {"i",   required_argument, NULL, 'i'},
                    {"o",   required_argument, NULL, 'o'},
                    {"c",   required_argument, NULL, 'c'},
                    {NULL,  0,                 NULL, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "i:o:c:", long_options, NULL)) != -1) {
        switch (arg) {
            case 'i':
                i_input_folder = optarg;
                break;
            case 'o':
                o_output_folder = optarg;
                break;
            case 'c':
                c_cores = optarg;
                break;
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for each parameter.
    for (string p: {i_input_folder, o_output_folder}) {
        if (p.size() < 1) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }
    int cores = stoi(SetOptionalInput(c_cores, "1"));

    vector <string> input_files;  // Prepare a vector with the file names.
    vector <string> input_suffix;  // Store the file suffix for the row indices.
    string run_prefix;
    string ls_out = GetStdoutFromCommand("ls " + i_input_folder + "/");
    istringstream read_ls_out(ls_out);
    for (string line; getline(read_ls_out, line);) {
        if (line.substr(line.size() - 3) == ".gz") {
            input_files.push_back(i_input_folder + "/" + line);
            string split_string = "TF_Gene_Affinities_";
            int suffix_start = line.find(split_string);
            int suffix_end = line.find(".txt.gz");
            input_suffix.push_back(line.substr(suffix_start + split_string.size(), suffix_end - (suffix_start + split_string.size())));
            run_prefix = line.substr(0, suffix_start);
        }
    }
    if (input_files.size() == 0) {
        cout << "No gzip-files to reshape found in the input folder." << endl;
        return 1;
    }
    // Now we can go through each line of the input files and write a separate file for each gene, the gene order in
    // the files must be the same.
    ofstream meta_out(o_output_folder + "/" + run_prefix + "StackedMetadata.txt");
    ofstream reshape_out(o_output_folder + "/" + run_prefix + "StackedAffinityMatrix.txt");
    for (int f=0; f < input_files.size(); f++) {
        string row;
        int row_counter = 0;
        istringstream file_read(GetStdoutFromCommand(string("zcat < ") + input_files[f]));
        while (!file_read.eof()) {  // Use the first file as index helper.
            getline(file_read, row);
            if (f == 0 and row_counter == 0) {
                meta_out << row.substr(row.find("\t") + 1) + "\n";
            }
            if (row_counter > 0) {
                if (f == 0) {  // Write the geneID into the metadata file.
                    meta_out << row.substr(0, row.find("\t")) + "\t";
                }
                if (row_counter == 1) {
                    reshape_out << input_suffix[f];
                }
                reshape_out << "\t" + row.substr(row.find("\t") + 1);
            }
            row_counter++;
        }
        reshape_out << "\n";
    }
    reshape_out.close();
    GzipFile(o_output_folder + "/" + run_prefix + "StackedAffinityMatrix.txt");
    auto stop0 = chrono::high_resolution_clock::now();
    auto duration0 = chrono::duration_cast<chrono::seconds>(stop0 - start0);
    cout << duration0.count() << "s stacked matrices" << endl;
}



