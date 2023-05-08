//
// Created by Dennis Hecker on 24.04.23.
//
# include <iostream>
# include <string>
# include <vector>
# include <limits>
# include <chrono>
# include <algorithm>
# include <unordered_map>
# include <unordered_set>
# include <getopt.h>
# include <stdio.h>
# include "STARE_MiscFunctions.h"

/*
 * Writes a bed-file with all unique promoter windows of size -w around the 5'TSS of all genes in a gtf-file.
 *
 * Compilation:
 * g++ Promoter_Windows.cpp STARE_MiscFunctions.cpp -std=c++11 -O3 -o Promoter_Windows}
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 */

class gene_information {
public:
    std::string chr;
    int tss;
};

int main(int argc, char **argv) {
    using namespace std;
    auto start0 = chrono::high_resolution_clock::now();
    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "\n-a gtf gene annotation"
                            "\n-u file with geneIDs/symbols to limit the output to"
                            "\n-w genewindow-size"
                            "\n-o prefix with which the files are written";

    vector <string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
    for (int i = 1; i < argc; ++i) {
        string current_arg;
        current_arg = argv[i];
        for (const string &h_flag : h_flags) {
            if (current_arg == h_flag) {
                cout << parameter_help << endl;
                return 1;
            }
        }
    }

    // Reading the input args:
    string a_promoterfile, u_genefile, w_genewindowsize, o_prefix;
    static struct option long_options[] =
            {
                    {"a",     required_argument, nullptr, 'a'},
                    {"u",     required_argument, nullptr, 'u'},
                    {"w",     required_argument, nullptr, 'w'},
                    {"o",     required_argument, nullptr, 'o'},
                    {nullptr, 0,                 nullptr, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "a:u:w:o:", long_options, nullptr)) != -1) {
        switch (arg) {
            case 'a':
                a_promoterfile = optarg;
                break;
            case 'u':
                u_genefile = optarg;
                break;
            case 'w':
                w_genewindowsize = optarg;
                break;
            case 'o':
                o_prefix = optarg;
                break;
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for each parameter.
    for (const string &p: {a_promoterfile, o_prefix}) {
        if (p.empty()) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    // If not given as input, set to default.
    // stod() works with scientific notation.
    int gene_windowsize = static_cast<int>(stod(SetOptionalInput(w_genewindowsize, "5000"))) / 2;

    // ____________________________________________________________
    // READ GENES TO FILTER FOR
    // ____________________________________________________________
    unordered_set <string> filter_genes;
    if (u_genefile.length() != 0 and u_genefile != "0") {
        ifstream read_genefile(u_genefile);
        if (!read_genefile) {
            cout << "ERROR could not open gene file\n" << u_genefile << endl;
            return 1;
        }
        string u_row;
        while (!read_genefile.eof()) {
            getline(read_genefile, u_row);
            filter_genes.insert(u_row);
        }
        read_genefile.close();
    }

    // ____________________________________________________________
    // PROCESS GTF GENE ANNOTATION
    // ____________________________________________________________
    cout << "Writing promoter windows" << endl;
    unordered_map<string, int> id_ensg_map;  // To be able to check if already saw a gene in the annotation.
    int gene_counter = 0;
    unordered_map<int, gene_information> promoter_map;  // Map of class to store id, name, chr and set of TSS.
    id_ensg_map.reserve(60000);
    promoter_map.reserve(60000);
    unordered_map <string, vector<int>> hic_boundaries;  // Store the min/max boundaries for each chromosome.
    FILE *Read_Gene_Annotation;
    if (a_promoterfile.substr(a_promoterfile.size() - 3, 3) == ".gz") {
        string promoter_cmd = string("zcat < ") + a_promoterfile + " 2>&1";
        Read_Gene_Annotation = popen(promoter_cmd.c_str(), "r");
    } else {
        Read_Gene_Annotation = fopen(a_promoterfile.c_str(), "rb");
    }
    int annot_len = FilePeek(a_promoterfile) * 5;
    char annot_buffer[annot_len];
    string row;
    bool first_line = true;  // To only look into the first line without the need to count lines.
    bool gene_chr_prefix;
    while (!feof(Read_Gene_Annotation)) {
        if (fgets(annot_buffer, annot_len, Read_Gene_Annotation) != nullptr) {
            row = annot_buffer;
            if ((row.size() > 1) and (row.rfind('#', 0) != 0)) {
                vector <string> columns = SplitTabLine(row);
                if (columns[2] == "gene") {
                    string id_delimiter = "gene_id ";
                    size_t id_start = columns[8].find(id_delimiter, 0) + id_delimiter.size();
                    int id_end = columns[8].find(';', id_start);
                    string gene_id = columns[8].substr(id_start + 1,
                                                       id_end - id_start - 2);  // Remove the quotation marks.
                    string name_delimiter = "gene_name ";
                    size_t name_start = columns[8].find(name_delimiter, 0) + name_delimiter.size();
                    int name_end = columns[8].find(';', name_start);
                    string gene_name = columns[8].substr(name_start + 1, name_end - name_start -
                                                                         2);  // Remove the quotation marks.

                    // Skip if it's not in -u.
                    if (u_genefile.length() != 0 and u_genefile != "0") {
                        if ((filter_genes.count(gene_id) == 0) and (filter_genes.count(gene_name) == 0)) {
                            continue;
                        }
                    }

                    string strand = columns[6];
                    string chr = columns[0];
                    if (first_line) {
                        if (chr.substr(0, 3) == "chr") {
                            gene_chr_prefix = true;
                        } else {
                            gene_chr_prefix = false;
                        }
                        first_line = false;
                    }
                    if (gene_chr_prefix) {
                        chr = chr.substr(3);
                    }
                    int gene_start;
                    if (strand == "+") {
                        gene_start = stoi(columns[3]);
                    } else {
                        gene_start = stoi(columns[4]);
                    }

                    auto counter_iter = id_ensg_map.find(gene_id);  // Check if we had this gene already.
                    if (counter_iter != id_ensg_map.end()) {
                    // Update if the new one is more in 5', must be on the same chromosome if non-unique.
                    // If the tss_identifier is gene, we change the set in the gene_information class.
                        if (((promoter_map[counter_iter->second].tss > gene_start) and (strand == "+") and
                             (promoter_map[counter_iter->second].chr == chr)) or
                            ((promoter_map[counter_iter->second].tss < gene_start) and (strand == "-") and
                             (promoter_map[counter_iter->second].chr == chr))) {
                            promoter_map[counter_iter->second].tss = gene_start;
                        }
                    } else {
                        id_ensg_map[gene_id] = gene_counter;
                        gene_information this_gene;
                        this_gene.chr = chr;
                        this_gene.tss = gene_start;
                        promoter_map[gene_counter] = this_gene;
                        gene_counter++;
                    }
                }
            }
        }
    }
    pclose(Read_Gene_Annotation);

    // Create a set of the strings to avoid identical windows for genes sharing a TSS.
    unordered_set <string> promoter_windows;
    for (auto const &gene : promoter_map) {  // gene_id, name, chr, TSS.
        // Gtf is 1-based, and bedtools 0-based and end-exclusive.
        int window_start = gene.second.tss - gene_windowsize - 1;
        int window_end = gene.second.tss + gene_windowsize;
        if (window_start < 0) {
            window_start = 0;
        }
        promoter_windows.insert(gene.second.chr + "\t" + to_string(window_start) + "\t" + to_string(window_end));
    }

    string temp_window_file = o_prefix + "_STARE_PromoterWindows.bed";
    ofstream window_out(temp_window_file);
    Test_outfile(window_out, temp_window_file);
    for (const string &promoter : promoter_windows) {
        window_out << promoter << "\n";
    }
    window_out.close();
}
