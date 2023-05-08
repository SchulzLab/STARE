//
// Created by Dennis Hecker on 28.07.21.
//

# include <iostream>
# include <string>
# include <vector>
# include <array>
# include <chrono>
# include <cmath>
# include <algorithm>
# include <unordered_map>
# include <unordered_set>
# include <map>
# include <set>
# include <regex>
# include <getopt.h>
# include <omp.h>
# include "STARE_MiscFunctions.h"

/* This script takes the TF-region affinities that was calculated earlier and returns the TF-gene affinities, based on
 * all the input parameters. The region-gene mapping is either done with overlapping the regions with a defined gene
 * window or is taken from the ABC-scoring (if done in advance).
 *
 * Compilation for MacOS:
 * g++ TF_Gene_Scorer.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -O3 -std=c++11 -o TF_Gene_Scorer
 * Compilation for Linux:
 * g++ TF_Gene_Scorer.cpp STARE_MiscFunctions.cpp -fopenmp -O3 -std=c++11 -o TF_Gene_Scorer
 *
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 */

// To check endianness with an example. Only necessary when reshape True.
enum hl_endianness : uint32_t {
    HL_LITTLE_ENDIAN   = 0x00000001,
    HL_BIG_ENDIAN      = 0x01000000,
    HL_PDP_ENDIAN      = 0x00010000,
    HL_UNKNOWN_ENDIAN  = 0xFFFFFFFF
};


int main(int argc, char **argv) {
    using namespace std;

    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-a gene annotation file"
                            "\n-u file with geneIDs/symbols to limit the output to"
                            "\n-b bed file with open chromatin regions"
                            "\n-n column(s) with the activity of regions in the -b file"
                            "\n-i TF-region affinity file from TRAP"
                            "\n-o prefix of the output files"
                            "\n-p PSEM file"
                            "\n-w window size around TSS (default 50kb)"
                            "\n-e decay (default true)"
                            "\n-c number of cores"
                            "\n-abc path to the ABC-scoring output, if available, otherwise set to 0"
                            "\n-z write output into one large binary file with separate metadata";

    vector <string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
    for (int i = 1; i < argc; ++i) {
        string current_arg;
        current_arg = argv[i];
        for (const string& h_flag : h_flags) {
            if (current_arg == h_flag) {
                cout << parameter_help;
                return 1;
            }
        }
    }

    // Reading the input args:
    string a_gene_annotation, u_genefile, b_scalingfile, n_scalingcol, i_tf_region_affinity, o_output_prefix, p_psem_file,
            w_window_size, e_decay, c_cores, abc_abc_output, z_reshape;
    static struct option long_options[] =
            {  // optional_argument throws a Segmentation fault, for whatever reason.
                    {"a",   required_argument, nullptr, 'a'},
                    {"u",   required_argument, nullptr, 'u'},
                    {"b",   required_argument, nullptr, 'b'},
                    {"n",   required_argument, nullptr, 'n'},
                    {"i",   required_argument, nullptr, 'i'},
                    {"o",   required_argument, nullptr, 'o'},
                    {"p",   required_argument, nullptr, 'p'},
                    {"w",   required_argument, nullptr, 'w'},
                    {"e",   required_argument, nullptr, 'e'},
                    {"c",   required_argument, nullptr, 'c'},
                    {"abc", required_argument, nullptr, 't'},
                    {"z",   required_argument, nullptr, 'z'},
                    {nullptr,  0,                 nullptr, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "a:u:b:n:i:o:p:w:e:c:t:z:", long_options, nullptr)) != -1) {
        switch (arg) {
            case 'a':
                a_gene_annotation = optarg;
                break;
            case 'u':
                u_genefile = optarg;
                break;
            case 'b':
                b_scalingfile = optarg;
                break;
            case 'n':
                n_scalingcol = optarg;
                break;
            case 'i':
                i_tf_region_affinity = optarg;
                break;
            case 'o':
                o_output_prefix = optarg;
                break;
            case 'p':
                p_psem_file = optarg;
                break;
            case 'w':
                w_window_size = optarg;
                break;
            case 'e':
                e_decay = optarg;
                break;
            case 'c':
                c_cores = optarg;
                break;
            case 't':
                abc_abc_output = optarg;
                break;
            case 'z':
                z_reshape = optarg;
                break;
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for the required parameters.
    for (const string& p: {a_gene_annotation, i_tf_region_affinity, o_output_prefix, p_psem_file}) {
        if (p.empty()) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    // If not given as input, set to default.
    int window_size = static_cast<int>(stod(SetOptionalInput(w_window_size, "50000"))) / 2;
    bool decay = SetBoolInput(e_decay, true);
    bool reshape = SetBoolInput(z_reshape, false);
    int cores = stoi(SetOptionalInput(c_cores, "1"));

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
        string row;
        while (!read_genefile.eof()) {
            getline(read_genefile, row);
            filter_genes.insert(row);
        }
        read_genefile.close();
    }

    // ____________________________________________________________
    // READ GENE ANNOTATION
    // ____________________________________________________________
    unordered_map <int, string> gene_id_map;  // To avoid storing the full gene id string for each data structure.
    unordered_map<string, int> gene_name_map;  // From the ABC-files we don't have the ID.
    int gene_counter = 0;
    unordered_map <int, vector<string>> tss_map;  // GeneID: [chr, 5'TSS]
    tss_map.reserve(60000);
    FILE *Read_Gene_Annotation;
    if (a_gene_annotation.substr(a_gene_annotation.size() - 3, 3) == ".gz") {
        string promoter_cmd = string("zcat < ") + a_gene_annotation + " 2>&1";
        Read_Gene_Annotation = popen(promoter_cmd.c_str(), "r");
    } else {
        Read_Gene_Annotation = fopen(a_gene_annotation.c_str(), "rb");
    }
    int annot_len = FilePeek(a_gene_annotation) * 5;
    char annot_buffer[annot_len];
    string row;
    while (!feof(Read_Gene_Annotation)) {
        if (fgets(annot_buffer, annot_len, Read_Gene_Annotation) != nullptr) {
            row = annot_buffer;
            if ((row.size() > 1) and (row.rfind('#', 0) != 0)) {
                vector <string> columns = SplitTabLine(row);
                if (columns[2] == "gene") {
                    string id_delimiter = "gene_id ";
                    size_t id_start = columns[8].find(id_delimiter, 0) + id_delimiter.size();
                    int id_end = columns[8].find(';', id_start);
                    string gene_string = columns[8].substr(id_start + 1,
                                                       id_end - id_start - 2);  // Removing the quotation marks.
                    string name_delimiter = "gene_name ";
                    size_t name_start = columns[8].find(name_delimiter, 0) + name_delimiter.size();
                    int name_end = columns[8].find(';', name_start);
                    string gene_name = columns[8].substr(name_start + 1, name_end - name_start - 2);
                    if (u_genefile.length() != 0 and u_genefile != "0") {
                        if ((filter_genes.count(gene_string) == 0) and (filter_genes.count(gene_name) == 0)) {
                            continue;
                        }
                    }
                    string strand = columns[6];
                    string gene_start;
                    if (strand == "+") {
                        gene_start = columns[3];
                    } else {
                        gene_start = columns[4];
                    }
                    auto tss_iter = tss_map.find(gene_counter);
                    if (tss_iter != tss_map.end()) {  // Replace if the new one is more in 5'.
                        if (((stoi(tss_iter->second[1]) > stoi(gene_start)) and (strand == "+")) or
                            ((stoi(tss_iter->second[1]) < stoi(gene_start)) and (strand == "-"))) {
                            tss_map[gene_counter][1] = gene_start;
                        }
                    } else {
                        gene_id_map[gene_counter] = gene_string;
                        gene_name_map[gene_string] = gene_counter;
                        tss_map[gene_counter] = vector < string > {columns[0], gene_start};
                        gene_counter++;
                    }
                }
            }
        }
    }
    pclose(Read_Gene_Annotation);

    // ____________________________________________________________
    // GET MOTIF LENGTHS
    // ____________________________________________________________
    unordered_map<string, int> motif_length_map;
    int length_helper = 0;
    string current_tf;
    ifstream Read_PSEM(p_psem_file);
    if (!Read_PSEM) {
        cout << "ERROR Could not open the PSEM file\n" + p_psem_file << endl;
        return 1;
    }
    while (!Read_PSEM.eof()) {
        getline(Read_PSEM, row);
        if ((row.rfind('>', 0) == 0)) {
            motif_length_map[current_tf] = length_helper;
            vector <string> columns = SplitTabLine(row);
            current_tf = columns[1];
            length_helper = 0;
        } else if (row.size() > 1) {
            length_helper++;
        }
    }
    motif_length_map[current_tf] = length_helper;  // For the last entry.

    ifstream Read_Affinity(i_tf_region_affinity);
    getline(Read_Affinity, row);  // Only need the first line.
    size_t num_tfs = count(row.begin(), row.end(), '\t'); // Already get to have it available in all blocks.
    Read_Affinity.close();
    Read_Affinity.clear();

    // ____________________________________________________________
    // GET MAPPING OF REGIONS TO GENES AND REGION-TF AFFINITIES
    // ____________________________________________________________
    unordered_map <int, set<int>> gene_region_map; // Storing the peak ids in its window.
    unordered_map <int, string> peak_id_map;  // Instead of storing the whole peak string.
    unordered_map <string, int> peak_string_map;  // Reverse of the former.
    vector <string> peak_locations;  // For the temporary bed file.
    unordered_map <int, string> kicked_genes;  // Track the genes that don't have any regions in their window.

    // Fetch the region-tf affinities and do the intersection with the gene windows if no ABC-scoring file was provided.
    // Already allocate two vectors in that datastructure to later store potential region scalings.
    unordered_map < int, array < vector < double >, 2 >> region_affinities;
    Read_Affinity.open(i_tf_region_affinity);
    int row_count = -1;
    vector<int> motif_lengths;  // Already write a vector with the motif lengths according to the file header.
    string tf_header;
    if (!Read_Affinity) {
        cout << "ERROR Could not open TRAP's affinity file\n" + i_tf_region_affinity << endl;
        return 1;
    }
    while (!Read_Affinity.eof()) {
        getline(Read_Affinity, row);
        if (row_count == -1) {
            tf_header = row;
            vector <string> tf_names = SplitTabLine(row);
            for (int i = 1; i < tf_names.size(); i++) {
                motif_lengths.push_back(motif_length_map[tf_names[i]]);
            }
        }
        if ((row.size() > 1) and (row_count != -1)) {
            int first_tab = row.find('\t', 0);
            string region_string = row.substr(0, first_tab);
            if (region_string.substr(0, 3) == "chr") {
                region_string = region_string.substr(3);
            }
            int colon_pos = region_string.find(':', 0);
            int hyphon_pos = region_string.find('-', 0);
            int region_start = stoi(region_string.substr(colon_pos + 1, hyphon_pos - colon_pos - 1));
            int region_end = stoi(region_string.substr(hyphon_pos + 1));

            peak_locations.push_back(
                    region_string.substr(0, colon_pos) + "\t" + to_string(region_start) + "\t" + to_string(region_end) + "\t" + to_string(row_count));
            string col_val;
            vector<double> affinity_vals(motif_lengths.size());
            stringstream row_read(row);
            int col_count = 0;
            while (getline(row_read, col_val, '\t')) {
                if (col_count > 0) {
                    affinity_vals[col_count - 1] = (stod(col_val));  // -1 because of the row index.
                }
                col_count++;
            }
            peak_id_map[row_count] = region_string;
            peak_string_map[region_string] = row_count;
            region_affinities[row_count][0] = affinity_vals;
        }
        row_count++;
    }
    cout << "TRAP file read" << endl;

    // ____________________________________________________________
    // GET SCALING COLS - if specified
    // ____________________________________________________________
    vector<int> col_indices;
    if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
        // Look into the peakfile to see with how many columns/files we are dealing with.
        ifstream peek_peak(b_scalingfile);
        string line_peek;
        while (!peek_peak.eof()) {
            getline(peek_peak, line_peek);
            if (line_peek.substr(0, 1) != "#") {
                int total_cols = SplitTabLine(line_peek).size();
                if (n_scalingcol.find('-') != string::npos) {
                    int hyphon_pos = n_scalingcol.find('-', 0);
                    for (int i = stoi(n_scalingcol.substr(0, hyphon_pos)) - 1;
                         i <= stoi(n_scalingcol.substr(hyphon_pos + 1)) - 1; i++) {
                        col_indices.push_back(i);
                    }
                } else if (n_scalingcol.find('+') != string::npos) {
                    for (int i = stoi(n_scalingcol.substr(0, n_scalingcol.find('+', 0))) - 1;
                         i <= total_cols - 1; i++) {
                        col_indices.push_back(i);
                    }
                } else if (n_scalingcol.find(',') != string::npos) {
                    stringstream find_cols(n_scalingcol);
                    string curr_val;
                    while (getline(find_cols, curr_val, ',')) {
                        if (stoi(curr_val) > 0 and stoi(curr_val) <= total_cols) {
                            col_indices.push_back(stoi(curr_val) - 1);
                        }
                    }
                } else {
                    col_indices.push_back(stoi(n_scalingcol) - 1);
                }
                break;
            }
        }
    }
    unordered_map<int, string> colname_map = FileHeaderMapping(b_scalingfile, col_indices);

    // Stores the intergenic score of region-gene interactions by a vector of {gene_id, peak_id, a_col}.
    map <vector<int>, double> interaction_scores;
    bool use_abc_scoring = false;
    if ((abc_abc_output != "False") and (abc_abc_output != "FALSE") and (abc_abc_output != "false") and
        (abc_abc_output != "0") and (abc_abc_output.size() > 1)) {
        cout << "Get the gene-region mapping from the ABC-scoring file(s)" << endl;
        use_abc_scoring = true;
        decay = false;
        window_size = 2500;  // To later add the peaks in promoter range.
        // ____________________________________________________________
        // GET INTERACTIONS FROM ABC OUT - if available
        // ____________________________________________________________
        // Instead of getting the affinities of all regions in a window around the gene, we use the regions that were
        // called to be regulating the gene by the ABC-score. We also store the ABC-score of the interaction to scale
        // the affinities, instead of the distance decay.
        // Fetch all ABC-scoring files.
        vector <string> abc_files;
        if (col_indices.size() == 1 and ifstream(abc_abc_output).is_open()) {
            abc_files = {abc_abc_output};
        } else {
            size_t col_suffix_pos;
            bool numeric_naming = false;
            bool header_naming = false;
            smatch re_col_match;
            for (int a_col : col_indices) {  // First check which naming convention was used.
                regex re_col("_c[0-9]+.txt.gz"); // First try to find cN.
                regex_search(abc_abc_output, re_col_match, re_col);
                if (!re_col_match.empty()) {
                    numeric_naming = true;
                    col_suffix_pos = re_col_match.position(re_col_match.size() - 1);
                } else {  // Then try to match to the column header.
                    regex re_header(colname_map[a_col] + ".txt.gz");  // Find the position of the column suffix.
                    regex_search(abc_abc_output, re_col_match, re_header);
                    if (!re_col_match.empty()) {
                        col_suffix_pos = re_col_match.position(re_col_match.size() - 1);
                        header_naming = true;
                    }
                }
            }
            for (int a_col : col_indices) {
                if (!numeric_naming and !header_naming) { // If there's no _cN.txt ending we cut it at the ending.
                    col_suffix_pos = abc_abc_output.size() - 7;  // - .txt.gz
                    header_naming = true;
                }
                if (header_naming) {
                    abc_files.push_back(abc_abc_output.substr(0, col_suffix_pos) + colname_map[a_col] + ".txt.gz");
                } else if (numeric_naming) {
                    abc_files.push_back(
                            abc_abc_output.substr(0, col_suffix_pos) + "_c" + to_string(a_col + 1) + ".txt.gz");
                }
            }
        }

        for (int abc = 0; abc < abc_files.size(); abc++) {
            unordered_map<string, int> abc_header_map;
            string abc_cmd = string("zcat < ") + abc_files[abc] + " 2>&1";
            FILE *abc_stream = popen(abc_cmd.c_str(), "r");
            char abc_buffer[512];  // Max in the usual format is 148.
            row_count = -1;
            string abc_row;
            if (abc_stream) {
                while (!feof(abc_stream)) {
                    if (fgets(abc_buffer, 512, abc_stream) != nullptr) {
                        abc_row = abc_buffer;
                        if (row_count == -1) {
                            vector <string> columns = SplitTabLine(abc_row);
                            for (int c = 0; c < columns.size(); c++) {
                                string non_comm = columns[c];  // To remove any potential # from the header.
                                non_comm.erase(remove(non_comm.begin(), non_comm.end(), '#'), non_comm.end());
                                abc_header_map[non_comm] = c;
                            }
                            for (const string &needed : {"Ensembl ID", "PeakID", "intergenicScore"}) {
                                if (abc_header_map.find(needed) == abc_header_map.end()) {
                                    cout << "ERROR missing column in the interaction file: " << needed << endl;
                                    cout << abc_files[abc] << endl;
                                    return 1;
                                }
                            }
                        } else if (!abc_row.empty()) {
                            // All peaks in the ABC files must exist in the TRAP file.
                            vector <string> columns = SplitTabLine(abc_row);
                            string gene_string = columns[abc_header_map["Ensembl ID"]];
                            int this_gene_id = gene_name_map[gene_string];
                            if (tss_map.count(this_gene_id) == 1) {  // TSS-map is filtered for -u gene file.
                                string peak_string = columns[abc_header_map["PeakID"]];
                                int this_peak_id = peak_string_map[peak_string];
                                interaction_scores[{this_gene_id, this_peak_id, abc}] = stod(
                                        columns[abc_header_map["intergenicScore"]]);
                                gene_region_map[this_gene_id].insert(this_peak_id);
                            }
                        }
                        row_count++;
                    }
                }
                pclose(abc_stream);
            } else {
                cout << "ERROR Could not open ABC-scored file, command was the following:\n" << abc_cmd << endl;
                return 1;
            }
        }
    }

    // ____________________________________________________________
    // WRITE TEMPORARY GENE WINDOW FILE
    // ____________________________________________________________
    // Is not written directly on gtf-file read, as we need to find the most 5'-TSS first.
    string temp_window_file = o_output_prefix + "_temp_windows.bed";
    ofstream window_out(temp_window_file);
    Test_outfile(window_out, temp_window_file);
    for (auto const &gene : tss_map) {
        int window_start = stoi(gene.second[1]) - window_size - 1;  // gtf is 1-based, bedtools 0-based.
        if (window_start < 0) {
            window_start = 0;
        }
        string chr = gene.second[0];
        if (chr.substr(0, 3) == "chr") {  // Remove potential chr prefix, we removed it for the peaks as well.
            chr = chr.substr(3);
        }
        window_out << chr << "\t" << to_string(window_start) << "\t"
                   << to_string(stoi(gene.second[1]) + window_size) << "\t" << gene.first << "\n";
    }
    window_out.close();

    // ____________________________________________________________
    // WRITE BED FILE FOR PEAK REGIONS
    // ____________________________________________________________
    string temp_region_file = o_output_prefix + "_temp_regions.bed";
    ofstream region_out(temp_region_file);
    Test_outfile(region_out, temp_region_file);
    for (const string& region: peak_locations) {
        region_out << region << "\n";
    }
    region_out.close();

    // ____________________________________________________________
    // INTERSECT REGIONS WITH GENE WINDOWS
    // ____________________________________________________________
    // If abc-scoring was done, still do the intersections to include the promoters.
    string bed_intersect_cmd;
    if (n_scalingcol.length() > 0) {
        bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + temp_region_file + " -wo 2>&1";
    }
    else {  // With promoter-windows we need a full intersection to not have redundancy.
        bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + temp_region_file + " -f 1 -F 1 -wo 2>&1";
    }
    char intersect_buffer[256];  // Only has the gene window and coordinates of the candidate enhancer.
    FILE *bed_intersect_stream = popen(bed_intersect_cmd.c_str(), "r");
    if (bed_intersect_stream) {
        while (!feof(bed_intersect_stream)) {
            if (fgets(intersect_buffer, 256, bed_intersect_stream) != nullptr) {
                string line = intersect_buffer;
                // Process the bed_output, by mapping the regions to the respective genes.
                vector <string> columns = SplitTabLine(line);
                gene_region_map[stoi(columns[3])].insert(stoi(columns[7]));
            }
        }
    } else {
        cout << "ERROR Intersection of gene windows and peaks was not possible, command was the following:\n"
             << bed_intersect_cmd << endl;
        return 1;
    }
    // Remove the temporary files again.
    GetStdoutFromCommand("rm " + temp_region_file);
    GetStdoutFromCommand("rm " + temp_window_file);

    // ____________________________________________________________
    // FETCH POTENTIAL AFFINITY SCALINGS
    // ____________________________________________________________
    // If one or more scaling columns were specified we go through the peak file and fetch them as vector for each region.
    if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
        // Find the total number of columns first.
        cout << "Fetch affinity scalings out of the region file." << endl;
        int row_counter = 0;
        ifstream scaling_reader(b_scalingfile);
        while (!scaling_reader.eof()) {
            getline(scaling_reader, row);
            vector <string> columns = SplitTabLine(row);
            if (!row.empty() and row.substr(0, 1) != "#") {
                string chr = columns[0];
                if (chr.substr(0, 3) == "chr") {
                    chr = chr.substr(3);
                }
                vector<double> scaling_vector;
                scaling_vector.reserve(col_indices.size());
                for (int c : col_indices) {
                    scaling_vector.push_back(stod(columns[c]));
                }
                region_affinities[peak_string_map[chr + ":" + columns[1] + "-" + columns[2]]][1] = scaling_vector;
            }
            row_counter++;
        }
    }
    cout << "affinities fetched" << endl;  // TODO
    // ____________________________________________________________
    // MAP REGION AFFINITIES TO GENES
    // ____________________________________________________________
    vector <int> candidate_genes;  // To have the order fixed. The tss_map was already filtered for -u if defined.
     for (auto const &gene : tss_map) {
        if (gene_region_map.count(gene.first) == 0) {    // Check which genes didn't have any regions in their window.
            kicked_genes[gene.first] = "No open region in gene window";
        }
        else {
            candidate_genes.push_back(gene.first);
        }
    }

    // Write one output file for each activity column.
    int out_iters = 1;  // At least once, even if no activity column was specified.
    if (!col_indices.empty()) {  // If we have activity columns, write as many outputs.
        out_iters = col_indices.size();
    }
    // ____________________________________________________________
    // PREPARE BINARY OUTPUT IF SPECIFIED
    // ____________________________________________________________
    // Already open the streams and create files to check the endianness of the system.
    if (reshape) {
        string endianness = "UNKNOWN_ENDIAN";
        if (1 == HL_LITTLE_ENDIAN) {
            endianness = "LITTLE_ENDIAN";
        } else if (1 == HL_BIG_ENDIAN) {
            endianness = "BIG_ENDIAN";
        }
        // Write the Endian metadata string file with numbers to check whether binary conversion worked.
        ofstream endian_meta(o_output_prefix + "_Endianness_check.txt");
        auto endian_binary = fstream(o_output_prefix + "_Endianness_check.bin", ios::out | ios::binary);
        endian_meta << endianness << "\n";
        vector<float> test_floats = {0.25, -0.75, 128.5, -8.5};
        for (float f : test_floats) {
            endian_meta << f << "\n";
            endian_binary.write((char *) &f, sizeof(f));
        }
        endian_meta.close();
        endian_binary.close();
        ofstream reshape_meta(o_output_prefix + "_Reshape_Meta.txt");  // Write it already, so we can close it.
        // First line is the total number of entries that will be written into the binary file.
        reshape_meta << out_iters * candidate_genes.size() * (num_tfs + 3) << "\n";  // +3 for Num- Distance and Size.
        // We can already fill the order of TFs, cells and genes, this will stay constant across all activity columns.
        reshape_meta
                << tf_header.substr(1, tf_header.size() - 1) + "\t" + "NumPeaks" + "\t" + "AvgPeakDistance" + "\t" +
                   "AvgPeakSize" + "\n";
        for (int c = 0; c < out_iters; ++c) {
            string column_suffix;
            if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
                column_suffix = colname_map[col_indices[c]];
                column_suffix = column_suffix.substr(1, column_suffix.size() - 1);  // To remove the leading _.
            }
            if (c < out_iters - 1) {
                reshape_meta << column_suffix << "\t";
            } else {
                reshape_meta << column_suffix;
            }
        }
        reshape_meta << "\n";
        for (int g = 0; g < candidate_genes.size() - 1; g++) {
            reshape_meta << gene_id_map[candidate_genes[g]] << "\t";
        }
        // Add the last one manually to avoid trailing tabs.
        reshape_meta << gene_id_map[candidate_genes[candidate_genes.size()-1]];
        reshape_meta.close();
    }

    // ____________________________________________________________
    // PREPARE OUTPUT FILES
    // ____________________________________________________________
    vector <string> out_files;
    fstream binary_stream;
    if (not reshape) {
        for (int c = 0; c < out_iters; ++c) {
            string column_suffix;
            if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
                column_suffix = colname_map[col_indices[c]];
            }
            // Open the output and write the header already.
            string out_path = o_output_prefix + "_TF_Gene_Affinities" + column_suffix + ".txt";
            ofstream gene_tf_out(out_path);
            Test_outfile(gene_tf_out, out_path);
            gene_tf_out << "geneID" + tf_header + "\t" + "NumPeaks" + "\t" + "AvgPeakDistance" + "\t" + "AvgPeakSize" +
                           "\n";
            gene_tf_out.close();
            out_files.push_back(out_path);
        }
    } else {
        binary_stream = fstream(o_output_prefix + "_Reshape_Binary.bin", ios::out | ios::binary);
    }

    // ____________________________________________________________
    // GET AND WRITE TF-GENE AFFINITIES
    // ____________________________________________________________
    auto start_out = chrono::high_resolution_clock::now();
    cout << "Summarize the region affinities on gene level." << endl;

    int batch_size = 1000;  // Collect data for a batch of genes and then write it to output in one go.
    int gene_out_counter = 0;
    map < int, vector < vector < double>>> gene_batch;
    for (const int& gene : candidate_genes) {
        gene_out_counter++;
        int gene_tss = stoi(tss_map[gene][1]);
        // Store affinities per out_iter.
        vector<vector<double>> these_gene_affinities(out_iters, vector<double>(num_tfs + 3));  // + Peak features.
        vector<int> distances(out_iters);
        vector<int> sizes(out_iters);
        vector<int> mapped_peaks(out_iters);

        for (const int &region: gene_region_map[gene]) {  // We already excluded genes without regions.
            auto region_iter = region_affinities.find(region);  // Should theoretically always be there.
            if (region_iter != region_affinities.end()) {
                const array<vector<double>, 2> &region_info = region_iter->second;  // Affs | Scalings
                if (!region_info[0].empty()) {  // Should only happen if someone ran ABC independently before.
                    string region_string = peak_id_map[region];
                    int colon = region_string.find(':', 0);
                    int hyphon = region_string.find('-', 0);
                    int start = stoi(region_string.substr(colon + 1, hyphon - colon - 1));
                    int end = stoi(region_string.substr(hyphon + 1));
                    int distance = GetDistance(start, end, gene_tss);
                    int peak_size = abs(end - start);

                    vector<double> affinity_scalers(out_iters, 1);  // Check the conditions and fill once.
                    if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
                        affinity_scalers = region_info[1];
                    }
                    // Decay is set to false if use_abs_scoring.
                    if (decay or (use_abc_scoring and distance <= 2500)) {
                        for (double &affinity_scaler : affinity_scalers) {
                            affinity_scaler *= exp(-(distance / 5000.0));
                        }
                    } else if (use_abc_scoring) {  // Scale with the intergenic score column.
                        for (int n = 0; n < affinity_scalers.size(); n++) {
                            auto inter_iter = interaction_scores.find({gene, region, n});
                            if (inter_iter != interaction_scores.end()) {
                                affinity_scalers[n] = inter_iter->second;
                            } else {  // If there was no score for that interaction, we don't consider that region.
                                affinity_scalers[n] = 0;
                            }
                        }
                    }
#pragma omp parallel for num_threads(cores)
                    for (int c = 0; c < out_iters; ++c) {
                        // If an interaction was not present in one of the abc-files, the interaction score is 0.
                        if (affinity_scalers[c] > 0) {
                            distances[c] += distance;
                            sizes[c] += peak_size;
                            mapped_peaks[c] += 1;
                        }
                        // Add the whole normalized affinity vector to the already existing one.
                        for (int i = 0; i < num_tfs; i++) {
                            these_gene_affinities[c][i] +=
                                    (region_info[0][i] / motif_lengths[i]) * affinity_scalers[c];
                        }
                    }
                }
            }
        }
        // After we're done with all regions, add the three additional peak features.
        for (int c = 0; c < out_iters; ++c) {
            float num_peaks = mapped_peaks[c];
            float avg_distance = 0;
            float avg_size = 0;
            if (num_peaks > 0) {
                avg_distance = distances[c] / num_peaks;
                avg_size = sizes[c] / num_peaks;
            }
            these_gene_affinities[c][num_tfs] = num_peaks;
            these_gene_affinities[c][num_tfs + 1] = avg_distance;
            these_gene_affinities[c][num_tfs + 2] = avg_size;
        }

        if (not reshape) {
            gene_batch[gene] = these_gene_affinities;
        }

        // ____________________________________________________________
        // WRITE TF-GENE AFFINITY OUTPUT
        // ____________________________________________________________
        // Write dense matrix output as strings when a batch is full or we reached the last gene.
        if (not reshape and (gene_batch.size() == batch_size or gene_out_counter == candidate_genes.size())) {
#pragma omp parallel for num_threads(cores)
            for (int c = 0; c < out_iters; ++c) {
                ofstream out_stream;
#pragma omp critical(out_write)
                {
                    out_stream.open(out_files[c].c_str(), ios_base::app);
                }
                out_stream.precision(10);
                for (auto const &curr_gene : gene_batch) {
                    out_stream << gene_id_map[curr_gene.first];
                    for (int i = 0; i < num_tfs + 3; i++) {
                        out_stream << "\t" << fixed << curr_gene.second[c][i];  // Fixed maintains the set precision.
                    }
                    out_stream << "\n";
                }
                out_stream.close();
            }
            gene_batch.clear();  // Reset to fill the next batch.
        }
        if (reshape) {  // Else write binary file, but for each gene, not batched.
            for (int c = 0; c < out_iters; ++c) {
                for (float i: these_gene_affinities[c]) {
                    binary_stream.write((char *) &i, sizeof(float));
                }
            }
        }
    }
    auto stop_out = chrono::high_resolution_clock::now();
    auto duration_out = chrono::duration_cast<chrono::milliseconds>(stop_out - start_out);
    cout << duration_out.count() << "ms; " << " Scores processed" << endl;

    if (reshape) {
        binary_stream.close();
    }

    auto start_comp = chrono::high_resolution_clock::now();
    cout << "Compressing output" << endl;
    if (not reshape) {
        #pragma omp parallel for num_threads(cores)
        for (int c = 0; c < out_iters; c++) {
            string column_suffix;
            if (n_scalingcol != "0" and n_scalingcol.length() > 0) {
                column_suffix = colname_map[col_indices[c]];
	        }
            GzipFile(o_output_prefix + "_TF_Gene_Affinities" + column_suffix + ".txt");
        }
    }
    else {
        binary_stream.close();  // Shouldn't be necessary but might cause issues with gzip.
        GzipFile(o_output_prefix + "_Reshape_Binary.bin");
    }
    auto stop_comp = chrono::high_resolution_clock::now();
    auto duration_comp = chrono::duration_cast<chrono::milliseconds>(stop_comp - start_comp);
    cout << duration_comp.count() << "ms; " << " Compression" << endl;

//     ____________________________________________________________
//     WRITE DISCARDED GENES FILE
//     ____________________________________________________________
    ofstream discarded_out(o_output_prefix + "_discarded_Genes.txt");
    Test_outfile(discarded_out, o_output_prefix + "_discarded_Genes.txt");
    for (auto const &gene : kicked_genes) {
        discarded_out << gene_id_map[gene.first] << "\t" << gene.second << "\n";
    }
    discarded_out.close();
}


