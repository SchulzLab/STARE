//
// Created by Dennis Hecker on 28.07.21.
//

# include <iostream>
# include <string>
# include <vector>
# include <array>
# include <numeric>
# include <fstream>
# include <math.h>
# include <sstream>
# include <algorithm>
# include <unordered_map>
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
 * MacOS:
 * g++ TF_Gene_Scorer.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -O3 -std=c++11 -o TF_Gene_Scorer
 * Linux:
 * g++ TF_Gene_Scorer.cpp STARE_MiscFunctions.cpp -fopenmp -O3 -std=c++11 -o TF_Gene_Scorer
 *
 * ./TF_Gene_Scorer -a ../Test/example_annotation.gtf -i region-TF-affinity-file -o output-dir -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM -w 50000 -e TRUE
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 */


int main(int argc, char **argv) {
    using namespace std;

    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-a gene annotation file"
                            "\n-b bed file with open chromatin regions"
                            "\n-n column(s) with the activity of regions in the -b file"
                            "\n-i TF-region affinity file from TRAP"
                            "\n-o prefix of the output files"
                            "\n-p PSEM file"
                            "\n-w window size around TSS (default 50kb)"
                            "\n-e decay (default true)"
                            "\n-c number of cores"
                            "\n-abc path to the ABC-scoring output, if available, otherwise set to 0\n";

    vector <string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
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

    // Reading the input args:
    string a_gene_annotation, b_scalingfile, n_scalingcol, i_tf_region_affinity, o_output_prefix, p_psem_file,
            w_window_size, e_decay, c_cores, abc_abc_output;
    static struct option long_options[] =
            {  // optional_argument throws a Segmentation fault, for whatever reason.
                    {"a",   required_argument, NULL, 'a'},
                    {"b",   required_argument, NULL, 'b'},
                    {"n",   required_argument, NULL, 'n'},
                    {"i",   required_argument, NULL, 'i'},
                    {"o",   required_argument, NULL, 'o'},
                    {"p",   required_argument, NULL, 'p'},
                    {"w",   required_argument, NULL, 'w'},
                    {"e",   required_argument, NULL, 'e'},
                    {"c",   required_argument, NULL, 'c'},
                    {"abc", required_argument, NULL, 't'},
                    {NULL,  0,                 NULL, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "a:b:n:i:o:p:w:e:c:t:", long_options, NULL)) != -1) {
        switch (arg) {
            case 'a':
                a_gene_annotation = optarg;
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
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for the required parameters.
    for (string p: {a_gene_annotation, i_tf_region_affinity, o_output_prefix, p_psem_file}) {
        if (p.size() < 1) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    // If not given as input, set to default.
    int window_size = stoi(SetOptionalInput(w_window_size, "50000")) / 2;
    bool decay = SetBoolInput(e_decay, true);
    int cores = stoi(SetOptionalInput(c_cores, "1"));

    // ____________________________________________________________
    // READ GENE ANNOTATION
    // ____________________________________________________________
    unordered_map <string, vector<string>> tss_map;  // GeneID: [chr, 5'TSS]
    string row;
    ifstream Read_Gene_Annotation(a_gene_annotation);
    if (!Read_Gene_Annotation) {
        cout << "ERROR Could not open the gene annotation file\n" + a_gene_annotation << endl;
        return 1;
    }
    cout << "reading gene annotation" << endl;
    int row_count = -1;
    while (!Read_Gene_Annotation.eof()) {
        getline(Read_Gene_Annotation, row);
        if ((row.size() > 1) and (row.rfind("#", 0) != 0)) {
            vector <string> columns = SplitTabLine(row);
            if (columns[2] == "gene") {
                string id_delimiter = "gene_id ";
                int id_start = columns[8].find(id_delimiter, 0) + id_delimiter.size();
                int id_end = columns[8].find(";", 0);
                string gene_id = columns[8].substr(id_start + 1,
                                                   id_end - id_start - 2);  // Removing the quotation marks.
                string strand = columns[6];
                string gene_start;
                if (strand == "+") {
                    gene_start = columns[3];
                } else {
                    gene_start = columns[4];
                }
                if (tss_map.find(gene_id) != tss_map.end()) {  // Replace if the new one is more in 5'.
                    if (((stoi(tss_map[gene_id][1]) > stoi(gene_start)) and (strand == "+")) or
                        ((stoi(tss_map[gene_id][1]) < stoi(gene_start)) and (strand == "-"))) {
                        tss_map[gene_id][1] = gene_start;
                    }
                } else {
                    tss_map[gene_id] = vector < string > {columns[0], gene_start};
                }
            }
        }
        row_count++;
    }
    Read_Gene_Annotation.close();

    // ____________________________________________________________
    // GET MOTIF LENGTHS
    // ____________________________________________________________
    unordered_map<string, int> motif_length_map;
    int length_helper = 0;
    string current_tf = "";
    ifstream Read_PSEM(p_psem_file);
    if (!Read_PSEM) {
        cout << "ERROR Could not open the PSEM file\n" + p_psem_file << endl;
        return 1;
    }
    while (!Read_PSEM.eof()) {
        getline(Read_PSEM, row);
        if ((row.rfind(">", 0) == 0)) {
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
    map <string, set<string>> gene_region_map; // Storing the names of peaks in its window.
    vector <string> peak_locations;  // For the temporary bed file.
    unordered_map <string, string> kicked_genes;  // Track the genes that don't have any regions in their window, non-ABC.

    // Fetch the region-tf affinities and do the intersection with the gene windows if no ABC-scoring file was provided.
    // Already allocate two vectors in that datastructure to later store potential region scalings.
    cout << "Read TRAP's affinity files" << endl;
    unordered_map < string, array < vector < double >, 2 >> region_affinities;
    Read_Affinity.open(i_tf_region_affinity);
    row_count = -1;
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
            for (int i = 1; i < tf_names.size(); i++) {  // No need to
                motif_lengths.push_back(motif_length_map[tf_names[i]]);
            }
        }
        if ((row.size() > 1) and (row_count != -1)) {
            int first_tab = row.find("\t", 0);
            string region_string = row.substr(0, first_tab);
            if (region_string.substr(0, 3) == "chr") {
                region_string = region_string.substr(3);
            }
            int colon_pos = region_string.find(":", 0);
            int hyphon_pos = region_string.find("-", 0);
            int region_start = stoi(region_string.substr(colon_pos + 1, hyphon_pos - colon_pos - 1));
            int region_end = stoi(region_string.substr(hyphon_pos + 1));

            peak_locations.push_back(
                    region_string.substr(0, colon_pos) + "\t" + to_string(region_start) + "\t" + to_string(region_end));
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
            region_affinities[region_string][0] = affinity_vals;
        }
        row_count++;
    }
    cout << "TRAP file read" << endl;
    // ____________________________________________________________
    // GET SCALING COLS - if specified
    // ____________________________________________________________
    vector<int> col_indices;
    if (n_scalingcol != "0") {
        // Look into the peakfile to see with how many columns/files we are dealing with.
        ifstream peek_peak(b_scalingfile);
        string line_peek;
        while (!peek_peak.eof()) {
            getline(peek_peak, line_peek);
            if (line_peek.substr(0, 1) != "#") {
                int total_cols = SplitTabLine(line_peek).size();
                if (n_scalingcol.find('-') != string::npos) {
                    int hyphon_pos = n_scalingcol.find("-", 0);
                    for (int i=stoi(n_scalingcol.substr(0, hyphon_pos)) - 1; i <= stoi(n_scalingcol.substr(hyphon_pos + 1)) - 1; i++) {
                        col_indices.push_back(i);
                    }
                } else if (n_scalingcol.find('+') != string::npos) {
                    for (int i=stoi(n_scalingcol.substr(0, n_scalingcol.find('+', 0))) - 1; i <= total_cols-1; i++) {
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

    // Stores the intergenic score of region-gene interactions.
    unordered_map <string, unordered_map<int, double>> interaction_scores;
    bool use_abc_scoring = false;
    if ((abc_abc_output != "False") and (abc_abc_output != "FALSE") and (abc_abc_output != "false") and
        (abc_abc_output != "0") and (abc_abc_output.size() > 1)) {
        cout << "Get the gene-region mapping from the ABC-scoring file(s)" << endl;
        use_abc_scoring = true;
        window_size = 2500;  // To later add the peaks in promoter range.
        // ____________________________________________________________
        // GET INTERACTIONS FROM ABC OUT - if available
        // ____________________________________________________________
        // Instead of getting the affinities of all regions in a window around the gene, we use the regions that were
        // called to be regulating the gene by the ABC-score. We also store the ABC-score of the interaction to scale
        // the affinities, instead of the distance decay.
        // Fetch all ABC-scoring files.
        vector<string> abc_files;
        if (col_indices.size() == 1 and ifstream(abc_abc_output).is_open()) {
            abc_files = {abc_abc_output};
        }
        else {
            int col_suffix_pos;
            smatch re_col_match;
            for (int a_col : col_indices) {
                regex re_col("_c[0-9]+.txt.gz"); // First try to find cN.
                regex_search(abc_abc_output, re_col_match, re_col);
                if (re_col_match.size() > 0) {
                    col_suffix_pos = re_col_match.position(re_col_match.size() - 1);
                    abc_files.push_back(abc_abc_output.substr(0, col_suffix_pos) + "_c" + to_string(a_col + 1) + ".txt.gz");
                }
                else {  // Then try to match to the column header.
                    regex re_header(colname_map[a_col] + ".txt.gz");  // Find the position of the column suffix.
                    regex_search(abc_abc_output, re_col_match, re_header);
                    if (re_col_match.size() > 0) {
                        col_suffix_pos = re_col_match.position(re_col_match.size() - 1);
                    } else {  // If there's no _cN.txt ending we cut it at the ending.
                        col_suffix_pos = abc_abc_output.size() - 7;  // - .txt.gz
                    }
                    abc_files.push_back(abc_abc_output.substr(0, col_suffix_pos) + colname_map[a_col] + ".txt.gz");
                }
            }
        }

        for (int abc = 0; abc < abc_files.size(); abc++) {
            unordered_map<string, int> abc_header_map;
            string abc_cmd = string("zcat < ") + abc_files[abc] + " 2>&1";
            FILE* abc_stream = popen(abc_cmd.c_str(), "r");
            char abc_buffer[512];  // Max in the usual format is 148.
            row_count = -1;
            if (abc_stream) {
                while (!feof(abc_stream)) {
                    if (fgets(abc_buffer, 512, abc_stream) != NULL) {
                        string row = abc_buffer;
                        if (row_count == -1) {
                            vector <string> columns = SplitTabLine(row);
                            for (int c = 0; c < columns.size(); c++) {
                                string non_comm = columns[c];  // To remove any potential # from the header.
                                non_comm.erase(remove(non_comm.begin(), non_comm.end(), '#'), non_comm.end());
                                abc_header_map[non_comm] = c;

                            }
                        } else if (row.size() > 0) {
                            vector <string> columns = SplitTabLine(row);
                            string gene_id = columns[abc_header_map["Ensembl ID"]];
                            string peak_id = columns[abc_header_map["PeakID"]];
                            gene_region_map[gene_id].insert(peak_id);
                            interaction_scores[gene_id + "-" + peak_id][abc] = stod(columns[abc_header_map["intergenicScore"]]);
                        }
                        row_count++;
                    }
                }
                pclose(abc_stream);
            }
            else {
                cout << "ERROR Could not open ABC-scored file, command was the following:\n" << abc_cmd << endl;
                return 1;
            }
        }
    }

    // ____________________________________________________________
    // WRITE TEMPORARY GENE WINDOW FILE
    // ____________________________________________________________
    // Is not written directly on gtf-file read, as we need to find the most 5'-TSS first.
    string temp_window_file = a_gene_annotation.substr(0, a_gene_annotation.size() - 4) + "_temp_windows.bed";
    ofstream window_out(temp_window_file);
    Test_outfile(window_out, temp_window_file);
    for (auto const &gene : tss_map) {
        int window_start = stoi(gene.second[1]) - window_size;
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
    string temp_region_file = i_tf_region_affinity.substr(0, i_tf_region_affinity.size() - 4) + "_temp_regions.bed";
    ofstream region_out(temp_region_file);
    Test_outfile(region_out, temp_region_file);
    for (string region: peak_locations) {
        region_out << region << "\n";
    }
    region_out.close();

    // ____________________________________________________________
    // INTERSECT REGIONS WITH GENE WINDOWS
    // ____________________________________________________________
    // If abc-scoring was done, still do the intersections to include the promoters.
    string bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + temp_region_file + " -wo 2>&1";
    char intersect_buffer[256];  // Only has the gene window and coordinates of the candidate enhancer.
    FILE* bed_intersect_stream = popen(bed_intersect_cmd.c_str(), "r");
    if (bed_intersect_stream) {
        while (!feof(bed_intersect_stream)) {
            if (fgets(intersect_buffer, 256, bed_intersect_stream) != NULL) {
                string line = intersect_buffer;
                // Process the bed_output, by mapping the regions to the respective genes.
                vector <string> columns = SplitTabLine(line);
                string region_name = columns[4] + ":" + columns[5] + "-" + columns[6];
                gene_region_map[columns[3]].insert(region_name);
            }
        }
    }
    else {
        cout << "ERROR Intersection of gene windows and peaks was not possible, command was the following:\n"
        << bed_intersect_cmd << endl;
        return 1;
    }
    // Remove the temporary files again.
    GetStdoutFromCommand("rm " + temp_region_file);
    GetStdoutFromCommand("rm " + temp_window_file);

    // Check which genes didn't have any regions in their window.
    for (auto const &gene : tss_map) {
        if (gene_region_map.find(gene.first) == gene_region_map.end()) {
            kicked_genes[gene.first] = "No open region in gene window";
        }
    }

    // ____________________________________________________________
    // FETCH POTENTIAL AFFINITY SCALINGS
    // ____________________________________________________________
    // If one or more scaling columns were specified we go through the peak file and fetch them as vector for each region.
    if (n_scalingcol != "0") {
        // Find the total number of columns first.
        cout << "Fetch affinity scalings out of the region file." << endl;
        int row_counter = 0;
        ifstream scaling_reader(b_scalingfile);
        while (!scaling_reader.eof()) {
            getline(scaling_reader, row);
            vector<string> columns = SplitTabLine(row);
            if (row.size() > 0 and row.substr(0, 1) != "#") {
                string chr = columns[0];
                if (chr.substr(0, 3) == "chr") {
                    chr = chr.substr(3);
                }
                vector<double> scaling_vector;
                for (int c : col_indices) {
                    scaling_vector.push_back(stod(columns[c]));
                }
                region_affinities[chr + ":" + columns[1] + "-" + columns[2]][1] = scaling_vector;
            }
            row_counter++;
        }
    }

    // ____________________________________________________________
    // MAP REGION AFFINITIES TO GENES
    // ____________________________________________________________
    vector<string> candidate_genes;  // Needed so that openmp can iterate through them.
    for (auto const& gene: gene_region_map) {
        candidate_genes.push_back(gene.first);
    }
    // Write one output file for each activity column.
    cout << "Summarize the region affinities on gene level." << endl;
    int out_iters = 1;  // At least once, even if no activity column was specified.
    if (col_indices.size() > 0) {  // If we have activity columns, write as many outputs.
        out_iters = col_indices.size();
    }
    for (int c = 0; c < out_iters; ++c) {
        string column_suffix = "";
        if (n_scalingcol != "0") {
            column_suffix = colname_map[col_indices[c]];
        }
        cout << "Starting activity column " << column_suffix << endl;
        // Allocate as vector so that each gene's affinities is at the same index as it's ID in candidate_genes.
        vector <vector<double>> gene_tf_affinities(candidate_genes.size(), vector<double>(num_tfs));
#pragma omp parallel for num_threads(cores)
        for (int g = 0; g < candidate_genes.size(); g++) {  // Runs once if zero or one column was given.
            string gene = candidate_genes[g];
            vector<double> these_gene_affinities(num_tfs + 2, 0.0);
            int distance_helper = 0;
            int gene_tss = stoi(tss_map[gene][1]);

            for (string region: gene_region_map[gene]) {
                double affinity_scaler = 1;  // Activity * Decay factor / ABC-score, or 1.
                array<vector<double>, 2> region_info = region_affinities[region];
                vector<double> affinity_vals = region_info[0];

                if (affinity_vals.size() > 0) {  // Should only happen if someone ran ABC-scoring independently before.
                    int colon = region.find(":", 0);
                    int hyphon = region.find("-", 0);
                    int start = stoi(region.substr(colon + 1, hyphon - colon - 1));
                    int end = stoi(region.substr(hyphon + 1));
                    double distance = min(abs(gene_tss - start), abs(gene_tss - end));
                    distance_helper += distance;
                    if (n_scalingcol != "0" and !use_abc_scoring) {
                        affinity_scaler = region_info[1][c];
                    }
                    if (decay) {
                        affinity_scaler *= exp(-(distance / 5000.0));
                    }
                    if (use_abc_scoring and distance <= 2500) {
                        affinity_scaler = region_info[1][c] * exp(-(distance / 5000.0));
                    }
                    else if (use_abc_scoring) {  // Scale with the intergenic score column.
                        affinity_scaler = interaction_scores[gene + "-" + region][c];  // 0 if not existent.
                    }
                    // Add the whole normalized affinity vector to the already existing one.
                    for (int i = 0; i < num_tfs; i++) {
                        these_gene_affinities[i] += (affinity_vals[i] / motif_lengths[i]) * affinity_scaler;
                    }
                }
            }
            these_gene_affinities[num_tfs] = gene_region_map[gene].size();
            these_gene_affinities[num_tfs + 1] = distance_helper / gene_region_map[gene].size();
            gene_tf_affinities[g] = these_gene_affinities;
        }

        // ____________________________________________________________
        // WRITE TF-GENE AFFINITY OUTPUT
        // ____________________________________________________________
        cout << "writing and compressing activity column" << endl;
        ofstream gene_tf_out(o_output_prefix + "_TF_Gene_Affinities" + column_suffix + ".txt");
        Test_outfile(gene_tf_out, o_output_prefix + "_TF_Gene_Affinities" + column_suffix + ".txt");
        gene_tf_out << "geneID" + tf_header + "\t" + "NumPeaks" + "\t" + "AvgPeakSize" + "\n";
        for (int g = 0; g < candidate_genes.size(); g++) {
            string gene = candidate_genes[g];
            vector<double> these_affinities = gene_tf_affinities[g];
            gene_tf_out << gene;
            for (int i = 0; i < num_tfs; i++) {
                gene_tf_out << "\t" << To_string_with_precision(these_affinities[i]);
            }
            gene_tf_out << "\t" << to_string(static_cast<double>(these_affinities[num_tfs]));  // NumPeaks
            gene_tf_out << "\t" << to_string(these_affinities[num_tfs + 1]);  // AvgPeakSize
            gene_tf_out << "\n";
        }
        gene_tf_out.close();  // Otherwise the return of the command would be written to output.
        GzipFile(o_output_prefix + "_TF_Gene_Affinities" + column_suffix + ".txt");
    }
    // ____________________________________________________________
    // WRITE DISCARDED GENES FILE
    // ____________________________________________________________
    ofstream discarded_out(o_output_prefix + "_discarded_Genes.txt");
    Test_outfile(discarded_out, o_output_prefix + "_discarded_Genes.txt");
    for (auto const &gene : kicked_genes) {
        discarded_out << gene.first << "\t" << gene.second << "\n";
    }
    discarded_out.close();
}

