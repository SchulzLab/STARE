//
// Created by Dennis Hecker on 10.12.21.
//

# include <iostream>
# include <iomanip>
# include <string>
# include <vector>
# include <math.h>
# include <limits>
# include <chrono>
# include <regex>
# include <algorithm>
# include <random>
# include <unordered_map>
# include <unordered_set>
# include <map>
# include <set>
# include <getopt.h>
# include <stdio.h>
# include <stdlib.h>
# include <omp.h>
# include <boost/numeric/ublas/matrix_sparse.hpp>
# include "STARE_MiscFunctions.h"

/*
 * When having chromatin contact data at hand, it is possible to get the regulatory interactions based on the
 * ABC-score. Either the regular version (-q False) or the adapted are possible (-q True).
 * Together with STARE for deriving affinities of TFs to genes, we can use these called interactions instead of all
 * the open regions in a window. This program also works independently of STARE.
 *
 * Compilation for MacOS:
 * g++ STARE_ABCpp.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o STARE_ABCpp
 * Compilation for or Linux:
 * g++ STARE_ABCpp.cpp STARE_MiscFunctions.cpp -fopenmp -std=c++11 -O3 -o STARE_ABCpp
 *
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 * C++ implementation of the ABC-scoring principle from Fulco et al.: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
 */


class gene_information {
public:
    std::string gene_id;
    std::string name;
    std::string chr;
    std::vector<int> tss;
};

class peak_information {
public:
    int start;
    int end;
    std::vector<double> signal;
    double contact_sum = 0;
};

class abc_hit {
public:
    int gene_id;
    int peak_id;
    double scaledActivity;
    double contact;
    double distance;
    double score;
};



int main(int argc, char **argv) {
    using namespace std;
    auto start0 = chrono::high_resolution_clock::now();
    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-b enhancer/peak-file"
                            "\n-n activity-column(s), start counting at 1"
                            "\n-a gtf gene annotation"
                            "\n-u file with geneIDs/symbols to limit the output to"
                            "\n-w genewindow-size"
                            "\n-f folder with the normalized hic-files per chromosome"
                            "\n-k binsize of the hic-files"
                            "\n-t cut-off for the ABC-score (default 0.02), set to 0 to get all scored interactions"
                            "\n-o prefix with which the files are written"
                            "\n-c number of cores used (default 1)"
                            "\n-x file with regions to be excluded"
                            "\n-i all_tss to average across TSS or 5_tss to use only the 5' TSS (default all_tss)"
                            "\n-d whether to use pseudocount for contact frequency (default True)"
                            "\n-q whether to adapt the activity for an enhancer's contacts (default True)"
                            "\n-m enhancer window size in which to consider contacts for adjustment";

    vector<string> h_flags = {"h", "-h", "--h", "help", "-help", "--help"};
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
    string b_peakfile, a_promoterfile, u_genefile, n_activitycol, w_genewindowsize, k_binsize, t_abc_cutoff, f_contactfolder,
            c_cores, o_prefix, x_exclude_regions, i_tss_type, d_pseudocount, q_adjusted_abc, m_enhwindowsize;
    static struct option long_options[] =
            {
                    {"b", required_argument, nullptr, 'b'},
                    {"a", required_argument, nullptr, 'a'},
                    {"u", required_argument, nullptr, 'u'},
                    {"n", required_argument, nullptr, 'n'},
                    {"w", required_argument, nullptr, 'w'},
                    {"k", required_argument, nullptr, 'k'},
                    {"t", required_argument, nullptr, 't'},
                    {"f", required_argument, nullptr, 'f'},
                    {"c", required_argument, nullptr, 'c'},
                    {"o", required_argument, nullptr, 'o'},
                    {"x", required_argument, nullptr, 'x'},
                    {"i", required_argument, nullptr, 'i'},
                    {"d", required_argument, nullptr, 'd'},
                    {"q", required_argument, nullptr, 'q'},
                    {"m", required_argument, nullptr, 'm'},
                    {nullptr, 0,             nullptr, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "b:a:u:n:w:k:t:c:f:o:x:i:d:h:q:m:", long_options, nullptr)) != -1) {
        switch (arg) {
            case 'b':
                b_peakfile = optarg;
                break;
            case 'a':
                a_promoterfile = optarg;
                break;
            case 'u':
                u_genefile = optarg;
                break;
            case 'n':
                n_activitycol = optarg;
                break;
            case 'w':
                w_genewindowsize = optarg;
                break;
            case 'k':
                k_binsize = optarg;
                break;
            case 't':
                t_abc_cutoff = optarg;
                break;
            case 'f':
                f_contactfolder = optarg;
                break;
            case 'c':
                c_cores = optarg;
                break;
            case 'o':
                o_prefix = optarg;
                break;
            case 'x':
                x_exclude_regions = optarg;
                break;
            case 'd':
                d_pseudocount = optarg;
                break;
            case 'i':
                i_tss_type = optarg;
                break;
            case 'q':
                q_adjusted_abc = optarg;
                break;
            case 'm':
                m_enhwindowsize = optarg;
                break;
            case '?':
                break;
            default:
                printf("Unknown flag\n");
        }
    }

    // Check if there is an argument for each parameter.
    for (const string &p: {b_peakfile, a_promoterfile, o_prefix}) {
        if (p.empty()) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    // If not given as input, set to default.
    // stod() works with scientific notation.
    int gene_windowsize = static_cast<int>(stod(SetOptionalInput(w_genewindowsize, "5000000"))) / 2;
    int bin_size = static_cast<int>(stod(SetOptionalInput(k_binsize, "5000")));
    // If contactfolder is False use the fractal function.
    string f_contactfolder_lower = f_contactfolder;
    transform(f_contactfolder_lower.begin(), f_contactfolder_lower.end(), f_contactfolder_lower.begin(),
                                       [](unsigned char c){ return tolower(c); });
    set<string> false_options = {"false", "0", "f"};
    if (f_contactfolder == " " or false_options.find(f_contactfolder_lower) != false_options.end()) {
        f_contactfolder = "";
    }
    string activity_col = SetOptionalInput(n_activitycol, "4");
    int cores = stoi(SetOptionalInput(c_cores, "1"));
    double abc_cutoff = stof(SetOptionalInput(t_abc_cutoff, "0.02"));
    string tss_type = SetOptionalInput(i_tss_type, "all_tss");
    bool do_pseudocount = SetBoolInput(d_pseudocount, true);
    int enh_windowsize = static_cast<int>(stod(SetOptionalInput(m_enhwindowsize, "5000000"))) / 2;
    bool do_adjusted_abc = SetBoolInput(q_adjusted_abc, true);

    string activity_header = "adaptedActivity";
    int intergenic_activity_col = 3;  // Vector-index where the contact-adjusted activity is stored.
    if (not do_adjusted_abc) {
        enh_windowsize = gene_windowsize;
        activity_header = "scaledActivity";
        intergenic_activity_col = 2;  // Vector-index for the unmodified region activity.
    }

    string tss_identifier;
    if (tss_type == "all_tss") {
        tss_identifier = "transcript";
    }
    else if (tss_type == "5_tss") {
        tss_identifier = "gene";
    }
    else {
        cout << "ERROR unknown option for -i  " << tss_identifier << "\noptions are all_tss or 5_tss" << endl;
        return 1;
    }

    // ____________________________________________________________
    // PROCESS GTF GENE ANNOTATION
    // ____________________________________________________________
    cout << "Reading annotation and intersecting peaks" << endl;
    map<string, unordered_set<int>> chr_gene_map;  // For each chromosome stores its genes.
    unordered_map<string, int> id_ensg_map;  // To be able to check if already saw a gene in the annotation.
    int gene_counter = 0;
    unordered_map<int, gene_information> promoter_map;  // Map of class to store id, name, chr and set of TSS.
    id_ensg_map.reserve(60000);
    promoter_map.reserve(60000);
    unordered_map<string, vector<int>> hic_boundaries;  // Store the min/max boundaries for each chromosome.
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
                vector<string> columns = SplitTabLine(row);
                if (columns[2] == tss_identifier) {
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
                        if (tss_identifier == "gene") {
                            // Update if the new one is more in 5', must be on the same chromosome if non-unique.
                            // If the tss_identifier is gene, we change the set in the gene_information class.
                            if (((promoter_map[counter_iter->second].tss[0] > gene_start) and (strand == "+") and
                                 (promoter_map[counter_iter->second].chr == chr)) or
                                ((promoter_map[counter_iter->second].tss[0] < gene_start) and (strand == "-") and
                                 (promoter_map[counter_iter->second].chr == chr))) {
                                promoter_map[counter_iter->second].tss[0] = {gene_start};
                            }
                        } else if (tss_identifier == "transcript") {  // Collect the unique start of each transcript.
                            if (promoter_map[counter_iter->second].chr == chr) {
                                promoter_map[counter_iter->second].tss.push_back(gene_start);
                            }
                        }
                    } else {
//                        gene_id_map[gene_counter] = gene_id;
                        id_ensg_map[gene_id] = gene_counter;
                        gene_information this_gene;
                        this_gene.gene_id = gene_id;
                        this_gene.name = gene_name;
                        this_gene.chr = chr;
                        this_gene.tss = {gene_start};
                        promoter_map[gene_counter] = this_gene;
                        chr_gene_map[chr].insert(gene_counter);
                        gene_counter++;
                    }
                    hic_boundaries[chr] = {0, 0};
                }
            }
        }
    }
    pclose(Read_Gene_Annotation);

    // ____________________________________________________________
    // WRITE TEMPORARY GENE WINDOW FILE
    // ____________________________________________________________
    // First check if the peak file has a chr-prefix. This can happen if the executable is called independent of the
    // TEPIC bash script which pre-processes the peak file.
    // Also directly check for the signal columns in case n++ was given.
    ifstream peek_peak(b_peakfile);
    if (!peek_peak) {
        cout << "ERROR could not open the peak file\n" << b_peakfile << endl;
        return 1;
    }
    string line_peek;
    string peak_chr_prefix;
    bool has_chr_prefix = false;
    vector<int> col_indices;  // Store the activity column indices.
    int peak_file_rowlen = FilePeek(b_peakfile) * 4;
    while (!peek_peak.eof()) {
        getline(peek_peak, line_peek);
        if (line_peek.substr(0, 1) != "#") {
            if (line_peek.substr(0, 3) == "chr") {
                peak_chr_prefix = "chr";
                has_chr_prefix = true;
            }
            int total_cols = SplitTabLine(line_peek).size();
            if (activity_col.find('-') != string::npos) {
                int hyphon_pos = activity_col.find('-', 0);
                for (int i = stoi(activity_col.substr(0, hyphon_pos)) - 1;
                     i <= stoi(activity_col.substr(hyphon_pos + 1)) - 1; i++) {
                    col_indices.push_back(i);
                }
            } else if (activity_col.find('+') != string::npos) {
                for (int i = stoi(activity_col.substr(0, activity_col.find('+', 0))) - 1; i <= total_cols - 1; i++) {
                    col_indices.push_back(i);
                }
            } else if (activity_col.find(',') != string::npos) {
                stringstream find_cols(activity_col);
                string curr_val;
                while (getline(find_cols, curr_val, ',')) {
                    if (stoi(curr_val) > 0 and stoi(curr_val) <= total_cols) {
                        col_indices.push_back(stoi(curr_val) - 1);
                    }
                }
            } else {
                col_indices.push_back(stoi(activity_col) - 1);
            }
            break;
        }
    }
    peek_peak.close();
    int col_num = col_indices.size();  // Total number of activity columns to iterate through.
    unordered_map<int, string> colname_map = FileHeaderMapping(b_peakfile, col_indices);

    // Is not written directly on gtf-file read, as we need to find the most 5'-TSS first.
    int intersect_window = gene_windowsize;
    if (enh_windowsize > gene_windowsize) {
        intersect_window = enh_windowsize;
    }

    string temp_window_file = o_prefix + "_ABCpp_Temp_GeneWindow.bed";
    ofstream window_out(temp_window_file);
    Test_outfile(window_out, temp_window_file);
    for (auto const &gene : promoter_map) {  // gene_id, name, chr, set<int>{TSS}.
        // Gtf is 1-based, and bedtools 0-based and end-exclusive.
        int window_start = *min_element(gene.second.tss.begin(), gene.second.tss.end()) - intersect_window - 1;
        int window_end = *max_element(gene.second.tss.begin(), gene.second.tss.end()) + intersect_window;
        if (window_start < 0) {
            window_start = 0;
        }
        if (window_start < hic_boundaries[gene.second.chr][0]) {
            hic_boundaries[gene.second.chr][0] = window_start;
        }
        if (window_end > hic_boundaries[gene.second.chr][1]) {
            hic_boundaries[gene.second.chr][1] = window_end;
        }
        window_out << peak_chr_prefix + gene.second.chr << "\t" << window_start << "\t" << window_end << "\t"
                   << gene.first << "\n";
    }
    window_out.close();

    // ____________________________________________________________
    // READ ACTIVITIES AND WRITE BED-FILE WITH LOCATIONS ONLY
    // ____________________________________________________________
    // With many activity columns it's more efficient to do the intersection with the gene window based on a bed-file
    // without the activity columns. So we read the bed-file, store the activity columns and write a temporary
    // reduced bed-file with a peak_id as 4th column.
    unordered_map <int, peak_information> peak_info_map;  // Storing an object for each peak, indexed by peakID.
    unordered_map <int, string> peak_id_map;  // So we can identify peaks with an int instead of the full string.
    peak_info_map.reserve(100000);
    string temp_bed_file = o_prefix + "_ABCpp_Temp_Regions.bed";
    ofstream bed_out(temp_bed_file);
    Test_outfile(bed_out, temp_bed_file);
    FILE *Read_Peak_file;
    Read_Peak_file = fopen(b_peakfile.c_str(), "rb");
    char peak_buffer[peak_file_rowlen];
    string bed_line;
    int peak_counter = 0;
    while (!feof(Read_Peak_file)) {
        if (fgets(peak_buffer, peak_file_rowlen, Read_Peak_file) != nullptr) {
            bed_line = peak_buffer;
            if (bed_line.substr(0, 1) != "#") {
                vector <string> columns = SplitTabLine(bed_line);
                bed_out << columns[0] << "\t" << columns[1] << "\t" << columns[2] << "\t" << peak_counter << "\n";
                string peak_id_string = columns[0] + ":" + columns[1] + "-" + columns[2];
                if (has_chr_prefix) {  // Internally we remove the chr-prefix, for the intersection above we need it.
                    peak_id_string = peak_id_string.substr(3);
                }
                peak_information this_peak;
                this_peak.start = stoi(columns[1]);
                this_peak.end = stoi(columns[2]);
                for (int a_col : col_indices) {
                    this_peak.signal.push_back(stof(columns[a_col]));
                }
                peak_info_map[peak_counter] = this_peak;
                peak_id_map[peak_counter] = peak_id_string;
                peak_counter++;
            }
        }
    }
    bed_out.close();

    // ____________________________________________________________
    // FIND EXCLUDED REGIONS
    // ____________________________________________________________
    // Find the peaks that overlap with an excluded region to later not add them to the peaks per gene.
    unordered_set <int> excluded_peaks;  // Collect the matching peak_ids.
    if (!x_exclude_regions.empty()) {
        string exclude_regions_intersect =
                "bedtools intersect -a " + temp_bed_file + " -b " + x_exclude_regions + " -u 2>&1";
        FILE *exclude_stream = popen(exclude_regions_intersect.c_str(), "r");
        char exclude_buffer[128];  // We only need chr, start and end of the peak.
        if (exclude_stream) {
            string line;
            while (!feof(exclude_stream)) {
                if (fgets(exclude_buffer, 128, exclude_stream) != nullptr) {
                    line = exclude_buffer;
                    size_t pos = 0;
                    int cnt = 0;
                    while (cnt != 3) {  // Find the position of the third tab.
                        pos++;
                        pos = line.find('\t', pos);
                        if (pos == string::npos)
                            break;
                        cnt++;
                    }
                    // If the chr-prefix is not the same for peak and exclude-file, a warning will be printed but the program
                    // continues.
                    if (line.substr(0, 13) == "***** WARNING") {
                        cout << line << "\n" << "excluding regions not possible" << endl;
                    }
                    excluded_peaks.insert(stoi(line.substr(pos+1, line.size()-pos-2)));  // \tpeak_id\n
                }
            }
            pclose(exclude_stream);
        } else {
            cout << "WARNING Excluding regions not possible" << endl;
        }
    }

    // ____________________________________________________________
    // INTERSECT GENE WINDOWS AND PEAKS
    // ____________________________________________________________
    // Call the bedtools intersect command via popen.
    string bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + temp_bed_file + " -wo 2>&1";
    // Process the whole bed_output, by mapping the peak_ids to the respective genes and also map the location and
    // activity of the peaks in a separate map, and gather the genes for each chromosome.
    unordered_set <string> peak_chromosomes; // Store on which chromosomes we have an intersection.
    unordered_map <int, set<int>> gene_peak_map; // Storing the ids of peaks in the window.
    gene_peak_map.reserve(promoter_map.size());

    // Using popen to iterate over the rows and fill both maps for genes and peaks.
    int bed_buffer_size = 256;  // Only the genes and the locations of the peaks.
    char bed_buffer[bed_buffer_size];
    FILE *bed_intersect_stream = popen(bed_intersect_cmd.c_str(), "r");
    if (bed_intersect_stream) {
        string line;
        while (!feof(bed_intersect_stream)) {
            if (fgets(bed_buffer, bed_buffer_size, bed_intersect_stream) != nullptr) {
                line = bed_buffer;
                // First only fetch the required columns and ignore the potential activity columns.
                vector <string> columns = SplitTabLine(line);
                stringstream row_read(line);
                if (excluded_peaks.count(stoi(columns[7])) == 1) {
                    continue;
                }
                string this_chr = columns[0];
                if (has_chr_prefix) {  // Already checked that when we peeked at the first line in the peak file.
                    this_chr = this_chr.substr(3);
                }
                peak_chromosomes.insert(this_chr);
                int intersect_gene = stoi(columns[3]);  // Fetch the gene_id.
                int peak_id = stoi(columns[7]);
                gene_peak_map[intersect_gene].insert(peak_id);
            }
        }
        pclose(bed_intersect_stream);
    } else {
        cout << "ERROR Could not intersect gene windows and peaks, command was the following:\n"
             << bed_intersect_cmd << endl;
        return 1;
    }
    GetStdoutFromCommand("rm " + temp_window_file);
    GetStdoutFromCommand("rm " + temp_bed_file);

    vector <string> chromosomes;  // New datastructure for better sorting and iteration.
    for (const string& chr : peak_chromosomes) {
        chromosomes.push_back(chr);
    }
    // Sort chromosomes to have the large ones first, potentially speeds up the parallel part.
    std::sort(chromosomes.begin(), chromosomes.end(), [](const std::string &first, const std::string &second) {
        return Chr_sorter(first, second);
    });
    auto stop_intersect = chrono::high_resolution_clock::now();
    auto duration_intersect = chrono::duration_cast<chrono::seconds>(stop_intersect - start0);
    cout << duration_intersect.count() << "s ABCpp: intersected genes and peaks" << endl;

    // ____________________________________________________________
    // PREPARE OUTPUT FILES
    // ____________________________________________________________
    // Open output files already, as well as the GeneInfo files.
    vector <string> out_files;
    vector <string> gene_info_files;
    for (int a_col = 0; a_col < col_num; a_col++) {
        string column_suffix = colname_map[col_indices[a_col]];
        // Open the output and write the header already. Only the first string has to be explicitly casted.
        ofstream abc_output(o_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        Test_outfile(abc_output, o_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        abc_output << string("#chr") + "\t" + "start" + "\t" + "end" + "\t" + "Ensembl ID" + "\t" + "Gene Name" +
                      "\t" + "PeakID" + "\t" + "signalValue" + "\t" + "Contact" + "\t" + activity_header + "\t" +
                      "scaledContact" + "\t" + "intergenicScore" + "\t" + "TSS-dist" + "\t" + "ABC-Score" + "\n";
        out_files.push_back(o_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        abc_output.close();
        // Also already write the GeneInfo header.
        ofstream gene_info_out(o_prefix + "_GeneInfo" + column_suffix + ".txt");
        gene_info_out << string("Ensembl ID") + "\t" + "Gene Name" + "\t" + "chr" + "\t" + "TSS" + "\t" + "#Enhancer"
                         + "\t" + "Avg_EnhancerActivity" + "\t" + "Avg_EnhancerContact" + "\t" +
                         "Avg_EnhancerDistance" + "\t" + "Failure" + "\n";
        gene_info_files.push_back(o_prefix + "_GeneInfo" + column_suffix + ".txt");
        gene_info_out.close();
    }

    // ____________________________________________________________
    // READ GENES TO FILTER FOR
    // ____________________________________________________________
    unordered_set<string> filter_genes;
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
    // ITERATE THROUGH THE CHROMOSOMES
    // ____________________________________________________________
    // First see which contact files are available.
    vector<string> contact_files;
    if (!f_contactfolder.empty()) {
        string ls_out = GetStdoutFromCommand("ls " + f_contactfolder + "/");
        istringstream read_ls_out(ls_out);
        if (!read_ls_out) {
            cout << "ERROR Could not access the folder with the contact files, command was the following:\nls " +
                    f_contactfolder + "/"
                 << endl;
            return 1;
        }
        for (string line; getline(read_ls_out, line);) {
            contact_files.push_back(line);
        }
    }

    // Prepare a structure to store genes that are discarded during processing.
    unordered_set <int> genes_wo_hic;  // Missing normalized contact file
    unordered_set <int> genes_wo_candidates;  // No candidate regions in gene window (ABC)
    unordered_set <int> written_genes;  // To track the genes we already wrote to the GeneInfo.
    written_genes.reserve(promoter_map.size());
    vector <omp_lock_t> file_locks;
    for (int i = 0; i < col_num; i++) {
        omp_lock_t new_lock;
        omp_init_lock(&new_lock);
        file_locks.push_back(new_lock);
    }
    // Iterate over the chromosomes, read the normalized matrix in, and iterate over the genes on the chr to calculate
    // the ABC-score.
#pragma omp parallel for num_threads(cores)
    for (unsigned int chr_idx = 0; chr_idx < chromosomes.size(); chr_idx++) {  // Can't be range-based due to openmp.
        string chr = chromosomes[chr_idx];
        if (chr_gene_map[chr].empty()) {  // In case we don't have any genes on that chromosome.
            continue;
        }
        // Already define the variables that are only filled with true Hi-C data, to have them in the outer scope.
        int min_bin, max_bin;
        double slope, intercept, close_pseudocount;
        boost::numeric::ublas::mapped_matrix<double> contact_matrix;
        if (!f_contactfolder.empty()) {
            // Find the matching contact file.
            string chr_file;
            regex re_chr(".*(C|c)hr" + chr + "[^0-9]{1}.*gz");
            for (const string &file : contact_files) {
                if (regex_match(file, re_chr)) {
                    chr_file = file;
                    cout << chr << " mapped to " << chr_file << endl;
                    continue;
                }
            }

            auto start_hic = chrono::high_resolution_clock::now();
            // Check in the first line of the file if there is an integer present. If the zcat returns an error, it will be
            // a character at the first position as well.
            string temp_row = "foo";
            // There were very odd crashes with a certain number of threads (≥16) in this section. The omp critical is
            // a workaround.
            FILE *peek_stream;
            string peek_cmd = string("zcat < ") + f_contactfolder + "/" + chr_file + " 2>&1 | head -n 1";
#pragma omp critical(peek_open)
            {
                peek_stream = popen(peek_cmd.c_str(), "r");
            }
            char dummy[1024];
            if (peek_stream != nullptr) {
                if (fgets(dummy, 1024, peek_stream) != nullptr) {
                    temp_row = dummy;
                }
                pclose(peek_stream);
            }
            regex int_match("[0-9]");
            if (temp_row == "foo" or chr_file.empty() or !regex_match(temp_row.substr(0, 1), int_match)) {
                for (const int &gene : chr_gene_map[chr]) {
                    genes_wo_hic.insert(gene);
                }
                cout << chr << " no contact matrix" << endl;
                continue;
            }

            // ____________________________________________________________
            // READ THE HI-C FILE
            // ____________________________________________________________
            // Find the highest bin and construct a max_bin*max_bin matrix, then filling the values into the upper
            // triangle.
            auto start_read = chrono::high_resolution_clock::now();
            min_bin = hic_boundaries[chr][0] / bin_size;
            max_bin = hic_boundaries[chr][1] / bin_size;
            // For the pseudocount we are not restricted to the hic_boundaries set by the genes.
            int min_pseudo_bin = numeric_limits<int>::max();
            int max_pseudo_bin = 0;
            // +1 To include max_bin.
            contact_matrix.resize(max_bin + 1 - min_bin, max_bin + 1 - min_bin, false);

            // Track the contacts and the number of non-zero for up to 1MB.
            vector<double> contact_sums(1000000 / bin_size);  // Store the summed contact per distance for the mean.
            FILE *hic_stream;
            string hic_cmd = string("zcat < ") + f_contactfolder + "/" + chr_file + " 2>&1";
#pragma omp critical(hic_open)
            {
                hic_stream = popen(hic_cmd.c_str(), "r");
            }
            char hic_buffer[256];
            string hic_row;
            string delimiter = "\t";
            while (!feof(hic_stream)) {
                if (fgets(hic_buffer, 256, hic_stream) != nullptr) {
                    hic_row = hic_buffer;
                    size_t first_tab = hic_row.find(delimiter);
                    size_t second_tab = hic_row.find(delimiter, first_tab + 1);
                    double contact = stof(hic_row.substr(second_tab + 1));
                    if (!isnan(contact) and contact > 0) {  // 0 Shouldn't appear if preprocessing done correctly.
                        int first_bin = stoi(hic_row.substr(0, first_tab)) / bin_size;
                        int second_bin =
                                stoi(hic_row.substr(first_tab + 1, hic_row.find(delimiter, second_tab))) / bin_size;
                        int bin1 = first_bin;
                        int bin2 = second_bin;
                        if (first_bin > second_bin) {  // Always index with the smaller bin first.
                            bin1 = second_bin;
                            bin2 = first_bin;
                        }
                        if (bin1 < min_pseudo_bin) {
                            min_pseudo_bin = bin1;
                        }
                        if (bin2 > max_pseudo_bin) {
                            max_pseudo_bin = bin2;
                        }
                        int distance = (bin2 - bin1) * bin_size;
                        if (bin_size <= distance and distance < 1000000) {
                            contact_sums[distance / bin_size - 1] += contact;
                        }
                        // The gene windows may not cover all available Hi-C bins.
                        if (bin1 >= min_bin and bin2 <= max_bin) {
                            contact_matrix(bin1 - min_bin, bin2 - min_bin) = contact;
                        }
                    }
                }
            }
#pragma omp critical(hic_close)
            {
                pclose(hic_stream);
            }


            // Fit a linear function to the mean contact at each distance within 1MB.
            double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, regression_count = 0;
            for (int c = 0; c < contact_sums.size(); c++) {
                if (contact_sums[c] > 0) {
                    regression_count++;
                    int possible_hits = (max_pseudo_bin - min_pseudo_bin) - (c + 1);
                    double log_distance = log((c + 1) * bin_size);  // Is base-e logarithm. Reverse is exp().
                    double log_contact = log(contact_sums[c] / possible_hits);
                    sumX += log_distance;
                    sumX2 += log_distance * log_distance;
                    sumY += log_contact;
                    sumXY += log_distance * log_contact;
                }
            }

            slope = (regression_count * sumXY - sumX * sumY) / (regression_count * sumX2 - sumX * sumX);
            intercept = (sumY - slope * sumX) / regression_count;
            if (do_pseudocount) {
                close_pseudocount = exp(slope * log(1000000) + intercept);  // For all distance <1 Mb.
                if (isnan(slope)) {  // In case a fit is not possible.
                    slope = 0;
                    intercept = 0;
                    close_pseudocount = 0;
                    cout << "Fit for pseudocount not possible, proceeding without pseudocounts." << endl;
                }
            } else {
                slope = 0;
                intercept = 0;
                close_pseudocount = 0;
            }

            // Correct the diagonal by replacing it with the maximum contact of its neighbours.
            for (unsigned i = 1;
                 i < contact_matrix.size1() - 1; ++i) {  // mapped_matrix returns 0 if no value was filled.
                double upper_neighbour = contact_matrix(i - 1, i);
                double right_neighbour = contact_matrix(i, i + 1);
                contact_matrix(i, i) = max(upper_neighbour, right_neighbour);
            }
            contact_matrix(0, 0) = contact_matrix(0, 1);  // Deal with the border cases manually.
            contact_matrix(max_bin - min_bin, max_bin - min_bin) = contact_matrix(max_bin - 1 - min_bin,
                                                                                  max_bin - min_bin);

            auto stop_hic = chrono::high_resolution_clock::now();
            auto duration_hic = chrono::duration_cast<chrono::milliseconds>(stop_hic - start_hic);
            cout << duration_hic.count() << "ms; " << chr << " contact matrix read" << endl;
        }
        // ____________________________________________________________
        // GET PEAK-GENE CONTACTS AND SUM THE CONTACT PER ENHANCER
        // ____________________________________________________________
        // Initialize an empty 2D vector to gather the interactions > threshold for each chromosomes.
        // Contact|scaledContact| + activity_cols x (signalValue|scaledActivity|ABC-Score)
        map<int, vector<abc_hit>> exceeders;  // Map with activity_cols, gathering all hits > threshold.
        double chr_max = 0;
        auto start_scores = chrono::high_resolution_clock::now();
        // Store contacts in same order as gene_peak_map.
        unordered_map<int, map<int, vector<double>>> gene_contact_map; //  {gene: {tss: [contacts]}}
        gene_contact_map.reserve(chr_gene_map[chr].size());
        for (const int &chr_gene : chr_gene_map[chr]) {
            if (gene_peak_map.count(chr_gene) == 1) {
                for (int gene_tss : promoter_map[chr_gene].tss) {
                    int gene_bin = gene_tss / bin_size;
                    vector<double> these_contacts;
                    // We collected all peaks that are in the window of any of the TSS of a gene.
                    for (auto peak = gene_peak_map[chr_gene].begin(); peak != gene_peak_map[chr_gene].end();) {
                        peak_information &matching_peak = peak_info_map[*peak];
                        double this_contact;
                        int peak_bin = ((matching_peak.end + matching_peak.start) / 2) / bin_size;
                        int distance_to_min = GetDistance(matching_peak.start, matching_peak.end, *min_element(promoter_map[chr_gene].tss.begin(), promoter_map[chr_gene].tss.end()));
                        int distance_to_max = GetDistance(matching_peak.start, matching_peak.end, *max_element(promoter_map[chr_gene].tss.begin(), promoter_map[chr_gene].tss.end()));
                        int distance = min(distance_to_min, distance_to_max);  // Closest TSS to that peak.

                        if (f_contactfolder.empty()) {
                            // Fractal module, approximate contact, but only for >5kb.
                            this_contact = pow(max(distance, 5000), -1);
                        } else {
                            double pseudocount = close_pseudocount;
                            if (abs(gene_bin - peak_bin) * bin_size > 1000000) {
                                pseudocount = exp(slope * log(abs(gene_bin - peak_bin) * bin_size) + intercept);
                            }
                            if ((peak_bin > max_bin) or (peak_bin < min_bin)) {
                                this_contact = pseudocount;
                            } else {
                                this_contact = pseudocount + contact_matrix(min(gene_bin, peak_bin) - min_bin,
                                                                            max(gene_bin, peak_bin) - min_bin);
                            }
                        }
                        matching_peak.contact_sum += this_contact;
                        if (distance <= gene_windowsize) {
                            if (this_contact > chr_max) {
                                chr_max = this_contact;
                            }
                            these_contacts.push_back(this_contact);
                            ++peak;
                        } else {  // Remove peaks that are only in the enh-window, only needed them for the contact_sum.
                            peak = gene_peak_map[chr_gene].erase(peak);
                        }
                    }
                    gene_contact_map[chr_gene][gene_tss] = these_contacts;
                }
            }
        }

        // ____________________________________________________________
        // SCORING FOR EACH GENE ON THE CHROMOSOME
        // ____________________________________________________________
        // Now that we mapped the peaks to the genes and fetched the respective information on the peaks, we can do
        // the gene-wise scoring.
        unordered_map<int, double> gene_contact_scaler; // No need to store it for each interaction of a gene.
        gene_contact_scaler.reserve(chr_gene_map[chr].size());
        for (const int &chr_gene : chr_gene_map[chr]) {
            gene_information &curr_gene = promoter_map[chr_gene];
            // Check if we can skip that gene, if it's not in the u_genefile with its ID or symbol.
            if (u_genefile.length() != 0 and u_genefile != "0") {
                if ((filter_genes.count(curr_gene.gene_id) == 0) and (filter_genes.count(curr_gene.name) == 0)) {
                    continue;
                }
            }
            auto gene_iter = gene_peak_map.find(chr_gene);
            if (gene_iter == gene_peak_map.end()) {  // Can also happen if all were overlapping with excluded regions.
                genes_wo_candidates.insert(chr_gene);
                continue;
            }
            int num_candidates = gene_iter->second.size();
            if (num_candidates == 0) {  // If all regions were only in the enhancer window, we still have an empty set.
                genes_wo_candidates.insert(chr_gene);
                continue;
            }
            // activity_cols x (scaled/adjustedActivity | Activity*Contact)
            vector<vector<double>> candidate_activities(num_candidates, vector<double>(2 * (col_num)));
            vector<double> candidate_contacts(num_candidates);
            vector<int> candidate_ids(num_candidates);  // Stores the peak_id.
            vector<double> max_signals(col_num);
            vector<double> abc_sums(col_num);
            double max_contact = 0;
            int peak_count_helper = 0;

            for (const int &peak_id : gene_iter->second) {
                peak_information &matching_peak = peak_info_map[peak_id];
                for (int tss : curr_gene.tss) {  // Sum up all information across TSS to later get the average.
                    // Should always have an entry if there were still num_candidates.
                    double current_contact = gene_contact_map[chr_gene][tss][peak_count_helper];
                    candidate_contacts[peak_count_helper] += current_contact;  // Sum over all TSS, A*C is further down.

                    for (int a_col = 0; a_col < col_num; a_col++) {
                        double current_activity = matching_peak.signal[a_col];
                        if (current_activity == 0) {
                            continue;
                        }
                        if (do_adjusted_abc) {
                            current_activity = current_activity * (current_contact / matching_peak.contact_sum);
                        }
                        candidate_activities[peak_count_helper][a_col * 2] += current_activity;
                        candidate_activities[peak_count_helper][a_col * 2 + 1] += current_activity * current_contact;
                        abc_sums[a_col] += current_activity * current_contact;
                    }
                    // Add the interaction so we can later trace back which peak and gene it was.
                    candidate_ids[peak_count_helper] = peak_id;
                }
                if (candidate_contacts[peak_count_helper] > max_contact) {
                    max_contact = candidate_contacts[peak_count_helper];
                }
                if (not do_adjusted_abc) {  // Only need to search for the maximum if we don't adjust the activity.
                    for (int a_col = 0; a_col < col_num; a_col++) {
                        if (candidate_activities[peak_count_helper][a_col * 2] > max_signals[a_col]) {
                            max_signals[a_col] =
                                    candidate_activities[peak_count_helper][a_col * 2] / curr_gene.tss.size();
                        }
                    }
                }
                peak_count_helper++;
            }
            // Rescaling the contact to a maximum of 100 only has to be done once. For the activity as many times as
            // there are activity columns. If none of the peaks has contact with the gene, the contact_rescaler and
            // thus the scores will be set to 0.
            gene_contact_scaler[chr_gene] = (max_contact != 0) ? 100 / (max_contact / curr_gene.tss.size()) : 0;
            // Write the relative ABC-score. And activity if not adjusting it.
            for (int r = 0; r < peak_count_helper; ++r) {
                double this_contact = candidate_contacts[r] / curr_gene.tss.size();  // Avg across TSS.
                for (int a_col = 0; a_col < col_num; a_col++) {
                    if (candidate_activities[r][a_col * 2] > 0) {  // W/o activity it will be 0 anyway.
                        double this_sum = abc_sums[a_col];
                        if (this_sum != 0) {
                            double this_score = candidate_activities[r][a_col * 2 + 1] / this_sum;
                            if (this_score >= abc_cutoff) {
                                abc_hit this_hit;
                                this_hit.gene_id = chr_gene;
                                // Before we stored this sequentially, not by peak-id, we save the lookup later.
                                this_hit.peak_id = candidate_ids[r];
                                if (not do_adjusted_abc) {
                                    this_hit.scaledActivity = (candidate_activities[r][a_col * 2] /
                                                               curr_gene.tss.size()) * (100 / max_signals[a_col]);
                                } else {
                                    this_hit.scaledActivity = candidate_activities[r][a_col * 2] / curr_gene.tss.size();
                                }
                                this_hit.contact = this_contact;
                                double distance = 0;
                                for (int tss : curr_gene.tss) {
                                    distance += GetDistance(peak_info_map[this_hit.peak_id].start, peak_info_map[this_hit.peak_id].end, tss);
                                }
                                this_hit.distance = distance / curr_gene.tss.size();
                                this_hit.score = this_score;
                                exceeders[a_col].push_back(this_hit);
                            }
                        }
                    }
                }
            }
        }

        auto stop_scores = chrono::high_resolution_clock::now();
        auto duration_scores = chrono::duration_cast<chrono::milliseconds>(stop_scores - start_scores);
        cout << duration_scores.count() << "ms; " << chr << " Scores processed" << endl;

        auto start_out = chrono::high_resolution_clock::now();
        // Shuffle the order in which the files are written to not block the locks by one chromosome.
        vector<int> a_col_shuffle;
        a_col_shuffle.reserve(col_num);
        for (int i = 0; i < col_num; i++) {
            a_col_shuffle.push_back(i);
        }
        random_device rd;
        mt19937 g(rd());
        shuffle(begin(a_col_shuffle), end(a_col_shuffle), g);

        for (int a_idx = 0; a_idx < col_num; a_idx++) {
            int a_col = a_col_shuffle[a_idx];
            if (!exceeders[a_col].empty()) {
                omp_set_lock(&file_locks[a_col]);  // Only one chromosome can write into the same file at a time.
                // Repeat the process for each activity column.
                ofstream out_stream(out_files[a_col], ios_base::app);
                ofstream gene_info_stream(gene_info_files[a_col], ios_base::app);
                gene_information gene_map = promoter_map[exceeders[a_col][0].gene_id];  // Get first entry, update for new genes.
                vector<double> gene_info(4, 0);
                int hit_helper = 0;
                for (auto hit : exceeders[a_col]) {
                    auto matching_peak = peak_info_map[hit.peak_id];
                    out_stream << gene_map.chr << "\t" << matching_peak.start << "\t" << matching_peak.end;
                    out_stream << "\t" << gene_map.gene_id << "\t" << gene_map.name;
                    out_stream << "\t" << peak_id_map[hit.peak_id];  // Get back the string version of the id.
                    out_stream << "\t" << setprecision(8) << peak_info_map[hit.peak_id].signal[a_col];  // signalValue
                    out_stream << "\t" << setprecision(8) << hit.contact;  // Contact
                    out_stream << "\t" << setprecision(8) << hit.scaledActivity;  // scaledActivity
                    out_stream << "\t" << setprecision(8) << hit.contact * gene_contact_scaler[hit.gene_id];  // scaledContact
                    // adjustedActivity OR chr-scaled Contact * signalValue
                    if (do_adjusted_abc) {
                      out_stream << "\t" << setprecision(8) << hit.scaledActivity;
                    } else {
                      out_stream << "\t" << setprecision(8) << hit.contact /
                      chr_max * peak_info_map[hit.peak_id].signal[a_col];
                    }
                    out_stream << "\t" << setprecision(8) << hit.distance;  // Distance
                    out_stream << "\t" << setprecision(8) << hit.score << "\n";  // ABC-Score
                    // Sum up the attributes of the interaction.
                    gene_info[0] += 1;
                    gene_info[1] += peak_info_map[hit.peak_id].signal[a_col];
                    gene_info[2] += hit.contact;
                    gene_info[3] += hit.distance;

                    // For the gene info keep track of the current gene, we filled the data ordered by genes.
                    if (hit_helper == exceeders[a_col].size() - 1 or exceeders[a_col][hit_helper+1].gene_id != hit.gene_id) {
                        written_genes.insert(hit.gene_id);
                        gene_info_stream << gene_map.gene_id << "\t" << gene_map.name << "\t" << gene_map.chr << "\t";
                        int comma_helper = 0;  // To make it csv without a trailing ",".
                        for (int tss : gene_map.tss) {
                            if (comma_helper > 0) {
                                gene_info_stream << ",";
                            }
                            gene_info_stream << to_string(tss);
                            comma_helper++;
                        }
                        int num_enhancer = static_cast<int>(gene_info[0]);
                        if (num_enhancer > 0) {
                            gene_info_stream << "\t" << num_enhancer << "\t" << gene_info[1] / num_enhancer <<
                                             "\t" << gene_info[2] / num_enhancer << "\t"
                                             << gene_info[3] / num_enhancer << "\t-\n";
                        } else {
                            gene_info_stream << "\t0\t0\t0\t0\tScore of candidates too low\n";
                        }
                        fill(gene_info.begin(), gene_info.end(), 0);  // Reset the vector.
                        if (hit_helper < exceeders[a_col].size() - 1) {  // Reduces the look up for each interaction.
                            gene_map = promoter_map[exceeders[a_col][hit_helper+1].gene_id];
                        }
                    }
                    hit_helper++;
                }
                out_stream.close();
                gene_info_stream.close();
                omp_unset_lock(&file_locks[a_col]);  // After writing, open the lock again.
            }
        }
        auto stop_out = chrono::high_resolution_clock::now();
        auto duration_out = chrono::duration_cast<chrono::milliseconds>(stop_out - start_out);
        cout << duration_out.count() << "ms; " << chr << " written to output" << endl;
    }
    for (int i = 0; i < col_num; i++) {
        omp_destroy_lock(&file_locks[i]);
    }

    // ____________________________________________________________
    // COMPLETE GENEINFO OUTPUT FILES AND COMPRESS OUTPUT
    // ____________________________________________________________
    auto start_compression = chrono::high_resolution_clock::now();
    cout << "Compressing files" << endl;
#pragma omp parallel for num_threads(cores)
    for (int a_col = 0; a_col < col_num; a_col++) {
        ofstream gene_info_stream(gene_info_files[a_col], ios_base::app);
        for (const auto& entry: promoter_map) {
            int gene_id = entry.first;
            auto gene_map = entry.second;
            if (u_genefile.length() != 0 and u_genefile != "0") {
                if ((filter_genes.count(gene_map.gene_id) == 0) and (filter_genes.count(gene_map.name) == 0)) {
                    continue;
                }
            }
            string fail_message;
            if (genes_wo_hic.count(gene_id) == 1) {
                fail_message = "Missing normalized contact file";
            } else if (genes_wo_candidates.count(gene_id) == 1) {
                fail_message = "No candidate regions in gene window";
            } else if (written_genes.count(gene_id) == 0) {
                if (peak_chromosomes.count(gene_map.chr) == 0) {
                    fail_message = "No candidate regions in gene window";
                } else {
                    fail_message = "Score of candidate regions too low";
                }
            }
            if (!fail_message.empty()) {
                gene_info_stream << gene_map.gene_id << "\t" << gene_map.name << "\t" << gene_map.chr << "\t";
                int comma_helper = 0;  // To make it csv without a trailing ",".
                for (int tss : gene_map.tss) {
                    if (comma_helper > 0) {
                        gene_info_stream << ",";
                    }
                    gene_info_stream << to_string(tss);
                    comma_helper++;
                }
                gene_info_stream << "\t0\t0\t0\t0\t" << fail_message << "\n";
            }
        }
        gene_info_stream.close();

        GzipFile(out_files[a_col]);
        GzipFile(gene_info_files[a_col]);
    }
    auto stop_compression = chrono::high_resolution_clock::now();
    auto duration_compression = chrono::duration_cast<chrono::milliseconds>(stop_compression - start_compression);
    cout << duration_compression.count() << "ms compression" << endl;

    auto stop0 = chrono::high_resolution_clock::now();
    auto duration0 = chrono::duration_cast<chrono::seconds>(stop0 - start0);
    cout << duration0.count() << "s ABCpp in total" << endl;
}



