//
// Created by Dennis Hecker on 10.12.21.
//

# include <iostream>
# include <string>
# include <vector>
# include <numeric>
# include <fstream>
# include <math.h>
# include <chrono>
# include <regex>
# include <sstream>
# include <algorithm>
# include <unordered_map>
# include <unordered_set>
# include <map>
# include <set>
# include <getopt.h>
# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <unistd.h>
# include <omp.h>
# include <boost/numeric/ublas/matrix_sparse.hpp>
# include "STARE_MiscFunctions.h"

/*
 * When having chromatin contact data at hand, it is possible to get the regulatory interactions based on the
 * ABC-score. Either the regular version (-q False) or the adapted are possible (-q True).
 * Together with STARE for deriving affinities of TFs to genes, we can use these called interactions instead of all
 * the open regions in a window. This program also works independently of STARE.
 *
 * For MacOS:
 * g++ STARE_ABCpp.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -lomp -std=c++11 -O3 -o STARE_ABCpp
 * For Linux:
 * g++ STARE_ABCpp.cpp STARE_MiscFunctions.cpp -fopenmp -std=c++11 -O3 -o STARE_ABCpp
 *
 * ./STARE_ABCpp -b ../Test/ABC_example_regions.bed -n 4 -a ../Test/ABC_example_annotation.gtf -gw 5000000 -bin 5000 -cf path-to-normalized-HiC-contacts -c 4 -t ABC-cut-off -d output-folder
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 * C++ implementation of the ABC-scoring principle from Fulco et al.: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
 */

class peak_information {
public:
    int start;
    int end;
    std::vector<double> signal;
    double contact_sum = 0;
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
                            "\n-gw genewindow-size"
                            "\n-cf folder with the normalized hic-files per chromosome"
                            "\n-bin binsize of the hic-files"
                            "\n-t cut-off for the ABC-score (default 0.02), set to 0 to get all scored interactions"
                            "\n-d prefix with which the files are written"
                            "\n-c number of cores used (default 1)"
                            "\n-x file with regions to be excluded"
                            "\n-p whether to use pseudocount for contact frequency (default True)"
                            "\n-q whether to adjust the activity for an enhancer's contacts (default True)"
                            "\n-m enhancer window size in which to consider contacts for adjustment";

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
    string b_peakfile, a_promoterfile, n_activitycol, gw_genewindowsize, bin_binsize, t_abc_cutoff, cf_contactfolder,
            c_cores, d_prefix, x_exclude_regions, p_pseudocount, q_adjusted_abc, m_enhwindowsize;
    static struct option long_options[] =
            {
                    {"b",   required_argument, NULL, 'b'},
                    {"a",   required_argument, NULL, 'a'},
                    {"n",   required_argument, NULL, 'n'},
                    {"gw",  required_argument, NULL, 'w'},
                    {"bin", required_argument, NULL, 'i'},
                    {"t",   required_argument, NULL, 't'},
                    {"cf",  required_argument, NULL, 'f'},
                    {"c",   required_argument, NULL, 'c'},
                    {"d",   required_argument, NULL, 'd'},
                    {"x",   required_argument, NULL, 'x'},
                    {"p",   required_argument, NULL, 'p'},
                    {"q",   required_argument, NULL, 'q'},
                    {"m",   required_argument, NULL, 'm'},
                    {NULL,  0,                 NULL, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "b:a:n:w:i:t:c:f:d:x:p:h:m:", long_options, NULL)) != -1) {
        switch (arg) {
            case 'b':
                b_peakfile = optarg;
                break;
            case 'a':
                a_promoterfile = optarg;
                break;
            case 'n':
                n_activitycol = optarg;
                break;
            case 'w':
                gw_genewindowsize = optarg;
                break;
            case 'i':
                bin_binsize = optarg;
                break;
            case 't':
                t_abc_cutoff = optarg;
                break;
            case 'f':
                cf_contactfolder = optarg;
                break;
            case 'c':
                c_cores = optarg;
                break;
            case 'd':
                d_prefix = optarg;
                break;
            case 'x':
                x_exclude_regions = optarg;
                break;
            case 'p':
                p_pseudocount = optarg;
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
    for (string p: {b_peakfile, a_promoterfile, cf_contactfolder, d_prefix}) {
        if (p.size() < 1) {
            cout << "Required parameter missing" << endl;
            cout << parameter_help << endl;
            return 1;
        }
    }

    // If not given as input, set to default.
    int gene_windowsize = stoi(SetOptionalInput(gw_genewindowsize, "5000000")) / 2;
    int bin_size = stoi(SetOptionalInput(bin_binsize, "5000"));
    string activity_col = SetOptionalInput(n_activitycol, "4");
    int cores = stoi(SetOptionalInput(c_cores, "1"));
    double abc_cutoff = stof(SetOptionalInput(t_abc_cutoff, "0.02"));
    bool do_pseudocount = SetBoolInput(p_pseudocount, true);
    int enh_windowsize = stoi(SetOptionalInput(m_enhwindowsize, "5000000")) / 2;
    bool do_adjusted_abc = SetBoolInput(q_adjusted_abc, true);

    string activity_header = "adjustedActivity";
    int intergenic_activity_col = 3;  // Vector-index where the contact-adjusted activity is stored.
    if (not do_adjusted_abc) {
        enh_windowsize = gene_windowsize;
        activity_header = "scaledActivity";
        intergenic_activity_col = 2;  // Vector-index for the unmodified region activity.
    }

    // ____________________________________________________________
    // PROCESS GTF GENE ANNOTATION
    // ____________________________________________________________
    map <string, unordered_set<string>> chr_gene_map;  // For each chromosome stores its genes.
    unordered_map <string, vector<string>> promoter_map;  // For ID, name, chr, TSS.
    promoter_map.reserve(10000);
    unordered_map <string, vector<int>> hic_boundaries;  // Store the min/max boundaries for each chromosome.
    cout << "Reading Gene annotation " << endl;
    int annot_len = FilePeek(a_promoterfile) * 3;
    FILE *Read_Gene_Annotation = fopen(a_promoterfile.c_str(), "rb");
    char annot_buffer[annot_len];
    string row;
    bool first_line = true;  // To only look into the first line without the need to count lines.
    bool gene_chr_prefix;
    while (!feof(Read_Gene_Annotation)) {
        if (fgets(annot_buffer, annot_len, Read_Gene_Annotation) != NULL) {
            row = annot_buffer;
            if ((row.size() > 1) and (row.rfind("#", 0) != 0)) {
                vector <string> columns = SplitTabLine(row);
                if (columns[2] == "gene") {
                    string id_delimiter = "gene_id ";
                    int id_start = columns[8].find(id_delimiter, 0) + id_delimiter.size();
                    int id_end = columns[8].find(";", 0);
                    string gene_id = columns[8].substr(id_start + 1,
                                                       id_end - id_start - 2);  // Remove the quotation marks.
                    string name_delimiter = "gene_name ";
                    int name_start = columns[8].find(name_delimiter, 0) + name_delimiter.size();
                    int name_end = columns[8].find(";", name_start);
                    string gene_name = columns[8].substr(name_start + 1, name_end - name_start -
                                                                         2);  // Remove the quotation marks.
                    string strand = columns[6];
                    string chr = columns[0];
                    if (first_line) {
                        if (chr.substr(0, 3) == "chr") {
                            gene_chr_prefix = true;
                        }
                        else {
                            gene_chr_prefix = false;
                        }
                        first_line = false;
                    }
                    if (gene_chr_prefix) {
                        chr = chr.substr(3);
                    }
                    string gene_start;

                    if (strand == "+") {
                        gene_start = columns[3];
                    } else {
                        gene_start = columns[4];
                    }
                    if (promoter_map.find(gene_id) != promoter_map.end()) {  // Replace if the new one is more in 5'.
                        if (((stoi(promoter_map[gene_id][3]) > stoi(gene_start)) and (strand == "+")) or
                            ((stoi(promoter_map[gene_id][3]) < stoi(gene_start)) and (strand == "-"))) {
                            promoter_map[gene_id][3] = gene_start;
                        }
                    } else {
                        promoter_map[gene_id] = vector < string > {gene_id, gene_name, chr, gene_start};
                    }
                    hic_boundaries[chr] = {0, 0};
                    chr_gene_map[chr].insert(gene_id);
                }
            }
        }
    }
    fclose(Read_Gene_Annotation);

    // ____________________________________________________________
    // WRITE TEMPORARY GENE WINDOW FILE
    // ____________________________________________________________
    // First check if the peak file has a chr-prefix. This can happen if the executable is called independent of the
    // TEPIC bash script which pre-processes the peak file.
    // Also directly check for the signal columns in case n++ was given.
    cout << "writing temp gene window file" << endl;
    ifstream peek_peak(b_peakfile);
    if (!peek_peak) {
        cout << "ERROR could not open the peak file\n" << b_peakfile << endl;
        return 1;
    }
    string line_peek;
    string peak_chr_prefix = "";
    bool has_chr_prefix = false;
    int start_col, last_col;
    int peak_file_rowlen;
    while (!peek_peak.eof()) {
        getline(peek_peak, line_peek);
        if (line_peek.substr(0, 1) != "#") {
            if (line_peek.substr(0, 3) == "chr") {
                peak_chr_prefix = "chr";
                has_chr_prefix = true;
            }
            int total_cols = SplitTabLine(line_peek).size();
            peak_file_rowlen = total_cols * 21;  // Spent 20 Characters for each column + the tabs.
            if (line_peek.size() > peak_file_rowlen) {  // Unless the line is even longer.
                peak_file_rowlen = line_peek.size() * 1.2;  // With additional buffer for longer lines.
            }
            if (activity_col.find('-') != string::npos) {
                int hyphon_pos = activity_col.find("-", 0);
                start_col = stoi(activity_col.substr(0, hyphon_pos)) - 1;
                last_col = stoi(activity_col.substr(hyphon_pos + 1)) - 1;
            } else if (activity_col.find('+') != string::npos) {
                start_col = stoi(activity_col.substr(0, activity_col.find('+', 0))) - 1;
                last_col = total_cols - 1;
            } else {
                start_col = stoi(activity_col) - 1;
                last_col = start_col;
            }
            break;
        }
    }
    int col_num = last_col - start_col;  // Total number of activity columns to iterate through.
    unordered_map<int, string> colname_map = FileHeaderMapping(b_peakfile, start_col, last_col);

    // Is not written directly on gtf-file read, as we need to find the most 5'-TSS first.
    int intersect_window = gene_windowsize;
    if (enh_windowsize > gene_windowsize) {
        intersect_window = enh_windowsize;
    }
    string temp_window_file = d_prefix + "_ABCpp_Temp_GeneWindow.bed";
    ofstream window_out(temp_window_file);
    Test_outfile(window_out, temp_window_file);
    for (auto const &gene : promoter_map) {  // [geneID]: ID, name, chr, TSS.
        int window_start = stoi(gene.second[3]) - intersect_window;
        int window_end = stoi(gene.second[3]) + intersect_window;
        if (window_start < 0) {
            window_start = 0;
        }
        if (window_start < hic_boundaries[gene.second[2]][0]) {
            hic_boundaries[gene.second[2]][0] = window_start;
        }
        if (window_end > hic_boundaries[gene.second[2]][1]) {
            hic_boundaries[gene.second[2]][1] = window_end;
        }
        window_out << peak_chr_prefix + gene.second[2] << "\t" << to_string(window_start) << "\t"
                   << to_string(window_end) << "\t" << gene.first << "\n";
    }
    window_out.close();

    // ____________________________________________________________
    // FIND EXCLUDED REGIONS
    // ____________________________________________________________
    // Find the peaks that overlap with an excluded region to later not add them to the peaks per gene.
    unordered_set <string> excluded_peaks;
    if (x_exclude_regions.size() > 0) {
        string exclude_regions_intersect =
                "bedtools intersect -a " + b_peakfile + " -b " + x_exclude_regions + " -u 2>&1";
        FILE *exclude_stream = popen(exclude_regions_intersect.c_str(), "r");
        char exclude_buffer[128];  // We only need chr, start and end of the peak.
        if (exclude_stream) {
            string line;
            while (!feof(exclude_stream)) {
                if (fgets(exclude_buffer, 128, exclude_stream) != NULL) {
                    line = exclude_buffer;
                    size_t pos = 0;
                    int cnt = 0;
                    while (cnt != 3) {  // Find the position of the third tab.
                        pos++;
                        pos = line.find("\t", pos);
                        if (pos == string::npos)
                            break;
                        cnt++;
                    }
                    // If the chr-prefix is not the same for peak and exclude-file, a warning will be printed but the program
                    // continues.
                    if (line.substr(0, 13) == "***** WARNING") {
                        cout << line << "\n" << "excluding regions not possible" << endl;
                    }
                    excluded_peaks.insert(line.substr(0, pos));  // Avoids using any ID, just using the coordinates.
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
    string bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + b_peakfile + " -wo 2>&1";
    cout << "INtersecting genes and peaks" << endl;
    // Process the whole bed_output, by mapping the peak_ids to the respective genes and also map the location and
    // activity of the peaks in a separate map, and gather the genes for each chromosome.
    unordered_map <string, set<string>> gene_peak_map; // Storing the names of peaks in the window.
    gene_peak_map.reserve(promoter_map.size());
    unordered_map <string, peak_information> peak_info_map;  // Storing an object for each peak, indexed by peakID.
    peak_info_map.reserve(100000);
    unordered_set <string> peak_chromosomes; // Store on which chromosomes the peaks are.
    // Using popen to iterate over the rows and fill both maps for genes and peaks.
    int bed_buffer_size = peak_file_rowlen * 2 + 100;  // +100 for the genes.
    char bed_buffer[bed_buffer_size];
    FILE *bed_intersect_stream = popen(bed_intersect_cmd.c_str(), "r");
    if (bed_intersect_stream) {
        string line;
        while (!feof(bed_intersect_stream)) {
            if (fgets(bed_buffer, bed_buffer_size, bed_intersect_stream) != NULL) {
                line = bed_buffer;
                vector <string> columns = SplitTabLine(line);
                if (excluded_peaks.find(columns[4] + "\t" + columns[5] + "\t" + columns[6]) != excluded_peaks.end()) {
                    continue;
                }
                string this_chr = columns[0];
                if (has_chr_prefix) {  // Already checked that when we peeked at the first line in the peak file.
                    this_chr = this_chr.substr(3);
                }
                peak_chromosomes.insert(this_chr);
                string intersect_gene = columns[3];
                int peak_start = stoi(columns[5]);
                int peak_end = stoi(columns[6]);
                string peak_id_string = this_chr + ":" + to_string(peak_start) + "-" + to_string(peak_end);
                peak_information this_peak;
                this_peak.start = peak_start;
                this_peak.end = peak_end;
                vector<double> activities;
                for (int a_col = start_col; a_col <= last_col; a_col++) {
                    activities.push_back(stof(columns[a_col + 4]));
                }
                this_peak.signal = activities;
                peak_info_map[peak_id_string] = this_peak;
                gene_peak_map[intersect_gene].insert(peak_id_string);
            }
        }
        pclose(bed_intersect_stream);
    } else {
        cout << "ERROR Could not intersect gene windows and peaks, command was the following:\n"
             << bed_intersect_cmd << endl;
        return 1;
    }
    GetStdoutFromCommand("rm " + temp_window_file);
    vector<string> chromosomes;  // New datastructure for better sorting and iteration.
    for (string chr : peak_chromosomes) {
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
    vector<ofstream> out_streams;
    vector<ofstream> gene_info_streams;
    for (int a_col = 0; a_col <= col_num; a_col++) {
        string column_suffix = colname_map[a_col + start_col];
        // Open the output and write the header already. Only the first string has to be explicitly casted.
        ofstream abc_output(d_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        Test_outfile(abc_output, d_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        abc_output << string("#chr") + "\t" + "Peak_Start" + "\t" + "Peak_End" + "\t" + "Ensembl ID" + "\t" + "Gene Name" +
                      "\t" + "PeakID" + "\t" + "signalValue" + "\t" + "Contact" + "\t" + activity_header + "\t" +
                      "scaledContact" + "\t" + "intergenicScore" + "\t" + "TSS-dist" + "\t" + "ABC-Score" + "\n";
        out_streams.push_back(move(abc_output));
        // Also already write the GeneInfo header.
        ofstream gene_info_out(d_prefix + "_GeneInfo" + column_suffix + ".txt");
        gene_info_out << string("Ensembl ID") + "\t" + "Gene Name" + "\t" + "chr" + "\t" + "TSS" + "\t" + "#Enhancer"
                         + "\t" + "Avg_EnhancerActivity" + "\t" + "Avg_EnhancerContact" + "\t" + "Avg_EnhancerDistance" + "\t" + "Failure" + "\n";
        gene_info_streams.push_back(move(gene_info_out));
    }

    // ____________________________________________________________
    // ITERATE THROUGH THE CHROMOSOMES
    // ____________________________________________________________
    // First see which contact files are available.
    string ls_out = GetStdoutFromCommand("ls " + cf_contactfolder + "/");
    vector <string> contact_files;
    istringstream read_ls_out(ls_out);
    if (!read_ls_out) {
        cout << "ERROR Could not access the folder with the contact files, command was the following:\nls " + cf_contactfolder + "/"
             << endl;
        return 1;
    }
    for (string line; getline(read_ls_out, line);) {
        contact_files.push_back(line);
    }

    // Prepare a structure to store genes that are discarded during processing.
    unordered_set <string> genes_wo_hic;  // Missing normalized contact file
    unordered_set <string> genes_wo_candidates;  // No candidate regions in gene window (ABC)
    unordered_set <string> written_genes;  // To track the genes we already wrote to the GeneInfo.
    written_genes.reserve(promoter_map.size());
    // Iterate over the chromosomes, read the normalized matrix in, and iterate over the genes on the chr to calculate
    // the ABC-score.
    cout << "Starting Chromosome-wise soring" << endl;
#pragma omp parallel for num_threads(cores)
    for (unsigned int chr_idx = 0; chr_idx < chromosomes.size(); chr_idx++) {
        string chr = chromosomes[chr_idx];
        if (chr_gene_map[chr].size() == 0) {  // In case we don't have any genes on that chromosome.
            continue;
        }
        // Find the matching contact file.
        string chr_file;
        regex re_chr(".*(C|c)hr" + chr + "[^0-9]{1}.*gz");
        for (string file : contact_files) {
            if (regex_match(file, re_chr)) {
                chr_file = file;
                cout << chr << " mapped to " << chr_file << endl;
                continue;
            }
        }
        auto start_hic = chrono::high_resolution_clock::now();
        // Check in the first line of the file if there is an integer present. If the zcat returns an error, it will be
        // a character at the first position as well.
        string peek_cmd = string("zcat < ") + cf_contactfolder + "/" + chr_file + " 2>&1 | head -n 1";
        FILE *peek_stream = popen(peek_cmd.c_str(), "r");
        char dummy[256];
        regex int_match("[0-9]");
        string temp_row = "foo";
        if (fgets(dummy, 256, peek_stream) != NULL) {
            temp_row = dummy;
        }
        pclose(peek_stream);
        if (!regex_match(temp_row.substr(0, 1), int_match) or chr_file.size() == 0) {
            for (string gene : chr_gene_map[chr]) {
                genes_wo_hic.insert(gene);
            }
            cout << chr << " no contact matrix" << endl;
            continue;
        }

        string hic_cmd = string("zcat < ") + cf_contactfolder + "/" + chr_file + " 2>&1";
        FILE *hic_stream = popen(hic_cmd.c_str(), "r");

        // ____________________________________________________________
        // READ THE HI-C FILE
        // ____________________________________________________________
        // Find the highest bin and construct a max_bin*max_bin matrix, then filling the values into the upper
        // triangle.
        auto start_read = chrono::high_resolution_clock::now();
        int min_bin = hic_boundaries[chr][0] / bin_size;
        int max_bin = hic_boundaries[chr][1] / bin_size;
        boost::numeric::ublas::mapped_matrix<double> contact_matrix(max_bin + 1 - min_bin,
                                                                    max_bin + 1 - min_bin);  // +1 To include max_bin.

        // Track the contacts and the number of non-zero for up to 1MB.
        vector<double> contact_sums(1000000 / bin_size);  // Store the summed contact per distance for the mean.
        char hic_buffer[256];
        string row;
        string delimiter = "\t";
        while (!feof(hic_stream)) {
            if (fgets(hic_buffer, 256, hic_stream) != NULL) {
                row = hic_buffer;
                size_t first_tab = row.find(delimiter);
                size_t second_tab = row.find(delimiter, first_tab + 1);
                double contact = stof(row.substr(second_tab + 1));
                if (!isnan(contact) and contact > 0) {  // 0 Shouldn't appear if preprocessing done correctly.
                    int first_bin = stoi(row.substr(0, first_tab)) / bin_size;
                    int second_bin = stoi(row.substr(first_tab + 1, row.find(delimiter, second_tab))) / bin_size;
                    int bin1 = first_bin;
                    int bin2 = second_bin;
                    if (first_bin > second_bin) {
                        bin1 = second_bin;
                        bin2 = first_bin;
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
        pclose(hic_stream);
        // Fit a linear function to the mean contact at each distance within 1MB.
        double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, regression_count = 0;
        for (int c = 0; c < contact_sums.size(); c++) {
            if (contact_sums[c] > 0) {
                regression_count++;
                int possible_hits = (max_bin - min_bin) - (c + 1);
                double log_distance = log((c + 1) * bin_size);  // Is base-e logarithm. Reverse is exp().
                double log_contact = log(contact_sums[c] / possible_hits);
                sumX += log_distance;
                sumX2 += log_distance * log_distance;
                sumY += log_contact;
                sumXY += log_distance * log_contact;
            }
        }

        double slope = (regression_count * sumXY - sumX * sumY) / (regression_count * sumX2 - sumX * sumX);
        double intercept = (sumY - slope * sumX) / regression_count;
        double close_pseudocount;
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
        contact_matrix(max_bin - min_bin, max_bin - min_bin) = contact_matrix(max_bin - 1 - min_bin, max_bin - min_bin);

        auto stop_hic = chrono::high_resolution_clock::now();
        auto duration_hic = chrono::duration_cast<chrono::seconds>(stop_hic - start_hic);
        cout << duration_hic.count() << "s; " << chr << " contact matrix read" << endl;


        // ____________________________________________________________
        // GET PEAK-GENE CONTACTS AND SUM THE CONTACT PER ENHANCER
        // ____________________________________________________________
        unordered_map <string, vector<double>> gene_contact_map; // Same order as peak_map, but lists the contacts.
        gene_contact_map.reserve(gene_peak_map.size());
        double chr_max = 0;
        for (string chr_gene : chr_gene_map[chr]) {
            int gene_tss = stoi(promoter_map[chr_gene][3]);
            int gene_bin = gene_tss / bin_size;
            vector<double> these_contacts;
            for (auto peak = gene_peak_map[chr_gene].begin(); peak != gene_peak_map[chr_gene].end();) {
                peak_information& matching_peak = peak_info_map[*peak];
                double this_contact;
                int peak_bin = ((matching_peak.end + matching_peak.start) / 2) / bin_size;
                int distance = min(gene_tss - matching_peak.start, gene_tss - matching_peak.end);

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
                matching_peak.contact_sum += this_contact;
                if (distance <= gene_windowsize) {
                    if (this_contact > chr_max) {
                        chr_max = this_contact;
                    }
                    these_contacts.push_back(this_contact);
                    ++peak;
                } else {  // Remove the peaks that are only in the enh-window, we only needed them for the contact_sum.
                    gene_peak_map[chr_gene].erase(*peak);
                }
            }
            gene_contact_map[chr_gene] = these_contacts;
        }

        // Initialize an empty 2D vector to gather the interactions > threshold for each chromosomes.
        // Contact|scaledContact| + activity_cols x (signalValue|scaledActivity|ABC-Score)
        vector <vector<double>> score_fetcher_floats;
        vector <string> interaction_tracker;
        
        // ____________________________________________________________
        // SCORING FOR EACH GENE ON THE CHROMOSOME
        // ____________________________________________________________
        // Now that we mapped the peaks to the genes and fetched the respective information on the peaks, we can do
        // the gene-wise scoring.
        for (string chr_gene : chr_gene_map[chr]) {
            int num_candidates = gene_peak_map[chr_gene].size();
            if (num_candidates == 0) {  // Can also happen if all were overlapping with excluded regions.
                genes_wo_candidates.insert(chr_gene);
                continue;
            }
            // Contact | scaledContact | + activity_cols x (SignalValue | scaled/adjustedActivity | ABC-Score)
            vector <vector<double>> candidate_scores(num_candidates,
                                                     vector<double>(2 + 3 * (col_num + 1)));
            vector <string> candidate_interactions(num_candidates);
            vector<double> max_signals(col_num + 1);
            double max_contact = 0;
            int peak_count_helper = 0;
            vector<double> abc_sums(col_num + 1);
            for (string peak_id : gene_peak_map[chr_gene]) {
                auto matching_peak = peak_info_map[peak_id];  // TODO reference not copy
                double current_contact = gene_contact_map[chr_gene][peak_count_helper];
                if (current_contact > max_contact) {
                    max_contact = current_contact;
                }
                for (int a_col = 0; a_col <= col_num; a_col++) {
                    double current_activity = matching_peak.signal[a_col];
                    candidate_scores[peak_count_helper][2 + a_col * 3] = current_activity;
                    if (do_adjusted_abc) {
                        current_activity = current_activity * (current_contact / matching_peak.contact_sum);
                    } else {  // Only need to search for the maximum if we don't adjust the activity.
                        if (current_activity > max_signals[a_col]) {
                            max_signals[a_col] = current_activity;
                        }
                    }
                    candidate_scores[peak_count_helper][3 + a_col * 3] = current_activity;
                    abc_sums[a_col] += current_activity * current_contact;
                }
                candidate_scores[peak_count_helper][0] = current_contact;
                // Add the interaction so we can later trace back which peak and gene it was.
                candidate_interactions[peak_count_helper] = chr_gene + "*" + peak_id;
                peak_count_helper++;
            }
            // Rescaling the contact to a maximum of 100 only has to be done once. For the activity as many times as
            // there are activity columns. If none of the peaks has contact with the gene, the contact_rescaler and
            // thus the scores will be set to 0.
            double contact_rescaler = (max_contact != 0) ? 100 / max_contact : 0;
            // Write the relative ABC-score. And activity if not adjusting it.
            for (int r = 0; r < peak_count_helper; ++r) {
                bool exceeds = false;
                candidate_scores[r][1] = candidate_scores[r][0] * contact_rescaler;  // Rescale the contact.
                for (int a_col = 0; a_col <= col_num; a_col++) {
                    double this_sum = abc_sums[a_col];
                    if (this_sum != 0) {
                        double this_score = (candidate_scores[r][0] * candidate_scores[r][3 + a_col * 3]) / this_sum;
                        candidate_scores[r][4 + a_col * 3] = this_score;  // Need to store as the sum is gene-specific.
                        if (this_score >= abc_cutoff) {
                            exceeds = true;  // If ≥ threshold with at least one activity column.
                        }
                    }
                    if (not do_adjusted_abc) {
                        candidate_scores[r][3 + a_col * 3] =
                                candidate_scores[r][2 + a_col * 3] * (100 / max_signals[a_col]);
                    }
                }
                if (exceeds) {
                    score_fetcher_floats.push_back(candidate_scores[r]);
                    interaction_tracker.push_back(candidate_interactions[r]);
                }
            }
        }

        auto start_out = chrono::high_resolution_clock::now();

        if (interaction_tracker.size() > 0) {
            interaction_tracker.push_back("Placeholder");  // We always peak one gene ahead.
#pragma omp critical(out_writer)  // Allows only one thread to write to output files at the same time.
            {
                // Rows: activity columns and columns: #Enhancer|EnhancerActivity|EnhancerContact|EnhancerDistance.
                vector<vector<double>> gene_info(col_num + 1, vector<double>(4));

                // Repeat the process for each activity column.
                for (int i = 0; i < score_fetcher_floats.size(); ++i) {
                    vector<double> row_floats = score_fetcher_floats[i];
                    string interaction = interaction_tracker[i];
                    // For the gene info keep track of the current gene, we filled the data ordered by genes.
                    string gene_id = interaction.substr(0, interaction.find("*"));
                    // Get the gene ID, gene name and chromosome out of the promoter_map.
                    auto gene_map = promoter_map[gene_id];

                    auto matching_peak = peak_info_map[interaction.substr(interaction.find("*") + 1)];
                    for (int a_col = 0; a_col <= col_num; a_col++) {
                        if (row_floats[4 + a_col * 3] >= abc_cutoff) {  // Check again, it's present if only one exceeds.
                            out_streams[a_col] << gene_map[2] + "\t" + to_string(matching_peak.start) + "\t" +
                                                  to_string(matching_peak.end);
                            out_streams[a_col] << "\t" + gene_map[0] + "\t" + gene_map[1];
                            out_streams[a_col] << "\t" + interaction.substr(interaction.find("*") + 1);
                            out_streams[a_col] << "\t" + to_string(row_floats[2 + a_col * 3]);  // signalValue
                            out_streams[a_col] << "\t" + to_string(row_floats[0]);  // Contact
                            out_streams[a_col] << "\t" + to_string(row_floats[3 + a_col * 3]);  // scaledActivity
                            out_streams[a_col] << "\t" + to_string(row_floats[1]);  // scaledContact
                            // chr-scaled Contact * signalValue OR adjustedActivity
                            if (do_adjusted_abc) {
                                out_streams[a_col] << "\t" + to_string(row_floats[intergenic_activity_col + a_col * 3]);
                            }
                            else {
                                out_streams[a_col] << "\t" + to_string(row_floats[0] / chr_max * row_floats[intergenic_activity_col + a_col * 3]);
                            }

                            int distance = min(abs(stoi(gene_map[3]) - matching_peak.start),
                                               abs(stoi(gene_map[3]) - matching_peak.end));

                            out_streams[a_col] << "\t" + to_string(distance);  // Distance
                            out_streams[a_col] << "\t" + to_string(row_floats[4 + a_col * 3])+ "\n";  // ABC-Score
                            // Sum up the attributes of the interaction.
                            gene_info[a_col][0] += 1;
                            gene_info[a_col][1] += row_floats[2 + a_col * 3];
                            gene_info[a_col][2] += row_floats[0];
                            gene_info[a_col][3] += static_cast<double>(distance);
                        }
                    }
                    // If there is a new gene next write to GeneInfo and reset.
                    if (interaction_tracker[i+1].substr(0, interaction_tracker[i+1].find("*")) != gene_id) {
                        written_genes.insert(gene_id);
                        for (int a_col = 0; a_col <= col_num; a_col++) {
                            gene_info_streams[a_col] << gene_id + "\t" + gene_map[1] + "\t" + gene_map[2] + "\t" + gene_map[3];
                            int num_enhancer = static_cast<int>(gene_info[a_col][0]);
                            if (num_enhancer > 0) {
                                gene_info_streams[a_col] << "\t" + to_string(num_enhancer) + "\t" + to_string(gene_info[a_col][1]/num_enhancer) +
                                                            "\t" + to_string(gene_info[a_col][2]/num_enhancer) + "\t" + to_string(gene_info[a_col][3]/num_enhancer) + "\t-\n";
                            }
                            else {
                                gene_info_streams[a_col] << "\t0\t0\t0\t0\tScore of candidates too low\n";
                            }
                            fill(gene_info[a_col].begin(), gene_info[a_col].end(), 0);  // Reset the vector.
                        }
                    }
                }
            }
        }
    }
    // ____________________________________________________________
    // COMPLETE GENEINFO OUTPUT FILES
    // ____________________________________________________________
    for (auto entry: promoter_map) {
        string gene_id = entry.first;
        auto gene_map = entry.second;
        string fail_message;
        if (genes_wo_hic.find(gene_id) != genes_wo_hic.end()) {
            fail_message = "Missing normalized contact file";
        }
        else if (genes_wo_candidates.find(gene_id) != genes_wo_candidates.end()) {
            fail_message = "No candidate regions in gene window";
        }
        else if (written_genes.find(gene_id) == written_genes.end()) {
            if (peak_chromosomes.find(gene_map[2]) == peak_chromosomes.find(gene_map[2])) {
                fail_message = "No candidate regions in gene window";
            }
            else {
                fail_message = "Score of candidate regions too low";
            }
        }
        if (fail_message.size() > 0) {
            for (int a_col = 0; a_col <= col_num; a_col++) {
                gene_info_streams[a_col] << gene_id + "\t" + gene_map[1] + "\t" + gene_map[2] + "\t" + gene_map[3] +
                                            "\t0\t0\t0\t0\t" + fail_message + "\n";
            }
        }
    }
    // Gzip all output files.
    cout << "Compressing files" << endl;
#pragma omp parallel for num_threads(cores)
    for (int a_col = 0; a_col <= col_num; a_col++) {  // Iterate through a_cols to have the correct suffix.
        string column_suffix = colname_map[a_col + start_col];
        out_streams[a_col].close();
        GzipFile(d_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        gene_info_streams[a_col].close();
        GzipFile(d_prefix + "_GeneInfo" + column_suffix + ".txt");
    }
    auto stop0 = chrono::high_resolution_clock::now();
    auto duration0 = chrono::duration_cast<chrono::seconds>(stop0 - start0);
    cout << duration0.count() << "s ABCpp in total" << endl;
}

