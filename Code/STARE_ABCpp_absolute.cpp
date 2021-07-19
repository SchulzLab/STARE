//
// Created by Dennis Hecker on 01.12.20.
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
 * ABC-score. So when examining the affinities of TFs to genes, we can use these called interactions instead of all
 * the open regions in a window. The scaling of TF affinity can then also be done with the ABC-score itself.
 * This program also works independently of STARE.
 *
 * For MacOS:
 * g++ STARE_ABCpp_absolute.cpp STARE_MiscFunctions.cpp -Xpreprocessor -fopenmp -std=c++11 -I path-to-boost -o STARE_ABCpp_absolute
 * For Linux:
 * g++ STARE_ABCpp_absolute.cpp STARE_MiscFunctions.cpp -fopenmp -std=c++11 -o STARE_ABCpp_absolute
 * ./STARE_ABCpp_absolute -b ../Test/ABC_example_regions.bed -n 4 -a ../Test/ABC_example_annotation.gtf -gw 5000000 -bin 5000 -cf path-to-normalized-HiC-contacts -c 4 -t ABC-cut-off -d output-folder
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 * C++ implementation of the ABC-scoring principle from Fulco et al.: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
 */

class peak_information {
    public:
        int start;
        int end;
        std::vector<double> norm_signal;
};


int main(int argc, char **argv) {
    using namespace std;
    auto start0 = chrono::high_resolution_clock::now();

    // ____________________________________________________________
    // FETCH AND CHECK INPUT ARGS
    // ____________________________________________________________
    string parameter_help = "-b enhancer/peak-file\n-n activity-column(s), start counting at 1\n-a gtf gene annotation\n-gw genewindow-size"
                            "\n-cf folder with the normalized hic-files per chromosome\n-bin binsize of the hic-files"
                            "\n-t cut-off for the ABC-score (default 0.02), set to 0 to get all scored interactions\n-d prefix with which the files are written"
                            "\n-c number of cores used (default 1)\n-x file with regions to be excluded";

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
            c_cores, d_prefix, x_exclude_regions;
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
                    {NULL,  0,                 NULL, 0}
            };

    int arg;
    while ((arg = getopt_long_only(argc, argv, "b:a:n:w:i:t:c:f:d:x:", long_options, NULL)) != -1) {
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
    int window_size = stoi(SetOptionalInput(gw_genewindowsize, "5000000")) / 2;
    int bin_size = stoi(SetOptionalInput(bin_binsize, "5000"));
    string activity_col = SetOptionalInput(n_activitycol, "4");
    int cores = stoi(SetOptionalInput(c_cores, "1"));
    double abc_cutoff = stof(SetOptionalInput(t_abc_cutoff, "0.02"));

    // ____________________________________________________________
    // PROCESS GTF GENE ANNOTATION
    // ____________________________________________________________
    vector <string> chromosomes;  // To store the unique chromosomes in and still be able to check their index.
    map <string, unordered_set<string>> chr_gene_map;  // For each chromosome stores its genes.
    unordered_map <string, vector<string>> promoter_map;  // For ID, name, chr, TSS.
    unordered_map <string, vector<int>> hic_boundaries;  // Store the min/max boundaries for each chromosome.

    string row;
    ifstream Read_Gene_Annotation(a_promoterfile);
    int row_count = -1;
    while (!Read_Gene_Annotation.eof()) {
        getline(Read_Gene_Annotation, row);
        if ((row.size() > 1) and (row.rfind("#", 0) != 0)) {
            vector <string> columns = SplitTabLine(row);
            if (columns[2] == "gene") {
                string id_delimiter = "gene_id ";
                int id_start = columns[8].find(id_delimiter, 0) + id_delimiter.size();
                int id_end = columns[8].find(";", 0);
                string gene_id = columns[8].substr(id_start + 1, id_end - id_start - 2);  // Remove the quotation marks.
                string name_delimiter = "gene_name ";
                int name_start = columns[8].find(name_delimiter, 0) + name_delimiter.size();
                int name_end = columns[8].find(";", name_start);
                string gene_name = columns[8].substr(name_start + 1, name_end - name_start -
                                                                     2);  // Remove the quotation marks.                string chr = columns[0];
                string strand = columns[6];
                string chr = columns[0];
                if (chr.substr(0, 3) == "chr") {
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
                if (find(chromosomes.begin(), chromosomes.end(), chr) != chromosomes.end() or
                    chr.size() > 5) {
                } else {  // Add the chromosome to the list to iterate over.
                    chromosomes.push_back(chr);
                    hic_boundaries[chr] = {0, 0};
                }
                chr_gene_map[chr].insert(gene_id);
            }
        }
        row_count++;
    }
    Read_Gene_Annotation.close();

    // ____________________________________________________________
    // WRITE TEMPORARY GENE WINDOW FILE
    // ____________________________________________________________
    // First check if the peak file has a chr-prefix. This can happen if the executable is called independent of the
    // TEPIC bash script which pre-processes the peak file.
    // Also directly check for the signal columns in case n++ was given.
    ifstream peek_peak(b_peakfile);
    string line_peek;
    string peak_chr_prefix = "";
    int start_col, last_col;
    while (!peek_peak.eof()) {
        getline(peek_peak, line_peek);
        if (line_peek.substr(0, 1) != "#") {
            if (line_peek.substr(0, 3) == "chr") {
                peak_chr_prefix = "chr";
            }
            int total_cols = SplitTabLine(line_peek).size();
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

    unordered_map<int, string> colname_map = FileHeaderMapping(b_peakfile, start_col, last_col);

    // Is not written directly on gtf-file read, as we need to find the most 5'-TSS first.
    string temp_window_file = d_prefix + "_ABCpp_Temp_GeneWindow.bed";
    ofstream window_out(temp_window_file);
    for (auto const &gene : promoter_map) {  // [geneID]: ID, name, chr, TSS.
        int window_start = stoi(gene.second[3]) - window_size;
        int window_end = stoi(gene.second[3]) + window_size;
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
    string exclude_regions_intersect = "bedtools intersect -a " + b_peakfile + " -b " + x_exclude_regions + " -u";
    string exclude_intersect_out = GetStdoutFromCommand(exclude_regions_intersect);
    unordered_set <string> excluded_peaks;
    istringstream excluded_line_reader(exclude_intersect_out);
    for (string line; getline(excluded_line_reader, line);) {
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

    // ____________________________________________________________
    // INTERSECT GENE WINDOWS AND PEAKS
    // ____________________________________________________________
    // Call the bedtools intersect command via popen.
    string bed_intersect_cmd = "bedtools intersect -a " + temp_window_file + " -b " + b_peakfile + " -wo";
    string bed_intersect_out = GetStdoutFromCommand(bed_intersect_cmd);
    GetStdoutFromCommand("rm " + temp_window_file);

    // Process the whole bed_output, by mapping the peak_ids to the respective genes and also map the location and
    // activity of the peaks in a separate map, and gather the genes for each chromosome.
    unordered_map <string, vector<string>> gene_peak_map; // Storing the names of peaks in its window.
    unordered_map <string, peak_information> peak_info_map;  // Storing an object for each peak, indexed by peakID.
    // Using istringstream to iterate over the rows and fill both maps for genes and peaks.
    istringstream bed_line_reader(bed_intersect_out);
    for (string line; getline(bed_line_reader, line);) {
        vector <string> columns = SplitTabLine(line);
        if (excluded_peaks.find(columns[4] + "\t" + columns[5] + "\t" + columns[6]) != excluded_peaks.end()) {
            continue;
        }
        string this_chr = columns[0];
        if (this_chr.substr(0, 3) == "chr") {
            this_chr = this_chr.substr(3);
        }
        string intersect_gene = columns[3];
        int peak_start = stoi(columns[5]);
        int peak_end = stoi(columns[6]);
        string peak_id_string = this_chr + ":" + to_string(peak_start) + "-" + to_string(peak_end);
        peak_information this_peak;
        this_peak.start = peak_start;
        this_peak.end = peak_end;
        vector<double> activities;
        for (int a_col = start_col; a_col <= last_col; a_col++) {
            activities.push_back(stof(columns[a_col + 4]) / (peak_end - peak_start));
        }
        this_peak.norm_signal = activities;
        peak_info_map[peak_id_string] = this_peak;
        gene_peak_map[intersect_gene].push_back(peak_id_string);
    }
    auto stop_intersect = chrono::high_resolution_clock::now();
    auto duration_intersect = chrono::duration_cast<chrono::seconds>(stop_intersect - start0);
    cout << duration_intersect.count() << "s ABCpp: intersected genes and peaks" << endl;

    // ____________________________________________________________
    // PREPARE OUTPUT FILES
    // ____________________________________________________________
    // Open output files already.
    vector<ofstream> out_streams;
    for (int a_col = 0; a_col <= last_col-start_col; a_col++) {
        string column_suffix = colname_map[a_col + start_col];
        // Open the output and write the header already. Only the first string has to be explicitly casted.
        ofstream abc_output(d_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
        abc_output << string("#chr") + "\t" + "Peak_Start" + "\t" + "Peak_End" + "\t" + "Ensembl ID" + "\t" + "Gene Name" +
                   "\t" + "PeakID" + "\t" + "norm_signalValue" + "\t" + "Contact" + "\t" + "scaledActivity" + "\t" +
                   "scaledContact" + "\t" + "ABC-Score" + "\t" + "TSS-dist" + "\n";
        out_streams.push_back(move(abc_output));
    }

    // ____________________________________________________________
    // ITERATE THROUGH THE CHROMOSOMES
    // ____________________________________________________________
    // First see which contact files are available.
    string ls_out = GetStdoutFromCommand("ls " + cf_contactfolder + "/");
    vector <string> contact_files;
    istringstream read_ls_out(ls_out);
    for (string line; getline(read_ls_out, line);) {
        contact_files.push_back(line);
    }
    // Prepare a structure to store genes that are discarded during processing.
    unordered_map <string, string> kicked_genes;
    unordered_set <string> genes_with_scores;  // To track the genes with interactions in the output.
    // Iterate over the chromosomes, read the normalized matrix in, and iterate over the genes on the chr to calculate
    // the ABC-score.
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
                cout << chr << " " << chr_file << endl;
                continue;
            }
        }
        auto start_hic = chrono::high_resolution_clock::now();
        stringstream ReadContactMatrix(GetStdoutFromCommand(string("zcat < ") + cf_contactfolder + "/" + chr_file));
        // If file does not exists or when it is empty.
        if (!ReadContactMatrix or ReadContactMatrix.peek() == ifstream::traits_type::eof() or chr_file.size() == 0) {
            for (string gene : chr_gene_map[chr]) {
                kicked_genes[gene] = "Missing normalized contact file";
            }
            cout << chr << " no contact matrix" << endl;
            continue;
        }
        // Initialize an empty 2D vector to gather the interactions > threshold for each chromosomes.
        // Contact|scaledContact| + activity_cols x (norm_SignalValue|scaledActivity|ABC-Score)
        vector <vector<double>> score_fetcher_floats;
        vector <string> interaction_tracker;

        // ____________________________________________________________
        // READ THE HI-C FILE
        // ____________________________________________________________
        // Find the highest bin and construct a max_bin*max_bin matrix, then filling the values into the upper
        // triangle.
        auto start_read = chrono::high_resolution_clock::now();
        string row;
        int min_bin = hic_boundaries[chr][0] / bin_size;
        int max_bin = hic_boundaries[chr][1] / bin_size;
        boost::numeric::ublas::mapped_matrix<double> contact_matrix(max_bin + 1 - min_bin,
                                                                    max_bin + 1 - min_bin);  // +1 To include max_bin.

        // Do a linear regression on ln(distance) ln(contact) to add pseudocounts to the interactions later.
        double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, regression_counter = 0;

        while (getline(ReadContactMatrix, row)) {
            string delimiter = "\t";
            size_t first_tab = row.find(delimiter);
            size_t second_tab = row.find(delimiter, first_tab + 1);
            double contact = stof(row.substr(second_tab + 1));
            if (!isnan(contact)) {
                int first_bin = stoi(row.substr(0, first_tab)) / bin_size;
                int second_bin = stoi(row.substr(first_tab + 1, row.find(delimiter, second_tab))) / bin_size;
                int bin1 = first_bin;
                int bin2 = second_bin;
                if (first_bin > second_bin) {
                    bin1 = second_bin;
                    bin2 = first_bin;
                }
                int distance = (bin2 - bin1) * bin_size;
                if (distance < 1000000 and distance > bin_size) {
                    double log_distance = log(distance);  // Is base-e logarithm. Reverse is exp().
                    double log_contact = log(contact);
                    regression_counter++;
                    sumX = sumX + log_distance;
                    sumX2 = sumX2 + log_distance * log_distance;
                    sumY = sumY + log_contact;
                    sumXY = sumXY + log_distance * log_contact;
                }
                // The gene windows may not cover all available Hi-C bins.
                if (bin1 >= min_bin and bin2 <= max_bin) {
                    contact_matrix(bin1 - min_bin, bin2 - min_bin) = contact;
                }
            }
        }
        double slope = (regression_counter * sumXY - sumX * sumY) / (regression_counter * sumX2 - sumX * sumX);
        double intercept = (sumY - slope * sumX) / regression_counter;
        double close_pseudocount = exp(slope * log(1000000) + intercept);  // For all distance <1 Mb.
        if (isnan(slope)) {  // In case a fit is not possible.
            slope = 0;
            intercept = 0;
            close_pseudocount = 0;
            cout << "Fit for pseudocount not possible, proceeding without pseudocounts." << endl;
        }
        auto duration_read = chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - start_read);
        cout << duration_read.count() << "ms " << chr << ", stringstream read and processed" << endl;

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
        auto duration_hic = chrono::duration_cast<chrono::milliseconds>(stop_hic - start_hic);
        cout << duration_hic.count() << "s " << chr << ", contact matrix read" << endl;

        // ____________________________________________________________
        // SCORING FOR EACH GENE ON THE CHROMOSOME
        // ____________________________________________________________
        // Now that we mapped the peaks to the genes and fetched the respective information on the peaks, we can do
        // the gene-wise scoring.
        for (string chr_gene : chr_gene_map[chr]) {
            int num_candidates = gene_peak_map[chr_gene].size();
            int gene_bin = stoi(promoter_map[chr_gene][3]) / bin_size;
            double self_contact = contact_matrix(gene_bin - min_bin, gene_bin - min_bin);
            if (self_contact == 0) {
                kicked_genes[chr_gene] = "No Hi-C contact at gene bin";
                continue;
            }
            if (num_candidates == 0) {  // Can also happen if all were overlapping with excluded regions.
                kicked_genes[chr_gene] = "No candidate regions in gene window (ABC)";
                continue;
            }

            // Contact|scaledContact| + activity_cols x (norm_SignalValue|scaledActivity|ABC-Score)
            vector <vector<double>> candidate_scores(num_candidates,
                                                     vector<double>(2 + 3 * (last_col - start_col + 1)));
            vector <string> candidate_interactions(num_candidates);
            vector<double> max_signals(last_col - start_col + 1);
            double max_contact = 0;
            int peak_count_helper = 0;
            for (string peak_id : gene_peak_map[chr_gene]) {
                auto matching_peak = peak_info_map[peak_id];
                int peak_bin = ((matching_peak.end + matching_peak.start) / 2) / bin_size;
                for (int a_col = 0; a_col <= last_col - start_col; a_col++) {
                    candidate_scores[peak_count_helper][2 + a_col * 3] = matching_peak.norm_signal[a_col];
                    if (matching_peak.norm_signal[a_col] > max_signals[a_col]) {
                        max_signals[a_col] = matching_peak.norm_signal[a_col];
                    }
                }
                double pseudocount = close_pseudocount;
                if (abs(gene_bin - peak_bin) * bin_size > 1000000) {
                    pseudocount = exp(slope * log(abs(gene_bin - peak_bin) * bin_size) + intercept);
                }
                if ((peak_bin > max_bin) or (peak_bin < min_bin)) {
                    candidate_scores[peak_count_helper][0] = pseudocount;
                } else {
                    candidate_scores[peak_count_helper][0] = pseudocount +
                                                             contact_matrix(min(gene_bin, peak_bin) - min_bin,
                                                                            max(gene_bin, peak_bin) - min_bin);
                }
                if (candidate_scores[peak_count_helper][0] > max_contact) {
                    max_contact = candidate_scores[peak_count_helper][0];
                }
                // Add the interaction so we can later trace back which peak and gene it was.
                candidate_interactions[peak_count_helper] = chr_gene + "*" + peak_id;
                peak_count_helper++;
            }
            // Rescaling the contact to a maximum of 100 only has to be done once. For the activity as many times as
            // there are activity columns. If non of the peaks has contact with the gene, the contact_rescaler and
            // thus the scores will be set to 0.
            double contact_rescaler = (max_contact != 0) ? 100 / max_contact : 0;
            // Get the summed ABC_score.
            vector<double> abc_sums(last_col - start_col + 1);
            // After fetching we can rescale contact and activity.
            for (int r = 0; r < peak_count_helper; ++r) {
                double rescaled_contact = candidate_scores[r][0] * contact_rescaler;
                candidate_scores[r][1] = rescaled_contact;
                for (int a_col = 0; a_col <= last_col - start_col; a_col++) {
                    double rescaled_activity = candidate_scores[r][2 + a_col * 3] * (100 / max_signals[a_col]);
                    candidate_scores[r][3 + a_col * 3] = rescaled_activity;
                    abc_sums[a_col] += rescaled_activity * rescaled_contact;
                }
            }

            // Write the relative ABC-score.
            for (int r = 0; r < peak_count_helper; ++r) {
                bool exceeds = false;
                for (int a_col = 0; a_col <= last_col - start_col; a_col++) {
                    double this_sum = abc_sums[a_col];
                    if (this_sum != 0) {
                        double this_score = (candidate_scores[r][1] * candidate_scores[r][3 + a_col * 3]) / this_sum;
                        candidate_scores[r][4 + a_col * 3] = this_score;  // Need to store as the sum is gene-specific.
                        if (this_score >= abc_cutoff) {
                            exceeds = true;  // If ≥ threshold with at least one activity column.
                        }
                    }
                }
                if (exceeds) {
                    score_fetcher_floats.push_back(candidate_scores[r]);
                    interaction_tracker.push_back(candidate_interactions[r]);
                }
            }
        }

#pragma omp critical(out_writer)  // Allows only one thread to write to output files at the same time.
        {
            // Repeat the process for each activity column.
            for (int i = 0; i < score_fetcher_floats.size(); ++i) {
                vector<double> row_floats = score_fetcher_floats[i];
                string interaction = interaction_tracker[i];
                // Get the gene ID, gene name and chromosome out of the promoter_map.
                auto matching_gene = promoter_map[interaction.substr(0, interaction.find("*"))];
                genes_with_scores.insert(matching_gene[0]);
                auto matching_peak = peak_info_map[interaction.substr(interaction.find("*") + 1)];
                for (int a_col = 0; a_col <= last_col - start_col; a_col++) {
                    if (row_floats[4 + a_col * 3] >= abc_cutoff) {  // Check again, it's present if only one exceeds.
                        out_streams[a_col] << matching_gene[2] + "\t" + to_string(matching_peak.start) + "\t" + to_string(matching_peak.end);
                        out_streams[a_col] << "\t" + matching_gene[0] + "\t" + matching_gene[1];
                        out_streams[a_col] << "\t" + interaction.substr(interaction.find("*") + 1);
                        out_streams[a_col] << "\t" + to_string(row_floats[2 + a_col * 3]);  // norm_signalValue
                        out_streams[a_col] << "\t" + to_string(row_floats[0]);  // Contact
                        out_streams[a_col] << "\t" + to_string(row_floats[3 + a_col * 3]);  // scaledActivity
                        out_streams[a_col] << "\t" + to_string(row_floats[1]);  // scaledContact
                        out_streams[a_col] << "\t" + to_string(row_floats[4 + a_col * 3]);  // ABC-Score
                        out_streams[a_col] << "\t" << to_string(min(abs(stoi(matching_gene[3]) - matching_peak.start),
                                                                    abs(stoi(matching_gene[3]) - matching_peak.end))) +
                                                      "\n";
                    }
                }
            }
        }
    }
    // Gzip all output files.
    for (int a_col = 0; a_col <= last_col-start_col; a_col++) {  // Iterate through a_cols to have the correct suffix.
        string column_suffix = colname_map[a_col + start_col];
        out_streams[a_col].close();
        GzipFile(d_prefix + "_ABCpp_scoredInteractions" + column_suffix + ".txt");
    }

    // ____________________________________________________________
    // WRITE DISCARDED GENES OUTPUT
    // ____________________________________________________________
    // Find the last genes that didn't make it to the output and write the file with discarded genes.
    unordered_set<string> already_discarded;
    ofstream discarded_out;
    discarded_out.open(d_prefix + "_discarded_Genes.txt", ios_base::app);
    for (auto const &gene : kicked_genes) {
        already_discarded.insert(gene.first);
        discarded_out << gene.first << "\t" << gene.second << "\n";
    }
    for (auto const &gene : promoter_map) {  // [geneID]: ID, name, chr, TSS.
        if (already_discarded.find(gene.first) == already_discarded.end()) {
            if (genes_with_scores.find(gene.first) == genes_with_scores.end()) {
                discarded_out << gene.first << "\t" << "Score of candidate enhancers too low" << "\n";
            }
        }
    }
    discarded_out.close();
    cout << "genes discarded" << endl;
    auto stop0 = chrono::high_resolution_clock::now();
    auto duration0 = chrono::duration_cast<chrono::seconds>(stop0 - start0);
    cout << duration0.count() << "s ABCpp in total" << endl;
}
