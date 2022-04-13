//
// Created by Dennis Hecker on 15.07.21.
//


#include "STARE_MiscFunctions.h"

/*
 * Collection of functions that are used by more or less all functions of STARE. Eases string formatting and file
 * handling.
 *
 * Part of STARE: https://github.com/SchulzLab/STARE
 */


std::string SetOptionalInput(std::string input_string, std::string default_value) {
    // If the input_string is empty, return the default value.
    std::string output_val;
    if (input_string.size() > 0) {
        output_val = input_string;
    } else {
        output_val = default_value;
    }
    return output_val;
}

bool SetBoolInput(std::string input_string, bool default_value) {
    // If the input_string is empty, return the default value.
    if (input_string.size() == 0) {
        return default_value;
    }
    else {
        if ((input_string == "True") or (input_string == "TRUE") or (input_string == "true") or (input_string == "T") or (input_string == "1")) {
            return true;
        }
        else if ((input_string == "False") or (input_string == "FALSE") or (input_string == "false") or (input_string == "F") or (input_string == "0")) {
            return false;
        }
        else {
            return default_value;
        }
    }
}

bool Chr_sorter(std::string first, std::string second) {
    // Sorts strings of potential ints and chars with the ints in increasing order and the chars at the end.
    std::regex re_chr("^[0-9]+$");
    if (std::regex_match(first, re_chr) and std::regex_match(second, re_chr)){
        return stoi(first) < stoi(second);
    }
    if (not std::regex_match(first, re_chr)) {
        return false;
    }
    else {
        return true;
    }
}

int FilePeek(std::string file_name) {
    // Looks at the first 5 non-# line in the file and returns its length. Exits when the file can't be opened or
    // when it's empty.
    std::ifstream peek_file(file_name.c_str());
    if (!peek_file) {
        std::cout << "ERROR Could not open the file:\n" << file_name << std::endl;
        return 1;
    }
    std::string line_peek;
    std::vector<int> peek_lengths;
    int row_count = 0;
    while (!peek_file.eof()) {
        getline(peek_file, line_peek);
        if (line_peek.substr(0, 1) != "#") {
            row_count++;
            peek_lengths.push_back(line_peek.size());
        }
        if (row_count >= 5) {
            break;
        }
    }
    int peek_len = *max_element(std::begin(peek_lengths), std::end(peek_lengths));
    if (peek_len == 0) {
        std::cout << "ERROR Empty file:\n" << file_name << std::endl;
        return 1;
    }
    return peek_len;
}

std::vector<std::string> SplitTabLine(std::string row) {
    // Splits a line by tabs and returns a vector of the strings.
    std::string col_val;
    std::vector<std::string> columns;
    std::stringstream row_read(row);
    while (std::getline(row_read, col_val, '\t')) {
        columns.push_back(col_val);
    }
    return columns;
}

std::string GetStdoutFromCommand(std::string cmd) {
// It creates a buffer, opens up a read-only stream, runs the command, and captures the output, stuffs it into the
// buffer, then returns it as a string.
    std::string data;
    FILE * stream;
    const int max_buffer = 128;
    char buffer[max_buffer];
    cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");

    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}

void Test_outfile(std::ofstream &out_file, std::string out_path) {
    // Test if the output path is viable.
    if (!out_file) {
        std::cout << "ERROR Can't open the output path, please check carefully:\n" + out_path << std::endl;
        exit(1);
    }
}

void GzipFile(std::string file) {
    // Call the gzip on a file. If it already exists it is skipped and the warning is printed
    // to stdout. Makes use of the function above for reading in the cout of the command call.
    std::string gzip_out = GetStdoutFromCommand("yes y | gzip " + file);
    if (gzip_out.size() > 1) {
        std::cout << gzip_out << std::endl;
    }
}

std::unordered_map<int, std::string> FileHeaderMapping(std::string file_name, std::vector<int> col_indices) {
// Reads the file header and creates a map based on the column names. Header has to start with #, otherwise columns
// will be named cN.
    std::unordered_map<int, std::string> header_map;
    std::ifstream header_read(file_name);
    std::string line_peek;
    getline(header_read, line_peek);
    if (line_peek.substr(0, 1) == "#") {
        std::vector<std::string> col_names = SplitTabLine(line_peek);
        for (int c : col_indices) {
            std::string non_whitespace = col_names[c];  // To remove any potential # from the header.
            std::replace(non_whitespace.begin(), non_whitespace.end(), ' ', '_');
            if (non_whitespace.size() > 0) {  // We can't work with empty strings.
                header_map[c] = "_" + non_whitespace;
            }
            else {
                header_map[c] = "_c" + std::to_string(c + 1);
            }
        }
    }
    else {
        for (int c : col_indices) {
            header_map[c] = "_c" + std::to_string(c + 1);  // 1-based counting for the suffixes.
        }
    }
    return header_map;
}


