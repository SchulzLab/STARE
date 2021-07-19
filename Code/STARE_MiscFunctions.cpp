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

std::string to_string_with_precision(double a_value, const int n)
{
    // to_string of doubles is restricted to 6 digits, this is a workaround.
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

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
        if ((input_string == "True") or (input_string == "TRUE") or (input_string == "true") or (input_string == "1")) {
            return true;
        }
        else if ((input_string == "False") or (input_string == "FALSE") or (input_string == "false") or (input_string == "0")) {
            return false;
        }
        else {
            return default_value;
        }
    }
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
    const int max_buffer = 256;
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

void GzipFile(std::string file) {
    // Call the gzip on a file. If it already exists it is skipped and the warning is printed
    // to stdout. Makes use of the function above for reading in the cout of the command call.
    std::string gzip_out = GetStdoutFromCommand("yes y | gzip " + file);
    if (gzip_out.size() > 1) {
        std::cout << gzip_out << std::endl;
    }
}

std::unordered_map<int, std::string> FileHeaderMapping(std::string file_name, int first_col, int last_col) {
// Reads the file header and creates a map based on the column names. Header has to start with #, otherwise columns
// will be named cN.
    std::unordered_map<int, std::string> header_map;
    std::ifstream header_read(file_name);
    std::string line_peek;
    getline(header_read, line_peek);
    if (line_peek.substr(0, 1) == "#") {
        std::vector<std::string> col_names = SplitTabLine(line_peek);
        for (int c = first_col; c <= last_col; c++) {
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
        for (int c = first_col; c <= last_col; c++) {
            header_map[c] = "_c" + std::to_string(c + 1);  // 1-based counting for the suffixes.
        }
    }
    return header_map;
}


