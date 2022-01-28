//
// Created by Dennis Hecker on 15.07.21.
//

#ifndef STARE_MISCFUNCTIONS_H
#define STARE_MISCFUNCTIONS_H

# include <string>
# include <iostream>
# include <string>
# include <vector>
# include <fstream>
# include <sstream>
# include <unordered_map>
# include <algorithm>
# include <regex>

class constants {
public:
    static const int precision = 10;
};

std::string To_string_with_precision(double a_value);

std::string SetOptionalInput(std::string input_string, std::string default_value);

bool SetBoolInput(std::string input_string, bool default_value);

bool Chr_sorter(std::string first, std::string second);

int FilePeek(std::string file_name);

std::vector <std::string> SplitTabLine(std::string row);

std::string GetStdoutFromCommand(std::string cmd);

void Test_outfile(std::ofstream &out_file, std::string out_path);

void GzipFile(std::string file);

std::unordered_map<int, std::string> FileHeaderMapping(std::string file_name, std::vector<int>);

#endif //STARE_MISCFUNCTIONS_H
