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

//class constants {
//public:
//    static const int precision = 10;
//};

int GetDistance(int start, int end, int other);

std::string SetOptionalInput(const std::string& input_string, std::string default_value);

bool SetBoolInput(const std::string& input_string, bool default_value);

bool Chr_sorter(const std::string& first, const std::string& second);

int FilePeek(const std::string& file_name);

std::vector <std::string> SplitTabLine(const std::string& row);

std::string GetStdoutFromCommand(std::string cmd);

void Test_outfile(std::ofstream &out_file, const std::string& out_path);

void GzipFile(const std::string& file);

std::unordered_map<int, std::string> FileHeaderMapping(const std::string& file_name, const std::vector<int>&);

#endif //STARE_MISCFUNCTIONS_H
