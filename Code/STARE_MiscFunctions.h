//
// Created by Dennis Hecker on 15.07.21.
//

#ifndef STARE_GAZE_STARE_MISCFUNCTIONS_H
#define STARE_GAZE_STARE_MISCFUNCTIONS_H

# include <string>
# include <iostream>
# include <string>
# include <vector>
# include <fstream>
# include <sstream>
# include <unordered_map>

std::string to_string_with_precision(double a_value, const int n = 10);

std::string SetOptionalInput(std::string input_string, std::string default_value);

bool SetBoolInput(std::string input_string, bool default_value);

std::vector <std::string> SplitTabLine(std::string row);

std::string GetStdoutFromCommand(std::string cmd);

void GzipFile(std::string file);

std::unordered_map<int, std::string> FileHeaderMapping(std::string file_name, int first_col, int last_col);

#endif //STARE_GAZE_STARE_MISCFUNCTIONS_H
