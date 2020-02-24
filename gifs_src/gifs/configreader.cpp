#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//
#include <iostream>
#include <regex>
//
#include "configreader.hpp"

inline
bool
is_str_substring(const char* substring, const char* string) 
{
    return strstr(string, substring) != NULL; 
};


inline
char* 
FileHandle::get_line()
{
    char* line = nullptr;
    size_t len;
    ssize_t read = getline(&line, &len, file.fh()); 
    if (read == -1) {
        return nullptr; 
    }
    return line;
};


inline
char* 
FileHandle::find_line(const std::string& key) {
    const char * ref = key.c_str();
    char * line = nullptr;

    while (true) {
        line = get_line();
        if (line == nullptr)  break; 
        if (is_str_substring(ref, line)) break; 
    }
    return line;
};

int
ConfigReader::parse(FileHandle& file) 
{
    // reset file pointer
    file.rewind_fh();
    // block
    std::regex block_expr("\\s*\\[(.*)\\]\\s*");
    std::smatch block_match;
    // line
    std::regex line_expr("\\s*(.*)\\s*=\\s*(.*)\\s*");
    std::smatch line_match;

    char* line = file.find_line(name);
    if (!line) {
        return -1;
    }
    // get next line
    line = file.get_line();
    while (line != nullptr) {
        std::string string(line);
        if (std::regex_match(string, block_match, block_expr)) {
            line = nullptr;
        } else {
            if(string.find_first_not_of("\t\n\v\f\r") != std::string::npos) {
                if (std::regex_match(string, line_match, line_expr)) {
                    parse_line(trim(line_match[1]), line_match[2]);
                } else {
                    std::cout << "Input Error: " << line << "\n";
                    return -1;
                }
            }
            line = file.get_line();
        }
    }
    return 0;

};


void 
ConfigReader::parse_line(const std::string& key, const std::string& value) 
{
    if (key_in_map(key, ints)) {
        parse_int(key, value);
    } else if (key_in_map(key, floats)) {
        parse_floats(key, value);
    } else if (key_in_map(key, strings)) {
        parse_string(key, value);
    } else {
    }
};
