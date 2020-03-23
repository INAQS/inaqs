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

Data::Data(const Data& rhs) {
    type = rhs.type;
    switch (type) {
        case types::INT:
            integer = rhs.integer;
            break;
        case types::DOUBLE:
            dbl = rhs.dbl;
            break;
        case types::STRING:
            str = rhs.str;
            break;
        default:
            type = types::NONE;
            break;
    }

}
void
Data::set_from_string(const std::string& value) 
{
    switch (type) {
        case types::INT:
            integer = std::stoi(value);
            break;
        case types::DOUBLE:
            dbl = std::stod(value);
            break;
        case types::STRING:
            str = value.c_str();
            break;
        default:
            break;
    }
};


bool
Data::get_data(double& out) {
    if (type != types::DOUBLE) {
        return false;
    } 
    out = dbl;
    return true;
    
};

std::string
Data::get_type() {
    switch (type) {
        case (types::INT):
            return "int";
        case (types::DOUBLE): 
            return "double";
        case (types::STRING): 
            return "break";
        default:
            return "unknown";
    }
};

bool
Data::get_data(int& out) {
    if (type != types::INT) {
        return false;
    }
    out = integer;
    return true;
}

bool
Data::get_data(std::string& out)
{
    if (type != types::STRING) {
        return false;
    }
    out = std::string(str);
    return true;
}

int
ConfigBlockReader::parse(FileHandle& file) 
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
ConfigBlockReader::parse_line(const std::string& key, const std::string& value)
{
    if (key_in_map(key, data)) {
        data[key].set_from_string(value);
    } else {
    }
};
