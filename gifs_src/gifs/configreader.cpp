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


std::vector<std::string>
split_string(const std::string& str, const std::string& delims)
{
    std::vector<std::string> outvec{};
    std::size_t current, previous = 0;
    current = str.find_first_of(delims);
    while (current != std::string::npos) {
        outvec.push_back(str.substr(previous, current-previous));
        previous = current + 1;
        current = str.find_first_of(delims, previous);
    };
    outvec.push_back(str.substr(previous, current-previous));
    return outvec;
};


std::vector<double>
string_to_dvec(const std::string& str)
{
    const auto value = split_string(str);
    std::vector<double> out{};
    for (auto& ele: value) {
        out.push_back(std::stod(ele));            
    }
    return out;
};


std::vector<int>
string_to_ivec(const std::string& str)
{
    const auto value = split_string(str);
    std::vector<int> out{};
    for (auto& ele: value) {
        out.push_back(std::stoi(ele));            
    }
    return out;
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
    isset = rhs.isset;
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
        case types::SVEC:
            svec = rhs.svec;
            break;
        case types::IVEC:
            ivec = rhs.ivec;
            break;
        case types::DVEC:
            dvec = rhs.dvec;
            break;
        default:
            type = types::NONE;
            break;
    }
}


void
Data::set_from_string(const std::string& value) 
{
    bool do_set = true;
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
        case types::IVEC:
            ivec = string_to_ivec(value);
            break;
        case types::DVEC:
            dvec = string_to_dvec(value);
            break;
        case types::SVEC:
            std::cout << "value = " << value << "\n";
            std::cout << "called split_string \n"; 
            svec = split_string(value);
            std::cout << "called split_string \n"; 
            break;
        default:
            do_set = false;
            break;
    }
    if (do_set) {
        isset = true;
    }
};


std::string
Data::get_type() {
    switch (type) {
        case (types::INT):
            return "int";
        case (types::DOUBLE): 
            return "double";
        case (types::STRING): 
            return "string";
        case (types::IVEC): 
            return "ivec";
        case (types::DVEC): 
            return "dvec";
        case (types::SVEC): 
            return "svec";
        default:
            return "unknown";
    }
};

bool
Data::get_data(double& out) {
    if (type != types::DOUBLE || !isset) {
        return false;
    } 
    out = dbl;
    return true;
};


bool
Data::get_data(int& out) {
    if (type != types::INT || !isset) {
        return false;
    }
    out = integer;
    return true;
}

bool
Data::get_data(std::string& out)
{
    if (type != types::STRING || !isset) {
        return false;
    }
    out = str;
    return true;
}

bool
Data::get_data(std::vector<std::string>& out) {
    if (type != types::DVEC || !isset) {
        return false;
    }
    out = svec;
    return true;
}

bool
Data::get_data(std::vector<int>& out) {
    if (type != types::IVEC || !isset) {
        return false;
    }
    out = ivec;
    return true;
}

bool
Data::get_data(std::vector<double>& out) {
    if (type != types::DVEC || !isset) {
        return false;
    }
    out = dvec;
    return true;
}

bool
Data::move_data(double& out) {
    if (type != types::DOUBLE || !isset) {
        return false;
    } 
    out = std::move(dbl);
    return true;
};


bool
Data::move_data(int& out) {
    if (type != types::INT || !isset) {
        return false;
    }
    out = std::move(integer);
    return true;
}

bool
Data::move_data(std::string& out)
{
    if (type != types::STRING || !isset) {
        return false;
    }
    out = std::move(str);
    return true;
}

bool
Data::move_data(std::vector<std::string>& out) {
    if (type != types::SVEC || !isset) {
        return false;
    }
    out = std::move(svec);
    return true;
}

bool
Data::move_data(std::vector<int>& out) {
    if (type != types::IVEC || !isset) {
        return false;
    }
    out = std::move(ivec);
    return true;
}

bool
Data::move_data(std::vector<double>& out) {
    if (type != types::DVEC || !isset) {
        return false;
    }
    out = std::move(dvec);
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
        std::cout << "Unkown field " << key << std::endl;
    }
};
