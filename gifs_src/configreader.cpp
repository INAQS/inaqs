#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//
#include <iostream>
#include <regex>
//
#include "configreader.hpp"
//
#define from_void(type, value) new type(*static_cast<type*>(value))
#define delete_void_as(type, value) delete static_cast<type*>(value)

template<typename T>
inline
void
copy_output(T& out, void* data) {
    out = *static_cast<T*>(data);
}


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
            data = from_void(int, rhs.data);
            break;
        case types::ULINT:
            data = from_void(size_t, rhs.data);
            break;
        case types::DOUBLE:
            data = from_void(double, rhs.data);
            break;
        case types::STRING:
            data = from_void(std::string, rhs.data);
            break;
        case types::SVEC:
            data = from_void(std::vector<std::string>, rhs.data);
            break;
        case types::IVEC:
            data = from_void(std::vector<int>, rhs.data);
            break;
        case types::DVEC:
            data = from_void(std::vector<double>, rhs.data);
            break;
        default:
            break;
    }
}

Data::~Data() {
    if (data == nullptr)
        return;
    switch (type) {
        case types::INT:
            delete_void_as(int, data);
            break;
	case types::ULINT:
            delete_void_as(size_t, data);
            break;
        case types::DOUBLE:
            delete_void_as(double, data);
            break;
        case types::STRING:
            delete_void_as(std::string, data);
            break;
        case types::SVEC:
            delete_void_as(std::vector<std::string>, data);
            break;
        case types::IVEC:
            delete_void_as(std::vector<int>, data);
            break;
        case types::DVEC:
            delete_void_as(std::vector<double>, data);
            break;
        default:
            break;
    }
}


void
Data::set_from_string(const std::string& value) 
{
    bool do_set = true;
    switch (type) {
        case types::INT:
            data = new int(std::stoi(value));
            break;
        case types::ULINT:
            data = new size_t(std::stoul(value));
            break;
        case types::DOUBLE:
            data = new double(std::stod(value));
            break;
        case types::STRING:
            data = new std::string(value);
            break;
        case types::IVEC:
            data = new std::vector<int>(string_to_ivec(value));
            break;
        case types::DVEC:
            data = new std::vector<double>(string_to_dvec(value));
            break;
        case types::SVEC:
            data = new std::vector<std::string>(split_string(value));
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
	case (types::ULINT):
            return "size_t";
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
    copy_output(out, data);
    return true;
};


bool
Data::get_data(int& out) {
    if (type != types::INT || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}

bool
Data::get_data(size_t& out) {
    if (type != types::ULINT || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}

bool
Data::get_data(std::string& out)
{
    if (type != types::STRING || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}

bool
Data::get_data(std::vector<std::string>& out) {
    if (type != types::DVEC || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}

bool
Data::get_data(std::vector<int>& out) {
    if (type != types::IVEC || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}

bool
Data::get_data(std::vector<double>& out) {
    if (type != types::DVEC || !isset) {
        return false;
    }
    copy_output(out, data);
    return true;
}


static const std::regex block_expr("\\s*\\[(.*)\\]\\s*");
static const std::regex line_expr("\\s*(.*)\\s*=\\s*(.*)\\s*");
static const std::regex comment_expr("^#|^\\s*$");

int
ConfigBlockReader::parse(FileHandle& file) 
{
    if (!file.is_open()) {
        return -2;
    }
    // reset file pointer
    file.rewind_fh();
    std::smatch block_match;
    std::smatch line_match;

    char* line = file.find_line("[" + name + "]");
    if (!line) {
      parsed = true;
      std::cerr << "Warning, no block found for "<< name << std::endl;
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
	        if (std::regex_search(string, comment_expr)){
                    // pass for comment
	        }
                else if (std::regex_match(string, line_match, line_expr)) {
                    parse_line(trim(line_match[1]), line_match[2]);
                } else {
		  std::cerr << "Input Error: " << line << std::endl;
                    return -1;
                }
            }
            line = file.get_line();
        }
    }
    parsed = true;
    return 0;
};


std::vector<std::string> ConfigBlockReader::enumerate_keys(void){
  std::vector<std::string> keys;
  keys.reserve(data.size());
  
  for (const auto &p: data){
    keys.push_back(p.first);
  }

  return keys;
}

void 
ConfigBlockReader::parse_line(const std::string& key, const std::string& value)
{
    if (key_in_map(key, data)) {
        data[key].set_from_string(value);
    } else {
        std::cerr << "Unknown field " << key << std::endl;
    }
};


FileHandle get_filehandle(const std::string filename) {
    try {
        return FileHandle(filename);
    } catch (std::exception& e) {
        return FileHandle{};
    }
};
