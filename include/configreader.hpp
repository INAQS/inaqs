#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <unordered_map>
#include <vector>
#include <string>

#include <exception>

std::vector<std::string> split_string(const std::string& str, const std::string& delims=",");

class CannotOpenFile: 
    public std::exception
{
public:
    CannotOpenFile(std::string msg) : _msg{msg} {

    }

    virtual const char* what() const throw() {
        return _msg.c_str();
    }
private:
    std::string _msg{};
};


inline std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
};


template<typename Tkey, typename Tmap>
bool
key_in_map(const Tkey& key, const Tmap& map) {
    auto itr = map.find(key);
    if (itr == map.end()) {
        return false;
    }
    return true;
}


class File
{
    public:
        using iterator = FILE*;
        // open new file
        explicit File(const std::string filename) : name{filename} {
            fh_ = fopen(name.c_str(), "r");
            if (fh_ == NULL) {
                throw CannotOpenFile("Cannot open file " + name);
            }
        }
        // init file with existing open file and position
        explicit File(const int fileno, const int pos) {
            fh_ =  fdopen(fileno, "r");
            if (fh_ == NULL) {
                throw CannotOpenFile("Cannot open file with number" + std::to_string(fileno));
            }
            seek(pos);
        }
        File() = default;
        File(File&& rhs) {
            name = rhs.name;    
            fh_ = rhs.fh_;
            rhs.fh_ = nullptr;
        }

        ~File() noexcept { if (fh_) fclose(fh_);}

        inline iterator fh() const noexcept {return fh_;}

        inline int tell() const noexcept {return ftell(fh_);}
        inline int seek(const int offset) {return fseek(fh_, offset, SEEK_SET);}
        inline int fnum() const noexcept {return fileno(fh_);}

    private:
        std::string name{};    
        iterator fh_{nullptr};
};


class FileHandle
{
    public:
        explicit FileHandle(const std::string filename): file{filename}, _is_open{true} { }
        explicit FileHandle(const int fileno, const int pos) : file{fileno, pos}, _is_open{true} { }
        FileHandle(){}

        // Helper
        char* get_line();
        char* find_line(const std::string& key); 
        //
        inline void rewind_fh() { rewind(file.fh()); };

        inline int tell() const noexcept {return file.tell();}
        inline int seek(const int offset) {return file.seek(offset);}

        File::iterator fh() const noexcept {return file.fh();}

        File& f() noexcept {return file;}

        bool is_open() const noexcept { return _is_open; }

    private:
        File file;
        bool _is_open{false};
};

enum class DATA_TYPE {
    INT,
    DOUBLE,
    STRING, 
    IVEC,
    DVEC,
    SVEC,
    NONE
};


class Data
{
public:
    using types = DATA_TYPE;
    //
    explicit Data(int input, bool in_isset=true)                      : isset{in_isset}, data{new int(input)}, type{DATA_TYPE::INT} {}
    explicit Data(double input, bool in_isset=true)                   : isset{in_isset}, data{new double(input)},     type{DATA_TYPE::DOUBLE} {}
    explicit Data(std::string input, bool in_isset=true)              : isset{in_isset}, data{new std::string{input}},     type{DATA_TYPE::STRING} {}
    explicit Data(std::vector<std::string> input, bool in_isset=true) : isset{in_isset}, data{new std::vector<std::string>{input}},    type{DATA_TYPE::SVEC} {}
    explicit Data(std::vector<int> input, bool in_isset=true)         : isset{in_isset}, data{new std::vector<int>{input}},    type{DATA_TYPE::IVEC} {}
    explicit Data(std::vector<double> input, bool in_isset=true)      : isset{in_isset}, data{new std::vector<double>{input}},    type{DATA_TYPE::DVEC} {}
    Data(): type{DATA_TYPE::NONE} { }
    ~Data(); 
    Data(const Data& rhs);
    // copy operation
    bool get_data(double& out); 
    bool get_data(int& out); 
    bool get_data(std::string& out); 
    bool get_data(std::vector<int>&);
    bool get_data(std::vector<double>&);
    bool get_data(std::vector<std::string>&);
    // move operation
    bool move_data(double& out); 
    bool move_data(int& out); 
    bool move_data(std::string& out); 
    bool move_data(std::vector<int>&);
    bool move_data(std::vector<double>&);
    bool move_data(std::vector<std::string>&);
    //
    void set_from_string(const std::string& validation);
    //
    //
    std::string get_type();
    //
private:
    bool isset{true};
    void* data{nullptr};

    DATA_TYPE type{};
};



class ConfigBlockReader
{
public:
    using types = Data::types;
    explicit ConfigBlockReader(std::string block): name{block}, data{} {}
    // add entry by type
    template<typename T>
    void add_entry(const std::string name, const T def) { data.emplace(name, Data{def}); }

    void add_entry(const std::string name, const Data::types type) { 
    switch (type) {
        case types::INT:
            data.emplace(name, Data{0, false}); 
            break;
        case types::DOUBLE:
            data.emplace(name, Data{0.0, false}); 
            break;
        case types::STRING:
            data.emplace(name, Data{std::string(""), false}); 
            break;
        case types::SVEC:
            data.emplace(name, Data{std::vector<std::string>{}, false}); 
            break;
        case types::IVEC:
            data.emplace(name, Data{std::vector<int>{}, false}); 
            break;
        case types::DVEC:
            data.emplace(name, Data{std::vector<double>{}, false}); 
            break;
        default:
            data.emplace(name, Data{}); 
            break;
    }
    }
    // parse file
    int parse(FileHandle& file);
    // get results
    template<typename T>
    inline bool get_data(const std::string& key, T& val) { return data.at(key).get_data(val); }
    template<typename T>
    inline bool move_data(const std::string& key, T& val) { return data.at(key).move_data(val); }
    //
  Data operator[](const std::string& key) {return data.at(key);}

private:
    void parse_line(const std::string& key, const std::string& value); 

    std::string name{};
    std::unordered_map<std::string, Data> data;
};

FileHandle get_filehandle(const std::string filename); 

#endif
