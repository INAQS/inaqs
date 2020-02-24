#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <unordered_map>
#include <vector>
#include <string>


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
            if (fh_ == NULL) exit(EXIT_FAILURE);
        }
        // init file with existing open file and position
        explicit File(const int fileno, const int pos) {
	  (void) pos;
            fh_ =  fdopen(fileno, "r");
            if (fh_ == NULL) exit(EXIT_FAILURE);
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
        explicit FileHandle(const std::string filename) : file{filename} { }
        explicit FileHandle(const int fileno, const int pos) : file{fileno, pos} { }
        // Helper
        char* get_line();
        char* find_line(const std::string& key); 
        //
        inline void rewind_fh() { rewind(file.fh()); };

        inline int tell() const noexcept {return file.tell();}
        inline int seek(const int offset) {return file.seek(offset);}

        File::iterator fh() const noexcept {return file.fh();}

        File f() noexcept {return file;}
    private:
        File file;
};

class Data
{
public:
    explicit Data(int input) : integer{input}, type{data_type::INT} {}
    explicit Data(double input) : dbl{input}, type{data_type::DOUBLE} {}
    explicit Data(std::string input) : str{input}, type{data_type::STRING} {}
    Data(): type{data_type::NONE} {}
    ~Data(){}
    Data(const Data& rhs);
    // copy operation
    bool get_data(double& out); 
    bool get_data(int& out); 
    bool get_data(std::string& out); 
    // move operation

    /*
    void move_data(double& out); 
    void move_data(int& out); 
    void move_data(std::vector<int>& out); 
    void move_data(std::vector<double>& out); 
    void move_data(std::string& out); 
    */
    void set_from_string(const std::string& validation);
private:
    std::string get_type();
    union
    {
        double dbl;
        int integer;
        std::string str{};
    };

    enum class data_type {
        INT,
        DOUBLE,
        STRING, 
        NONE
    };
    data_type type{};
    using types = Data::data_type;
};


class ConfigBlockReader
{
public:
    explicit ConfigBlockReader(std::string block): name{block}, data{} {}
    // add entry by type
    void add_entry(const std::string name, const int def) { data.emplace(name, Data{def}); }
    void add_entry(const std::string name, const double def) { data.emplace(name, Data{def}); }
    void add_entry(const std::string name, const std::string def) { data.emplace(name, Data{def}); }
    // parse file
    int parse(FileHandle& file);
    // get results
    inline bool get_data(const std::string& key, int& val) { return data[key].get_data(val); }
    inline bool get_data(const std::string& key, double& val) { return data[key].get_data(val); }
    inline bool get_data(const std::string& key, std::string& val) { return data[key].get_data(val); }

    Data operator[](const std::string& key) {return data[key];}

private:
    void parse_line(const std::string& key, const std::string& value); 

    std::string name{};
    std::unordered_map<std::string, Data> data;
};


#endif
