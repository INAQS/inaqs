#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <unordered_map>
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


class ConfigReader
{
public:
    explicit ConfigReader(std::string block): ints{}, name{block} {}
    // add entry by type
    void add_entry(std::string name, int def) { ints.emplace(name, def); }
    void add_entry(std::string name, double def) { floats.emplace(name, def); }
    void add_entry(std::string name, std::string def) { strings.emplace(name, def); }
    // parse file
    int parse(FileHandle& file);
    // get results
    inline int get_int(const std::string& key) { return ints[key]; }
    inline double get_double(const std::string& key) { return floats[key]; }
    inline std::string get_string(const std::string& key) { return strings[key]; }

private:
    void parse_line(const std::string& key, const std::string& value); 
    void parse_string(const std::string& key, const std::string& value) { strings[key] = trim(value); }
    void parse_floats(const std::string& key, const std::string& value) { floats[key] = std::stod(value); }
    void parse_int(const std::string& key, const std::string& value) { ints[key] = std::stoi(value); }

    std::string name{};
    std::unordered_map<std::string, int> ints{};
    std::unordered_map<std::string, double> floats{};
    std::unordered_map<std::string, std::string> strings{};
};


#endif
