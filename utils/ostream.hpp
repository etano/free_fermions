#ifndef UTILS_OSTREAM_HPP
#define UTILS_OSTREAM_HPP

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vs){
    for(unsigned i=0; i<vs.size(); i++)
        out << i << " " << vs[i] << std::endl;
    return out;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& out, const std::vector<std::pair<T,U>>& vs){
    for(unsigned i=0; i<vs.size(); i++)
        out << i << " " << vs[i].first << " " << vs[i].second << std::endl;
    return out;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& out, const std::map<T,U>& vs){
    for(const auto& v : vs)
        out << v.first << " " << v.second << std::endl;
    return out;
}

template<class T, class U, class V>
std::ostream& operator<<(std::ostream& out, const std::map<T,std::pair<U,V>>& vs){
    for(const auto& v : vs)
        out << v.first << " " << v.second.first << " " << v.second.second << std::endl;
    return out;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& out, const std::unordered_map<T,U>& vs){
    out << std::map<T,U>(vs.begin(),vs.end());
    return out;
}

template<class T, class U, class V>
std::ostream& operator<<(std::ostream& out, const std::unordered_map<T,std::pair<U,V>>& vs){
    out << std::map<T,std::pair<U,V>>(vs.begin(),vs.end());
    return out;
}

namespace utils{

template<class T>
void write(const std::string& f, const T val, const unsigned digits=16, const bool print=false){
    std::ofstream f_stream;
    f_stream.precision(digits);
    f_stream.open(f);
    f_stream << val << "\n";
    f_stream.close();
    if(print)
        std::cout << val << std::endl;
}

void split(const std::string& s, char delim, std::vector<std::string>& elem) {
    std::istringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        if(!item.empty())
            elem.push_back(item);
    }
}

std::vector<std::vector<std::string>> read_file(const std::string& file_name, char delim, const unsigned n_columns){
    std::vector<std::vector<std::string>> data(n_columns);
    std::ifstream infile(file_name);
    for(std::string line; std::getline(infile,line);){
        std::vector<std::string> elements;
        utils::split(line,delim,elements);
        if(!elements.empty()){
            for(unsigned i=0; i<elements.size(); i++)
                data[i].push_back(elements[i]);
        }
    }
    return std::move(data);
}

}

#endif
