#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include "stringutil.hpp"
#include "logging.hpp"

using namespace std;

struct csv_line {
    size_t user_id;
    size_t url_id;
    short reaction;

    csv_line(size_t uid, size_t rid, short react) : user_id(uid), url_id(rid), reaction(react) {}
    csv_line() : user_id(0), url_id(0), reaction(0) {}
};

csv_line parse_csv_line(const string & line, char input_delim, vector<string> & input_features) {
    size_t user_id, url_id;
    split(line, input_delim, input_features);
    if (input_features.size() < 2) return csv_line(0, 0, 0);
    istringstream(input_features[0]) >> user_id;
    istringstream(input_features[1]) >> url_id;
    return csv_line(user_id, url_id, 0);
}

csv_line parse_csv_line_reaction(const string & line, char input_delim, vector<string> & input_features) {
    csv_line csvl = parse_csv_line(line, input_delim, input_features);
    if (csvl.user_id > 0) {
        istringstream(input_features[input_features.size() - 1]) >> csvl.reaction;
    }
    return csvl;
}

vector<string> readFileToMemory(string input, bool do_log = true) {
    ifstream ins(input.c_str());
    if (!ins.is_open())
    {
        cerr << "Cannot open " << input << " for reading!" << endl;
        exit(-1);
    }

    if (do_log) LOG("Reading file " << input << "...");
    vector<string> myLines;
    copy(istream_iterator<string>(ins),
         istream_iterator<string>(),
         back_inserter(myLines));
    if (do_log) LOG("  ...file read into memory. " << myLines.size() << " lines.");
    ins.close();
    return myLines;
}

vector<string> readLinesFromFile(istream & ins, size_t num_lines) {
    vector<string> res;
    if (!ins.good()) return res;
    string line;
    for (size_t i=0; i<num_lines; ++i) {
        if (!ins.good()) break;
        getline(ins, line);
        res.push_back(line);
    }
    return res;
}

bool open_for_reading(const string & fname, ifstream & ifs) {
    ifs.open(fname.c_str());
    if (!ifs.is_open()) {
        LOG("Cannot open " << fname << " for reading!");
        return false;
    }
    return true;
}

bool open_for_writing(const string & fname, ofstream & ofs) {
    ofs.open(fname.c_str());
    if (!ofs.is_open()) {
        LOG("Cannot open " << fname << " for writing!");
        return false;
    }
    return true;
}

void write_to_file(const string & filename, std::function<void (ofstream &)> do_write, bool do_log = false) {
    ofstream ofs;
    if (open_for_writing(filename, ofs)) {
        do_write(ofs);
        ofs.close();
    } else {
        if (do_log) LOG("ERROR: Cannot open file " << filename << " for writing!");
    }
}

void process_file_buffered(const string & filename, std::function<void (const size_t i, const string &)> do_with_line, bool do_log = false, size_t buffer_size = 5000000) {
    ifstream ins;
    size_t line_counter = 0;
    if (open_for_reading(filename, ins)) {
        while (ins.good()) {
            vector<string> myLines = readLinesFromFile(ins, buffer_size);
            for (size_t i = 0; i < myLines.size(); ++i) {
                string l = trim(myLines[i]);
                if (l.size() > 0) do_with_line(i, l);
            }
            line_counter += myLines.size();
            if (do_log) LOG("\t\t..." << line_counter << " lines processed...");
        }
        ins.close();
    }
}

void process_csv_file_buffered(const string & filename, const string & separators,
        std::function<void (const size_t i, const vector<string> &)> do_with_split_line, bool do_log = false, size_t buffer_size = 5000000) {
    vector<string> parts;
    process_file_buffered(filename,  [&] (const size_t i, const string & line) {
        boost::split(parts, line, boost::is_any_of(separators));
        do_with_split_line(i, parts);
    });
}

vector<boost::filesystem::path> get_list_of_files(const string & input) {
    boost::filesystem::path p(input);
    vector<boost::filesystem::path> infiles;
    if (boost::filesystem::exists(p)) {
        if (boost::filesystem::is_directory(p)) {
            LOG(p << " is a directory. Processing all .csv and .txt files in it...");
            boost::filesystem::recursive_directory_iterator endit;
            for (boost::filesystem::recursive_directory_iterator it(p); it != endit; ++it) {
                boost::filesystem::path cur_file = *it;
                if( !boost::filesystem::is_regular_file( cur_file ) ) continue;
                if( cur_file.extension() != ".csv" && cur_file.extension() != ".txt" ) continue;
                infiles.push_back(cur_file);
            }
        } else if (boost::filesystem::is_regular_file(p)) {
            LOG(p << " is a file. Processing...");
            infiles.push_back(p);
        } else {
            cerr << p << " is neither a dir nor a file. Exiting." << endl;
            exit(0);
        }
    } else {
        cerr << p << " does not exist." << endl;
        exit(0);
    }
    return infiles;
}

template<typename T>
void print_vec_to_cout_tabbed(size_t size, T* v) {
    for (size_t j = 0 ; j < size; ++j) {
        cout << "\t" << v[j];
    }
    cout << endl;
}

template<typename T>
void print_1d_to_cout(const string & name, size_t size, T* v) {
    cout << name << endl;
    for (size_t i = 0 ; i < size; ++i) {
        cout << i << "\t" << v[i] << endl;
    }
}

template<typename T>
void print_2d_to_cout(const string & name, size_t dim1, size_t dim2, T** v) {
    cout << name << endl;
    for (size_t i = 0 ; i < dim1; i ++) {
        cout << i << "\t";
        for (size_t j = 0 ; j < dim2; ++j) {
            cout << "\t" << v[i][j];
        }
        cout << endl;
    }
}
