//
// Created by Dominik on 23.09.26.
//

#include <fstream>
#include <iosfwd>
#include <iostream>
using std::ifstream;
using std::string;

bool readFile(const char* filename, string& outstr) {


    ifstream f(filename);

    bool success = false;

    if (f.is_open()) {
        string line;
        while (getline(f, line)) {
            outstr.append(line);
            outstr.append("\n");
        }
        f.close();

        success = true;
    }
    else {
        std::cerr << "File not found" << std::endl;
    }

    return success;
}