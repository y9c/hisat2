#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <libgen.h>
#include <stdlib.h>
#include <limits.h>

std::string get_executable_dir() {
    char result[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    if (count != -1) {
        std::string path(result, count);
        size_t last_slash = path.find_last_of('/');
        if (last_slash != std::string::npos) {
            return path.substr(0, last_slash);
        }
    }
    return ".";
}

int main(int argc, char *argv[]) {
    std::string base_bin = "hisat2-align-s";
    std::string fake_name = "hisat-3n";
    bool shift = false;

    if (argc >= 2) {
        std::string cmd = argv[1];
        if (cmd == "build") {
            base_bin = "hisat2-build-s";
            fake_name = "hisat-3n-build";
            shift = true;
        } else if (cmd == "inspect") {
            base_bin = "hisat2-inspect-s";
            fake_name = "hisat-3n-inspect";
            shift = true;
        } else if (cmd == "align") {
            base_bin = "hisat2-align-s";
            fake_name = "hisat-3n";
            shift = true;
        }
    }

    std::string bin_path = get_executable_dir() + "/" + base_bin;
    setenv("HISAT2_WRAPPER_NAME", "hisat3n", 1);

    std::vector<char*> args;
    args.push_back((char*)fake_name.c_str());
    
    int start_idx = shift ? 2 : 1;
    for (int i = start_idx; i < argc; ++i) {
        args.push_back(argv[i]);
    }
    args.push_back(NULL);

    execv(bin_path.c_str(), args.data());

    std::cerr << "Error: Could not execute " << bin_path << std::endl;
    return 1;
}
