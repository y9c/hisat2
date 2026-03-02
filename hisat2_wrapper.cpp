#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <libgen.h>
#include <stdlib.h>

void print_usage() {
    std::cerr << "HISAT3N (Optimized HISAT2-3N Unified Tool)" << std::endl;
    std::cerr << "Usage: hisat3n <command> [options]" << std::endl;
    std::cerr << "Commands:" << std::endl;
    std::cerr << "  align    Align reads to reference (default)" << std::endl;
    std::cerr << "  build    Build index from reference" << std::endl;
    std::cerr << "  inspect  Inspect index properties" << std::endl;
}

int main(int argc, char *argv[]) {
    std::string bin_path = "hisat2-align-s";
    bool shift = false;

    // Suppress the "run directly" warning
    setenv("HISAT2_WRAPPER_NAME", "hisat3n", 1);

    if (argc >= 2) {
        std::string cmd = argv[1];
        if (cmd == "build") {
            bin_path = "hisat2-build-s";
            shift = true;
        } else if (cmd == "inspect") {
            bin_path = "hisat2-inspect-s";
            shift = true;
        } else if (cmd == "align") {
            bin_path = "hisat2-align-s";
            shift = true;
        }
    }

    std::vector<char*> args;
    args.push_back((char*)bin_path.c_str());
    
    int start_idx = shift ? 2 : 1;
    for (int i = start_idx; i < argc; ++i) {
        args.push_back(argv[i]);
    }
    args.push_back(NULL);

    execvp(bin_path.c_str(), args.data());

    std::cerr << "Error: Could not execute " << bin_path << ". Ensure it is in your PATH." << std::endl;
    return 1;
}
