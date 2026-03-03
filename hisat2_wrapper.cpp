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

bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

int main(int argc, char *argv[]) {
    std::string base_bin = "hisat2-align-s";
    std::string fake_name = "hisat-3n";
    bool shift = false;
    bool is_align = true;

    if (argc >= 2) {
        std::string cmd = argv[1];
        if (cmd == "build") {
            base_bin = "hisat2-build-s";
            fake_name = "hisat-3n-build";
            shift = true;
            is_align = false;
        } else if (cmd == "inspect") {
            base_bin = "hisat2-inspect-s";
            fake_name = "hisat-3n-inspect";
            shift = true;
            is_align = false;
        } else if (cmd == "align") {
            base_bin = "hisat2-align-s";
            fake_name = "hisat-3n";
            shift = true;
        }
    }

    std::string bin_path = get_executable_dir() + "/" + base_bin;
    setenv("HISAT2_WRAPPER_NAME", "hisat3n", 1);

    std::vector<std::string> new_args_storage;
    new_args_storage.push_back(fake_name);
    
    int start_idx = shift ? 2 : 1;
    for (int i = start_idx; i < argc; ++i) {
        std::string arg = argv[i];
        
        // Map --index to -x
        if (arg == "--index") {
            new_args_storage.push_back("-x");
            continue;
        }

        // Handle compressed files for alignment
        if (is_align && (arg == "-U" || arg == "-1" || arg == "-2" || arg == "--12")) {
            new_args_storage.push_back(arg);
            if (i + 1 < argc) {
                std::string files = argv[++i];
                // Split comma-separated files
                size_t pos = 0;
                std::string token;
                std::string processed_files = "";
                while (true) {
                    size_t next_pos = files.find(',', pos);
                    token = files.substr(pos, next_pos - pos);
                    
                    if (ends_with(token, ".gz")) {
                        processed_files += "<(gzip -dc " + token + ")";
                    } else if (ends_with(token, ".bz2")) {
                        processed_files += "<(bzip2 -dc " + token + ")";
                    } else {
                        processed_files += token;
                    }
                    
                    if (next_pos == std::string::npos) break;
                    processed_files += ",";
                    pos = next_pos + 1;
                }
                new_args_storage.push_back(processed_files);
            }
            continue;
        }
        
        new_args_storage.push_back(arg);
    }

    std::vector<char*> final_args;
    for (const auto& s : new_args_storage) {
        final_args.push_back(const_cast<char*>(s.c_str()));
    }
    final_args.push_back(NULL);

    // Note: process substitution <() is a shell feature. 
    // To use it, we must execute via shell.
    if (is_align) {
        std::string shell_cmd = bin_path;
        for (size_t i = 0; i < new_args_storage.size(); ++i) {
            // Skip the fake name which is argv[0]
            if (i == 0) continue;
            shell_cmd += " " + new_args_storage[i];
        }
        
        // We use bash -c to support process substitution
        execl("/bin/bash", "bash", "-c", shell_cmd.c_str(), (char*)NULL);
    } else {
        execv(bin_path.c_str(), final_args.data());
    }

    std::cerr << "Error: Could not execute " << bin_path << std::endl;
    return 1;
}
