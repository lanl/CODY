///
// A generic main program for an OpenCL mini-app
///

#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "app.hpp"

void usage(const char *name)
{
    std::cerr << "Usage: " << name << "[options] [args]\n"
              << "       " << name << "-h | --help\n";
    exit(1);
}


void help(const char *name)
{
    std::cout << "Usage: " << name << "[options]\n"
              << "\n"
              << "Options:\n"
              << "  -h | --help         print help (this message)\n"
              << "  -D | --device-list  comma-separated list of devices to use\n"
              << "  -d | --debug        enable debugging\n"
              << "  -p | --profile      enable profiling\n"
              << "  -v | --verbose      verbose output\n"
              << "\n";
    exit(0);
}


int csv_to_list(const char *string, std::vector<int> & list)
{
    std::string sep(",");
    std::string s(string);
    size_t start = 0;
    size_t end = 0;

    do {
        int val;
        end = s.find(sep, start);
        if (!(std::istringstream(s.substr(start, end - start)) >> val)) {
            return -1;
        }
        list.push_back(val);
        start = end + sep.size();
    } while (end != std::string::npos);
    return 0;
}


int main (int argc, char *argv[])
{
    // default options
    int debug = 0;
    int profile = 0;
    int verbose = 0;
    std::vector<int> device_list;

    // process options
    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"debug", no_argument, NULL, 'd'},
            {"profile", no_argument, NULL, 'p'},
            {"verbose", no_argument, NULL, 'v'},
            {"device-list", required_argument, NULL, 'D'},
            {NULL, 0, NULL, 0}
        };
        int option_index = 0;
        int c = getopt_long(argc, argv, "dhpvD:", long_options, &option_index);
        if (c == -1) {
            break;
        } else if (c == '?') {
            usage(argv[0]);
            return 1;
        } else if ('h' == c) {
            help(argv[0]);
            return 0;
        } else if ('d' == c) {
            debug = 1;
        } else if ('p' == c) {
            profile = 1;
        } else if ('v' == c) {
            verbose = 1;
        } else if ('D' == c) {
            if (csv_to_list(optarg, device_list) < 0) {
                fprintf(stderr, "Invalid device list: %s\n", optarg);
                usage(argv[0]);
                return 1;
            }
        } else {
            return 1;
        }
    }
    if (debug) {
        std::cerr << "[Options]\n"
                  << "  debug = " << debug << "\n"
                  << "  profile = " << profile << "\n"
                  << "  verbose = " << verbose << "\n"
                  << "  device list = ";
        for (size_t i = 0; i < device_list.size(); ++i) {
            std::cerr << " " << device_list[i];
        }
        std::cerr << "\n";
        if (optind < argc) {
            std::cerr << "Non-option arguments: ";
            while (optind < argc) {
                std::cerr << argv[optind++];
            }
            std::cerr << "\n";
        }
    }

    // start work
    try {
        App app(debug,
                profile,
                verbose,
                device_list);
		app.build_program();
        app.host_run();
        
    }
    catch (cl::Error const & e) {
        std::cerr << "ERROR: OpenCL: "
                  << e.what()
                  << "("
                  << App::opencl_error_string(e.err())
                  << ")\n";
        return 1;
    }

    return 0;
}
