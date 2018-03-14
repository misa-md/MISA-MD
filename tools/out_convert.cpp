#include <iostream>
#include <fstream>

#define BUF_SIZE 1024
#define HEADER_SIZE 1024
#define IN_FILENAME_DEF "md.out"
#define OUT_FILENAME_DEF "md.out.txt"

typedef char byte;

void parseArgv(std::string &inFilename, std::string &outFilename, int argc, char *argv[]);

/**
 * this file convert binary output file of CrystalMD program to readable txt format.
 * usage: ./conv -f [input_filename] -o [output_filename]
 * @param argc count of argv string.
 * @param argv example "-f md.out -o md.out.txt"
 * @return returns 0 if conversion is ok, otherwise returns non-zero.
 */
int main(int argc, char *argv[]) {
    std::string inFilename, outFilename;
    parseArgv(inFilename, outFilename, argc, argv);

    std::ifstream infile(inFilename, std::ios::in | std::ios::binary); // td::ios::ate
    if (!infile.good()) {
        std::cerr << "can not access the input file" << std::endl;
        return 1;
    }

    std::ofstream outfile(outFilename, std::ios::out);
    if (!outfile.good()) {
        std::cerr << "can not access the output file" << std::endl;
        return 1;
    }

    infile.seekg(HEADER_SIZE, std::ios::beg); // head.
    byte *buffer = new byte[BUF_SIZE];
    double *x;

    while (!infile.eof()) {
        infile.read(buffer, BUF_SIZE);
        x = (double *) buffer; // convert byte to double.
        for (int i = 0; i < BUF_SIZE / sizeof(double); i++) {
            if (i % 4 == 3) {
                outfile << x[i] << std::endl;
            } else {
                outfile << x[i] << "\t";
            }
        }
    }

    infile.close();
    outfile.close();
    return 0;
}

void parseArgv(std::string &inFilename, std::string &outFilename, int argc, char *argv[]) {
    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "-f" && i + 1 < argc) { // has "-f" and not the last param.
            inFilename = std::string(argv[i + 1]);
            i++; // skip value of inFilename.
            continue;
        }
        if (std::string(argv[i]) == "-o" && i + 1 < argc) { // has "-o" and not the last param.
            outFilename = std::string(argv[i + 1]);
            i++; // skip value of outFilename.
            continue;
        }
    }

    // set default values if is empty string.
    if (inFilename == "") {
        inFilename = IN_FILENAME_DEF;
    }
    if (outFilename == "") {
        outFilename = OUT_FILENAME_DEF;
    }
    std::clog << "input file will be " << inFilename << std::endl;
    std::clog << "output file will be " << outFilename << std::endl;
}