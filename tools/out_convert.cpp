#include <iostream>
#include <fstream>

#define BUF_SIZE (1024) // can't be more than twice of block_size.
#define HEADER_SIZE 128
#define _LOCAL_HEADER_SIZE 128
#define IN_FILENAME_DEF "md.out"
#define OUT_FILENAME_DEF "md.out.txt"
//#define DEBUG

typedef char byte;

typedef struct {
    size_t atoms_count;
} localHeader;

typedef struct {
    unsigned long id;
    size_t step;
    int type;
    double atom_location[3]; // atom location
    double atom_velocity[3]; // atom velocity
} type_atom;

const long BLOCK_SIZE = (1024 * sizeof(type_atom));

void parseArgv(int &n, std::string &inFilename, std::string &outFilename, int argc, char *argv[]);

void readBytes(std::ifstream &infile, byte *buff, size_t buff_size, int rank_n, int rank);

void checkoutHeaders(std::ifstream &infile, int rank);

void checkoutToAtoms(std::ifstream &infile, int n_ranks, int rank);

std::string getNameByEleName(int type);

/**
 * this file convert binary output file of CrystalMD program to readable txt format.
 * usage: ./conv -f [input_filename] -o [output_filename]
 * @param argc count of argv string.
 * @param argv example "-f md.out -o md.out.txt"
 * @return returns 0 if conversion is ok, otherwise returns non-zero.
 */
int main(int argc, char *argv[]) {
    std::string inFilename, outFilename;
    int n_rank; // count of ranks
    parseArgv(n_rank, inFilename, outFilename, argc, argv);

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
    auto buffer = new type_atom[BUF_SIZE];

    outfile << "id \tstep \ttype \tlocate.x \tlocate.y \tlocate.z \tv.x \tv.y \tv.z" << std::endl;
//    size_t atom_read = 10;
//    TEST(infile, (byte *) buffer, atom_read * sizeof(type_atom), n_rank, 0);
    for (int rank = 0; rank < n_rank; rank++) {
        checkoutHeaders(infile, rank);
        localHeader local_header;
        infile.read((byte *) &local_header, sizeof(localHeader));

        checkoutToAtoms(infile, n_rank, rank);

#ifdef DEBUG
        std::cout << "start rank " << rank << "\n";
#endif
        while (local_header.atoms_count > 0) { // left atoms count
            size_t atom_read = local_header.atoms_count > BUF_SIZE ? BUF_SIZE : local_header.atoms_count;
            readBytes(infile, (byte *) buffer, atom_read * sizeof(type_atom), n_rank, rank);
            for (int i = 0; i < atom_read; i++) {
                outfile << buffer[i].id << "\t" << buffer[i].step << "\t"
                        << getNameByEleName(buffer[i].type) << "\t"
                        << buffer[i].atom_location[0] << "\t" << buffer[i].atom_location[1] << "\t"
                        << buffer[i].atom_location[2] << "\t"
                        << buffer[i].atom_velocity[0] << "\t" << buffer[i].atom_velocity[1] << "\t"
                        << buffer[i].atom_velocity[2] << "\t" << std::endl;
            }
            local_header.atoms_count -= atom_read;
        }
    }

    infile.close();
    outfile.close();
    return 0;
}

void parseArgv(int &n, std::string &inFilename, std::string &outFilename, int argc, char *argv[]) {
    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "-f" && i + 1 < argc) { // has "-f" and not the last param.
            inFilename = std::string(argv[i + 1]);
            i++; // skip value of inFilename.
            continue;
        }
        if (std::string(argv[i]) == "-n" && i + 1 < argc) { // has "-f" and not the last param.
            n = std::stoi(argv[i + 1]);
            i++;
            continue;
        }
        if (std::string(argv[i]) == "-o" && i + 1 < argc) { // has "-o" and not the last param.
            outFilename = std::string(argv[i + 1]);
            i++; // skip value of outFilename.
            continue;
        }
    }

    // set default values if is empty string.
    if (inFilename.empty()) {
        inFilename = IN_FILENAME_DEF;
    }
    if (outFilename.empty()) {
        outFilename = OUT_FILENAME_DEF;
    }
    std::clog << "input file will be " << inFilename << std::endl;
    std::clog << "output file will be " << outFilename << std::endl;
}

void readBytes(std::ifstream &infile, byte *buff, size_t buff_size, int rank_n, int rank) {
    long begin = infile.tellg();
    long H = HEADER_SIZE + _LOCAL_HEADER_SIZE * rank_n;
    long k = (begin - H) / BLOCK_SIZE / rank_n; // now it is on the (k+1)th blocks.
    long next_rank_block = H + (rank_n * (k + 1) + rank) * BLOCK_SIZE;
    long next_block = H + (rank_n * k + rank + 1) * BLOCK_SIZE;

    long left_this_block = next_block - begin > buff_size ? buff_size : next_block - begin;
#ifdef DEBUG
    std::cout << "k:" << k << " next_rank_block: " << next_rank_block << " next_block: " << next_block
              << " left_this_block: " << left_this_block << ".\n";

    std::cout << "(1) read from " << infile.tellg() << " of " << left_this_block << " bytes. \n";
#endif

    infile.read(buff, left_this_block);
    int gap = sizeof(type_atom) - (BLOCK_SIZE * (rank + 1)) % sizeof(type_atom);
    if (left_this_block < buff_size) { // read more
        infile.seekg(next_rank_block, std::ios::beg);
#ifdef DEBUG
        std::cout << "(2) read from " << infile.tellg() << " of " << buff_size - left_this_block
                  << " bytes " << +".\n";
#endif
        infile.read(buff + left_this_block, gap); // todo fixme bug: can not combine.
        infile.read(buff + left_this_block + gap, buff_size - left_this_block - gap);
        infile.read(buff + left_this_block, buff_size - left_this_block);
    }
    if (infile.tellg() == next_block) { // automatic goes to next rank block.
        infile.seekg(next_rank_block, std::ios::beg);
    }
#ifdef DEBUG
    std::cout << "\n";
#endif
}

void checkoutHeaders(std::ifstream &infile, int rank) {
    infile.seekg(HEADER_SIZE + _LOCAL_HEADER_SIZE * rank, std::ios::beg);
}

void checkoutToAtoms(std::ifstream &infile, int n_ranks, int rank) {
    infile.seekg(HEADER_SIZE + _LOCAL_HEADER_SIZE * n_ranks + BLOCK_SIZE * rank, std::ios::beg);
}

std::string getNameByEleName(int type) {
    switch (type) {
        case 0:
            return "Fe";
        case 1:
            return "Cu";
        case 2:
            return "Ni";
        default:
            return "Unknown";
    }
}