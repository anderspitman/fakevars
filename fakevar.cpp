#include <stdint.h>
#include <cstddef>
#include <iostream>
#include <fstream>

struct File {
    struct Locus* loci;
    size_t num_loci;
    uint32_t depth;
    size_t num_samples;
};

struct Locus {
    struct Sample* samples;
};

struct Sample {
    uint8_t true_genotype;
    struct Base* bases;
};

struct Base {
    uint32_t qual;
    char c[3];
};

int main() {

    std::ofstream file("file.bin");

    Base base;
    //base.c = 'H';
    base.qual = 42;

    std::cout << sizeof(Base) << std::endl;

    file.write((char*)&base, sizeof(Base));

    return 0;
}
