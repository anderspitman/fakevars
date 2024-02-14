#include <stdint.h>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <random>

//enum class Genotype : uint8_t {
//    AA = 0x30, AC, AG, AT,
//    CC, CG, CT,
//    GG, GT,
//    TT,
//};

//const int NUM_GENOTYPES = 10;
//const char* GENOTYPE_TABLE[NUM_GENOTYPES] = {
//    "A/A", "A/C", "A/G", "A/T",
//    "C/C", "C/G", "C/T",
//    "G/G", "G/T",
//    "T/T",
//};

const uint64_t NUM_BASES = 4;
const char BASES[NUM_BASES] = { 'A', 'C', 'G', 'T' };

const uint64_t GENOTYPE_SIZE = 2;

const uint64_t BASE_SIZE = 2;
struct Base {
    char val;
    uint8_t qual;
};

struct Sample {
    char true_genotype[GENOTYPE_SIZE];
    Base* bases;
};

const uint64_t ALLELES_SIZE = 2;
struct Locus {
    char alleles[ALLELES_SIZE];
    Sample* samples;
};

struct AllData {
    Locus* loci;
};

struct Fakevar {
    uint64_t num_loci;
    uint64_t num_samples;
    uint32_t depth;
    uint8_t* data;
    uint64_t locus_offset;
    uint64_t size;
};

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> base_dist(0, NUM_BASES-1);
std::uniform_int_distribution<std::mt19937::result_type> geno_dist(0, GENOTYPE_SIZE-1);
std::uniform_int_distribution<std::mt19937::result_type> allele_dist(0, ALLELES_SIZE-1);
const uint64_t MIN_PHRED_QUAL = 0x21;
const uint64_t MAX_PHRED_QUAL = 0x7e;
std::uniform_int_distribution<std::mt19937::result_type> qual_dist(MIN_PHRED_QUAL, MAX_PHRED_QUAL);
std::uniform_real_distribution<> error_dist(0.0, 1.0);

const double ERROR_RATE = 0.02;


void gen_base_array(uint8_t* data, uint64_t depth, char* true_genotype) {

    uint8_t* ptr;

    for (uint64_t i = 0; i < depth; i++) {

        ptr = (uint8_t*)&(data[i*BASE_SIZE]);

        double error = error_dist(rng);
        if (error < ERROR_RATE) {
            // Simulate error by choosing random base
            auto idx = base_dist(rng);
            ptr[0] = BASES[idx];
        }
        else {
            // Sample randomly from true_genotype
            auto idx = geno_dist(rng);
            ptr[0] = true_genotype[idx];
        }

        auto qual = qual_dist(rng);
        ptr[1] = qual;
    }
}

void gen_sample(uint8_t* data, uint64_t depth, char* locus_alleles) {

    auto base_1 = allele_dist(rng);
    auto base_2 = allele_dist(rng);

    // true_genotype
    data[0] = locus_alleles[base_1];
    data[1] = locus_alleles[base_2];
    gen_base_array(data + GENOTYPE_SIZE, depth, (char*)data);
}

void gen_locus(uint8_t* data, uint64_t num_samples, uint64_t depth) {

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;

    auto base_1 = base_dist(rng);
    auto base_2 = base_dist(rng);

    // alleles
    data[0] = BASES[base_1];
    data[1] = BASES[base_2];

    uint8_t* sample_ptr = &data[ALLELES_SIZE];

    for (uint64_t i = 0; i < num_samples; i++) {
        gen_sample(sample_ptr, depth, (char*)data);
        sample_ptr += sample_size;
    }
}

Fakevar* fakevar_create(uint64_t num_loci, uint64_t num_samples, uint32_t depth) {

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
    const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;
    const uint64_t total_size = locus_size*num_loci;

    std::cout << "base_array_size: " << base_array_size << std::endl;
    std::cout << "sample_size: " << sample_size << std::endl;
    std::cout << "locus_size: " << locus_size << std::endl;
    std::cout << "total_size: " << total_size << std::endl;

    Fakevar* fv = (Fakevar*)malloc(sizeof(Fakevar));
    fv->num_loci = num_loci;
    fv->num_samples = num_samples;
    fv->depth = depth;
    fv->data = (uint8_t*)malloc(total_size);
    fv->locus_offset = 0;
    fv->size = total_size;

    for (uint64_t i = 0; i < total_size; i++) {
        fv->data[i] = 0x58;
    }

    uint8_t* ptr = fv->data;
    for (uint64_t i = 0; i < num_loci; i++) {
        gen_locus(ptr, num_samples, depth);
        ptr += locus_size;
    }

    return fv;
}

void do_indent(uint64_t indent) {
    for (uint64_t i = 0; i < indent; i++) {
        std::cout << " ";
    }
}

void print_base(uint8_t* data, uint64_t indent) {
    do_indent(indent);
    std::cout << "Base:" << std::endl;
    do_indent(indent + 2);
    std::cout << "value: " << data[0] << std::endl;
    do_indent(indent + 2);
    std::cout << "qual: " << data[1] << std::endl;
}

void print_sample(uint8_t* data, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "Sample:" << std::endl;
    do_indent(indent + 2);
    std::cout << "true_genotype: " << data[0] << "/" << data[1] << std::endl;

    do_indent(indent + 2);
    std::cout << "bases:" << std::endl;
    do_indent(indent + 4);
    for (uint64_t i = 0; i < depth; i++) {
        uint8_t* ptr = &(data[i*BASE_SIZE + GENOTYPE_SIZE]);
        std::cout << ptr[0];
        //print_base(ptr, indent + 4);
    }
    std::cout << std::endl;

    do_indent(indent + 2);
    std::cout << "quals:" << std::endl;
    do_indent(indent + 4);
    for (uint64_t i = 0; i < depth; i++) {
        uint8_t* ptr = &(data[i*BASE_SIZE + GENOTYPE_SIZE]);
        std::cout << ptr[1];
    }
    std::cout << std::endl;
}

void print_locus(uint8_t* data, uint64_t num_samples, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "Locus:" << std::endl;
    do_indent(indent + 2);
    std::cout << "alleles: " << data[0] << "/" << data[1] << std::endl;
    do_indent(indent + 2);
    std::cout << "samples: " << std::endl;

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;

    for (uint64_t i = 0; i < num_samples; i++) {
        uint8_t* ptr = &(data[ALLELES_SIZE + i*sample_size]);
        print_sample(ptr, depth, indent + 4);
    }
}

void print(uint8_t* data, uint64_t num_loci, uint64_t num_samples, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "loci:" << std::endl;

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
    const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;

    for (uint64_t i = 0; i < num_loci; i++) {
        uint8_t* ptr = &(data[i*locus_size]);
        print_locus(ptr, num_samples, depth, indent + 2);
    }
}


int main() {

    std::ofstream file("file.bin");

    const uint64_t num_loci = 1;
    const uint64_t num_samples = 2;
    const uint64_t depth = 60;
    auto fv = fakevar_create(num_loci, num_samples, depth);

    print(fv->data, num_loci, num_samples, depth, 0);

    file.write((char*)fv->data, fv->size);

    return 0;
}
