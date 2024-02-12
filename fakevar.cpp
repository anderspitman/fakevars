#include <stdint.h>
#include <cstddef>
#include <iostream>
#include <fstream>

const uint64_t BASE_SIZE = 2;
struct Base {
    uint8_t genotype;
    uint8_t qual;
};

const uint64_t TRUE_GENO_SIZE = 1;
struct Sample {
    uint8_t true_genotype;
    Base* bases;
};

struct Locus {
    Sample* samples;
};

struct Fakevar {
    uint64_t num_loci;
    uint64_t num_samples;
    uint32_t depth;
    uint8_t* data;
    uint64_t locus_offset;
    uint64_t size;
};


void gen_base_array(uint8_t* data, uint64_t depth) {

    uint8_t* ptr;
    for (uint64_t i = 0; i < depth; i++) {
        ptr = (uint8_t*)&(data[i*BASE_SIZE]);
        // genotype
        ptr[0] = 0x41;
        // qual
        ptr[1] = 0x42;
    }
}

void gen_sample(uint8_t* data, uint64_t depth) {
    data[0] = 0x43;
    gen_base_array(data + TRUE_GENO_SIZE, depth);
}

void gen_locus(uint8_t* data, uint64_t num_samples, uint64_t depth) {

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = TRUE_GENO_SIZE + base_array_size;

    uint8_t* ptr = data;

    for (uint64_t i = 0; i < num_samples; i++) {
        gen_sample(ptr, depth);
        ptr += sample_size;
    }
}

Fakevar* fakevar_create(uint64_t num_loci, uint64_t num_samples, uint32_t depth) {

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = TRUE_GENO_SIZE + base_array_size;
    const uint64_t locus_size = sample_size*num_samples;
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
        fv->data[i] = 0x55;
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
    std::cout << "genotype: " << data[0] << std::endl;
    do_indent(indent + 2);
    std::cout << "qual: " << data[1] << std::endl;
}

void print_sample(uint8_t* data, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "Sample:" << std::endl;
    do_indent(indent + 2);
    std::cout << "true_genotype: " << data[0] << std::endl;
    do_indent(indent + 2);
    std::cout << "bases:" << std::endl;

    for (uint64_t i = 0; i < depth; i++) {
        uint8_t* ptr = &(data[i*BASE_SIZE]);
        print_base(ptr, indent + 4);
    }
}

void print_locus(uint8_t* data, uint64_t num_samples, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "Locus:" << std::endl;
    do_indent(indent + 2);
    std::cout << "samples: " << std::endl;

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = TRUE_GENO_SIZE + base_array_size;

    for (uint64_t i = 0; i < num_samples; i++) {
        uint8_t* ptr = &(data[i*sample_size + TRUE_GENO_SIZE]);
        print_sample(ptr, depth, indent + 4);
    }
}

void print(uint8_t* data, uint64_t num_loci, uint64_t num_samples, uint64_t depth, uint64_t indent) {
    do_indent(indent);
    std::cout << "loci:" << std::endl;

    const uint64_t base_array_size = depth*BASE_SIZE;
    const uint64_t sample_size = TRUE_GENO_SIZE + base_array_size;
    const uint64_t locus_size = sample_size*num_samples;

    for (uint64_t i = 0; i < num_loci; i++) {
        uint8_t* ptr = &(data[i*locus_size]);
        print_locus(ptr, num_samples, depth, indent + 2);
    }
}


int main() {

    std::ofstream file("file.bin");

    const uint64_t num_loci = 2;
    const uint64_t num_samples = 2;
    const uint64_t depth = 2;
    auto fv = fakevar_create(num_loci, num_samples, depth);

    print(fv->data, num_loci, num_samples, depth, 0);

    file.write((char*)fv->data, fv->size);

    return 0;
}
