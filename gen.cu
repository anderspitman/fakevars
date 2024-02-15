#include <stdint.h>
#include <stdio.h>
#include <curand_kernel.h>
#include <assert.h>
#include <inttypes.h>

#define cudaCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}


const uint64_t ALLELES_SIZE = 2;
const uint64_t BASE_SIZE = 2;
const uint64_t MIN_PHRED_QUAL = 0x21;
const uint64_t MAX_PHRED_QUAL = 0x7e;
const double ERROR_RATE = 0.02;
const uint64_t NUM_BASES = 4;
__device__ const char BASES[NUM_BASES] = { 'A', 'C', 'G', 'T' };

const uint64_t GENOTYPE_SIZE = 2;

__device__ uint64_t rand_in_range(uint64_t min, uint64_t max, curandState_t* state) {
        float rand = curand_uniform(state);
        return (rand * (max - min + 0.999999)) + min;
}

__global__ void init_rand_states(curandState* state, uint64_t seed) {
            int tid = threadIdx.x + blockIdx.x * blockDim.x;
            curand_init(seed, tid, 0, &state[tid]);
}

__device__ void gen_base_array_kernel(
        uint8_t* data,
        uint64_t depth,
        char* true_genotype,
        curandState* local_rand_state)
{
        //int tid = threadIdx.x + blockIdx.x * blockDim.x;
        //curandState local_rand_state = rand_state[tid];

        uint8_t* ptr;
        for (uint64_t i = 0; i < depth; i++) {
                ptr = (uint8_t*)&(data[i*BASE_SIZE]);
                float error = curand_uniform(local_rand_state);

                if (error < ERROR_RATE) {
                        // Simulate error by choosing random base
                        uint64_t idx = rand_in_range(0, NUM_BASES-1, local_rand_state);
                        ptr[0] = BASES[idx];
                }
                else {
                        // Sample randomly from true_genotype
                        uint64_t idx = rand_in_range(0, GENOTYPE_SIZE-1, local_rand_state);
                        //printf("%" PRIu64 "\n", idx);
                        ptr[0] = true_genotype[idx];
                }

                uint64_t qual = rand_in_range(MIN_PHRED_QUAL, MAX_PHRED_QUAL, local_rand_state);
                ptr[1] = qual;
        }
}

__device__ void gen_sample_kernel(
        uint8_t* data,
        uint64_t depth,
        char* locus_alleles,
        curandState* local_rand_state)
{

        uint64_t base_1_idx = rand_in_range(0, ALLELES_SIZE-1, local_rand_state);
        uint64_t base_2_idx = rand_in_range(0, ALLELES_SIZE-1, local_rand_state);

        // true_genotype
        data[0] = locus_alleles[base_1_idx];
        data[1] = locus_alleles[base_2_idx];

        gen_base_array_kernel(data + GENOTYPE_SIZE, depth, (char*)data, local_rand_state);
}

__device__ void gen_locus_kernel(
        uint8_t* data,
        uint64_t num_samples,
        uint64_t depth,
        curandState* local_rand_state)
{

        const uint64_t base_array_size = depth*BASE_SIZE;
        const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;

        auto base_1 = rand_in_range(0, NUM_BASES-1, local_rand_state);
        auto base_2 = rand_in_range(0, NUM_BASES-1, local_rand_state);

        // alleles
        data[0] = BASES[base_1];
        data[1] = BASES[base_2];

        uint8_t* sample_ptr = &data[ALLELES_SIZE];

        for (uint64_t i = 0; i < num_samples; i++) {
                gen_sample_kernel(sample_ptr, depth, (char*)data, local_rand_state);
                sample_ptr += sample_size;
        }
}

__global__ void gen_loci_kernel(
        uint8_t* data,
        uint64_t num_loci,
        uint64_t num_samples,
        uint64_t depth,
        curandState* rand_state)
{
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        curandState local_rand_state = rand_state[tid];

        const uint64_t base_array_size = depth*BASE_SIZE;
        const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
        const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;

        uint8_t* ptr = data;
        for (uint64_t i = 0; i < num_loci; i++) {
                gen_locus_kernel(ptr, num_samples, depth, &local_rand_state);
                ptr += locus_size;
        }
}


uint8_t* gen_data(uint64_t num_loci, uint64_t num_samples, uint64_t depth) {
        const uint64_t n_blocks = 256;
        const uint64_t n_threads = 256;
        const uint64_t seed = 9999;

        curandState* d_rand_states;
        cudaError_t err;

        err = cudaMalloc(&d_rand_states, n_blocks*n_threads*sizeof(curandState));
        assert(err == cudaSuccess);

        init_rand_states<<<n_blocks, n_threads>>>(d_rand_states, seed);
        cudaDeviceSynchronize();


        uint8_t* h_data;

        const uint64_t base_array_size = depth*BASE_SIZE;
        const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
        const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;
        const uint64_t total_size = locus_size*num_loci;

        const uint64_t size = total_size;

        err = cudaMallocHost(&h_data, size);
        assert(err == cudaSuccess);

        uint8_t* d_data;
        err = cudaMalloc(&d_data, size);
        assert(err == cudaSuccess);

        char h_locus_alleles[2] = { 'A', 'T' };

        char* d_locus_alleles;
        cudaMalloc(&d_locus_alleles, ALLELES_SIZE);
        cudaCheckError();
        cudaMemcpy(d_locus_alleles, h_locus_alleles, GENOTYPE_SIZE, cudaMemcpyHostToDevice);
        cudaCheckError();

        gen_loci_kernel<<<n_blocks, n_threads>>>(d_data, num_loci, num_samples, depth, d_rand_states);
        cudaDeviceSynchronize();

        err = cudaMemcpy(h_data, d_data, size, cudaMemcpyDeviceToHost);
        cudaCheckError();

        return h_data;
}


