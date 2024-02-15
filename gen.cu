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
        uint8_t* all_data,
        uint64_t num_loci,
        uint64_t num_samples,
        uint64_t depth,
        uint64_t rand_seed,
        curandState* rand_state)
{
        uint64_t tid = threadIdx.x + blockIdx.x * blockDim.x;


        if (tid > num_loci - 1) {
                //printf("ditching tid: %llu\n", tid);
                return;
        }

        curand_init(rand_seed, tid, 0, &rand_state[tid]);

        curandState local_rand_state = rand_state[tid];

        const uint64_t base_array_size = depth*BASE_SIZE;
        const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
        const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;

        const uint64_t data_idx = tid*locus_size;
        uint8_t* data = &all_data[data_idx];

        for (uint64_t i = 0; i < locus_size; i++) {
                data[i] = 232;
        }

        auto base_1 = rand_in_range(0, NUM_BASES-1, &local_rand_state);
        auto base_2 = rand_in_range(0, NUM_BASES-1, &local_rand_state);

        // alleles
        data[0] = BASES[base_1];
        data[1] = BASES[base_2];

        uint8_t* sample_ptr = &data[ALLELES_SIZE];

        for (uint64_t i = 0; i < num_samples; i++) {
                gen_sample_kernel(sample_ptr, depth, (char*)data, &local_rand_state);
                sample_ptr += sample_size;
        }
}

void gen_data(uint64_t num_loci, uint64_t num_samples, uint64_t depth, uint8_t** d_data, uint64_t *size) {

        const uint64_t n_blocks = ceil(num_loci/256.0);
        const uint64_t n_threads = 256;
        const uint64_t seed = time(NULL);

        curandState* d_rand_states;

        cudaMalloc(&d_rand_states, num_loci*sizeof(curandState));
        cudaCheckError();

        const uint64_t base_array_size = depth*BASE_SIZE;
        const uint64_t sample_size = GENOTYPE_SIZE + base_array_size;
        const uint64_t locus_size = ALLELES_SIZE + sample_size*num_samples;
        const uint64_t total_size = locus_size*num_loci;

        *size = total_size;

        cudaMalloc(d_data, *size);
        cudaCheckError();

        gen_loci_kernel<<<n_blocks, n_threads>>>(*d_data, num_loci, num_samples, depth, seed, d_rand_states);
        cudaDeviceSynchronize();
}


void gen_data_gpu(uint64_t num_loci, uint64_t num_samples, uint64_t depth, uint8_t** h_data, uint64_t *size) {

        uint8_t* d_data;

        gen_data(num_loci, num_samples, depth, &d_data, size);

        cudaMallocHost(h_data, *size);
        cudaCheckError();

        cudaMemcpy(*h_data, d_data, *size, cudaMemcpyDeviceToHost);
        cudaCheckError();
}
