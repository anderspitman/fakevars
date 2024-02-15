NVCC=nvcc

fakevar: fakevar.cpp
	g++ -o fakevar fakevar.cpp

#run: fakevar
#	./fakevar && xxd file.bin && ls -l

fakevarnv: gen.cu fakevar.cpp
	NVCC_PREPEND_FLAGS='-ccbin g++-12' $(NVCC) -Xcompiler -fopenmp -g -o $@ $^

run: fakevarnv
	./fakevarnv
