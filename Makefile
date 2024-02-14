fakevar: fakevar.cpp
	g++ -o fakevar fakevar.cpp

run: fakevar
	./fakevar && xxd file.bin && ls -l
