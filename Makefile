CC=g++
CFLAGS=-std=c++11 -O3 -fcilkplus -fopenmp
DEPS = params.h partition.h main.cc partition.cc

main: $(DEPS)
	$(CC) -o $@ $^ $(CFLAGS)
clean:
	rm -f main main_cpy $(OBJ)

