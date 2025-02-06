OBJS=kmeans.o
EXE=mykmeans

CFLAGS=-g -O0
LIB = -lm

all: $(EXE)

mykmeans: $(OBJS) main.o
	$(CC) $(CFLAGS) $(LIB) $^ -o $@

clean:
	rm -f $(EXE) *.mod *.o *.txt

