
redfir: redfir.o
	$(CC) -o $@ redfir.o -lm

redfir.o: redfir.h

clean:
	rm -f redfir *.o

