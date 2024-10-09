
redfir: redfir.o
	$(CC) -o $@ redfir.o

redfir.o: redfir.h

clean:
	rm -f redfir *.o

