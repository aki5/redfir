
redfir: redfir.o
	$(CC) -o $@ redfir.o

clean:
	rm -f redfir *.o

