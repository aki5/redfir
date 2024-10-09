
redfir_test: redfir_test.o redfir.o
	$(CC) -o $@ redfir_test.o redfir.o -lm

redfir_test.o redfir.o: redfir.h

clean:
	rm -f redfir_test *.o *.txt

