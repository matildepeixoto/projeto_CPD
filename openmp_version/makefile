#make file - this is a comment section

CC=gcc  #compiler
TARGET=ballAlg #target file name

all:    ballAlg.o 
	$(CC) -g ballAlg.c -o $(TARGET) -fopenmp -lm 

clean:
	rm *.o $(TARGET)