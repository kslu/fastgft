CC=gcc
CFLAGS=-O2

all: gft_speed

gft_speed: txfm.o gsptools.o
	@echo "txfm.o"
	$(CC) -o $@ $?

txfm.o: txfm.c
	$(CC) -c $<

gsptools.o: gsptools.c
	$(CC) -c $<

clean:
		rm -rf *.o gft_speed

