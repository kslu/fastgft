CC=gcc
CFLAGS=-O2

all: speed_sk15 speed_star10 speed_bd4x4 speed_dct4x4

speed_sk15: txfm.o speed_sk15.o
	$(CC) -o $@ $?

speed_star10: txfm.o speed_star10.o
	$(CC) -o $@ $?

speed_bd4x4: txfm.o speed_bd4x4.o
	$(CC) -o $@ $?

speed_dct4x4: txfm.o speed_dct4x4.o
	$(CC) -o $@ $?

txfm.o: txfm.c
	@echo "txfm.o"
	$(CC) -c $<

speed_sk15.o: speed_sk15.c
	@echo "speed_sk15.o"
	$(CC) -c $<

speed_star10.o: speed_star10.c
	@echo "speed_star10.o"
	$(CC) -c $<

speed_bd4x4.o: speed_bd4x4.c
	@echo "speed_bd4x4.o"
	$(CC) -c $<

speed_dct4x4.o: speed_dct4x4.c
	@echo "speed_dct4x4.o"
	$(CC) -c $<

clean:
		rm -rf *.o gft_speed

