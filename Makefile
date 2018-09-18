CC=gcc
CFLAGS=-O2

all: speed_sk15 speed_sk25 speed_star10 speed_star100 speed_bd4x4 speed_bd8x8 speed_dct4x4 speed_dct8x8 speed_z4x4 speed_z8x8

speed_sk15: txfm.o speed_sk15.o
	$(CC) -o $@ $?

speed_sk25: txfm.o speed_sk25.o
	$(CC) -o $@ $?

speed_star10: txfm.o speed_star10.o
	$(CC) -o $@ $?

speed_star100: txfm.o speed_star100.o
	$(CC) -o $@ $?

speed_bd4x4: txfm.o speed_bd4x4.o
	$(CC) -o $@ $?

speed_bd8x8: txfm.o speed_bd8x8.o
	$(CC) -o $@ $?

speed_dct4x4: txfm.o speed_dct4x4.o
	$(CC) -o $@ $?

speed_dct8x8: txfm.o speed_dct8x8.o
	$(CC) -o $@ $?

speed_z4x4: txfm.o speed_z4x4.o
	$(CC) -o $@ $?

speed_z8x8: txfm.o speed_z8x8.o
	$(CC) -o $@ $?

txfm.o: txfm.c
	@echo "txfm.o"
	$(CC) -c $<

speed_sk15.o: speed_sk15.c
	@echo "speed_sk15.o"
	$(CC) -c $<

speed_sk25.o: speed_sk25.c
	@echo "speed_sk25.o"
	$(CC) -c $<

speed_star10.o: speed_star10.c
	@echo "speed_star10.o"
	$(CC) -c $<

speed_star100.o: speed_star100.c
	@echo "speed_star100.o"
	$(CC) -c $<

speed_bd4x4.o: speed_bd4x4.c
	@echo "speed_bd4x4.o"
	$(CC) -c $<

speed_bd8x8.o: speed_bd8x8.c
	@echo "speed_bd8x8.o"
	$(CC) -c $<

speed_dct4x4.o: speed_dct4x4.c
	@echo "speed_dct4x4.o"
	$(CC) -c $<

speed_dct8x8.o: speed_dct8x8.c
	@echo "speed_dct8x8.o"
	$(CC) -c $<

speed_z4x4.o: speed_z4x4.c
	@echo "speed_z4x4.o"
	$(CC) -c $<

speed_z8x8.o: speed_z8x8.c
	@echo "speed_z8x8.o"
	$(CC) -c $<

clean:
		rm -rf *.o gft_speed

