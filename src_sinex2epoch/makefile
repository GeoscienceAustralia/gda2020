OBJS = sinex2epoch.o geodetic.o matrix.o time.o
CC = gcc
#CFLAGS = -Wall -Ofast
CFLAGS = -Wall -static -B/g/data/ms9/.envbin/glibc-static/lib -I${ENV_PATH}/glibc-static/include -L${ENV_PATH}/glibc-static/lib

sinex2epoch: $(OBJS) sinex2epoch.h
	$(CC) $(CFLAGS) $(OBJS) -o sinex2epoch -lm

sinex2epoch.o: sinex2epoch.c sinex2epoch.h
	$(CC) $(CFLAGS) -c sinex2epoch.c

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c matrix.c

geodetic.o: geodetic.c geodetic.h
	$(CC) $(CFLAGS) -c geodetic.c

time.o: time.c time.h
	$(CC) $(CFLAGS) -c time.c

