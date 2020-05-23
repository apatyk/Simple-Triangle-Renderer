CC = gcc
CFLAGS = -Wall -O4 -lm

BINS = render

all: $(BINS)

$(BINS):  $(BINS).c
	$(CC) $(BINS).c $(CFLAGS) -o $(BINS)

style:
	astyle --style=java --break-blocks --pad-oper --pad-header --align-pointer=name --delete-empty-lines *.c

clean:
	rm $(BINS)

cleano:
	rm *.orig
