CC = cc
LD     = $(CC)
FFLAGS = -g 
CFLAGS = -g -O2
DEFS= #-D__INLINE__=inline $(BVT) $(PCHECK)
LDFLAGS = $(CFLAGS)

LIBS=${MAINLIBS}

TARGET = aggregate
OBJS =  main.o getranks.o

all: grid


%.o : %.c
	$(CC) $(CFLAGS) $(DEFS) -c $<

grid: $(OBJS)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJS) $(MPITRACE_LIB) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS)

