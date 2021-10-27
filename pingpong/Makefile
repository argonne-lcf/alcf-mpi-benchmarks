CC  = cc
LD     = $(CC)

# Do not change -O2!
FFLAGS = -g -O2
CFLAGS = -g -O2
DEFS= -D__INLINE__=inline $(BVT) $(PCHECK)
LDFLAGS = $(CFLAGS)

TARGET = pingpong
OBJS =  main.o getranks.o 

all: grid


%.o : %.c
	$(CC) $(CFLAGS) $(DEFS) -c $<

grid: $(OBJS)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJS) $(MPITRACE_LIB)

clean:
	rm -f $(TARGET) $(OBJS)

