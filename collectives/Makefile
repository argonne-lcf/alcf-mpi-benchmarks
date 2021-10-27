CC=cc
LD = $(CC)
FFLAGS = -g 
CFLAGS = -g  -O2
DEFS= #-D__INLINE__=inline $(BVT) $(PCHECK)
LDFLAGS = $(CFLAGS) -shared-intel
LIBS=${MAINLIBS} -L/home/chunduri/variability/sst-dumpi/install/lib -ldumpi

TARGET = collectives
OBJS =  main.o split.o

all: grid


%.o : %.c
	$(CC) $(CFLAGS) $(DEFS) -c $<

grid: $(OBJS)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJS) $(MPITRACE_LIB) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS)


