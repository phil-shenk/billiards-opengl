CC       = gcc 
CFLAGS   = -O3 -Wall 
LIBS     = -lXi -lXmu -lglut -lGLEW -lGLU -lm -lGL
WLIBS    = -lXi -lXmu -lglut -lGLEW -lglu32 -lm -lopengl32
OBJDIR   = ../linalglib
OBJS     = initShader.o $(OBJDIR)/linalglib.o

billiards: billiards.c $(OBJS)
	$(CC) -o billiards billiards.c $(OBJS) $(CFLAGS) $(LIBS)

$(OBJDIR)/%.o: %.c %.h
	$(CC) -c @< -o $@ $(CFLAGS)

