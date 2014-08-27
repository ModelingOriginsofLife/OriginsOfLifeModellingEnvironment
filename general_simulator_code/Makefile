DISPLAYDEFINES = -DFB -DSCRNSHOT -DLINUX -DCURS -DMIKMOD 
#DISPLAYDEFINES = -DX11 -DSCRNSHOT -DLINUX
#DISPLAYDEFINES = -DSVGA -DSCRNSHOT -DMIKMOD
#DISPLAYDEFINES = -DSDL -DSCRNSHOT -DLINUX -DMIKMOD

MIKMODLIB = `libmikmod-config --libs`
MIKMODFLAGS = `libmikmod-config --cflags`
#MIKMODFLAGS = 
#MIKMODLIB = 

#LIBS = -L/usr/X11R6/lib/ -lX11 -lm -lncurses -lSDL -lpthread -lpng  $(MIKMODLIB)
#LIBS = -L/usr/X11R6/lib/ -lX11 -lm -lncurses -lpng  $(MIKMODLIB)
LIBS = -lm -lncurses -lpng -lz $(MIKMODLIB) 
#LIBS = -lvga -lm -lncurses -lpng   $(MIKMODLIB)

# Note- I've had some problems with gcc 4.0.4, so try with an earlier version
# if things don't work right
CFLAGS = -finline-functions -Wno-main $(MIKMODFLAGS) -O6  -g
CC = gcc $(DISPLAYDEFINES) $(CFLAGS)
ELDRITCHFLAGS = -o eldritch $(LIBS) 
MAPEDITFLAGS = -o mapedit $(LIBS)  
#CMAPEDITFLAGS = -o cmapedit $(LIBS) -O6 

EXTRAOBJ = 
include Makefile.core

clean: ; rm *.o ; rm eldritch ; rm mapedit 
