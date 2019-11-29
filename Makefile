
# Compiler and linker options
CC        = g++
CFLAGS    = -c -v
LDFLAGS   =

# User directories to search
BUILDDIR  = ./build
INCDIR    = ./include
SRCDIR    = ./src

# 3rd Party library directories to search
LIBDIR    = /usr/local/lib
LIBINCDIR = 

# List sources, object files and libraries to be used
SOURCES   = $(wildcard $(SRCDIR)/*.cpp )
OBJECTS   = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))
LIBRARIES = 

# Create command line arguments
LIBCMD    = $(addprefix -l,$(LIBRARIES))
LIBDIRCMD = $(addprefix -L,$(LIBDIR))
INCCMD    = $(addprefix -I,$(INCDIR) $(LIBINCDIR)) 

# The final executable name
TARGET    = ellipse

# The rules

all: $(SOURCES) $(TARGET)
    
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCCMD) $< -o $@

$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBDIRCMD) $(LIBCMD) -o $@


