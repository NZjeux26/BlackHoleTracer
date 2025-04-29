# Compiler and flags
CC = clang
CFLAGS = -g -Wall -I/opt/homebrew/opt/sdl2/include -I/opt/homebrew/Cellar/sdl2/2.32.4/include
LDFLAGS = -L/opt/homebrew/opt/sdl2/lib -lSDL2

# Source files and target
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
TARGET = BlackHoleTracer

# Build rule
$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Compile .c to .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)
