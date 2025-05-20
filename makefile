# Compiler and flags
CC = clang
CFLAGS = -g -Wall \
 -I/opt/homebrew/opt/sdl2/include/SDL2 \
 -I/opt/homebrew/Cellar/sdl2_image/2.8.8/include/SDL2

LDFLAGS = -L/opt/homebrew/opt/sdl2/lib \
 -L/opt/homebrew/Cellar/sdl2_image/2.8.8/lib \
-lSDL2 -lSDL2_image

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
