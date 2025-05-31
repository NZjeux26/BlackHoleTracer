SRC_DIR = src
OBJDIR = build
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJDIR)/%.o, $(SRC))
TARGET = BlackHoleTracer

CFLAGS = -g -Wall -Iinclude \
 -I/opt/homebrew/opt/sdl2/include/SDL2 \
 -I/opt/homebrew/Cellar/sdl2_image/2.8.8/include/SDL2 \
 -I/opt/homebrew/Cellar/glew/2.2.0_1/include

LDFLAGS = -L/opt/homebrew/opt/sdl2/lib \
 -L/opt/homebrew/Cellar/sdl2_image/2.8.8/lib \
 -L/opt/homebrew/Cellar/glew/2.2.0_1/lib \
 -lSDL2 -lSDL2_image -lGLEW -framework OpenGL

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
