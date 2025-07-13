SRC_DIR = src
OBJDIR = build
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJDIR)/%.o, $(SRC))
TARGET = BlackHoleTracer

# Homebrew libomp paths
OMP_INC = /opt/homebrew/opt/libomp/include
OMP_LIB = /opt/homebrew/opt/libomp/lib

CFLAGS = -g -Wall -Iinclude \
 -I/opt/homebrew/opt/sdl2/include/SDL2 \
 -I/opt/homebrew/Cellar/sdl2_image/2.8.8/include/SDL2 \
 -I/opt/homebrew/Cellar/glew/2.2.0_1/include \
 -Xclang -fopenmp -I$(OMP_INC)

LDFLAGS = -L/opt/homebrew/opt/sdl2/lib \
 -L/opt/homebrew/Cellar/sdl2_image/2.8.8/lib \
 -L/opt/homebrew/Cellar/glew/2.2.0_1/lib \
 -L$(OMP_LIB) -lomp \
 -lSDL2 -lSDL2_image -lGLEW -framework OpenGL

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
