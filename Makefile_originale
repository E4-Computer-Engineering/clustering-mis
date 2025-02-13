# Directories
SRC_DIR = .
KMEANS_LIB_DIR = vendor/kmeans
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

# Files
OBJS = $(OBJ_DIR)/kmeans.o $(OBJ_DIR)/main.o
EXE = $(BIN_DIR)/mykmeans

# Compilation flags
CFLAGS = -g -O0 -I$(KMEANS_LIB_DIR)
LIB = -lm

# Create build directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# Default target
all: $(EXE)

# Link executable
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

# Compile source files into build directory
$(OBJ_DIR)/kmeans.o: $(KMEANS_LIB_DIR)/kmeans.c $(KMEANS_LIB_DIR)/kmeans.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
