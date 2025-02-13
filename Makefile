# Directories
SRC_DIR = .
INC_DIR = ./include
KMEANS_LIB_DIR = vendor/kmeans
DBSCAN_LIB_DIR = vendor/dbscan
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

LIB = -lmpi
# Files
OBJS = $(OBJ_DIR)/kmeans.o $(OBJ_DIR)/dbscan.o $(OBJ_DIR)/clustering.o
EXE = $(BIN_DIR)/clustering

# Compilation flags
CFLAGS = -O3 -I$(INC_DIR) -I$(MPI_INC) -L$(MPI_LIB) -I$(KMEANS_LIB_DIR) -I$(DBSCAN_LIB_DIR)

# Create build directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# Default target
all: $(EXE)

# Link executable
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

#Compile source files into build directory
$(OBJ_DIR)/kmeans.o: $(KMEANS_LIB_DIR)/kmeans.c $(KMEANS_LIB_DIR)/kmeans.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/dbscan.o: $(DBSCAN_LIB_DIR)/dbscan.c $(DBSCAN_LIB_DIR)/dbscan.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/clustering.o: $(SRC_DIR)/clustering.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
