# Directories
SRC_DIR = .
INC_DIR = ./include
KMEANS_LIB_DIR = vendor/kmeans
DBSCAN_LIB_DIR = vendor/dbscan
HIERARCHICAL_LIB_DIR = vendor/hclust
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

LIB = -lmpi -lm
# Files
CLUS_OBJS = $(OBJ_DIR)/common.o $(OBJ_DIR)/kmeans_cl.o $(OBJ_DIR)/kmeans.o $(OBJ_DIR)/dbscan.o $(OBJ_DIR)/clustering.o $(OBJ_DIR)/fastcluster.o $(OBJ_DIR)/hclust.o
CLUS_EXE = $(BIN_DIR)/clustering
SIL_OBJS = $(OBJ_DIR)/common.o $(OBJ_DIR)/silhouette.o
SIL_EXE = $(BIN_DIR)/silhouette

# Compilation flags
CFLAGS = -O3 -I$(INC_DIR) -I$(MPI_INC) -L$(MPI_LIB) -I$(KMEANS_LIB_DIR) -I$(DBSCAN_LIB_DIR)
CXXFLAGS = -O3 -I$(INC_DIR) -I$(MPI_INC) -L$(MPI_LIB) -I$(KMEANS_LIB_DIR) -I$(DBSCAN_LIB_DIR) -I$(HIERARCHICAL_LIB_DIR) -std=c++20

# Create build directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# Default target
all: $(CLUS_EXE) $(SIL_EXE)

# Link executables
$(CLUS_EXE): $(CLUS_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

$(SIL_EXE): $(SIL_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

#Compile source files into build directory
$(OBJ_DIR)/common.o: $(SRC_DIR)/common.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/kmeans_cl.o: $(KMEANS_LIB_DIR)/kmeans_cl.c $(KMEANS_LIB_DIR)/kmeans_cl.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/kmeans.o: $(KMEANS_LIB_DIR)/kmeans.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/dbscan.o: $(DBSCAN_LIB_DIR)/dbscan.cpp $(DBSCAN_LIB_DIR)/dbscan.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/fastcluster.o: $(HIERARCHICAL_LIB_DIR)/fastcluster.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/hclust.o: $(HIERARCHICAL_LIB_DIR)/hclust.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/clustering.o: $(SRC_DIR)/clustering.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/silhouette.o: $(SRC_DIR)/silhouette.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
