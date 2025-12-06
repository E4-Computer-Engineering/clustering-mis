# Directories
SRC_DIR = .
INC_DIR = ./include
KMEANS_LIB_DIR = vendor/kmeans
DBSCAN_LIB_DIR = vendor/dbscan
HIERARCHICAL_LIB_DIR = vendor/hclust
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

#LIB = -lmpi -lm
CXX = mpicc # no mpicxx
CC = mpicc
PAPI = -I/home/agarcia/local/include/ -L/home/agarcia/local/lib/ -lpapi
GLPK = -I/home/agarcia/local/include/ -L/home/agarcia/local/lib/ -lglpk
IC = -I/home/agarcia/local/include -L/home/agarcia/local/lib -licc
FMPI = -I/home/agarcia/local/include/empi -L/home/agarcia/local/lib -lempi
#FMPI = -I/home/agarcia/FlexMPI/include -L/home/agarcia/FlexMPI/lib -lempi
ADMDEPS = $(PAPI) $(GLPK) $(IC) $(FMPI)
LIB = $(ADMDEPS) -lm -lstdc++

# Files
OBJS = $(OBJ_DIR)/points.o $(OBJ_DIR)/kmeans_cl.o $(OBJ_DIR)/kmeans.o $(OBJ_DIR)/dbscan.o $(OBJ_DIR)/clustering.o $(OBJ_DIR)/fastcluster.o $(OBJ_DIR)/hclust.o
EXE = $(BIN_DIR)/clustering

# Compilation flags
CFLAGS = -O3 -I$(INC_DIR) -I$(MPI_INC) -L$(MPI_LIB) -I$(KMEANS_LIB_DIR) -I$(DBSCAN_LIB_DIR) 
CXXFLAGS = -O3 -I$(INC_DIR) -I$(MPI_INC) -L$(MPI_LIB) -I$(KMEANS_LIB_DIR) -I$(DBSCAN_LIB_DIR) -I$(HIERARCHICAL_LIB_DIR) -std=c++20 

# Create build directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# Default target
all: $(EXE)

# Link executable
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIB) $^ -o $@ 

#Compile source files into build directory
$(OBJ_DIR)/points.o: $(SRC_DIR)/points.cpp
	$(CXX) $(CXXFLAGS) $(LIB) -c $< -o $@

$(OBJ_DIR)/kmeans_cl.o: $(KMEANS_LIB_DIR)/kmeans_cl.c $(KMEANS_LIB_DIR)/kmeans_cl.h
	$(CC) $(CFLAGS)  $(LIB) -c $< -o $@

$(OBJ_DIR)/kmeans.o: $(KMEANS_LIB_DIR)/kmeans.c
	$(CC) $(CFLAGS) $(LIB) -c $< -o $@

$(OBJ_DIR)/dbscan.o: $(DBSCAN_LIB_DIR)/dbscan.cpp $(DBSCAN_LIB_DIR)/dbscan.hpp
	$(CXX) $(CXXFLAGS) $(LIB) -c $< -o $@

$(OBJ_DIR)/fastcluster.o: $(HIERARCHICAL_LIB_DIR)/fastcluster.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/hclust.o: $(HIERARCHICAL_LIB_DIR)/hclust.cpp
	$(CXX) $(CXXFLAGS) $(LIB) -c $< -o $@

$(OBJ_DIR)/clustering.o: $(SRC_DIR)/clustering.cpp
	$(CXX) $(CXXFLAGS) $(LIB) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
