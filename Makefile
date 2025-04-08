CXX = g++
CXXFLAGS = -std=c++14 -Wall -Wextra -O2
LDFLAGS =

SRC_DIR = src
LIB_DIR = lib
BRKGA_DIR = $(LIB_DIR)/brkgaAPI
OBJ_DIR = obj
BIN_DIR = bin

$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

BRKGA_SRCS = $(BRKGA_DIR)/Population.cpp
BRKGA_OBJS = $(patsubst $(BRKGA_DIR)/%.cpp,$(OBJ_DIR)/brkga_%.o,$(notdir $(BRKGA_SRCS)))

TARGET = $(BIN_DIR)/program

INCLUDES = -I$(SRC_DIR) -I$(BRKGA_DIR)

all: $(TARGET)

$(TARGET): $(OBJS) $(BRKGA_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/brkga_%.o: $(BRKGA_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)/* $(BIN_DIR)/*

format:
	find $(SRC_DIR) -name "*.cpp" -o -name "*.h" | xargs clang-format -style=file -i
	find $(BRKGA_DIR) -name "*.cpp" -o -name "*.h" | xargs clang-format -style=file -i

tidy:
	find $(SRC_DIR) -name "*.cpp" | xargs clang-tidy -p=. -checks=*,-fuchsia-*,-google-*,-llvm-header-guard,-llvm-include-order,-hicpp-*,-misc-unused-parameters,-modernize-use-trailing-return-type,-cppcoreguidelines-avoid-magic-numbers,-readability-magic-numbers,-readability-else-after-return

check-syntax:
	$(CXX) -fsyntax-only $(CXXFLAGS) $(INCLUDES) $(CHK_SOURCES)

.PHONY: all clean format tidy check-syntax
