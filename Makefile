# Compilador e flags
CXX = g++
CXXFLAGS = -std=c++14 -Wall -Wextra -O3 -fopenmp

# Diretórios
SRC_DIR = src
LIB_DIR = lib
BRKGA_DIR = $(LIB_DIR)/brkgaAPI
OBJ_DIR = obj
BIN_DIR = bin

# Criar diretórios se não existirem
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# Arquivos fonte
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
BRKGA_SRCS = $(BRKGA_DIR)/Population.cpp

# Objetos gerados
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
BRKGA_OBJS = $(OBJ_DIR)/Population.o

# Programa final
TARGET = $(BIN_DIR)/program

# Flags de inclusão
INCLUDES = -I$(SRC_DIR) -I$(BRKGA_DIR)

# Regra padrão
all: $(TARGET)

# Compilar o programa
$(TARGET): $(OBJS) $(BRKGA_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilar arquivos do src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compilar arquivos do BRKGA
$(OBJ_DIR)/Population.o: $(BRKGA_DIR)/Population.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Limpar arquivos temporários
clean:
	rm -rf $(OBJ_DIR)/* $(BIN_DIR)/*

# Formatar código
format:
	find $(SRC_DIR) -name "*.cpp" -o -name "*.h" -o -name "*.hpp" | xargs clang-format -style=file -i
	find $(BRKGA_DIR) -name "*.cpp" -o -name "*.h" | xargs clang-format -style=file -i

# Verificar código com clang-tidy
tidy:
	find $(SRC_DIR) -name "*.cpp" | xargs clang-tidy -p=. -checks=*,-fuchsia-*,-google-*,-llvm-header-guard,-llvm-include-order,-hicpp-*,-misc-unused-parameters,-modernize-use-trailing-return-type,-cppcoreguidelines-avoid-magic-numbers,-readability-magic-numbers,-readability-else-after-return

# Verificar sintaxe
check-syntax:
	$(CXX) -fsyntax-only $(CXXFLAGS) $(INCLUDES) $(CHK_SOURCES)

.PHONY: all clean format tidy check-syntax
