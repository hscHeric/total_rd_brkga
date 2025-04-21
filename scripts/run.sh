#!/bin/bash

# Diretórios
INPUT_DIR="data/edges/"
OUTPUT_DIR="data/results/"

# Parâmetros do algoritmo BRKGA
POP=2.5427
ELITE_SET=0.4835
MUTANTS=0.1465
INHERIT_ELITE=0.6712
POPULATIONS=30
EXCHANGE_INTERVAL=256
EXCHANGE_NUMBER=1
MAX_GENERATIONS=1573
STAGNATION_LIMIT=413
TRIALS=10
PARALLEL=8

# Criar o diretório de resultados, se não existir
mkdir -p "$OUTPUT_DIR"

# Estatísticas
total_graphs=0
skipped_graphs=0
processed_graphs=0
start_total_time=$(date +%s)

# Coletar arquivos e ordená-los por número de linhas (tamanho do grafo)
declare -a files
while IFS= read -r graph_file; do
  num_lines=$(wc -l <"$graph_file")
  files+=("$num_lines:$graph_file")
done < <(find "$INPUT_DIR" -type f -name "*.txt")

IFS=$'\n' sorted_files=($(sort -n <<<"${files[*]}"))
unset IFS

# Processar os grafos
for entry in "${sorted_files[@]}"; do
  num_lines=${entry%%:*}
  graph_file=${entry#*:}

  relative_path="${graph_file#$INPUT_DIR}"
  graph_name=$(basename "$graph_file" .txt)
  subdir=$(dirname "$relative_path")
  output_subdir="$OUTPUT_DIR$subdir"
  mkdir -p "$output_subdir"

  output_file="$output_subdir/${graph_name}.csv"

  if [ -f "$output_file" ]; then
    echo "[INFO] Resultado já existe para: $graph_name. Pulando..."
    skipped_graphs=$((skipped_graphs + 1))
    continue
  fi

  echo "[INFO] Processando: $graph_name"
  echo "       Linhas: $num_lines | Subpasta: $subdir"

  start_time=$(date +%s)

  # Executar seu programa
  ./bin/cl_total_rdbrkga "$graph_file" \
    --population "$POP" \
    --elite-set "$ELITE_SET" \
    --mutants "$MUTANTS" \
    --inherit-elite-probability "$INHERIT_ELITE" \
    --populations "$POPULATIONS" \
    --exchange-interval "$EXCHANGE_INTERVAL" \
    --exchange-number "$EXCHANGE_NUMBER" \
    --max-generations "$MAX_GENERATIONS" \
    --stagnation-limit "$STAGNATION_LIMIT" \
    --trials "$TRIALS" \
    --output "$output_file" \
    --parallel "$PARALLEL"

  end_time=$(date +%s)
  elapsed=$((end_time - start_time))

  echo "[INFO] Tempo: ${elapsed}s | Resultado: $output_file"
  echo "----------------------------------------"

  total_graphs=$((total_graphs + 1))
  processed_graphs=$((processed_graphs + 1))
done

# Estatísticas finais
end_total_time=$(date +%s)
total_time=$((end_total_time - start_total_time))

echo "----------------------------------------"
echo "[INFO] FIM DA EXECUÇÃO"
echo "       Processados: $processed_graphs"
echo "       Pulados: $skipped_graphs"
echo "       Tempo total: ${total_time}s"
echo "Resultados estão em: $OUTPUT_DIR"
echo "----------------------------------------"
