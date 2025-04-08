#include "../lib/brkgaAPI/BRKGA.h"
#include "../lib/brkgaAPI/MTRand.h"
#include "Decoder.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"

#include <cfloat>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>

#define DEBUG

#if defined( DEBUG )
  #define DEBUG_PRINT( msg )     std::cout << "DEBUG: " << msg << std::endl;
  #define DEBUG_PRINT_VAR( var ) std::cout << "DEBUG: " #var " = " << var << std::endl;
#else
  #define DEBUG_PRINT( msg )
  #define DEBUG_PRINT_VAR( var )
#endif

struct AlgorithmParams {
  unsigned      n           = 0;      // Size of chromosomes
  unsigned      p           = 1000;   // Size of population
  double        pe          = 0.20;   // Fraction of population to be the elite-set
  double        pm          = 0.057;  // Fraction of population to be replaced by mutants
  double        rhoe        = 0.70;   // Probability offspring inherits elite parent allele
  unsigned      K           = 3;      // Number of independent populations
  unsigned      MAXT        = 2;      // Number of threads for parallel decoding
  unsigned      x_intvl     = 100;    // Exchange best individuals every x_intvl generations
  unsigned      x_number    = 2;      // Number of best individuals to exchange
  unsigned      max_gens    = 586;    // Maximum number of generations
  unsigned      stag_limit  = 363;    // Stagnation limit (number of generations without improvement)
  long unsigned rng_seed    = 0;      // Random number generator seed (0 means generate random seed)
  bool          random_seed = true;   // Flag to indicate if we should use a random seed
  std::string   file_path;            // File path for the input file
  std::string   output_file = "results.csv";
  unsigned      trials      = 1;
};

struct TrialResult {
  std::string graph_name;
  size_t      node_count;
  size_t      edge_count;
  float       graph_density;
  int         fitness;
  uint64_t    elapsed_micros;
  bool        is_dense;
};

AlgorithmParams parse_args( int argc, char * argv[] );
void            load_and_normalize_graph( const std::string & filename, bool & graph_is_matrix, ListGraph & graph_l, MatrixGraph & graph_m );
long unsigned   generate_random_seed();
void            write_to_csv( const std::string & filename, const std::vector<TrialResult> & results );
std::string     get_filename_from_path( const std::string & path );

long unsigned generate_random_seed() {
  std::random_device                           rd;
  std::mt19937_64                              gen( rd() );
  std::uniform_int_distribution<long unsigned> dis;
  return dis( gen );
}

int main( int argc, char * argv[] ) {
  auto params = parse_args( argc, argv );

  double BEST = DBL_MAX;

  ListGraph   graph_l( 1 );
  MatrixGraph graph_m( 1 );

  bool graph_is_matrix = false;
  load_and_normalize_graph( params.file_path, graph_is_matrix, graph_l, graph_m );

  // Determinar a semente RNG a ser usada
  long unsigned seed_to_use = params.random_seed ? generate_random_seed() : params.rng_seed;
  DEBUG_PRINT( "Usando semente RNG: " << seed_to_use );

  MTRand rng( seed_to_use );

  // Vetor para armazenar os resultados de cada trial
  std::vector<TrialResult> all_results;

  // Laço para os diferentes trials
  for ( unsigned trial = 0; trial < params.trials; ++trial ) {
    TrialResult result;

    if ( graph_is_matrix ) {
      params.n = graph_m.order();

      TRDDecoder<MatrixGraph> decoder( graph_m );
      DEBUG_PRINT( "Executando BRKGA com MatrixGraph" );

      BRKGA<TRDDecoder<MatrixGraph>, MTRand> algorithm( params.n, params.p, params.pe, params.pm, params.rhoe, decoder, rng, params.K, params.MAXT );

      unsigned generation       = 1;
      unsigned stagnation_count = 0;
      double   best_fitness     = std::numeric_limits<double>::infinity();

      auto start_time = std::chrono::high_resolution_clock::now();
      do {
        algorithm.evolve();

        double current_best = algorithm.getBestFitness();
        if ( current_best < best_fitness ) {
          best_fitness     = current_best;
          stagnation_count = 0;
          DEBUG_PRINT( "Nova melhor solução encontrada: " << best_fitness );
        } else {
          stagnation_count++;
          if ( stagnation_count >= params.stag_limit ) {
            DEBUG_PRINT( "Critério de estagnação atingido após " << stagnation_count << " gerações sem melhoria" );
            break;
          }
        }

        if ( ( ++generation ) % params.x_intvl == 0 ) {
          algorithm.exchangeElite( params.x_number );
          DEBUG_PRINT( "Troca de elites realizada na geração " << generation );
        }
      } while ( generation < params.max_gens );

      auto end_time = std::chrono::high_resolution_clock::now();
      auto elapsed  = std::chrono::duration_cast<std::chrono::microseconds>( end_time - start_time );

      DEBUG_PRINT( "Melhor solução encontrada tem valor objetivo = " << algorithm.getBestFitness() );
      DEBUG_PRINT( "Número total de gerações executadas: " << generation );

      BEST = algorithm.getBestFitness();

      result.graph_name     = get_filename_from_path( params.file_path );
      result.node_count     = graph_m.order();
      result.edge_count     = graph_m.size();
      result.fitness        = static_cast<int>( BEST );
      result.elapsed_micros = elapsed.count();
      result.is_dense       = ( graph_m.size() / ( graph_m.order() * ( graph_m.order() - 1 ) / 2.0 ) ) > 0.5;

    } else {
      params.n = graph_l.order();

      TRDDecoder<ListGraph> decoder( graph_l );
      DEBUG_PRINT( "Executando BRKGA com ListGraph" );

      BRKGA<TRDDecoder<ListGraph>, MTRand> algorithm( params.n, params.p, params.pe, params.pm, params.rhoe, decoder, rng, params.K, params.MAXT );

      unsigned generation       = 1;
      unsigned stagnation_count = 0;
      double   best_fitness     = std::numeric_limits<double>::infinity();

      auto start_time = std::chrono::high_resolution_clock::now();
      do {
        algorithm.evolve();

        double current_best = algorithm.getBestFitness();
        if ( current_best < best_fitness ) {
          best_fitness     = current_best;
          stagnation_count = 0;
          DEBUG_PRINT( "Nova melhor solução encontrada: " << best_fitness );
        } else {
          stagnation_count++;
          if ( stagnation_count >= params.stag_limit ) {
            DEBUG_PRINT( "Critério de estagnação atingido após " << stagnation_count << " gerações sem melhoria" );
            break;
          }
        }

        if ( ( ++generation ) % params.x_intvl == 0 ) {
          algorithm.exchangeElite( params.x_number );
          DEBUG_PRINT( "Troca de elites realizada na geração " << generation );
        }
      } while ( generation < params.max_gens );

      auto end_time = std::chrono::high_resolution_clock::now();
      auto elapsed  = std::chrono::duration_cast<std::chrono::microseconds>( end_time - start_time );

      DEBUG_PRINT( "Melhor solução encontrada tem valor objetivo = " << algorithm.getBestFitness() );
      DEBUG_PRINT( "Número total de gerações executadas: " << generation );

      BEST = algorithm.getBestFitness();

      result.graph_name     = get_filename_from_path( params.file_path );
      result.node_count     = graph_l.order();
      result.edge_count     = graph_l.size();
      result.fitness        = static_cast<int>( BEST );
      result.elapsed_micros = elapsed.count();
      result.is_dense       = ( graph_l.size() / ( graph_l.order() * ( graph_l.order() - 1 ) / 2.0 ) ) > 0.5;
    }

    // Armazenar o resultado do trial no vetor
    all_results.push_back( result );
  }

  write_to_csv( params.output_file, all_results );

  std::cout << "Melhor fitness encontrado: " << BEST << "\n";
  return 0;
}

AlgorithmParams parse_args( int argc, char * argv[] ) {
  AlgorithmParams params;

  if ( argc < 2 ) {
    std::cerr << "Usage: " << argv[0] << " <input_file> [options]\n"
              << "Options:\n"
              << "  --population SIZE\n"
              << "  --elite-set FRACTION\n"
              << "  --mutants FRACTION\n"
              << "  --inherit-elite-probability VALUE\n"
              << "  --populations NUMBER\n"
              << "  --parallel NUMBER\n"
              << "  --exchange-interval GENERATIONS\n"
              << "  --exchange-number COUNT\n"
              << "  --max-generations NUMBER\n"
              << "  --stagnation-limit GENERATIONS\n"
              << "  --rng-seed VALUE  (0 ou omitido para semente aleatória)\n"
              << "  --output FILE\n"
              << "  --trials NUMBER\n";
    exit( 1 );
  }

  params.file_path = argv[1];

  for ( int i = 2; i < argc; i++ ) {
    std::string arg = argv[i];

    if ( arg == "--population" && i + 1 < argc ) {
      params.p = std::stoul( argv[++i] );
    } else if ( arg == "--elite-set" && i + 1 < argc ) {
      params.pe = std::stod( argv[++i] );
    } else if ( arg == "--mutants" && i + 1 < argc ) {
      params.pm = std::stod( argv[++i] );
    } else if ( arg == "--inherit-elite-probability" && i + 1 < argc ) {
      params.rhoe = std::stod( argv[++i] );
    } else if ( arg == "--populations" && i + 1 < argc ) {
      params.K = std::stoul( argv[++i] );
    } else if ( arg == "--parallel" && i + 1 < argc ) {
      params.MAXT = std::stoul( argv[++i] );
    } else if ( arg == "--exchange-interval" && i + 1 < argc ) {
      params.x_intvl = std::stoul( argv[++i] );
    } else if ( arg == "--exchange-number" && i + 1 < argc ) {
      params.x_number = std::stoul( argv[++i] );
    } else if ( arg == "--max-generations" && i + 1 < argc ) {
      params.max_gens = std::stoul( argv[++i] );
    } else if ( arg == "--stagnation-limit" && i + 1 < argc ) {
      params.stag_limit = std::stoul( argv[++i] );
    } else if ( arg == "--rng-seed" && i + 1 < argc ) {
      params.rng_seed    = std::stoul( argv[++i] );
      params.random_seed = ( params.rng_seed == 0 );  // Se for 0, ainda gera aleatório
    } else if ( arg == "--output" && i + 1 < argc ) {
      params.output_file = argv[++i];
    } else if ( arg == "--trials" && i + 1 < argc ) {
      params.trials = std::stoul( argv[++i] );
    } else {
      std::cerr << "Usage: " << argv[0] << " <input_file> [options]\n"
                << "Options:\n"
                << "  --population SIZE\n"
                << "  --elite-set FRACTION\n"
                << "  --mutants FRACTION\n"
                << "  --inherit-elite-probability VALUE\n"
                << "  --populations NUMBER\n"
                << "  --parallel NUMBER\n"
                << "  --exchange-interval GENERATIONS\n"
                << "  --exchange-number COUNT\n"
                << "  --max-generations NUMBER\n"
                << "  --stagnation-limit GENERATIONS\n"
                << "  --rng-seed VALUE  (0 ou omitido para semente aleatória)\n"
                << "  --output FILE\n"
                << "  --trials NUMBER\n";
      exit( 1 );
    }
  }

  return params;
}

void load_and_normalize_graph( const std::string & filename, bool & graph_is_matrix, ListGraph & graph_l, MatrixGraph & graph_m ) {
  std::ifstream file( filename );
  if ( !file.is_open() ) {
    throw std::runtime_error( "Não foi possível abrir o arquivo: " + filename );
  }

  std::set<int>                    vertices;
  std::vector<std::pair<int, int>> edges;
  std::string                      line;

  while ( std::getline( file, line ) ) {
    std::istringstream iss( line );
    int                u, v;
    if ( !( iss >> u >> v ) ) {
      continue;
    }
    vertices.insert( u );
    vertices.insert( v );
    if ( u != v ) {
      edges.emplace_back( u, v );
    }
  }

  if ( vertices.empty() ) {
    throw std::runtime_error( "Nenhum vértice encontrado na entrada" );
  }

  std::vector<int> missing_vertices;
  if ( !vertices.empty() ) {
    int min_vertex = *vertices.begin();
    int max_vertex = *vertices.rbegin();

    for ( int i = min_vertex; i <= max_vertex; i++ ) {
      if ( vertices.find( i ) == vertices.end() ) {
        missing_vertices.push_back( i );
      }
    }
  }

  std::unordered_map<int, int> id_map;
  int                          new_id = 0;
  for ( int original_id : vertices ) {
    id_map[original_id] = new_id++;
  }

  int                           num_vertices       = vertices.size();
  double                        max_possible_edges = num_vertices * ( num_vertices - 1 ) / 2.0;
  std::set<std::pair<int, int>> unique_edges;

  for ( const auto & edge : edges ) {
    int norm_u = id_map[edge.first];
    int norm_v = id_map[edge.second];
    if ( norm_u > norm_v ) {
      std::swap( norm_u, norm_v );
    }
    unique_edges.insert( { norm_u, norm_v } );
  }

  double density = unique_edges.size() / max_possible_edges;

  if ( density > 0.5 ) {
    graph_is_matrix = true;

    for ( const auto & edge : unique_edges ) {
      graph_m.add_edge( edge.first, edge.second );
    }
  } else {
    graph_is_matrix = false;
    graph_l         = ListGraph( num_vertices );

    for ( const auto & edge : unique_edges ) {
      graph_l.add_edge( edge.first, edge.second );
    }
  }
}

void write_to_csv( const std::string & filename, const std::vector<TrialResult> & results ) {
  std::ofstream file( filename, std::ios::app );
  if ( file.tellp() == 0 ) {
    file << "graph_name,graph_order,graph_size,fitness_value,elapsed_time(microseconds),is_dense\n";  // Cabeçalho
  }
  for ( const auto & result : results ) {
    file << result.graph_name << "," << result.node_count << "," << result.edge_count << "," << result.fitness << "," << result.elapsed_micros << ","
         << ( result.is_dense ? "yes" : "no" ) << "\n";
  }
  file.close();
}

std::string get_filename_from_path( const std::string & path ) {
  size_t pos = path.find_last_of( "/\\" );
  if ( pos == std::string::npos ) {
    return path;
  }

  std::string filename = path.substr( pos + 1 );

  size_t ext_pos = filename.find_last_of( '.' );
  if ( ext_pos != std::string::npos ) {
    filename = filename.substr( 0, ext_pos );
  }

  return filename;
}
