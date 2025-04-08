#include "ListGraph.hpp"
#include "MatrixGraph.hpp"

#include <fstream>
#include <iostream>
#include <set>

#define DEBUG

#if defined( DEBUG )
  #define DEBUG_PRINT( msg )     std::cout << "DEBUG: " << msg << std::endl;
  #define DEBUG_PRINT_VAR( var ) std::cout << "DEBUG: " #var " = " << var << std::endl;
#else
  #define DEBUG_PRINT( msg )
  #define DEBUG_PRINT_VAR( var )
#endif

struct AlgorithmParams {
  unsigned    n    = 100;   // Size of chromosomes
  unsigned    p    = 1000;  // Size of population
  double      pe   = 0.20;  // Fraction of population to be the elite-set
  double      pm   = 0.10;  // Fraction of population to be replaced by mutants
  double      rhoe = 0.70;  // Probability offspring inherits elite parent allele
  unsigned    K    = 3;     // Number of independent populations
  unsigned    MAXT = 2;     // Number of threads for parallel decoding
  std::string file_path;    // File path for the input file
  std::string output_file = "results.csv";
};

AlgorithmParams parse_args( int argc, char * argv[] );
void            load_and_normalize_graph( const std::string & filename, bool & graph_is_matrix, ListGraph & graph_l, MatrixGraph & graph_m );

int main( int argc, char * argv[] ) {
  auto params = parse_args( argc, argv );

  ListGraph   graph_l( 1 );
  MatrixGraph graph_m( 1 );

  bool graph_is_matrix = false;
  load_and_normalize_graph( params.file_path, graph_is_matrix, graph_l, graph_m );

  return 0;
}

AlgorithmParams parse_args( int argc, char * argv[] ) {
  AlgorithmParams params;

  if ( argc < 2 ) {
    std::cerr << "Usage: " << argv[0] << " <input_file> [options]\n"
              << "Options:\n"
              << "  --chromosomes SIZE\n"
              << "  --population SIZE\n"
              << "  --elite-set FRACTION\n"
              << "  --mutants FRACTION\n"
              << "  --inherit-elite-probability VALUE\n"
              << "  --populations NUMBER\n"
              << "  --threads NUMBER\n"
              << "  --output FILE\n";
    exit( 1 );
  }

  params.file_path = argv[1];

  for ( int i = 2; i < argc; i++ ) {
    std::string arg = argv[i];

    if ( arg == "--chromosomes" && i + 1 < argc ) {
      params.n = std::stoul( argv[++i] );
    } else if ( arg == "--population" && i + 1 < argc ) {
      params.p = std::stoul( argv[++i] );
    } else if ( arg == "--elite-set" && i + 1 < argc ) {
      params.pe = std::stod( argv[++i] );
    } else if ( arg == "--mutants" && i + 1 < argc ) {
      params.pm = std::stod( argv[++i] );
    } else if ( arg == "--inherit-elite-probability" && i + 1 < argc ) {
      params.rhoe = std::stod( argv[++i] );
    } else if ( arg == "--populations" && i + 1 < argc ) {
      params.K = std::stoul( argv[++i] );
    } else if ( arg == "--threads" && i + 1 < argc ) {
      params.MAXT = std::stoul( argv[++i] );
    } else if ( arg == "--output" && i + 1 < argc ) {
      params.output_file = argv[++i];
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

  // Identificar lacunas nos IDs de vértices
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
    graph_m         = MatrixGraph( num_vertices );

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
