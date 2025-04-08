#include <iostream>

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

int main() {
  std::cout << "Hello World!\n";

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
