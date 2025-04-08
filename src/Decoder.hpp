#pragma once
#include "Chromosome.hpp"

#include <type_traits>
#include <vector>

template <typename GraphType>
class TRDDecoder {
private:
  const GraphType & graph;

  template <typename T = GraphType>
  typename std::enable_if<std::is_same<T, ListGraph>::value, void>::type fix_chromosome( Chromosome & chr ) const {
    chr.fix_l( graph );
  }

  template <typename T = GraphType>
  typename std::enable_if<std::is_same<T, MatrixGraph>::value, void>::type fix_chromosome( Chromosome & chr ) const {
    chr.fix_m( graph );
  }

public:
  explicit TRDDecoder( const GraphType & g ) : graph( g ) {}

  double decode( const std::vector<double> & chromosome ) const {
    const unsigned size = graph.order();
    Chromosome     chr( size );

    for ( size_t i = 0; i < size && i < chromosome.size(); ++i ) {
      double gene = chromosome[i];
      int    label;
      if ( gene < 0.3333 ) {
        label = LABEL_ZERO;
      } else if ( gene < 0.6666 ) {
        label = LABEL_ONE;
      } else {
        label = LABEL_TWO;
      }
      chr.set_value( i, label );
    }

    fix_chromosome<GraphType>( chr );

    return chr.calculate_fitness();
  }
};
