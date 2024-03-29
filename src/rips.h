#include <string>
#include <vector>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>


// for Dionysus
#include <tdautils/dionysusUtils.h>

// for Rips
#include <tdautils/ripsL2.h>
#include <tdautils/ripsArbit.h>



// ripsFiltration
/** \brief Interface for R code, construct the rips filtration on the input
  *         set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  dist           "euclidean" for Euclidean distance,
  *                            "arbitrary" for an arbitrary distance
  * @param[in]  library        "Dionysus"
  * @param[in]  printProgress  Is progress printed?
  * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
  *                            diagram. Diagram must point to enough memory space for
  *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
  *                            write nothing after.
  */
template< typename IntVector, typename RealMatrix, typename VectorList,
          typename RealVector, typename Print >
inline void ripsFiltration(
    const RealMatrix  & X,
    const unsigned      nSample,
    const unsigned      nDim,
    const int           maxdimension,
    const double        maxscale,
    const std::string & dist,
    const std::string & library,
    const bool          printProgress,
    const Print       & print,
    VectorList        & cmplx,
    RealVector        & values,
    VectorList        & boundary
) {
  if (dist[0] == 'e') {
    // RipsDiag for L2 distance
    filtrationDionysusToTda< IntVector >(RipsFiltrationDionysus< PairDistances, Generator, FltrR >(X, nSample, nDim, false, maxdimension, maxscale, printProgress, print), cmplx, values, boundary);
  }
  else {
    // RipsDiag for arbitrary distance
    filtrationDionysusToTda< IntVector >(RipsFiltrationDionysus< PairDistancesA, GeneratorA, FltrRA >(X, nSample, nDim, true, maxdimension, maxscale, printProgress, print), cmplx, values, boundary);
  }
}



// ripsDiag
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  dist           "euclidean" for Euclidean distance,
  *                            "arbitrary" for an arbitrary distance
  * @param[in]  libraryFiltration  "Dionysus"
  * @param[in]  libraryDiag        "Dionysus"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
  *                            diagram. Diagram must point to enough memory space for
  *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
  *                            write nothing after.
  */
template< typename RealMatrix, typename Print >
inline void ripsDiag(
    const RealMatrix  & X,
    const unsigned      nSample,
    const unsigned      nDim,
    const int           maxdimension,
    const double        maxscale,
    const std::string & dist,
    const std::string & libraryFiltration,
    const std::string & libraryDiag,
    const bool          location,
    const bool          printProgress,
    const Print       & print,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {
  
  if (dist[0] == 'e') {
    // RipsDiag for L2 distance
    FltrR filtration =
      RipsFiltrationDionysus< PairDistances, Generator, FltrR >(X, nSample, nDim, false, maxdimension, maxscale, printProgress, print);
    
    if (libraryDiag[0] == 'D') {
      FiltrationDiagDionysus< Persistence >(filtration, maxdimension, location, printProgress, persDgm,persLoc, persCycle);
    }
  }
  else {
    // RipsDiag for arbitrary distance
    FltrRA filtration =
      RipsFiltrationDionysus< PairDistancesA, GeneratorA, FltrRA >(X, nSample, nDim, true, maxdimension, maxscale, printProgress, print);
    
    if (libraryDiag[0] == 'D') {
      FiltrationDiagDionysus< Persistence >(filtration, maxdimension, location, printProgress, persDgm, persLoc, persCycle);
    }
  }
}

