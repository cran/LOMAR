#ifndef __TYPECASTUTILS_H__
#define __TYPECASTUTILS_H__

#include <vector>
#include <map>
#include <algorithm>



template<typename PersistenceDiagram, typename RcppMatrix>
inline PersistenceDiagram RcppToDionysus(const RcppMatrix& rcppMatrix) {
	PersistenceDiagram dionysusDiagram;
	const unsigned rowNum = rcppMatrix.nrow();
	for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx)
	{
		dionysusDiagram.push_back(typename PersistenceDiagram::Point(
				rcppMatrix[rowIdx + 0 * rowNum], rcppMatrix[rowIdx + 1 * rowNum]));
	}
	return dionysusDiagram;
}



template< typename StlMatrix, typename RealMatrix >
inline StlMatrix TdaToStl(const RealMatrix & rcppMatrix,
    const unsigned nRow, const unsigned nCol, bool is_row_names = false) {

  if (is_row_names) {
    StlMatrix stlMatrix(nRow, typename StlMatrix::value_type(nCol + 1));
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      stlMatrix[rowIdx][0] = rowIdx + 1;
    }
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      for (unsigned colIdx = 0; colIdx < nCol; ++colIdx) {
        stlMatrix[rowIdx][colIdx + 1] = rcppMatrix[rowIdx + colIdx * nRow];
      }
    }
    return stlMatrix;
  }
  else {
    StlMatrix stlMatrix(nRow, typename StlMatrix::value_type(nCol));
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      for (unsigned colIdx = 0; colIdx < nCol; ++colIdx) {
        stlMatrix[rowIdx][colIdx] = rcppMatrix[rowIdx + colIdx * nRow];
      }
    }
    return stlMatrix;
  }
}



template< typename StlMatrix, typename RcppMatrix >
inline StlMatrix RcppToStl(const RcppMatrix& rcppMatrix,
		bool is_row_names = false) {

	const unsigned rowNum = rcppMatrix.nrow();
	const unsigned colNum = rcppMatrix.ncol();
	if (is_row_names) {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum + 1));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			stlMatrix[rowIdx][0] = rowIdx + 1;
		}
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx + 1] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	else {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	
}








template< typename VertexVector, typename RcppVector, typename RcppList >
inline std::vector< VertexVector > RcppCmplxToStl(
    const RcppList & rcppCmplx, const int idxShift) {

  const unsigned nCmplx = rcppCmplx.size();
  std::vector< VertexVector > stlCmplx(nCmplx);

  typename RcppList::const_iterator iRcppVec = rcppCmplx.begin();
  typename std::vector< VertexVector >::iterator iStlVec = stlCmplx.begin();
  for (; iRcppVec != rcppCmplx.end(); ++iRcppVec, ++iStlVec) {
    RcppVector cmplxVec(*iRcppVec);
    *iStlVec = VertexVector(cmplxVec.size());

    typename RcppVector::const_iterator iRcpp = cmplxVec.begin();
    typename VertexVector::iterator iStl = iStlVec->begin();
    for (; iRcpp != cmplxVec.end(); ++iRcpp, ++iStl) {
      *iStl = *iRcpp - idxShift;
    }
  }

  return stlCmplx;
}



template< typename RcppVector, typename RcppList, typename VectorList >
inline RcppList StlCmplxToRcpp(
    const VectorList & stlCmplx, const int idxShift) {

  const unsigned nCmplx = stlCmplx.size();
  RcppList rcppCmplx(nCmplx);

  typename VectorList::const_iterator iStlVec = stlCmplx.begin();
  typename RcppList::iterator iRcppVec = rcppCmplx.begin();
  for (; iStlVec != stlCmplx.end(); ++iStlVec, ++iRcppVec) {
    RcppVector cmplxVec(iStlVec->size());

    typename VectorList::value_type::const_iterator iStl = iStlVec->begin();
    typename RcppVector::iterator iRcpp = cmplxVec.begin();
    for (; iStl != iStlVec->end(); ++iStl, ++iRcpp) {
      *iRcpp = *iStl + idxShift;
    }
    *iRcppVec = cmplxVec;
  }

  return rcppCmplx;
}



template<typename RcppMatrix, typename StlMatrix>
inline RcppMatrix concatStlToRcpp(const std::vector< StlMatrix >& stlMatrices,
		bool includeIndex, unsigned colNum) {
	unsigned rowNum = 0;

	typename std::vector< StlMatrix >::const_iterator vecItr;
	for (vecItr = stlMatrices.begin(); vecItr != stlMatrices.end(); ++vecItr) {
		rowNum += vecItr->size();
	}
	RcppMatrix rcppMatrix(rowNum, colNum);

	unsigned vecIdx, rowIdx, colIdx;
	for (vecIdx = 0, rowIdx = 0; vecIdx < stlMatrices.size(); ++vecIdx) {
		typename StlMatrix::const_iterator matItr;
		for (matItr = stlMatrices[vecIdx].begin();
				matItr != stlMatrices[vecIdx].end(); ++matItr, ++rowIdx) {
			if (includeIndex) {
				rcppMatrix[rowIdx] = vecIdx;
				for (colIdx = 0; colIdx < colNum - 1; ++colIdx) {
					rcppMatrix[rowIdx + (colIdx + 1) * rowNum] = (*matItr)[colIdx];
				}
			}
			else {
				for (colIdx = 0; colIdx < colNum; ++colIdx) {
					rcppMatrix[rowIdx + colIdx * rowNum] = (*matItr)[colIdx];
				}
			}
		}
	}

	return rcppMatrix;
}



template<typename RcppList, typename RcppVector, typename StlSet>
inline RcppList StlToRcppList(
		const std::vector< std::vector< StlSet > >& stlSets) {
	unsigned rowNum = 0;

	typename std::vector< std::vector< StlSet > >::const_iterator vecsItr;
	for (vecsItr = stlSets.begin(); vecsItr != stlSets.end(); ++vecsItr) {
		rowNum += vecsItr->size();
	}
	RcppList rcppList(rowNum);

	typename RcppList::iterator listItr;
	typename std::vector< StlSet >::const_iterator setVecItr;
	typename StlSet::const_iterator setItr;
	unsigned setIdx;
	for (vecsItr = stlSets.begin(), listItr = rcppList.begin();
			vecsItr != stlSets.end(); ++vecsItr) {

		for (setVecItr = vecsItr->begin(); setVecItr != vecsItr->end(); ++setVecItr, ++listItr) {
			RcppVector rcppVec(setVecItr->size());
			for (setIdx = 0, setItr = setVecItr->begin(); setItr != setVecItr->end(); ++setItr, ++setIdx) {
				rcppVec[setIdx] = *setItr;
			}
			*listItr = rcppVec;
		}
	}

	return rcppList;
}



template< typename RcppList, typename RcppMatrix, typename StlVector >
inline RcppList StlToRcppMatrixList(
	const std::vector< std::vector< std::vector< StlVector > > >& stlArrays) {
	unsigned listNum = 0;

	typename std::vector< std::vector< std::vector< StlVector > > >::const_iterator vecsItr;
	for (vecsItr = stlArrays.begin(); vecsItr != stlArrays.end(); ++vecsItr) {
		listNum += vecsItr->size();
	}
	RcppList rcppList(listNum);

	typename RcppList::iterator listItr;
	typename std::vector< std::vector< StlVector > >::const_iterator matrixItr;
	typename std::vector< StlVector >::const_iterator rowItr;
	typename StlVector::const_iterator colItr;
	unsigned rowIdx, colIdx, rowNum;
	for (vecsItr = stlArrays.begin(), listItr = rcppList.begin();
	vecsItr != stlArrays.end(); ++vecsItr) {

		for (matrixItr = vecsItr->begin(); matrixItr != vecsItr->end();
		++matrixItr, ++listItr) {
			rowNum = matrixItr->size();
			if (rowNum != 0) {
				RcppMatrix rcppMatrix(rowNum, (*matrixItr)[0].size());
				for (rowIdx = 0, rowItr = matrixItr->begin();
				rowItr != matrixItr->end(); ++rowIdx, ++rowItr) {
					for (colIdx = 0, colItr = rowItr->begin(); colItr != rowItr->end();
					++colIdx, ++colItr) {
						rcppMatrix[rowIdx + colIdx * rowNum] = *colItr;
					}
				}
				*listItr = rcppMatrix;
			}
			else {
				RcppMatrix rcppMatrix(0, 0);
				*listItr = rcppMatrix;
			}
		}
	}

	return rcppList;
}

template< typename Simplex, typename SimplexMap, typename RealVector >
inline void filtrationDionysusOne(
  const Simplex & c, const SimplexMap & simplex_map, const int idxShift,
  RealVector & cmplxVec, double & value, RealVector & boundaryVec) {

  const unsigned nVtx = c.dimension() + 1;

  cmplxVec = RealVector(nVtx);
  typename RealVector::iterator iCmplxVec = cmplxVec.begin();
  for (typename Simplex::VertexContainer::const_iterator vit =
       c.vertices().begin(); vit != c.vertices().end(); ++vit, ++iCmplxVec) {
    // R is 1-base, while C++ is 0-base
    *iCmplxVec = *vit + idxShift;
  }

  value = c.data();

  // might need to change for cubical complex
  if (nVtx > 1) {
    boundaryVec = RealVector(nVtx);
  }
  typename RealVector::iterator iBdyVec = boundaryVec.begin();
  for (typename Simplex::BoundaryIterator bit = c.boundary_begin();
       bit != c.boundary_end(); ++bit, ++iBdyVec) {
    // R is 1-base, while C++ is 0-base
    *iBdyVec = simplex_map.find(*bit)->second + idxShift;
  }
}



template< typename IntegerVector, typename Filtration, typename VectorList,
          typename RealVector >
inline void filtrationDionysusToTda(
    const Filtration & filtration, VectorList & cmplx, RealVector & values,
    VectorList & boundary) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, unsigned,
      typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary = VectorList(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();
  typename VectorList::iterator iBdy = boundary.begin();

  for (typename Filtration::Index it = filtration.begin();
      it != filtration.end(); ++it, ++iCmplx, ++iValue, ++iBdy) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    IntegerVector cmplxVec;
    IntegerVector boundaryVec;
    filtrationDionysusOne(c, simplex_map, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;

    simplex_map.insert(typename
        std::map< typename Filtration::Simplex, unsigned >::value_type(
        c, size_of_simplex_map++));
  }
}



template< typename RcppList, typename RcppVector, typename Filtration >
inline RcppList filtrationDionysusToRcpp(const Filtration & filtration) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, unsigned,
    typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;

  RcppList cmplx(nFltr);
  RcppVector values(nFltr);
  RcppList boundary(nFltr);
  typename RcppList::iterator iCmplx = cmplx.begin();
  typename RcppVector::iterator iValue = values.begin();
  typename RcppList::iterator iBdy = boundary.begin();

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it, ++iCmplx, ++iValue, ++iBdy) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    RcppVector cmplxVec;
    RcppVector boundaryVec;
    filtrationDionysusOne(c, simplex_map, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;

    simplex_map.insert(typename
        std::map< typename Filtration::Simplex, unsigned >::value_type(
            c, size_of_simplex_map++));
  }

  return RcppList::create(cmplx, values, boundary);
}



template< typename IntegerVector, typename Filtration, typename VectorList,
          typename RealVector >
inline Filtration filtrationTdaToDionysus(
    const VectorList & cmplx, const RealVector & values,
    const unsigned idxShift) {

  Filtration filtration;

  typename VectorList::const_iterator iCmplx = cmplx.begin();
  typename RealVector::const_iterator iValue = values.begin();
  for (; iCmplx != cmplx.end(); ++iCmplx, ++iValue) {
    const IntegerVector tdaVec(*iCmplx);
    IntegerVector dionysusVec(tdaVec.size());
    typename IntegerVector::const_iterator iTda = tdaVec.begin();
    typename IntegerVector::iterator iDionysus = dionysusVec.begin();
    for (; iTda != tdaVec.end(); ++iTda, ++iDionysus) {
      // R is 1-base, while C++ is 0-base
      *iDionysus = *iTda - idxShift;
    }
    filtration.push_back(typename Filtration::Simplex(
        dionysusVec.begin(), dionysusVec.end(), *iValue));
  }

  return filtration;
}



template< typename Filtration, typename RcppVector, typename RcppList >
inline Filtration filtrationRcppToDionysus(const RcppList & rcppList) {

  const RcppList rcppComplex(rcppList[0]);
  const RcppVector rcppValue(rcppList[1]);
  Filtration filtration;

  typename RcppList::const_iterator iCmplx = rcppComplex.begin();
  typename RcppVector::const_iterator iValue = rcppValue.begin();
  for (; iCmplx != rcppComplex.end(); ++iCmplx, ++iValue) {
    const RcppVector rcppVec(*iCmplx);
    RcppVector dionysusVec(rcppVec.size());
    typename RcppVector::const_iterator iRcpp = rcppVec.begin();
    typename RcppVector::iterator iDionysus = dionysusVec.begin();
    for (; iRcpp != rcppVec.end(); ++iRcpp, ++iDionysus) {
      // R is 1-base, while C++ is 0-base
      *iDionysus = *iRcpp - 1;
    }
    filtration.push_back(typename Filtration::Simplex(
        dionysusVec.begin(), dionysusVec.end(), *iValue));
  }

  return filtration;
}


# endif // __TYPECASTUTILS_H__
