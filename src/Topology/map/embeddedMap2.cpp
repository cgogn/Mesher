/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* version 0.1                                                                  *
* Copyright (C) 2009-2012, IGG Team, LSIIT, University of Strasbourg           *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#include <vector>
#include <algorithm>

#include "Topology/map/embeddedMap2.h"
#include "Topology/generic/traversor2.h"

namespace CGoGN
{

Dart EmbeddedMap2::newPolyLine(unsigned int nbEdges)
{
	Dart d = Map2::newPolyLine(nbEdges) ;

	if (isOrbitEmbedded<VERTEX>())
	{
		Dart e = d ;
		for (unsigned int i = 0 ; i <= nbEdges ; ++i)
		{
			initOrbitEmbeddingNewCell<VERTEX>(e) ;
			e = this->phi1(e) ;
		}
	}

	if (isOrbitEmbedded<EDGE>())
	{
		Dart e = d ;
		for (unsigned int i = 0 ; i < nbEdges ; ++i)
		{
			initOrbitEmbeddingNewCell<EDGE>(e) ;
			e = this->phi1(e) ;
		}
	}

	if (isOrbitEmbedded<FACE>())
	{
		initOrbitEmbeddingNewCell<FACE>(d) ;
	}
	return d ;
}

Dart EmbeddedMap2::newFace(unsigned int nbEdges, bool withBoundary)
{
	Dart d = Map2::newFace(nbEdges, withBoundary);

	if(withBoundary)
	{
		if (isOrbitEmbedded<VERTEX>())
		{
/*
			Traversor2FV<EmbeddedMap2> t(*this, d);
			for(Dart it = t.begin(); it != t.end(); it = t.next())
				initOrbitEmbeddingNewCell<VERTEX>(it) ;
*/
			Dart e = d;
			do
			{
				initOrbitEmbeddingNewCell<VERTEX>(e) ;
				e = this->phi1(e);
			} while (d!=e);
		}

		if(isOrbitEmbedded<EDGE>())
		{
/*
			Traversor2FE<EmbeddedMap2> t(*this, d);
			for(Dart it = t.begin(); it != t.end(); it = t.next())
				initOrbitEmbeddingNewCell<EDGE>(it) ;
*/
			Dart e = d;
			do
			{
				initOrbitEmbeddingNewCell<EDGE>(e) ;
				e = this->phi1(e);
			} while (d!=e);
		}

		if(isOrbitEmbedded<FACE>())
		{
			initOrbitEmbeddingNewCell<FACE>(d) ;
			initOrbitEmbeddingNewCell<FACE>(phi2(d)) ;
		}
	}
//	else
//	{
//		if (isOrbitEmbedded<VERTEX>())
//		{
///*
//			Traversor2FV<EmbeddedMap2> t(*this, d);
//			for(Dart it = t.begin(); it != t.end(); it = t.next())
//				initOrbitEmbeddingNewCell<VERTEX>(it) ;
//*/
//			Dart e = d;
//			do
//			{
//				initDartEmbedding<VERTEX>(e,newCell<VERTEX>());
//				e = this->phi1(e);
//			} while (d!=e);
//		}
//
//		if(isOrbitEmbedded<FACE>())
//			initOrbitEmbeddingNewCell<FACE>(d) ;
//
//	}
	return d ;
}

void EmbeddedMap2::splitVertex(Dart d, Dart e)
{
	Dart dd = phi2(d) ;
	Dart ee = phi2(e) ;

	Map2::splitVertex(d, e) ;

	if (isOrbitEmbedded<VERTEX>())
	{
		initDartEmbedding<VERTEX>(phi1(dd), getEmbedding<VERTEX>(d)) ;
		setOrbitEmbeddingOnNewCell<VERTEX>(e) ;
		copyCell<VERTEX>(e, d) ;
	}

	if(isOrbitEmbedded<EDGE>())
	{
		initOrbitEmbeddingNewCell<EDGE>(phi1(dd)) ;
	}

	if(isOrbitEmbedded<FACE>())
	{
		initDartEmbedding<FACE>(phi1(dd), getEmbedding<FACE>(dd)) ;
		initDartEmbedding<FACE>(phi1(ee), getEmbedding<FACE>(ee)) ;
	}
}

Dart EmbeddedMap2::deleteVertex(Dart d)
{
	Dart f = Map2::deleteVertex(d) ;
	if(f != NIL)
	{
		if (isOrbitEmbedded<FACE>())
		{
			setOrbitEmbedding<FACE>(f, getEmbedding<FACE>(f)) ;
		}
	}
	return f ;
}

Dart EmbeddedMap2::cutEdge(Dart d)
{
	Dart nd = Map2::cutEdge(d) ;

	if(isOrbitEmbedded<VERTEX>())
	{
		initOrbitEmbeddingNewCell<VERTEX>(nd) ;
	}

	if (isOrbitEmbedded<EDGE>())
	{
		initDartEmbedding<EDGE>(phi2(d), getEmbedding<EDGE>(d)) ;
		setOrbitEmbeddingOnNewCell<EDGE>(nd) ;
		copyCell<EDGE>(nd, d) ;
	}

	if(isOrbitEmbedded<FACE>())
	{
		initDartEmbedding<FACE>(nd, getEmbedding<FACE>(d)) ;
		Dart e = phi2(nd) ;
		initDartEmbedding<FACE>(phi1(e), getEmbedding<FACE>(e)) ;
	}

	return nd;
}

bool EmbeddedMap2::uncutEdge(Dart d)
{
	if(Map2::uncutEdge(d))
	{
		if(isOrbitEmbedded<EDGE>())
		{
			copyDartEmbedding<EDGE>(phi2(d), d) ;
		}
		return true ;
	}
	return false ;
}

bool EmbeddedMap2::edgeCanCollapse(Dart d)
{
	if(isBoundaryVertex(d) || isBoundaryVertex(phi1(d)))
		return false ;

	unsigned int val_v1 = vertexDegree(d) ;
	unsigned int val_v2 = vertexDegree(phi1(d)) ;

	if(val_v1 + val_v2 < 8 || val_v1 + val_v2 > 14)
		return false ;

	if(faceDegree(d) == 3)
	{
		if(vertexDegree(phi_1(d)) < 4)
			return false ;
	}

	Dart dd = phi2(d) ;
	if(faceDegree(dd) == 3)
	{
		if(vertexDegree(phi_1(dd)) < 4)
			return false ;
	}

	// Check vertex sharing condition
	std::vector<unsigned int> vu1 ;
	vu1.reserve(32) ;
	Dart vit1 = alpha1(alpha1(d)) ;
	Dart end = phi1(dd) ;
	do
	{
		unsigned int ve = getEmbedding<VERTEX>(phi2(vit1)) ;
		vu1.push_back(ve) ;
		vit1 = alpha1(vit1) ;
	} while(vit1 != end) ;
	end = phi1(d) ;
	Dart vit2 = alpha1(alpha1(dd)) ;
	do
	{
		unsigned int ve = getEmbedding<VERTEX>(phi2(vit2)) ;
		std::vector<unsigned int>::iterator it = std::find(vu1.begin(), vu1.end(), ve) ;
		if(it != vu1.end())
			return false ;
		vit2 = alpha1(vit2) ;
	} while(vit2 != end) ;

	return true ;
}

Dart EmbeddedMap2::collapseEdge(Dart d, bool delDegenerateFaces)
{
	unsigned int vEmb = EMBNULL ;
	if (isOrbitEmbedded<VERTEX>())
	{
		vEmb = getEmbedding<VERTEX>(d) ;
	}

	Dart dV = Map2::collapseEdge(d, delDegenerateFaces);

	if (isOrbitEmbedded<VERTEX>())
	{
		setOrbitEmbedding<VERTEX>(dV, vEmb) ;
	}
	
	return dV ;
}

bool EmbeddedMap2::flipEdge(Dart d)
{
	if(Map2::flipEdge(d))
	{
		Dart e = phi2(d) ;

		if (isOrbitEmbedded<VERTEX>())
		{
			copyDartEmbedding<VERTEX>(d, phi1(e)) ;
			copyDartEmbedding<VERTEX>(e, phi1(d)) ;
		}

		if (isOrbitEmbedded<FACE>())
		{
			copyDartEmbedding<FACE>(phi_1(d), d) ;
			copyDartEmbedding<FACE>(phi_1(e), e) ;
		}

		return true ;
	}
	return false ;
}

bool EmbeddedMap2::flipBackEdge(Dart d)
{
	if(Map2::flipBackEdge(d))
	{
		Dart e = phi2(d) ;

		if (isOrbitEmbedded<VERTEX>())
		{
			copyDartEmbedding<VERTEX>(d, phi1(e)) ;
			copyDartEmbedding<VERTEX>(e, phi1(d)) ;
		}

		if (isOrbitEmbedded<FACE>())
		{
			copyDartEmbedding<FACE>(phi1(d), d) ;
			copyDartEmbedding<FACE>(phi1(e), e) ;
		}

		return true ;
	}
	return false ;
}

void EmbeddedMap2::swapEdges(Dart d, Dart e)
{
	Dart d2 = phi2(d);
	Dart e2 = phi2(e);
	Map2::swapEdges(d,e);

	if(isOrbitEmbedded<VERTEX>())
	{
		copyDartEmbedding<VERTEX>(d, phi2(phi_1(d)));
		copyDartEmbedding<VERTEX>(e, phi2(phi_1(e)));
		copyDartEmbedding<VERTEX>(d2, phi2(phi_1(d2)));
		copyDartEmbedding<VERTEX>(e2, phi2(phi_1(e2)));
	}

	if(isOrbitEmbedded<EDGE>())
	{

	}

	if(isOrbitEmbedded<VOLUME>())
		setOrbitEmbeddingOnNewCell<VOLUME>(d);
}

void EmbeddedMap2::insertEdgeInVertex(Dart d, Dart e)
{
	Map2::insertEdgeInVertex(d, e);

	if (isOrbitEmbedded<VERTEX>())
	{
		copyDartEmbedding<VERTEX>(e, d) ;
	}

	if (isOrbitEmbedded<FACE>())
	{
		if(!sameFace(d,e))
		{
			setOrbitEmbeddingOnNewCell<FACE>(e);
			copyCell<FACE>(e, d);
		}
		else
		{
			setOrbitEmbedding<FACE>(d, getEmbedding<FACE>(d)) ;
		}
	}
}

bool EmbeddedMap2::removeEdgeFromVertex(Dart d)
{
	Dart dPrev = alpha_1(d);

	if (dPrev == d) return false ;

	bool b = Map2::removeEdgeFromVertex(d);

	if (isOrbitEmbedded<VERTEX>())
	{
		setOrbitEmbeddingOnNewCell<VERTEX>(d);
		copyCell<VERTEX>(d, dPrev);
	}

	if (isOrbitEmbedded<FACE>())
	{
		if(!sameFace(d, dPrev))
		{
			setOrbitEmbeddingOnNewCell<FACE>(d);
			copyCell<FACE>(d, dPrev);
		}
		else
		{
			setDartEmbedding<FACE>(d, getEmbedding<FACE>(d)) ;
		}
	}
	return b ;
}

void EmbeddedMap2::sewFaces(Dart d, Dart e, bool withBoundary)
{
	if (!withBoundary)
	{
		Map2::sewFaces(d, e, false) ;

//		if(isOrbitEmbedded<EDGE>())
//		{
///*
//			initOrbitEmbeddingNewCell<EDGE>(d) ;
//*/
//			unsigned int emb = newCell<EDGE>();
//			initDartEmbedding<EDGE>(d,emb);
//			initDartEmbedding<EDGE>(e,emb);
//		}

		return ;
	}

	Map2::sewFaces(d, e, withBoundary) ;

	if (isOrbitEmbedded<VERTEX>())
	{
		setOrbitEmbedding<VERTEX>(d, getEmbedding<VERTEX>(d)) ;
		setOrbitEmbedding<VERTEX>(e, getEmbedding<VERTEX>(phi1(d))) ;
	}

	if (isOrbitEmbedded<EDGE>())
	{
		copyDartEmbedding<EDGE>(e, d) ;
	}
}

void EmbeddedMap2::unsewFaces(Dart d, bool withBoundary)
{
	if (!withBoundary)
	{
		Map2::unsewFaces(d, false) ;
		return ;
	}

	Dart e = phi2(d) ;
	Map2::unsewFaces(d) ;

	if (isOrbitEmbedded<VERTEX>())
	{
		copyDartEmbedding<VERTEX>(phi2(e), d) ;
		copyDartEmbedding<VERTEX>(phi2(d), e) ;

		Dart ee = phi1(e) ;
		if(!sameVertex(d, ee))
		{
			setOrbitEmbeddingOnNewCell<VERTEX>(ee);
			copyCell<VERTEX>(ee, d);
		}

		Dart dd = phi1(d) ;
		if(!sameVertex(e, dd))
		{
			setOrbitEmbeddingOnNewCell<VERTEX>(dd);
			copyCell<VERTEX>(dd, e);
		}
	}

	if (isOrbitEmbedded<EDGE>())
	{
		setOrbitEmbeddingOnNewCell<EDGE>(e);
		copyCell<EDGE>(e, d);
	}
}

bool EmbeddedMap2::collapseDegeneratedFace(Dart d)
{
	Dart e = phi2(d) ;
	bool updateEdgeEmb = false ;
	if(phi1(d) != d)
		updateEdgeEmb = true ;

	if(Map2::collapseDegeneratedFace(d))
	{
		if (isOrbitEmbedded<EDGE>() && updateEdgeEmb)
		{
			copyDartEmbedding<EDGE>(phi2(e), e) ;
		}
		return true ;
	}
	return false ;
}

void EmbeddedMap2::splitFace(Dart d, Dart e)
{
	Map2::splitFace(d, e) ;

	if (isOrbitEmbedded<VERTEX>())
	{
		initDartEmbedding<VERTEX>(phi_1(e), getEmbedding<VERTEX>(d)) ;
		initDartEmbedding<VERTEX>(phi_1(d), getEmbedding<VERTEX>(e)) ;
	}

	if(isOrbitEmbedded<EDGE>())
	{
		initOrbitEmbeddingNewCell<EDGE>(phi_1(d)) ;
	}

	if (isOrbitEmbedded<FACE>())
	{
		initDartEmbedding<FACE>(phi_1(d), getEmbedding<FACE>(d)) ;
		setOrbitEmbeddingOnNewCell<FACE>(e) ;
		copyCell<FACE>(e, d) ;
	}
}

bool EmbeddedMap2::mergeFaces(Dart d)
{
	Dart dNext = phi1(d) ;

	if(Map2::mergeFaces(d))
	{
		if (isOrbitEmbedded<FACE>())
		{
			setOrbitEmbedding<FACE>(dNext, getEmbedding<FACE>(dNext)) ;
		}
		return true ;
	}
	return false ;
}

bool EmbeddedMap2::mergeVolumes(Dart d, Dart e, bool deleteFace)
{
	std::vector<Dart> darts ;
	std::vector<unsigned int> vEmb ;
	vEmb.reserve(32) ;
	std::vector<unsigned int> eEmb ;
	eEmb.reserve(32) ;
	Dart fit = d ;
	do
	{
		darts.push_back(phi2(fit)) ;

		if (isOrbitEmbedded<VERTEX>())
		{
			vEmb.push_back(getEmbedding<VERTEX>(phi2(fit))) ;
		}

		if (isOrbitEmbedded<EDGE>())
		{
			eEmb.push_back(getEmbedding<EDGE>(fit)) ;
		}
		
		fit = phi1(fit) ;
	} while (fit != d) ;

	if(Map2::mergeVolumes(d, e, deleteFace))
	{
		for(unsigned int i = 0; i < darts.size(); ++i)
		{
			if (isOrbitEmbedded<VERTEX>())
			{
				setOrbitEmbedding<VERTEX>(darts[i], vEmb[i]) ;
			}

			if (isOrbitEmbedded<EDGE>())
			{
				setOrbitEmbedding<EDGE>(darts[i], eEmb[i]) ;
			}
		}
		return true ;
	}
	return false ;
}

void EmbeddedMap2::splitSurface(std::vector<Dart>& vd, bool firstSideClosed, bool secondSideClosed)
{
	std::vector<Dart> darts ;
	darts.reserve(vd.size());

	// save the edge neighbors darts
	for(std::vector<Dart>::iterator it = vd.begin() ; it != vd.end() ; ++it)
	{
		darts.push_back(phi2(*it));
	}

	assert(darts.size() == vd.size());

	Map2::splitSurface(vd, firstSideClosed, secondSideClosed);

	// follow the edge path a second time to embed the vertex, edge and volume orbits
	for(unsigned int i = 0; i < vd.size(); ++i)
	{
		Dart dit = vd[i];
		Dart dit2 = darts[i];

		// embed the vertex embedded from the origin volume to the new darts
		if(isOrbitEmbedded<VERTEX>())
		{
			initDartEmbedding<VERTEX>(phi2(dit), getEmbedding<VERTEX>(phi1(dit)));
			initDartEmbedding<VERTEX>(phi2(dit2), getEmbedding<VERTEX>(phi1(dit2)));
		}

		// embed the edge embedded from the origin volume to the new darts
		if(isOrbitEmbedded<EDGE>())
		{
			initDartEmbedding<EDGE>(phi2(dit), getEmbedding<EDGE>(dit));
			setOrbitEmbeddingOnNewCell<EDGE>(phi2(dit2));
			copyCell<EDGE>(dit2, dit);
		}

		// embed the volume embedded from the origin volume to the new darts
		if(isOrbitEmbedded<VOLUME>())
		{
			initDartEmbedding<VOLUME>(phi2(dit), getEmbedding<VOLUME>(dit));
			initDartEmbedding<VOLUME>(phi2(dit2), getEmbedding<VOLUME>(dit2));
		}
	}
}

unsigned int EmbeddedMap2::closeHole(Dart d, bool forboundary)
{
	unsigned int nbE = Map2::closeHole(d, forboundary) ;
	Dart dd = phi2(d) ;
	Dart f = dd ;
	do
	{
		if (isOrbitEmbedded<VERTEX>())
		{
			unsigned int emb = getEmbedding<VERTEX>(phi1(phi2(f)));
			if (emb == EMBNULL)
				initOrbitEmbeddingNewCell<VERTEX>(f) ;
			else
				initDartEmbedding<VERTEX>(f, emb) ;
		}

		if (isOrbitEmbedded<EDGE>())
		{
			unsigned int emb = getEmbedding<EDGE>(phi2(f));
			if (emb == EMBNULL)
				initOrbitEmbeddingNewCell<EDGE>(f) ;
			else
				initDartEmbedding<EDGE>(f, emb) ;
		}

		f = phi1(f) ;
	} while(dd != f) ;

	if(isOrbitEmbedded<FACE>())
	{
		initOrbitEmbeddingNewCell<FACE>(dd) ;
	}

	return nbE ;
}

bool EmbeddedMap2::check()
{
	bool topo = Map2::check() ;
	if (!topo)
		return false ;

	CGoGNout << "Check: embedding begin" << CGoGNendl ;
	for(Dart d = begin(); d != end(); next(d))
	{
		if (isOrbitEmbedded<VERTEX>())
		{
			if (getEmbedding<VERTEX>(d) != getEmbedding<VERTEX>(alpha1(d)))
			{
				CGoGNout << "Check: different embeddings on vertex" << CGoGNendl ;
				return false ;
			}
		}

		if (isOrbitEmbedded<EDGE>())
		{
			if (getEmbedding<EDGE>(d) != getEmbedding<EDGE>(phi2(d)))
			{
				CGoGNout << "Check: different embeddings on edge" << CGoGNendl ;
				return false ;
			}
		}

//		if (isOrbitEmbedded(ORIENTED_FACE))
//		{
//			if (getEmbedding(ORIENTED_FACE, d) != getEmbedding(ORIENTED_FACE, phi1(d)))
//		{
//				CGoGNout << "Check: different embeddings on oriented face" << CGoGNendl ;
//				return false ;
//			}
//		}

		if (isOrbitEmbedded<FACE>())
		{
			if (getEmbedding<FACE>(d) != getEmbedding<FACE>(phi1(d)))
			{
				CGoGNout << "Check: different embeddings on face" << CGoGNendl ;
				return false ;
			}
		}
	}

	CGoGNout << "Check: embedding ok" << CGoGNendl ;

    std::cout << "nb vertex orbits : " << getNbOrbits<VERTEX>() << std::endl ;
    std::cout << "nb vertex cells : " << m_attribs[VERTEX].size() << std::endl ;

    std::cout << "nb edge orbits : " << getNbOrbits<EDGE>() << std::endl ;
    std::cout << "nb edge cells : " << m_attribs[EDGE].size() << std::endl ;

    std::cout << "nb face orbits : " << getNbOrbits<FACE>() << std::endl ;
    std::cout << "nb face cells : " << m_attribs[FACE].size() << std::endl ;

	return true ;
}

} // namespace CGoGN
