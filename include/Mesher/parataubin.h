/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* version 0.1                                                                  *
* Copyright (C) 2009, IGG Team, LSIIT, University of Strasbourg                *
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
* Web site: https://iggservis.u-strasbg.fr/CGoGN/                              *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef __PARA_TAUBIN_H_
#define __PARA_TAUBIN_H_

#include "Algo/Parallel/parallel_foreach.h"
#include "Algo/Selection/collector.h"
#include "Algo/Filtering/functors.h"


namespace CGoGN
{
namespace Mesher
{


template<typename PFP>
class PassOne: public FunctorMapThreaded<typename PFP::MAP>
{
protected:
	const VertexAttribute<typename PFP::VEC3>& m_position1;
	VertexAttribute<typename PFP::VEC3>& m_position2;
	Algo::Surface::Selection::Collector_OneRing<PFP> m_collect;
	Algo::Surface::Filtering::FunctorAverage<typename PFP::VEC3, VERTEX> m_fa;

public:
	PassOne( typename PFP::MAP& map, const VertexAttribute<typename PFP::VEC3>& position1, VertexAttribute<typename PFP::VEC3>& position2):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_position1(position1), m_position2(position2), m_collect(map), m_fa(position1)
	{}

	FunctorMapThreaded<typename PFP::MAP>* duplicate() const
	{
		PassOne<PFP>* copy = new PassOne<PFP>(this->m_map ,m_position1,m_position2);
		FunctorMapThreaded<typename PFP::MAP>* ptr = reinterpret_cast< FunctorMapThreaded<typename PFP::MAP>* >(copy);
		return ptr;
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		const float lambda = 0.6307 ;
		if(!this->m_map.isBoundaryVertex(d))
		{
			m_collect.collectBorder(d) ;
			m_fa.reset() ;
			m_collect.applyOnBorder(m_fa) ;
			typename PFP::VEC3 p = m_position1[d] ;
			typename PFP::VEC3 displ = m_fa.getAverage() - p ;
			displ *= lambda ;
			m_position2[d] = p + displ ;
		}
		else
			m_position2[d] = m_position1[d] ;
	}
};


template<typename PFP>
class PassTwo: public  FunctorMapThreaded<typename PFP::MAP>
{
protected:
	VertexAttribute<typename PFP::VEC3>& m_position1;
	const VertexAttribute<typename PFP::VEC3>& m_position2;
	Algo::Surface::Selection::Collector_OneRing<PFP> m_collect;
	Algo::Surface::Filtering::FunctorAverage<typename PFP::VEC3, VERTEX> m_fa;


public:
	PassTwo( typename PFP::MAP& map, VertexAttribute<typename PFP::VEC3>& position1, const VertexAttribute<typename PFP::VEC3>& position2):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_position1(position1), m_position2(position2), m_collect(map), m_fa(position2)
	{}

	FunctorMapThreaded<typename PFP::MAP>* duplicate() const
	{
		PassTwo<PFP>* copy = new PassTwo<PFP>(this->m_map ,m_position1,m_position2);
		return reinterpret_cast< FunctorMapThreaded<typename PFP::MAP>* >(copy);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		const float mu = -0.6732 ;
		if(!this->m_map.isBoundaryVertex(d))
		{
			m_collect.collectBorder(d) ;
			m_fa.reset() ;
			m_collect.applyOnBorder(m_fa) ;
			typename PFP::VEC3 p = m_position2[d] ;
			typename PFP::VEC3 displ = m_fa.getAverage() - p ;
			displ *= mu ;
			m_position1[d] = p + displ ;
		}
		else
			m_position1[d] = m_position2[d] ;
	}

};



template <typename PFP>
void ParallelfilterTaubin(typename PFP::MAP& map, VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<typename PFP::VEC3>& position2, unsigned int nbloops, unsigned int nbth)
{
	PassOne<PFP> funcFront(map, position, position2);
	PassTwo<PFP> funcBack(map, position, position2);
	/*
	for (unsigned int i=0; i< nbloops; ++i)
	{
		Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funcFront, nbth, true, select);
		Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funcBack, nbth, true, select);
	}
	*/

	Algo::Parallel::foreach_cell2Pass<typename PFP::MAP,VERTEX>(map, funcFront, funcBack, nbloops, nbth, true);
}

}
}

#endif /* TAUBIN_H_ */
