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
#ifndef SIMPLIFVOXMESH_H_
#define SIMPLIFVOXMESH_H_


#include <set>

#include "Topology/generic/parameters.h"
#include "Topology/map/map2.h"
#include "Topology/generic/attributeHandler.h"
#include "Topology/generic/functor.h"
#include "Topology/map/embeddedMap2.h"
#include "Algo/MC/marchingcube.h"
#include "Algo/MC/marchingcubeGen.h"

#include "dartQueueBlock.h"


namespace CGoGN
{

namespace Mesher
{

template <typename PFP>
class SimplifVoxMesh
{
protected:
	static const unsigned int NB_QUEUES = 256;

	Utils::DartQueue m_criterias[NB_QUEUES];
	unsigned int m_first_non_empty;

	std::vector<Dart> m_edgeBuffer;

	typename PFP::MAP& m_map;
	VertexAttribute<typename PFP::VEC3> m_positions;
	EdgeAttribute<unsigned long long> m_crits;
	VertexAttribute<float> m_curvatures;
	VertexAttribute<unsigned char> m_valences;
	VertexAttribute<unsigned char> m_bounds;

	float m_coef[2];
	CellMarker<EDGE> m_marker;

	unsigned int m_nbFaces;

	float m_length_base;
	
	float m_coef_base;

	/// compute the queue index from the length
	unsigned char queueIndex(float length, float curv);

	/// insert an edge in the queue table
	void insertEdgeQueue(Dart d);

	/// remove an edge from the queue table
	void removeEdgeQueue(Dart d);

	/// find the first edge to collapse in the queue table.
	Dart getFirstEdgeToCollapse();

public:
	SimplifVoxMesh(typename PFP::MAP& m);

	unsigned int init(float v1, float v2);
	
	void setInitialLength(float l);

	~SimplifVoxMesh();

	/**
	 * compute new position of vertex of edge to collapse
	 */
	typename PFP::VEC3 computeNewPosition(Dart d, float& curvature);

	/**
	 * collapse edge and update all attributes and multiset.
	 */
	void collapseEdge();

	/**
	 * simplifcation until a nb of face is reached
	 */
	void until(unsigned int nbFaces);
	
	VertexAttribute<unsigned char>& boundAttrib() { return m_bounds;}

	VertexAttribute<float>& curvatureAttrib() { return m_curvatures;}

	VertexAttribute<unsigned char>& valenceAttrib() { return m_valences;}

};

}
}

#include "simplifvoxmesh.hpp"

#endif /* SIMPLIFVOXMESH_H_ */
