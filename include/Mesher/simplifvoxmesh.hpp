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
#include "Algo/Parallel/parallel_foreach.h"

namespace CGoGN
{

namespace Mesher
{

template <typename PFP>
SimplifVoxMesh<PFP>::SimplifVoxMesh(typename PFP::MAP& m):
m_map(m),
m_marker(m)
{
	m_crits      = m.template addAttribute<unsigned long long, EDGE>("crits");
	m_curvatures = m.template addAttribute<float, VERTEX>("curvatures");
	m_bounds     = m.template addAttribute<unsigned char, VERTEX>("bounds");
//	m_valences   = m.template addAttribute<unsigned char, VERTEX>("valences");

	m_coef[0]=0.0f;
	m_coef[1]=0.0f;

	// init table of queues
	for (unsigned int i=1; i<NB_QUEUES; ++i)
		m_criterias[i].setShareBlockWith(&(m_criterias[0]));
	m_first_non_empty = 0;


	m_length_base = 1.0f;
}


template <typename PFP>
SimplifVoxMesh<PFP>::~SimplifVoxMesh()
{
	m_map.template removeAttribute<unsigned char, VERTEX>(m_bounds);
	m_map.template removeAttribute<float, VERTEX>(m_curvatures);
//	m_map.template removeAttribute<unsigned char, VERTEX>(m_valences);
	m_map.template removeAttribute< unsigned long long, EDGE>(m_crits);
}

template <typename PFP>
void SimplifVoxMesh<PFP>::setInitialLength(float l)
{
	m_length_base = l;
}


template <typename PFP>
unsigned int SimplifVoxMesh<PFP>::init(float v1, float v2)
{
	// now we need position
	m_positions  = m_map.template getAttribute<typename PFP::VEC3, VERTEX>("position");

	m_coef[0] = v1;
	m_coef[1] = v2;
	
	m_coef_base = (30.0f - 15.0f*v1) / m_length_base;

	m_criterias[0].clearBlocks();
	for (unsigned int i=0; i<NB_QUEUES;++i)
		m_criterias[i].clear();

	TraversorE<typename PFP::MAP> trav(m_map);
	unsigned int nb=0;
	for (Dart d = trav.begin(); d!= trav.end(); d= trav.next())
	{
		insertEdgeQueue(d);	
		if (m_map.isBoundaryEdge(d))
			nb++;
		else 
			nb+=2;
	}
	m_nbFaces = nb/3; 

	return m_nbFaces;
}




template <typename PFP>
unsigned char SimplifVoxMesh<PFP>::queueIndex(float length, float curv)
{
	float le = m_coef_base * length * (1.0f  + 5.0f*m_coef[0]*curv);
	unsigned int idx = (unsigned int)(round(le));
	// clamp to max number of queue
	if (idx >= NB_QUEUES)
		idx = NB_QUEUES-1;

	return idx;
}

template <typename PFP>
void SimplifVoxMesh<PFP>::insertEdgeQueue(Dart d)
{
	Dart d1 = m_map.phi1(d);

//	int ve = m_valences[d] + m_valences[d1];
	Dart e = m_map.phi2(d);
//	int vo = m_valences[m_map.phi_1(d)];
//	int voo= m_valences[m_map.phi_1(e)];

	typename PFP::VEC3 V = m_positions[d1] - m_positions[d];
	float length = V.norm();

/*
	if ((ve>12) || (vo<6) || (voo<6))
	{
		length *= 2.0f;
	}
*/
	float curv = (m_curvatures[d] + m_curvatures[m_map.phi1(d)])/2.0f;
	unsigned char idq = queueIndex(length, curv);

	if (idq < m_first_non_empty)
		m_first_non_empty = idq;

	unsigned int id = m_criterias[idq].pushBack(d);

	m_marker.unmark(d);
	m_crits[d] = ((unsigned long long)(idq) << 32) | id;
}

template <typename PFP>
void SimplifVoxMesh<PFP>::removeEdgeQueue(Dart d)
{
	if (m_marker.isMarked(d))
			return;

	unsigned long long ll = m_crits[d];
	unsigned int id = ll & 0xffffffff;
	unsigned char q = ll >> 32;

	m_criterias[q].remove(id);
	m_marker.mark(d);

}




template <typename PFP>
inline typename PFP::VEC3 SimplifVoxMesh<PFP>::computeNewPosition(Dart d, float& curvature)
{
	Dart e = m_map.phi1(d);

	if ((m_bounds[d] == 0) && (m_bounds[e] != 0))
	{
		curvature = m_curvatures[e];
		return m_positions[e];
	}

	if ((m_bounds[d] != 0) && (m_bounds[e] == 0))
	{
		curvature = m_curvatures[d];
		return m_positions[d];
	}

	curvature = (m_curvatures[d] + m_curvatures[e])/2.0f;

	if (fabs(m_curvatures[d]-m_curvatures[e]) < 0.1f)
		return  (m_positions[d]+m_positions[e])/2.0f;

//	return (m_curvatures[d]*m_positions[d] + m_curvatures[e]*m_positions[e])/(m_curvatures[d]+m_curvatures[e]);
	return (m_curvatures[d]*m_positions[d] + m_curvatures[e]*m_positions[e])/(2.0f*curvature);
}


template <typename PFP>
Dart SimplifVoxMesh<PFP>::getFirstEdgeToCollapse()
{
	while ((m_first_non_empty<NB_QUEUES)  && (m_criterias[m_first_non_empty].empty()))
		m_first_non_empty++;

	Dart d = m_criterias[m_first_non_empty].popFront();
	m_marker.mark(d);
	return d;

}

template <typename PFP>
void SimplifVoxMesh<PFP>::collapseEdge()
{
	Dart d = getFirstEdgeToCollapse();

	if (! m_map.edgeCanCollapse(d))
	{
		m_marker.mark(d);
		removeEdgeQueue(d);
		return;
	}

	Dart e = m_map.phi1(d);

	float curv;
	typename PFP::VEC3 P = computeNewPosition(d, curv);

//	unsigned short val = m_valences[d] + m_valences[e] - 4;

	e = d;
	do
	{
		removeEdgeQueue(e);
		e = m_map.alpha1(e);
	} while (e!=d);

	d = m_map.phi2(d);
	e = m_map.alpha1(d);
	do
	{
		removeEdgeQueue(e);
		e = m_map.alpha1(e);
	} while (e!=d);

	Dart f = m_map.collapseEdge(d,true);
	m_positions[f] = P;
	m_curvatures[f] = curv;
//	m_valences[f] = (unsigned char)(val);

	e = f;
	do
	{
		insertEdgeQueue(e);
		e = m_map.alpha1(e);
	} while (e!=f);

	m_nbFaces -= 2;
}


template <typename PFP>
void SimplifVoxMesh<PFP>::until(unsigned int nbFaces)
{
	while ((m_criterias[0].nbTotal() > 0) && (m_nbFaces > nbFaces))
	{
		collapseEdge();
		
	}
}

}
}

