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

#ifndef PREP_VOXMESH_H_
#define PREP_VOXMESH_H_

#include "Topology/generic/parameters.h"
#include "Topology/map/map2.h"
#include "Topology/generic/attributeHandler.h"
#include "Algo/MC/marchingcube.h"
#include "Algo/Parallel/parallel_foreach.h"

namespace CGoGN
{
namespace Mesher
{
template <typename PFP>
void smoothNormals(typename PFP::MAP& map, VertexAttribute<Geom::Vec3f>& normals);


template <typename PFP>
void  smoothCurvature(typename PFP::MAP& map, const VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors);

template <typename PFP>
void  smoothCurvature(typename PFP::MAP& map, const VertexAttribute<float>& curvatures);


template <typename PFP>
void  parallelsmoothCurvature(typename PFP::MAP& map, VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors, int nbth);

				
template <typename PFP>
void  parallelsmoothCurvature(typename PFP::MAP& map, VertexAttribute<float>& curvatures, int nbth);



template <typename PFP,typename DATATYPE>
void parallelComputeCurvatureValence( typename PFP::MAP& map,
										const VertexAttribute<typename PFP::VEC3>& position,
										VertexAttribute<float>& curvature,
										VertexAttribute<Geom::Vec3f>& colors,
										VertexAttribute<unsigned char>& valences,
										Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val,
										float vx, float vy, float vz, int nbth);

template <typename PFP,typename DATATYPE>
void parallelComputeCurvature( typename PFP::MAP& map,
								const VertexAttribute<typename PFP::VEC3>& position,
								VertexAttribute<float>& curvature,
								VertexAttribute<Geom::Vec3f>& colors,
								Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val,
								float vx, float vy, float vz, int nbth);


template <typename PFP>
void diffuseCurvature(typename PFP::MAP& map, Dart d, VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors);


template <typename PFP>
void diffuseCurvature(typename PFP::MAP& map, Dart d, VertexAttribute<float>& curvatures);

}
}

#include "prep_mesh.hpp"

#endif

