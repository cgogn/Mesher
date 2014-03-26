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

namespace CGoGN
{
namespace Mesher
{

template <typename PFP>
void smoothNormals(typename PFP::MAP& map, VertexAttribute<Geom::Vec3f>& normals)
{
	VertexAttribute<Geom::Vec3f> normals2 = map.template addAttribute<Geom::Vec3f,VERTEX>("normals2");
	CellMarker<VERTEX> m(map);
	for (Dart d=map.begin(); d!=map.end();map.next(d))
	{
		if (!m.isMarked(d))
		{
			// store for omp optimized loop bellow
			Geom::Vec3f N = normals[d];
			Dart dd = d;
			do
			{
				m.mark(dd);
				dd = map.alpha1(dd);
				N += normals[map.phi1(dd)];
			} while (dd!=d);
			N.normalize();
			normals2[d] = N;
			m.mark(d);
		}
	}
	map.swapAttributes(normals,normals2);
	map.template removeAttribute<Geom::Vec3f,VERTEX>(normals2);
}


template <typename PFP>
void  smoothCurvature(typename PFP::MAP& map, const VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors)
{
	VertexAttribute<float>  curv2 = map.template addAttribute<float,VERTEX>("curv2");

	CellMarker<VERTEX> m(map);
	for (Dart d=map.begin(); d!=map.end();map.next(d))
	{
		if (!m.isMarked(d))
		{
			float c = curvatures[d];
			int nb = 1;
			Dart dd = d;
			do
			{
				m.mark(dd);
				dd = map.alpha1(dd);
				c += curvatures[map.phi1(dd)];
				nb++;
			} while (dd!=d);
			curv2[d] = c/float(nb);
			colors[d] = Geom::Vec3f(curv2[d],0.0f,1.0f-curv2[d]);
			m.mark(d);
		}
	}
	map.swapAttributes(curvatures,curv2);
	map.template removeAttribute<float,VERTEX>(curv2);
}

template <typename PFP>
void  smoothCurvature(typename PFP::MAP& map, const VertexAttribute<float>& curvatures)
{
	VertexAttribute<float>  curv2 = map.template addAttribute<float,VERTEX>("curv2");

	CellMarker<VERTEX> m(map);
	for (Dart d=map.begin(); d!=map.end();map.next(d))
	{
		if (!m.isMarked(d))
		{
			float c = curvatures[d];
			int nb = 1;
			Dart dd = d;
			do
			{
				m.mark(dd);
				dd = map.alpha1(dd);
				c += curvatures[map.phi1(dd)];
				nb++;
			} while (dd!=d);
			curv2[d] = c/float(nb);
			m.mark(d);
		}
	}
	map.swapAttributes(curvatures,curv2);
	map.template removeAttribute<float,VERTEX>(curv2);
}




template<typename PFP>
class ParFuncSmoothCurvCol: public  FunctorMapThreaded<typename PFP::MAP>
{
protected:
	const VertexAttribute<float>& m_curv1;
	VertexAttribute<float>& m_curv2;
	VertexAttribute<typename PFP::VEC3>& m_colors;

public:
	ParFuncSmoothCurvCol( typename PFP::MAP& map, const VertexAttribute<float>& curv1, VertexAttribute<float>& curv2, VertexAttribute<typename PFP::VEC3>& colors):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_curv1(curv1), m_curv2(curv2), m_colors(colors)
	{}

	FunctorMapThreaded<typename PFP::MAP>* duplicate()
	{
		ParFuncSmoothCurvCol<PFP>* copy = new ParFuncSmoothCurvCol<PFP>(this->m_map ,m_curv1,m_curv2,m_colors);
		return reinterpret_cast< FunctorMapThreaded<typename PFP::MAP>* >(copy);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		float c = m_curv1[d];
		int nb = 1;
		Dart dd = d;
		do
		{
			dd = this->m_map.alpha1(dd);
			c += m_curv1[this->m_map.phi1(dd)];
			nb++;
		} while (dd!=d);
		m_curv2[d] = c/float(nb);
		m_colors[d] = typename PFP::VEC3(m_curv2[d],0.0f,1.0f-m_curv2[d]);
	}
};

template <typename PFP>
void  parallelsmoothCurvature(typename PFP::MAP& map, VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors, int nbth)
{
	VertexAttribute<float>  curv2 = map.template addAttribute<float,VERTEX>("curv2");

	ParFuncSmoothCurvCol<PFP> funct(map, curvatures, curv2, colors);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth);
	map.swapAttributes(curvatures,curv2);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth);
	map.swapAttributes(curvatures,curv2);
	
	map.template removeAttribute<float,VERTEX>(curv2);
}





template<typename PFP>
class ParFuncSmoothCurv: public  FunctorMapThreaded<typename PFP::MAP>
{
protected:
	const VertexAttribute<float>& m_curv1;
	VertexAttribute<float>& m_curv2;

public:
	ParFuncSmoothCurv( typename PFP::MAP& map, const VertexAttribute<float>& curv1, VertexAttribute<float>& curv2):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_curv1(curv1), m_curv2(curv2)
	{}

	FunctorMapThreaded<typename PFP::MAP>* duplicate()
	{
		ParFuncSmoothCurv<PFP>* copy = new ParFuncSmoothCurv<PFP>(this->m_map ,m_curv1,m_curv2);
		return reinterpret_cast< FunctorMapThreaded<typename PFP::MAP>* >(copy);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		float c = m_curv1[d];
		int nb = 1;
		Dart dd = d;
		do
		{
			dd = this->m_map.alpha1(dd);
			c += m_curv1[this->m_map.phi1(dd)];
			nb++;
		} while (dd!=d);
		m_curv2[d] = c/float(nb);
	}
};

template <typename PFP>
void  parallelsmoothCurvature(typename PFP::MAP& map, VertexAttribute<float>& curvatures, int nbth)
{
	VertexAttribute<float>  curv2 = map.template addAttribute<float,VERTEX>("curv2");

	ParFuncSmoothCurv<PFP> funct(map, curvatures, curv2);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth);
	map.swapAttributes(curvatures,curv2);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth);
	map.swapAttributes(curvatures,curv2);
	
	map.template removeAttribute<float,VERTEX>(curv2);
}






template<typename PFP, typename DATATYPE>
class ParallelFunctorCurvScaleCol: public FunctorMapThreaded<typename PFP::MAP>
{
	typedef typename PFP::MAP MAP;

protected:
	MAP& m_map;
	const VertexAttribute<typename PFP::VEC3>& m_position;
	VertexAttribute<float>& m_curvature;
	VertexAttribute<Geom::Vec3f>& m_colors;
	VertexAttribute<unsigned char>& m_valences;
	Algo::Surface::MC::Image<DATATYPE>* m_img;
	std::vector<int> m_sphere;
	int m_rad;
	DATATYPE m_val;
	float m_vx;
	float m_vy;
	float m_vz;

public:
	ParallelFunctorCurvScaleCol(MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								VertexAttribute<Geom::Vec3f>& colors, VertexAttribute<unsigned char>& valences,
								Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_map(map), m_position(position), m_curvature(curvature), m_colors(colors), m_valences(valences),
		m_img(im), m_rad(rad), m_val(val),m_vx(vx),m_vy(vy),m_vz(vz)
	{
		if (rad == 0)
			rad = 1;
		im->createMaskOffsetSphere(m_sphere,rad);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		const typename PFP::VEC3& pos = m_position[d];
//		DATATYPE* ptr = m_img->getVoxelPtr(int(round(pos[0]/m_vx)), int(round(pos[1]/m_vy)), int(round(pos[2]/m_vz)));

		float  curv = m_img->computeCurvatureCount(pos[0],pos[1],pos[2], m_sphere, m_val);
		m_curvature[d] = curv;
		m_colors[d] = Geom::Vec3f(curv,0.0f,1.0f-curv);
		unsigned char val=0;
		Dart dd = d;
		do
		{
			val++;
			dd = m_map.alpha1(dd);
		} while (dd!=d);
		m_valences[d] = val;
	}
};


template <typename PFP,typename DATATYPE>
void parallelComputeCurvatureValence( typename PFP::MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								VertexAttribute<Geom::Vec3f>& colors, VertexAttribute<unsigned char>& valences, Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz, int nbth)
	{
	ParallelFunctorCurvScaleCol<PFP,DATATYPE> funct(map, position, curvature, colors, valences, im, rad, val, vx, vy, vz);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth, false);

}





template<typename PFP, typename DATATYPE>
class ParallelFunctorCurvScale: public FunctorMapThreaded<typename PFP::MAP>
{
	typedef typename PFP::MAP MAP;

protected:
	MAP& m_map;
	const VertexAttribute<typename PFP::VEC3>& m_position;
	VertexAttribute<float>& m_curvature;
	VertexAttribute<unsigned char>& m_valences;
	Algo::Surface::MC::Image<DATATYPE>* m_img;
	std::vector<int> m_sphere;
	int m_rad;
	DATATYPE m_val;
	float m_vx;
	float m_vy;
	float m_vz;

public:
	ParallelFunctorCurvScale(MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								VertexAttribute<unsigned char>& valences,
								Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_map(map), m_position(position), m_curvature(curvature), m_valences(valences),
		m_img(im), m_rad(rad), m_val(val),m_vx(vx),m_vy(vy),m_vz(vz)
	{
		if (rad == 0)
			rad = 1;
		im->createMaskOffsetSphere(m_sphere,rad);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		const typename PFP::VEC3& pos = m_position[d];
//		DATATYPE* ptr = m_img->getVoxelPtr(int(round(pos[0]/m_vx)), int(round(pos[1]/m_vy)), int(round(pos[2]/m_vz)));

		float  curv = m_img->computeCurvatureCount(pos[0],pos[1],pos[2], m_sphere, m_val);
		m_curvature[d] = curv;
		unsigned char val=0;
		Dart dd = d;
		do
		{
			val++;
			dd = m_map.alpha1(dd);
		} while (dd!=d);
		m_valences[d] = val;
	}
};


template <typename PFP,typename DATATYPE>
void parallelComputeCurvatureValence( typename PFP::MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								VertexAttribute<unsigned char>& valences, Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz, int nbth)
	{
	ParallelFunctorCurvScale<PFP,DATATYPE> funct(map, position, curvature, valences, im, rad, val, vx, vy, vz);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth, false);

}







template<typename PFP, typename DATATYPE>
class ParallelFunctorCurv: public FunctorMapThreaded<typename PFP::MAP>
{
	typedef typename PFP::MAP MAP;

protected:
	MAP& m_map;
	const VertexAttribute<typename PFP::VEC3>& m_position;
	VertexAttribute<float>& m_curvature;
	Algo::Surface::MC::Image<DATATYPE>* m_img;
	std::vector<int> m_sphere;
	int m_rad;
	DATATYPE m_val;
	float m_vx;
	float m_vy;
	float m_vz;

public:
	ParallelFunctorCurv(MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz):
		FunctorMapThreaded<typename PFP::MAP>(map),
		m_map(map), m_position(position), m_curvature(curvature),
		m_img(im), m_rad(rad), m_val(val),m_vx(vx),m_vy(vy),m_vz(vz)
	{
		if (rad == 0)
			rad = 1;
		im->createMaskOffsetSphere(m_sphere,rad);
	}

	void run(Dart d, unsigned int /*threadID*/)
	{
		const typename PFP::VEC3& pos = m_position[d];
//		DATATYPE* ptr = m_img->getVoxelPtr(int(round(pos[0]/m_vx)), int(round(pos[1]/m_vy)), int(round(pos[2]/m_vz)));

		float  curv = m_img->computeCurvatureCount(pos[0],pos[1],pos[2], m_sphere, m_val);
		m_curvature[d] = curv;
	}
};


template <typename PFP,typename DATATYPE>
void parallelComputeCurvature( typename PFP::MAP& map, const VertexAttribute<typename PFP::VEC3>& position, VertexAttribute<float>& curvature,
								 Algo::Surface::MC::Image<DATATYPE>* im, int rad, DATATYPE val, float vx, float vy, float vz, int nbth)
{
	ParallelFunctorCurv<PFP,DATATYPE> funct(map, position, curvature, im, rad, val, vx, vy, vz);
	Algo::Parallel::foreach_cell<typename PFP::MAP,VERTEX>(map, funct, nbth, false);

}










template <typename PFP>
void diffuseCurvature(typename PFP::MAP& map, Dart d, VertexAttribute<float>& curvatures, VertexAttribute<Geom::Vec3f>& colors)
{
	float c = curvatures[d];
	colors[d] = Geom::Vec3f(curvatures[d],0.0f,1.0f-curvatures[d]);
	Dart dd = d;
	do
	{
		dd = map.alpha1(dd);
		Dart e = map.phi1(dd);
		curvatures[e] = (c+curvatures[e])/2.0f;
		colors[e] = Geom::Vec3f(curvatures[e],0.0f,1.0f-curvatures[e]);
	} while (dd!=d);
}


template <typename PFP>
void diffuseCurvature(typename PFP::MAP& map, Dart d, VertexAttribute<float>& curvatures)
{
	float c = curvatures[d];
	Dart dd = d;
	do
	{
		dd = map.alpha1(dd);
		Dart e = map.phi1(dd);
		curvatures[e] = (c+curvatures[e])/2.0f;
	} while (dd!=d);
}


}
}


