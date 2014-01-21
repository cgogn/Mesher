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
#ifndef DART_BLOCK_H_
#define DART_BLOCK_H_

#include <cassert>
#include <vector>
#include <list>
#include "Topology/generic/dart.h"



namespace CGoGN
{
namespace Utils
{


class DartBlock
{
public:
	static const unsigned int NBBITS = 13;
	static const unsigned int SZBLK = 8192; // 4096 = 2^12  1024 = 2^10

protected:
	Dart m_data[SZBLK];
	unsigned int m_first;
	unsigned int m_last;
	int m_nb;

public:

	DartBlock():m_first(0),m_last(0),m_nb(0) {}
	
	inline void reset()
	{
		m_first = 0;
		m_last = 0;
		m_nb = 0;
	}	

	inline static unsigned int createIdx(unsigned int idb, unsigned int idd) { return (idb << NBBITS) | idd ;}

	inline static unsigned int idxBlock(unsigned int idx) { return idx >> NBBITS; }

	inline static unsigned int idxInBlock(unsigned int idx) { return idx & (SZBLK-1); }

	inline bool full() const { return m_last >= SZBLK;}

	inline bool empty() const { return m_first == m_last;}

	inline unsigned int size() const { return m_last - m_first;}

//	inline bool emptyNb() const { return m_nb == 0;}

	inline unsigned int pushBack(Dart d)
	{
		assert(!full());
		m_data[m_last] = d;
		return m_last++;
	}

	inline bool remove(unsigned int idx)
	{
		if (m_data[idx] == EMBNULL)
			return false;
		m_data[idx] = EMBNULL;
		return true;
	}

	inline Dart popFront()
	{
		return m_data[m_first++];
	}

	inline bool checkBeforePopFront()
	{
		while (!empty() && (m_data[m_first] == EMBNULL))	// skip removed values
		{
			m_first++;
		}
		return empty();
	}

};



//
//class DartQueuesTable
//{
//protected:
//	std::vector<DartBlock*> m_tableOfBlocks;
//	std::vector<unsigned int> m_freeBlocks;
//
//public:
//	DartQueuesTable()
//	{
//		m_tableOfBlocks.reserve(512);
//		m_freeBlocks.reserve(16);
//	}
//
//	~DartQueuesTable()
//	{
//		std::cout << "Size at destruction: "<< m_tableOfBlocks.size() << " / "<< m_freeBlocks.size() << std::endl;
//		for (std::vector<DartBlock*>::iterator it = m_tableOfBlocks.begin(); it !=  m_tableOfBlocks.end(); ++it)
//		{
//			delete *it;
//		}
//	}
//
//	inline DartBlock* newBlock(unsigned int& id)
//	{
//		if (!m_freeBlocks.empty())
//		{
//			id =  m_freeBlocks.back();
//			m_freeBlocks.pop_back();
//			m_tableOfBlocks[id]->reset();
//			return m_tableOfBlocks[id];
//		}
//		else
//		{
//			DartBlock* db = new DartBlock();
//			id = m_tableOfBlocks.size();
//			m_tableOfBlocks.push_back(db);
//			return db;
//		}
//	}
//
//	inline void freeBlock(unsigned int id)
//	{
//		m_freeBlocks.push_back(id);
//	}
//
//	inline DartBlock* getBlock(unsigned int id)
//	{
//		return m_tableOfBlocks[id];
//	}
//
//	void remove(unsigned int id)
//	{
//		unsigned int idb = DartBlock::idxBlock(id);
//		unsigned int idd = DartBlock::idxInBlock(id);
//		m_tableOfBlocks[idb]->remove(idd);
//	}
//
//	void dump() const
//	{
//		std::cout << "============= DartQueueTable =====================" << std::endl;
//		std::cout << "NB BLOCKS " << m_tableOfBlocks.size()<< std::endl;
//		std::cout << "NB FREE BLOCKS " << m_freeBlocks.size() << std::endl;
//		std::cout << "=============================================" << std::endl;
//	}
//};
//

class DartQueue
{
protected:
	std::vector<DartBlock*> m_tableOfBlocks;
	std::vector<unsigned int> m_freeBlocks;

protected:
	DartQueue* m_shared;

	std::list<unsigned int> m_listofBlocks;

	DartBlock* m_firstBlock;
	
	DartBlock* m_lastBlock;

	int m_nbElts;

	int m_nbTotal;

	int id;

public:
	DartQueue():
		m_shared(this),m_firstBlock(NULL),m_lastBlock(NULL),m_nbElts(0),m_nbTotal(0)
	{
		m_tableOfBlocks.reserve(512);
		m_freeBlocks.reserve(16);
	}

	DartQueue(DartQueue* dq):
		m_shared(dq),m_firstBlock(NULL),m_lastBlock(NULL),m_nbElts(0),m_nbTotal(0)
	{
	}


	~DartQueue()
	{
		if (m_shared == this)
			for (std::vector<DartBlock*>::iterator it = m_tableOfBlocks.begin(); it !=  m_tableOfBlocks.end(); ++it)
				delete *it;
	}

//	void clearBlocks()
//	{
//		for (std::vector<DartBlock*>::iterator it = m_shared->m_tableOfBlocks.begin(); it !=  m_shared->m_tableOfBlocks.end(); ++it)
//			delete *it;
//		m_tableOfBlocks.clear();
//		m_freeBlocks.clear();
	void clearBlocks()
	{
		unsigned int nb = m_shared->m_tableOfBlocks.size();
		m_shared->m_freeBlocks.resize(nb);
		for (unsigned int i=0; i<nb; ++i)
			m_shared->m_freeBlocks[i] = i;
	}

	void clear()
	{
		m_firstBlock	= NULL;
		m_lastBlock		= NULL;
		m_nbElts		= 0;
		m_nbTotal		= 0;
		m_listofBlocks.clear();
	}

	/// set sharing after construction (for table of DartQueue for example)
	inline void setShareBlockWith(DartQueue* dq)
	{
		m_shared = dq;
	}

	unsigned int pushBack(Dart d)
	{
		if (m_firstBlock == NULL)
		{
			unsigned int id;
			m_firstBlock = m_shared->newBlock(id);
			m_lastBlock = m_firstBlock;
			m_listofBlocks.push_back(id);
		}
		
		m_nbElts++;
		m_shared->m_nbTotal++;
		if (m_lastBlock->full())
		{
			unsigned int idb;
			m_lastBlock = m_shared->newBlock(idb);
			m_listofBlocks.push_back(idb);
			unsigned int idd  = m_lastBlock->pushBack(d);
			return DartBlock::createIdx(idb,idd);
		}
		// else
		unsigned int idd  = m_lastBlock->pushBack(d);
		return DartBlock::createIdx(m_listofBlocks.back(),idd);
	}

	Dart popFront()
	{
//		std::cout << long(this)<< "_PopFront ->";
		assert(!empty());
//		while (m_firstBlock->empty())
		while (m_firstBlock->checkBeforePopFront())
		{
			m_shared->freeBlock(m_listofBlocks.front());
			m_listofBlocks.pop_front();  // remove block from list of block of queue
			if (!m_listofBlocks.empty())
				m_firstBlock = m_shared->getBlock(m_listofBlocks.front());
			else 
				m_firstBlock = NULL;
		}

		Dart d = m_firstBlock->popFront();
		--m_nbElts;
		m_shared->m_nbTotal--;
		if (m_firstBlock->empty())
		{
			m_shared->freeBlock(m_listofBlocks.front());
			m_listofBlocks.pop_front();  // remove block from list of block of queue
			if (!m_listofBlocks.empty())
				m_firstBlock = m_shared->getBlock(m_listofBlocks.front());
			else 
				m_firstBlock = NULL;
		}
		return d;
	}

	inline int size()
	{
		return m_nbElts;
	}

	inline bool empty()
	{
		return m_nbElts <= 0;
	}

	void remove(unsigned int id)
	{
		unsigned int idb = DartBlock::idxBlock(id);
		unsigned int idd = DartBlock::idxInBlock(id);
		if (m_shared->m_tableOfBlocks[idb]->remove(idd))
		{
			m_nbElts--;
			m_shared->m_nbTotal--;
//			std::cout << long(this)<< "_Remove "<< idb << " / "<<idd<< std::endl;
		}
//		else
//		{
//			std::cout << long(this)<< "_Not_Remove "<< idb << " / "<<idd<< std::endl;
//		}


	}


	inline int nbTotal() const
	{
		return m_shared->m_nbTotal;
	}


	// TO CALL ONLY FROM m_SHARED !!
protected:
	inline DartBlock* newBlock(unsigned int& id)
	{
		if (!m_freeBlocks.empty())
		{
			id =  m_freeBlocks.back();
			m_freeBlocks.pop_back();
			m_tableOfBlocks[id]->reset();
			return m_tableOfBlocks[id];
		}
		else
		{
			DartBlock* db = new DartBlock();
			id = m_tableOfBlocks.size();
			m_tableOfBlocks.push_back(db);
			return db;
		}
	}

	inline void freeBlock(unsigned int id)
	{
		m_freeBlocks.push_back(id);
	}

	inline DartBlock* getBlock(unsigned int id)
	{
		return m_tableOfBlocks[id];
	}

	
};


}
}

#endif 
