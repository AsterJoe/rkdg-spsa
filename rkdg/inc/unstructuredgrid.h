/**
 * ��ά�ǽṹ������
 * ����Ϊ�����ⷨ���ṩ����ӿ�
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 */
 
#pragma once

#include "cppstdheaders.h"
#include "vertice2d.h"
#include "edge.h"
#include "config.h"
#include "triangleinfo.h"
#include "defines.h"
#include<set>
using namespace std;

/**
 * ��ά�ǽṹ������
 * ����Ϊ�����ⷨ���ṩ����ӿڣ�ͬʱ�ṩ����Ԫ����ϸ��Ϣ
 * ������Ҫ��������Ϣ�� ���񶥵㣬����ıߣ�����������Լ���ϸ����������Ϣ���ھӣ��������ȣ�
 */
class CUnstructuredGrid {
	// ����
	public:
		string config_file;				/**< ���������ļ� */
		string output_filename;			/**< ��������ļ� */

		vector<CVertice2D> vertice;		/**< ���񶥵� */
		vector<CEdge> edge;				/**< ����� */
		vector<int> tri_vertice;		/**< �����εĶ������ */
		vector<int> tri_edge;			/**< �����εı���� */

		vector<int> farfiled;
		vector<int> wall;
		vector<int> upper_node;
		vector<int> lower_node;
		vector<int> wall_elem;
		int *ver_flag;
		
		// ��������Ҫ�����Ƶ�GPU��
		int *tri_neighbour;				/**< �������ھ���Ϣ */
		int *tri_sharedEdge;			/**< �����ι������Ϣ */
		int *tri_flag;					/**< �����α�ǣ��߽���������Ϣ�� */

		int *airfoil;                   /**< �������ĵ�Ԫ>*/
		

		CTriangleInfo triangle_infos;	/**< �����ζ�����Ϣ */
		
	protected:
		int _vertice_num;				/**< ������Ŀ */
		int _edge_num;					/**< ����Ŀ */
		int _triangle_num;				/**< ��������Ŀ */
		int _ghost_triangle_num;		/**< �������ǵ�Ԫ��Ŀ */
		int _cell_num;					/**< �������ܵ�Ԫ��Ŀ�������������� */
		
	private:
	
	// ����
	public:
		/** Ĭ�Ϲ��캯�� */
		CUnstructuredGrid();
		
		/**
		 * ����������������Ŀ, ���ĳ����Ԫ�ھӱ��Ϊ-1������Ҫ������������
		 * @param[in] trineigh_filename string �������ھӱ���ļ�
		 */
		void parseGhostNum(const string& trineigh_filename);

		/** ��ʼ���������� */
		void initializeGrid(void);
		
		/** ��ʼ������������ */
		void initializeTriangleInfos(void);
		
		/** 
		 * ���������ζ�����
		 * @param[in] trivertice_file string �����ζ������ļ�
		 */
		void readTriangleVertice(const string& trivertice_file);
	
		/** 
		* �������������α��
		* @param[in] trineigh_file string �������ھӱ���ļ�
		*/
		void readTriangleNeighbour(const string& trineigh_file);
	
		/**
		* ���������α߱��
		* @param triedge_file string �����α�����ļ�
		*/
		void readTriangleEdge(const string& triedge_file);
	
		/** 
		* ���붥������
		* @param[in] vertice_file string ���������ļ�
		*/
		void readVertice(const string& vertice_file);
	
		/** 
		* ����ߵĶ�����
		* @param[in] edge_file �ߵĶ������ļ�
		*/
		void readEdgeVertice(const string& edge_file);
	
		/** 
		* ��ʼ�����ⵥԪ
		*/
		void initGhostGrid(void);
		
		/** @return integer ��������Ŀ */
		int getTriangleNumber(void) const;

		/** @return integer ����Ŀ */
		int getEdgeNumber(void) const;
		
		/** @return integer ������Ŀ */
		int getVerticeNumber(void) const;
		
		/** @return integer ������������Ŀ */
		int getGhostTriangleNumber(void) const;

		/** @return integer �����ܵ�Ԫ��Ŀ */
		int getCellNumber(void) const;

		/** ��ʼ����Ա֮�乲��� */
		void initSharedEdge(void);

		/** �������߽����� */
		void markBoundaryTriangles(void);
		
		/** ������� */
		void outputGrid(void) const;

		/** ��������Ԫ�ǲ��Ƕ�����ʱ���� */
		void testTrianglesAntiwise(void) const;
		
		/**
		 * �����������(������������)
		 * @param[in] filename string ������������ļ�
		 */
		void outputGridWithGhostCells(const string& filename) const;

		void getUpperNode2D(void);

		void sortUpperNode2D(vector<int>* vec);

		double calArea(void);
		
		/** �������� */
		~CUnstructuredGrid();
	protected:
	
	private:
};