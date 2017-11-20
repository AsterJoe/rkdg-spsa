/**
 * 二维非结构网格类
 * 该类为流场解法器提供网格接口
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
 * 二维非结构网格类
 * 该类为流场解法器提供网格接口，同时提供网格单元的详细信息
 * 该类主要包含的信息： 网格顶点，网格的边，三角形组成以及详细的三角形信息（邻居，基函数等）
 */
class CUnstructuredGrid {
	// 属性
	public:
		string config_file;				/**< 网格配置文件 */
		string output_filename;			/**< 网格输出文件 */

		vector<CVertice2D> vertice;		/**< 网格顶点 */
		vector<CEdge> edge;				/**< 网格边 */
		vector<int> tri_vertice;		/**< 三角形的顶点组成 */
		vector<int> tri_edge;			/**< 三角形的边组成 */

		vector<int> farfiled;
		vector<int> wall;
		vector<int> upper_node;
		vector<int> lower_node;
		vector<int> wall_elem;
		int *ver_flag;
		
		// 此三项需要被复制到GPU上
		int *tri_neighbour;				/**< 三角形邻居信息 */
		int *tri_sharedEdge;			/**< 三角形共享边信息 */
		int *tri_flag;					/**< 三角形标记（边界条件等信息） */

		int *airfoil;                   /**< 机翼表面的单元>*/
		

		CTriangleInfo triangle_infos;	/**< 三角形额外信息 */
		
	protected:
		int _vertice_num;				/**< 顶点数目 */
		int _edge_num;					/**< 边数目 */
		int _triangle_num;				/**< 三角形数目 */
		int _ghost_triangle_num;		/**< 虚拟三角单元数目 */
		int _cell_num;					/**< 网格总总单元数目（包括虚拟网格） */
		
	private:
	
	// 方法
	public:
		/** 默认构造函数 */
		CUnstructuredGrid();
		
		/**
		 * 解析虚拟三角形数目, 如果某个单元邻居编号为-1，则需要生成虚拟网格
		 * @param[in] trineigh_filename string 三角形邻居编号文件
		 */
		void parseGhostNum(const string& trineigh_filename);

		/** 初始化网格数据 */
		void initializeGrid(void);
		
		/** 初始化三角形数据 */
		void initializeTriangleInfos(void);
		
		/** 
		 * 导入三角形顶点编号
		 * @param[in] trivertice_file string 三角形顶点编号文件
		 */
		void readTriangleVertice(const string& trivertice_file);
	
		/** 
		* 导入相邻三角形编号
		* @param[in] trineigh_file string 三角形邻居编号文件
		*/
		void readTriangleNeighbour(const string& trineigh_file);
	
		/**
		* 导入三角形边编号
		* @param triedge_file string 三角形边组成文件
		*/
		void readTriangleEdge(const string& triedge_file);
	
		/** 
		* 导入顶点坐标
		* @param[in] vertice_file string 顶点坐标文件
		*/
		void readVertice(const string& vertice_file);
	
		/** 
		* 导入边的顶点编号
		* @param[in] edge_file 边的顶点编号文件
		*/
		void readEdgeVertice(const string& edge_file);
	
		/** 
		* 初始化虚拟单元
		*/
		void initGhostGrid(void);
		
		/** @return integer 三角形数目 */
		int getTriangleNumber(void) const;

		/** @return integer 边数目 */
		int getEdgeNumber(void) const;
		
		/** @return integer 顶点数目 */
		int getVerticeNumber(void) const;
		
		/** @return integer 虚拟三角形数目 */
		int getGhostTriangleNumber(void) const;

		/** @return integer 返回总单元数目 */
		int getCellNumber(void) const;

		/** 初始化党员之间共享边 */
		void initSharedEdge(void);

		/** 标记网格边界条件 */
		void markBoundaryTriangles(void);
		
		/** 输出网格 */
		void outputGrid(void) const;

		/** 测试网格单元是不是都是逆时针编号 */
		void testTrianglesAntiwise(void) const;
		
		/**
		 * 输出计算网格(包括虚拟网格)
		 * @param[in] filename string 虚拟网格输出文件
		 */
		void outputGridWithGhostCells(const string& filename) const;

		void getUpperNode2D(void);

		void sortUpperNode2D(vector<int>* vec);

		double calArea(void);
		
		/** 析构函数 */
		~CUnstructuredGrid();
	protected:
	
	private:
};