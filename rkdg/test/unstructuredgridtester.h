/**
 * 网格测试类，该类测试网格的信息是否正确
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#pragma once

#include "../inc/cppstdheaders.h"
#include "../inc/unstructuredgrid.h"

using namespace std;

/** 
 * 网格测试类，负责测试网格的信息是否正确
 */
class CUnstructuredGridTester: public CUnstructuredGrid
{
	// 属性
	public:

	protected:

	// 方法
	public:
		/** 构造函数 */
		CUnstructuredGridTester();

		/**
		 * 测试三角形邻居
		 * @param[in] output string 邻居信息的输出文件
		 */
		void testNeighbour(const string& output="output/neighbour.dat");

		/**
		 * 测试三角形共享边
		 * @param[in] output string 共享边信息的输出文件
		 */
		void testSharedEdge(const string& output="output/sharededge.dat");

		/**
		 * 测试三角形标记
		 * @param[in] output string 标记信息的输出文件
		 */
		void testTriangleFlag(const string& output="output/triangleflag.dat");

		/**
		 * 测试三角形面积
		 * @param[in] output string 三角形面积输出文件
		 */
		void testArea(const string& output="output/area.dat");

		/**
		 * 测试三角形周长
		 * @param[in] output string 三角形周长输出文件
		 */
		void testPerimeter(const string& output="output/perimeter.dat");


		/**
		 * 测试三角形边外法向量
		 * @param[in] output string 三角形外法向量输出文件
		 */
		void testOuterNormalVector(const string& output="output/outernormalvector.dat");

		/**
		 * 测试三角形基函数在体高斯积分节点上的值
		 */
		void testVolumeBfValue(void);


		/**
		 * 测试三角形基函数的质量系数
		 * @param[in] output string 三角形基函数质量系数输出文件
		 */
		void testMassCoeff(const string& output="output/masscoeff.dat");

		/**
		 * 测试三角形基函数导数在体高斯积分节点上的值
		 */
		void testVolumeBdfValue(void);

		/**
		 * 测试三角形基函数在边高斯积分节点上的值
		 */
		void testEdgeBfValue(void);

		/**
		 * 测试三角形体积分的权重
		 * @param[in] output string 三角形体积分的高斯权重输出文件
		 */
		void testVolumeGaussWeight(const string& output="output/volumegaussweight.dat");

		/**
		 * 测试三角边体积分的权重
		 * @param[in] output string 三角形边积分的高斯权重输出文件
		 */
		void testEdgeGaussWeight(const string& output="output/edgegaussweight.dat");


		/**
		 * 输出二维变量到文件查看
		 * @param[in] rnum 变量行个数 
		 * @param[in] cnum 变量列个数
		 * @param[in] var double[] 需要输出的变量
		 * @param[in] output string 输出的结果文件
		 */
		template<typename T>
		void output2DVar(int rnum, int cnum, const T *var, const string &output);

		/** 析构函数 */
		~CUnstructuredGridTester();

	protected:

};