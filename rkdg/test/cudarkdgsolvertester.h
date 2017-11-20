/**
 * 类CCUDARkdgSolver的测试类
 * 该类主要负责测试CCUDARkdgSolver类函数的正确与否
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */

#pragma once

#include "../inc/cudarkdgsolver.h"
#include "unstructuredgridtester.h"

class CCUDARkdgSolverTester: public CCUDARkdgSolver
{
	// 属性
	public:

	//	CUnstructuredGridTester grid;  /**< 使用测试网格类代替 */
	
	protected:
		cudaError_t _cuerror;	/** CUDA错误 */
	
	// 方法
	public:
		/** 默认构造函数 */
		CCUDARkdgSolverTester();
		
		/** 准备好测试前的环境 */
		void initEnvironment(void);

		/**
		 * 输出二维变量到文件查看
		 * @param[in] rnum 变量行个数 
		 * @param[in] cnum 变量列个数
		 * @param[in] var double[] 需要输出的变量
		 * @param[in] output string 输出的结果文件
		 */
		void output2DVar(int rnum, int cnum, const double *var, const string &output);

		/** 执行测试 */
		void runTest(void);

		/** 析构函数 */
		~CCUDARkdgSolverTester();
	protected:
	
};