/**
 * ��CCUDARkdgSolver�Ĳ�����
 * ������Ҫ�������CCUDARkdgSolver�ຯ������ȷ���
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#pragma once

#include "../inc/cudarkdgsolver.h"
#include "unstructuredgridtester.h"

class CCUDARkdgSolverTester: public CCUDARkdgSolver
{
	// ����
	public:

	//	CUnstructuredGridTester grid;  /**< ʹ�ò������������ */
	
	protected:
		cudaError_t _cuerror;	/** CUDA���� */
	
	// ����
	public:
		/** Ĭ�Ϲ��캯�� */
		CCUDARkdgSolverTester();
		
		/** ׼���ò���ǰ�Ļ��� */
		void initEnvironment(void);

		/**
		 * �����ά�������ļ��鿴
		 * @param[in] rnum �����и��� 
		 * @param[in] cnum �����и���
		 * @param[in] var double[] ��Ҫ����ı���
		 * @param[in] output string ����Ľ���ļ�
		 */
		void output2DVar(int rnum, int cnum, const double *var, const string &output);

		/** ִ�в��� */
		void runTest(void);

		/** �������� */
		~CCUDARkdgSolverTester();
	protected:
	
};