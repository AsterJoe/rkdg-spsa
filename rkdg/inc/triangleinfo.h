/**
 * ��������Ϣ�࣬���ౣ���˳����������֮��Ķ�����Ϣ��
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once

#include "cppstdheaders.h"
#include "gauss.h"
#include "vertice2d.h"
#include "defines.h"
using namespace std;

/**
 * ��������Ϣ�࣬���ౣ������㡢��֮��Ķ�����Ϣ
 * ���ģ��ܳ���������ⷨ��������Ϣ���������ڸ�����
 */
class CTriangleInfo{
	// ����
	public:
		
		vector<CVertice2D> barycenter;				/**< ���������� */
		
		vector<CVertice2D> edge_middle_vertice;		/**< ���е� */
		
		vector<double> radius_of_inscribed_circle;	/**< ����Բ�뾶 */
		
		
		vector<CVertice2D> vol_gauss_vertice; 		/**< ����ֵĸ�˹ϵ�� */
		vector<CVertice2D> edge_gauss_vertice;		/**< ���ϵĸ�˹���ֵ� */
		
		vector<double> basis_fun_coeff;				/**< ������ϵ�� */
	
	protected:
		int _cell_num;								/**< �ܵ�Ԫ��Ŀ */


	// ����������Ҫ�����Ƶ�GPU��	
	public:
		double *area;								/**< ��� */
		double *perimeter;							/**< �ܳ� */
		double *outer_normal_vector;				/**< ��Ԫ�ߵ��ⷨ���� */
		
		double *mass_coeff;							/**< ������������ϵ�� */
		
		double *vol_bf_value;						/**< �����������˹���ֽڵ��ϵ�ֵ */
		
		double *vol_bdf_value;						/**< �������ĵ��������˹���ֽڵ��ֵ */
		
		double *edge_bf_value; 						/**< �������ڱ��ϸ�˹���ֽڵ��ֵ */
	
		double *vol_gauss_weight;					/**< ��˹����ֵ�Ȩ�� */

		double *edge_gauss_weight;					/**< ��˹�߻��ֵ�Ȩ�� */
		
	// ��Ϊ
	protected:
		/** 
		 * ��ʼ��ϵ�����󣬱��������Գ̽�
		 * @param[in] index integer �����α��
		 * @param[in] x double[] �����ζ��������
		 * @param[in] y double[] �����ζ��������
		 * @param[in|out] s double[] δ֪
		 */
		void initXishuMatrix(int index, double x[], double y[], double s[]);


		/** 
		 * ��������ϵ���� �ú������Գ̽�
		 * @param[in] index integer �����α��
		 * @param[in|out] s double[] δ֪
		 */
		void solveCoefficient(int index, double s[]);

		/** 
		 * ��ʼ��ϵ������
		 * @param[in] index integer �����α��
		 */
		void initMassMatrix(int index);
	
	public:
		/** ���캯�� */
		CTriangleInfo();

		// getter
		/** 
		 * ��ȡָ�������ε�Ԫ������
		 * @param[in] index integer �����α��
		 * @return CVertice2D ����������
 		 */
		CVertice2D& getBarycenter(int index);

		/**
		 * ���ر��е�
		 * @param[in] tindex integer �����α��
		 * @param[in] eindex integer �ڼ�����
		 * @return CVertice2D �ñ��е�
		 */
		CVertice2D& getEdgeMiddleVertice(int tindex, int eindex);

		/**
		 * @param[in] index integer �����α��
		 * @return double ����������Բ�뾶
		 */
		double getRadiusOfInscribedCircle(int index);
		
		/**
		 * ������
		 * @param[in] tindex integer �����α��
		 * @param[in] findex integer �ڼ���������
		 * @param[in] x	double ������
		 * @param[in] y double ������
		 * @return double �������ڸõ��ֵ
		 */
		double basisFunction(int tindex, int findex, double x, double y);
		
		/**
		 * �������ĵ�����ָ�����ֵ
		 * @param[in] tindex integer �����α��
		 * @param[in] findex integer �ڼ����������������ĵ�����
		 * @param[in] x	double ������
		 * @param[in] y double ������
		 * @return double �������ڸõ��ֵ
		 */
		double derivativeFunction(int tindex, int findex, double x, double y);
		

		/**
		 * ��������������Բ�뾶
		 * @param[in] tindex integer ���Ǳ��
		 * @param[in] radius double ����������Բ�뾶
		 */
		void setRadiusOfInscribedCircle(int tindex, double radius);


		/** 
		 * ���ݶ��������ʼ��������Ϣ
		 * @param[in] index integer �����α��
		 * @param[in] x double[] ��������ĺ�����
		 * @param[in] y double[] ���������������
		 */
		void initInformation(int index, double x[3], double y[3]);

		
		/**
		 * ��������ڴ�
		 * @param[in] num integer �����ܵ�Ԫ��Ŀ
		 */
		void allocateMemory(int num);

		/** �������� */
		~CTriangleInfo();
	
};