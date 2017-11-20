/**
 * ��������࣬��������������Ϣ�Ƿ���ȷ
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#pragma once

#include "../inc/cppstdheaders.h"
#include "../inc/unstructuredgrid.h"

using namespace std;

/** 
 * ��������࣬��������������Ϣ�Ƿ���ȷ
 */
class CUnstructuredGridTester: public CUnstructuredGrid
{
	// ����
	public:

	protected:

	// ����
	public:
		/** ���캯�� */
		CUnstructuredGridTester();

		/**
		 * �����������ھ�
		 * @param[in] output string �ھ���Ϣ������ļ�
		 */
		void testNeighbour(const string& output="output/neighbour.dat");

		/**
		 * ���������ι����
		 * @param[in] output string �������Ϣ������ļ�
		 */
		void testSharedEdge(const string& output="output/sharededge.dat");

		/**
		 * ���������α��
		 * @param[in] output string �����Ϣ������ļ�
		 */
		void testTriangleFlag(const string& output="output/triangleflag.dat");

		/**
		 * �������������
		 * @param[in] output string �������������ļ�
		 */
		void testArea(const string& output="output/area.dat");

		/**
		 * �����������ܳ�
		 * @param[in] output string �������ܳ�����ļ�
		 */
		void testPerimeter(const string& output="output/perimeter.dat");


		/**
		 * ���������α��ⷨ����
		 * @param[in] output string �������ⷨ��������ļ�
		 */
		void testOuterNormalVector(const string& output="output/outernormalvector.dat");

		/**
		 * ���������λ����������˹���ֽڵ��ϵ�ֵ
		 */
		void testVolumeBfValue(void);


		/**
		 * ���������λ�����������ϵ��
		 * @param[in] output string �����λ���������ϵ������ļ�
		 */
		void testMassCoeff(const string& output="output/masscoeff.dat");

		/**
		 * ���������λ��������������˹���ֽڵ��ϵ�ֵ
		 */
		void testVolumeBdfValue(void);

		/**
		 * ���������λ������ڱ߸�˹���ֽڵ��ϵ�ֵ
		 */
		void testEdgeBfValue(void);

		/**
		 * ��������������ֵ�Ȩ��
		 * @param[in] output string ����������ֵĸ�˹Ȩ������ļ�
		 */
		void testVolumeGaussWeight(const string& output="output/volumegaussweight.dat");

		/**
		 * �������Ǳ�����ֵ�Ȩ��
		 * @param[in] output string �����α߻��ֵĸ�˹Ȩ������ļ�
		 */
		void testEdgeGaussWeight(const string& output="output/edgegaussweight.dat");


		/**
		 * �����ά�������ļ��鿴
		 * @param[in] rnum �����и��� 
		 * @param[in] cnum �����и���
		 * @param[in] var double[] ��Ҫ����ı���
		 * @param[in] output string ����Ľ���ļ�
		 */
		template<typename T>
		void output2DVar(int rnum, int cnum, const T *var, const string &output);

		/** �������� */
		~CUnstructuredGridTester();

	protected:

};