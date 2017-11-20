/**
 * ���ļ������ά������
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once
#include "cppstdheaders.h"

using namespace std;

/**
 * ��ά�����࣬������Ϊ�����x, y����
 * ��ͨ���ӿڷ���������
 */
class CVertice2D
{
	// ����
	public:
	
	protected:
		double _x; /**< ����ĺ����� */
		double _y; /**< ����������� */
	
	// ����
	public:
		/**
		 * Ĭ�Ϲ��캯��
		 * @param[in] x double ����ĺ�����
		 * @param[in] y double �����������
		 */
		CVertice2D(double x=0, double y=0);

		// setter
		/**
		 * ���ö�������
		 * @param[in] x double ������º�����
		 * @param[in] y double �������������
		 */
		void setVertice(double x, double y);
		
		/**
		 * ���ö���ĺ�����
		 * @param[in] x double ���������
		 */
		void setX(double x);
		
		/**
		 * ���ö����������
		 * @param[in] y double ����������
		 */
		void setY(double y);
		
		// getter
		
		/**
		 * ��ȡ����ĺ�����
		 * @return double ����ĺ�����
		 */
		double getX(void) const;
		
		/**
		 * ��ȡ�����������
		 * @return double �����������
		 */
		double getY(void) const;
		
		
		/** �������� */
		~CVertice2D();
	private:
	
};