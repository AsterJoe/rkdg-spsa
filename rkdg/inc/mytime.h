/** 
 * �Զ����ʱ����
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#pragma once

#ifndef MAX_LENGTH
#define MAX_LENGTH 256
#endif

#include "cppstdheaders.h"
using namespace std;

/**
 * �Զ����ʱ���࣬��ʵ��CPU��ʱ��ǽ��ʱ���ʱ�Լ���ʽ������ʱ��
 */
class CMyTime {
	
	// ����
	protected:
		clock_t _cpu_start;			/**< ����ʼ�Ĵ�����ʱ�� */
		time_t _wall_start; 		/**< ����ʼ��ǽ��ʱ�� */
		clock_t _cpu_end;			/**< ��������Ĵ�����ʱ�� */
		time_t _wall_end;			/**< ���������ǽ��ʱ�� */
	
	public:
	
	// ����
	protected:
	
	public:
		// ���캯��
		CMyTime();
		
		
		// ��ʱ��ʼ
		void beginTimer( void );
		
		// ��ʱ����
		void endTimer( void );
		
		/**
		 * ��ȡCPUʱ��
		 * @return ������ʽ���ص�CPUʱ��
		 */
		double getCPUElapsedTime( void );
		
		/**
		 * ��ȡǽ��ʱ��
		 * @return ������ʽ���ص�ǽ��ʱ��
		 */
		double getWallElapsedTime( void );
		
		
		/**
		 * ��ȡ��ǰʱ�̣����ַ�����ʾ
		 */
		string getCurrentTime( void );
		
		// ��������
		~CMyTime();
};