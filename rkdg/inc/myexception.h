/**
 * �Զ�����쳣��
 * ���ಶ׽���˱�׼�쳣������쳣
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 

#pragma once

#include "cppstdheaders.h"
using namespace std;

/** 
 * �Զ�����쳣��
 */
class CMyException: public exception
{
	// ����
	public:

	protected:
		string _msg;		/**< ������Ϣ���� */

	// ����
	public:
		// ���캯��
		CMyException(const string &error="Unknow exception.");
		
		/**< ��д��������what���� */
		const char *what()const throw();

		/** �������� */
		~CMyException() throw();
	protected:
};