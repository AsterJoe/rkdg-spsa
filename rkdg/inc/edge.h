/**
 * �õ�Ԫ�����
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
 * ���嵥Ԫ�ıߣ����巽ʽΪָ���ߵ������յ���
 * @pre Ϊ���ܹ�����ʹ�ø��࣬���ȳ�ʼ����������
 */
class CEdge{
	// ����
	public:
	
	protected:
		int _start; 		/**< ����� */
		int _terminal; 		/**< �յ��� */
	
	// ��Ϊ
	public:
		/**
		 * Ĭ�Ϲ��캯��
		 * @param[in] start integer ���
		 * @param[in] terminal integer �յ�
		 */
		CEdge(int start=0, int terminal=1);
		
		// setter
		/**
		 * �������ñ�
		 * @param[in] start integer ���
		 * @param[in] terminal integer �յ�
		 */
		void setEdge(int start, int terminal);
		
		/**
		 * ���ñߵ����
		 * @param[in] start integer ���
		 */
		void setStart(int start);
		
		/**
		 * ���ñߵ��յ�
		 * @param[in] terminal integer �յ�
		 */
		void setTerminal(int terminal);
		
		// getter
		
		/**
		 * ��ȡ���
		 */
		int getStart(void) const;
		
		/**
		 * ��ȡ�յ�
		 */
		int getTerminal(void) const;
		
		/** �������� */
		~CEdge();
	private:
};