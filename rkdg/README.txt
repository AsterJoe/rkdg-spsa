在GPU上使用RKDG方法计算二维非结构网格上NACA0012机翼扰流，特点如下：
 1). 在每个三角单元上使用二次正交多项式作为基函数
 2). 三阶龙格库塔法进行时间推进
作者： tlanyan<tag.yuan@gmail.com>
主要参考文献： Jun Zhu, Jianxian Qiu, Chi-Wang Shu, Michael Dumbser, Runge-Kutta discontinuous Galerkin method using WENO limiters II： unstructured meshes.
感谢： CJBUAA 提供参考源代码

[项目说明]
	创建日期： 2013-5-3
	代码编写完成日期： 2013-5-5
	代码测试完成日期： 2013-5-11
	代码整理日期： 2013-5-13
	

[程序结构]
	inc : 该文件夹存放程序的头文件
	input ： 该文件夹存放程序的输入文件（控制文件）
	src ： 该文件夹包含程序的源代码文件
	output ： 该文件夹用来存放程序的输出文件
	test:  该文件夹存放测试函数
	
