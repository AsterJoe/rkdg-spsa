# rkdg-cuda上程序的makefile
# @author tlanyan<tag.yuan@gmail.com>
# @link www.tlanyan.me

# 目标名称
EXECUTABLE = rkdg-cuda

# 定义编译器和编译选项
LD = icpc
CXX = icpc
CXFLAG = -O2  -std=c++0x
NC = nvcc
NCFLAG = -O2 -m64 -arch=compute_30 -code=sm_30 -maxrregcount=0


CUDA_HOME =/soft/cuda/6.0.37

CUDA_INC = -I$(CUDA_HOME)/include
CUDA_LIB = $(CUDA_HOME)/lib64/libcudart.so.6.0.37

# 定义源文件和目标文件
CXSRCS = $(wildcard src/*.cpp)
NCSRCS = $(wildcard src/*.cu)

CXOBJS = $(CXSRCS:.cpp=.o)
CUOBJS = $(NCSRCS:.cu=.o)

OBJS = $(CXOBJS) $(CUOBJS)
	
	
# 定义目标规则
$(EXECUTABLE): $(OBJS)
	$(LD) $(CUDA_LIB) -o $(EXECUTABLE) $(OBJS) 
	
%.o: %.cpp
	$(CXX) $(CXFLAG) -c $< -o $@ $(CUDA_INC)
	
%.o: %.cu
	$(NC)  $(NCFLAG) $(CUDA_INC) -c $< -o $@
	
.PHONY:clean
clean:
	rm -f $(EXECUTABLE) src/*.o` 
	
	
