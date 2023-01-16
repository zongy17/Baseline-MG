SOFT_HOME = /storage/hpcauser/zongyi/HUAWEI/software

HYPRE_F64_VER = 2.25.0
HYPRE_F64_HOME = $(SOFT_HOME)/hypre/$(HYPRE_F64_VER)
HYPRE_F64_LIB  = $(HYPRE_F64_HOME)/lib
HYPRE_F64_INCLUDE = $(HYPRE_F64_HOME)/include
L_F64_HYPRE = -L$(HYPRE_F64_LIB) -lHYPRE

HYPRE_F32_VER = 2.25.0-f32
HYPRE_F32_HOME = $(SOFT_HOME)/hypre/$(HYPRE_F32_VER)
HYPRE_F32_LIB  = $(HYPRE_F32_HOME)/lib
HYPRE_F32_INCLUDE = $(HYPRE_F32_HOME)/include
L_F32_HYPRE = -L$(HYPRE_F32_LIB) -lHYPRE

HYPRE_F32_I64_VER = 2.25.0-f32-i64
HYPRE_F32_I64_HOME = $(SOFT_HOME)/hypre/$(HYPRE_F32_I64_VER)
HYPRE_F32_I64_LIB  = $(HYPRE_F32_I64_HOME)/lib
HYPRE_F32_I64_INCLUDE = $(HYPRE_F32_I64_HOME)/include
L_F32_I64_HYPRE = -L$(HYPRE_F32_I64_LIB) -lHYPRE

CXX= mpicxx
CFLAGS = -O3 -std=c++11 -fopenmp -g

CL = mpicxx
LFLAGS = -lm -fopenmp -ldl

all:	base_sstruct_f64.exe\
	base_sstruct_f32.exe\
	base_struct_f64.exe\
	base_struct_f32.exe\
	base_IJ_f32.exe\
	base_IJ_f64.exe\
	base_IJ_f32_i64.exe

base_struct_f64.exe : base_struct.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F64_INCLUDE) $^ $(L_F64_HYPRE) $(LFLAGS) -o $@

base_struct_f32.exe : base_struct.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F32_INCLUDE) $^ $(L_F32_HYPRE) $(LFLAGS) -o $@

base_sstruct_f64.exe : base_sstruct.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F64_INCLUDE) $^ $(L_F64_HYPRE) $(LFLAGS) -o $@

base_sstruct_f32.exe : base_sstruct_grapes.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F32_INCLUDE) $^ $(L_F32_HYPRE) $(LFLAGS) -o $@

base_IJ_f64.exe : base_IJ.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F64_INCLUDE) $^ $(L_F64_HYPRE) $(LFLAGS) -o $@

base_IJ_f32.exe : base_IJ.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F32_INCLUDE) $^ $(L_F32_HYPRE) $(LFLAGS) -o $@

base_IJ_f32_i64.exe : base_IJ.c
	$(CXX) $(CFLAGS) -I$(HYPRE_F32_I64_INCLUDE) $^ $(L_F32_I64_HYPRE) $(LFLAGS) -o $@


default: all
.PHONY: all clean

clean:
	rm *.o *.exe

