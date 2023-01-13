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

# SUPERLU_HOME = $(SOFT_HOME)/superlu/5.3.0
# LSUPERLU = -L$(SUPERLU_HOME)/lib64 -lsuperlu

# LBLAS = -L${BLAS_LIB} -lopenblas

# METIS_HOME = $(SOFT_HOME)/metis/5.1.0
# PAR_METIS_HOME = $(SOFT_HOME)/parmetis/4.0.3
# LMETIS = -L$(PAR_METIS_HOME)/lib -lparmetis -L$(METIS_HOME)/lib -lmetis

CXX= mpicxx
CFLAGS = -O3 -std=c++11 -fopenmp -g

CL = mpicxx
LFLAGS = -lm -fopenmp -ldl

all: base_struct.exe base_sstruct.exe base_hypre_f64.exe base_hypre_f32.exe base_hypre_f32_i64.exe

ALL_F64_LINK = $(L_F64_HYPRE) #$(LSUPERLU) $(LBLAS) $(LMETIS)
ALL_F32_LINK = $(L_F32_HYPRE)
ALL_F32_I64_LINK = $(L_F32_I64_HYPRE)

base_struct.exe : base_struct.o
	$(CL) $^ $(ALL_F64_LINK) $(LFLAGS) -o $@
base_struct.o : base_struct.c
	$(CXX) -c $^ $(CFLAGS) -I$(HYPRE_F64_INCLUDE) -o $@

base_sstruct.exe : base_sstruct.o
	$(CL) $^ $(ALL_F64_LINK) $(LFLAGS) -o $@
base_sstruct.o : base_sstruct.c
	$(CXX) -c $^ $(CFLAGS) -I$(HYPRE_F64_INCLUDE) -o $@

base_hypre_f64.exe : base_hypre_f64.o
	$(CL) $^ $(ALL_F64_LINK) $(LFLAGS) -o $@
base_hypre_f64.o : base_hypre_f64.c
	$(CXX) -c $^ $(CFLAGS) -I$(HYPRE_F64_INCLUDE) -o $@

base_hypre_f32.exe : base_hypre_f32.o
	$(CL) $^ $(ALL_F32_LINK) $(LFLAGS) -o $@
base_hypre_f32.o : base_hypre_f32.c
	$(CXX) -c $^ $(CFLAGS) -I$(HYPRE_F32_INCLUDE) -o $@

base_hypre_f32_i64.exe : base_hypre_f32_i64.o
	$(CL) $^ $(ALL_F32_I64_LINK) $(LFLAGS) -o $@
base_hypre_f32_i64.o : base_hypre_f32.c
	$(CXX) -c $^ $(CFLAGS) -I$(HYPRE_F32_I64_INCLUDE) -o $@

default: all
.PHONY: all clean

clean:
	rm *.o *.exe

