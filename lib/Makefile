INC= -I../include -I/homec/ecy00/ecy011/qcd_kyriakos/lime-1.3.2/include

MPICXX = mpicxx
CXXFLAGS = -Wall -O3 -Wwrite-strings -qopenmp $(INC)
QKXTM = libqkxTM_MPI.a

QKXTM_OBJS = lattice_util.o util.o read_conf_lime.o blas.o solvers.o blas_qkxTM.o read_prop_lime.o write_conf_lime.o contractions.o read_eigenVector.o

QKXTM_HDRS = lattice_util.h qkxTM.h 


all: $(QKXTM)

$(QKXTM): $(QKXTM_OBJS)
	ar cru $@ $(QKXTM_OBJS)

clean:
	-rm -f *.o $(QKXTM)

contractions.o: contractions.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

read_eigenVector.o: read_eigenVector.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

write_conf_lime.o: write_conf_lime.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

lattice_util.o: lattice_util.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

util.o: util.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

read_conf_lime.o: read_conf_lime.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

read_prop_lime.o: read_prop_lime.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

blas.o: blas.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

solvers.o: solvers.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

blas_qkxTM.o: blas_qkxTM.cpp
	$(MPICXX) $(INC) $(CXXFLAGS) $< -c -o $@

