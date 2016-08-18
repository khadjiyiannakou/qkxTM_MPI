QKXTM = ../lib/libqkxTM_MPI.a

INC= -I../include
MPICXX = mpicxx
CXXFLAGS = -Wall -O3 -qopenmp $(INC) 
LIBS = /homec/ecy00/ecy011/qcd_kyriakos/lime-1.3.2/lib/liblime.a

#EXECUTABLES = test_geometry test_comm test_read_gauge test_ape_smearing test_gauss_smearing test_applyTM test_invert test_invert_prop test_invert_end_trick test_invert_end_trick_1 test_invert_tDil test_HPE test_HPE_phi test_seq_source test_gluonLoop test_read_prop createTwopWtSLProp createTwopWtSSProp write_smeared_gauge write_prop plaquette twop-hadrons

#EXECUTABLES = test_read_prop test_gauss_smearing createTwopWtSLProp createTwopWtSSProp
#EXECUTABLES = twop-hadrons write_smeared_gauge createTwopWtSLProp test_applyTM_PC test_applyTM test_applyTM_chiral write_unit_gauge  test_invert_end_trick test_invert_volumeSource

EXECUTABLES = test_invert_end_trick test_field_strength

all: $(EXECUTABLES)

write_unit_gauge: write_unit_gauge.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_applyTM_PC: test_applyTM_PC.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_applyTM_chiral: test_applyTM_chiral.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

plaquette: plaquette.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

write_prop: write_prop.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

createTwopWtSLProp: createTwopWtSLProp.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

write_smeared_gauge: write_smeared_gauge.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

createTwopWtSSProp: createTwopWtSSProp.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_read_prop: test_read_prop.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_geometry: test_geometry.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_comm: test_comm.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_read_gauge: test_read_gauge.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_ape_smearing: test_ape_smearing.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_gauss_smearing: test_gauss_smearing.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_applyTM: test_applyTM.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert: test_invert.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert_prop: test_invert_prop.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert_end_trick: test_invert_end_trick.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert_volumeSource: test_invert_volumeSource.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert_end_trick_1: test_invert_end_trick_1.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_invert_tDil: test_invert_tDil.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_HPE: test_HPE.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_HPE_phi: test_HPE_phi.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_seq_source: test_seq_source.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_gluonLoop: test_gluonLoop.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

twop-hadrons: twop-hadrons.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

twop-hadrons_cache: twop-hadrons_cache.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

test_field_strength: test_field_strength.cpp $(QKXTM)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	-rm -f $(EXECUTABLES) 

%.o: %.cpp 
	$(MPICXX) $(CXXFLAGS) $< -c -o $@
