
CXX           = g++ -Wno-write-strings -Wno-pragmas
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINCLUDE  := -I$(shell root-config --incdir)


all:
	make model

model: main.$(ObjSuf) anti_rot.$(ObjSuf)  global.$(ObjSuf)  res_param.$(ObjSuf) corr_centrifug_bar_rho.$(ObjSuf) rho_ls_to_hel.$(ObjSuf) Clebsch.$(ObjSuf) amp_elmag.$(ObjSuf) Wigner_d.$(ObjSuf) amp_rho_pipi.$(ObjSuf) prop_rho_pipi.$(ObjSuf)
	$(CXX) -g -o $@ $^ $(ROOTGLIBS) -lEG -lMathMore 

%.$(ObjSuf): %.$(SrcSuf)
	$(CXX) -g -c $(ROOTINCLUDE) -c $<

clean:
	rm -f *.o
#	rm -f *Dict.*
#	rm -f G__*
	rm model
