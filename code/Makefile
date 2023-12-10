##
## Makefile for running_coupling.x
##

CXXFLAGS += -Wall `root-config --cflags`
INCLUDES += -I
LIBS     += `root-config --glibs`

SRCS1 := analytic_coupling.cpp BetaFn.cpp
OBJS1 := analytic_coupling.o BetaFn.o
EXEC1 := analytic_coupling.x

SRCS2 := numerical_coupling.cpp BetaFn.cpp
OBJS2 := numerical_coupling.o BetaFn.o
EXEC2 := numerical_coupling.x

SRCS3 := fit_nf0.cpp BetaFn.cpp
OBJS3 := fit_nf0.o BetaFn.o
EXEC3 := fit_nf0.x

SRCS4 := fit_nf2.cpp BetaFn.cpp
OBJS4 := fit_nf2.o BetaFn.o
EXEC4 := fit_nf2.x

SRCS5 := improved_convergence.cpp BetaFn.cpp
OBJS5 := improved_convergence.o BetaFn.o
EXEC5 := improved_convergence.x

SRCS6 := find_lambda3.cpp BetaFn.cpp
OBJS6 := find_lambda3.o BetaFn.o
EXEC6 := find_lambda3.x

SRCS7 := global_fit_for_r0a.cpp
OBJS7 := global_fit_for_r0a.o
EXEC7 := global_fit_for_r0a.x

SRCS8 := nonglobal_fit_for_plaq.cpp 
OBJS8 := nonglobal_fit_for_plaq.o 
EXEC8 := nonglobal_fit_for_plaq.x

SRCS9 := replicate_fit_nf2.cpp BetaFn.cpp
OBJS9 := replicate_fit_nf2.o BetaFn.o
EXEC9 := replicate_fit_nf2.x

SRCS10 := running_coupling.cpp BetaFn.cpp
OBJS10 := running_coupling.o BetaFn.o
EXEC10 := running_coupling.x

SRCS11 := compare_results.cpp
OBJS11 := compare_results.o
EXEC11 := compare_results.x
    
all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11)

$(EXEC1): $(OBJS1) $(SRCS1)
	$(CXX) -o analytic_coupling.x $(LIBS) ${OBJS1}

analytic_coupling.o : $(SRCS1)
	$(CXX) -c $(CXXFLAGS) analytic_coupling.cpp

$(EXEC2): $(OBJS2) $(SRCS2)
	$(CXX) -o numerical_coupling.x $(LIBS) ${OBJS2}

numerical_coupling.o : $(SRCS2)
	$(CXX) -c $(CXXFLAGS) numerical_coupling.cpp

$(EXEC3): $(OBJS3) $(SRCS3)
	$(CXX) -o fit_nf0.x $(LIBS) ${OBJS3}

fit_nf0.o : $(SRCS3)
	$(CXX) -c $(CXXFLAGS) fit_nf0.cpp

$(EXEC4): $(OBJS4) $(SRCS4)
	$(CXX) -o fit_nf2.x $(LIBS) ${OBJS4}

fit_nf2.o : $(SRCS4)
	$(CXX) -c $(CXXFLAGS) fit_nf2.cpp

$(EXEC5): $(OBJS5) $(SRCS5)
	$(CXX) -o improved_convergence.x $(LIBS) ${OBJS5}

improved_convergence.o : $(SRCS5)
	$(CXX) -c $(CXXFLAGS) improved_convergence.cpp

$(EXEC6): $(OBJS6) $(SRCS6)
	$(CXX) -o find_lambda3.x $(LIBS) ${OBJS6}

find_lambda3.o : $(SRCS6)
	$(CXX) -c $(CXXFLAGS) find_lambda3.cpp

$(EXEC7): $(OBJS7) $(SRCS7)
	$(CXX) -o global_fit_for_r0a.x $(LIBS) ${OBJS7}

global_fit_for_r0a.o : $(SRCS7)
	$(CXX) -c $(CXXFLAGS) global_fit_for_r0a.cpp

$(EXEC8): $(OBJS8) $(SRCS8)
	$(CXX) -o nonglobal_fit_for_plaq.x $(LIBS) ${OBJS8}

nonglobal_fit_for_plaq.o : $(SRCS8)
	$(CXX) -c $(CXXFLAGS) nonglobal_fit_for_plaq.cpp

$(EXEC9): $(OBJS9) $(SRCS9)
	$(CXX) -o replicate_fit_nf2.x $(LIBS) ${OBJS9}

replicate_fit_nf2.o : $(SRCS9)
	$(CXX) -c $(CXXFLAGS) replicate_fit_nf2.cpp

$(EXEC10): $(OBJS10) $(SRCS10)
	$(CXX) -o running_coupling.x $(LIBS) ${OBJS10}

running_coupling.o : $(SRCS10)
	$(CXX) -c $(CXXFLAGS) running_coupling.cpp

$(EXEC11): $(OBJS11) $(SRCS11)
	$(CXX) -o compare_results.x $(LIBS) ${OBJS11}

compare_results.o : $(SRCS11)
	$(CXX) -c $(CXXFLAGS) compare_results.cpp

BetaFn.o : $(SRCS1) BetaFn.h
	$(CXX) -c $(CXXFLAGS) BetaFn.cpp

clean:
	rm -f ${OBJS1} $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7) $(OBJS8) $(OBJS9) $(OBJS10) $(OBJS11)

distclean: clean
