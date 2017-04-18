CXX := g++
OUT := tpfn
$(OUT) : main.o filtre.o bmp.o quantlm.o
	$(CXX) $^ -o $(OUT)
bmp.o : bmp.cpp
	$(CXX) -c $^ -o $@
quantlm.o : quantlm.cpp
	$(CXX) -c $^ -o $@
filtre.o : filtre.cpp
	$(CXX) -c $^ -o $@
main.o : main.cpp
	$(CXX) -c $^ -o $@

.PHONY : clean
clean ::
	 rm -f *.o
