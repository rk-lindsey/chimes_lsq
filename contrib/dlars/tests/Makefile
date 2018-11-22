COMPARE=perl ../../compare/compare.pl
all:
	../src/dlars Xcpp.txt Ycpp.txt > dlars.cpp.txt
	$(COMPARE) dlars.cpp.txt correct_output/dlars.cpp.txt
	../src/dlars X.64.32.txt Y.64.32.txt > dlars.64.32.txt
	$(COMPARE) dlars.64.32.txt correct_output/dlars.64.32.txt
	../src/dlars A.txt b.txt > dlars.A.txt
	$(COMPARE) dlars.A.txt correct_output/dlars.A.txt

clean:
	rm -f dlars.*.txt
