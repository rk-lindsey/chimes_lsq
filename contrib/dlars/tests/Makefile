COMPARE=perl ../../compare/compare.pl
all: lars lasso

lars:
	../src/dlars Xcpp.txt Ycpp.txt Xcpp.dim lars > dlars.cpp.txt
	-$(COMPARE) dlars.cpp.txt correct_output/dlars.cpp.txt
	../src/dlars X.64.32.txt Y.64.32.txt X.64.32.dim lars > dlars.64.32.txt
	-$(COMPARE) dlars.64.32.txt correct_output/dlars.64.32.txt
	../src/dlars A.txt b.txt A.dim lars > dlars.A.txt
	-$(COMPARE) dlars.A.txt correct_output/dlars.A.txt


lasso:
	../src/dlars Xcpp.txt Ycpp.txt Xcpp.dim lasso > dlasso.cpp.txt
	-$(COMPARE) dlasso.cpp.txt correct_output/dlasso.cpp.txt
	../src/dlars X.64.32.txt Y.64.32.txt X.64.32.dim lasso > dlasso.64.32.txt
	-$(COMPARE) dlasso.64.32.txt correct_output/dlasso.64.32.txt
	../src/dlars A.txt b.txt A.dim lasso > dlasso.A.txt
	-$(COMPARE) dlasso.A.txt correct_output/dlasso.A.txt

clean:
	rm -f dlars.*.txt dlasso.*.txt

