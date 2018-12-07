RUN=srun -n 7 ../src/dlars
COMPARE=perl ../../compare/compare.pl
all: lars lasso

lars:
	$(RUN) Xcpp.txt Ycpp.txt Xcpp.dim --algorithm=lars > dlars.cpp.txt
	-$(COMPARE) dlars.cpp.txt correct_output/dlars.cpp.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --algorithm=lars > dlars.64.32.txt
	-$(COMPARE) dlars.64.32.txt correct_output/dlars.64.32.txt
	$(RUN) A.txt b.txt A.dim --algorithm=lars > dlars.A.txt
	-$(COMPARE) dlars.A.txt correct_output/dlars.A.txt

lasso:
	$(RUN) Xcpp.txt Ycpp.txt Xcpp.dim > dlasso.cpp.txt
	-$(COMPARE) dlasso.cpp.txt correct_output/dlasso.cpp.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim > dlasso.64.32.txt
	-$(COMPARE) dlasso.64.32.txt correct_output/dlasso.64.32.txt
	$(RUN) A.txt b.txt A.dim > dlasso.A.txt
	-$(COMPARE) dlasso.A.txt correct_output/dlasso.A.txt

stopping:
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --max_norm=20 > dlasso.64.32.mn.txt
	-$(COMPARE) dlasso.64.32.txt correct_output/dlasso.64.32.mn.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --lambda=5.0e-03 > dlasso.64.32.lam.txt
	-$(COMPARE) dlasso.64.32.txt correct_output/dlasso.64.32.lam.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --iterations=10 > dlasso.64.32.its.txt
	-$(COMPARE) dlasso.64.32.txt correct_output/dlasso.64.32.its.txt


clean:
	rm -f dlars.*.txt dlasso.*.txt

