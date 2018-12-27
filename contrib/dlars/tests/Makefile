RUN=srun -n 7 ../src/dlars
COMPARE=perl ../../compare/compare.pl
all: lars lasso stopping split weights

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
	-$(COMPARE) dlasso.64.32.mn.txt correct_output/dlasso.64.32.mn.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --lambda=5.0e-03 > dlasso.64.32.lam.txt
	-$(COMPARE) dlasso.64.32.lam.txt correct_output/dlasso.64.32.lam.txt
	$(RUN) X.64.32.txt Y.64.32.txt X.64.32.dim --iterations=10 > dlasso.64.32.its.txt
	-$(COMPARE) dlasso.64.32.its.txt correct_output/dlasso.64.32.its.txt

split:
	srun -n 2 ../src/dlars A.txt b.txt dim.txt --split_files >& dlasso.A.split.txt
	-$(COMPARE) dlasso.A.split.txt correct_output/dlasso.A.split.txt
	srun -n 4 ../src/dlars A.txt b.txt dim.txt --split_files >& dlasso.A.split.4.txt
	-$(COMPARE) dlasso.A.split.4.txt correct_output/dlasso.A.split.4.txt
	srun -n 7 ../src/dlars A.txt b.txt dim.txt --split_files >& dlasso.A.split.7.txt
	-$(COMPARE) dlasso.A.split.7.txt correct_output/dlasso.A.split.7.txt

weights:
	$(RUN) Xcpp.txt Ycpp.txt Xcpp.dim --weights=Ycpp.weights > dlasso.weights.txt
	-$(COMPARE) dlasso.weights.txt correct_output/dlasso.weights.txt

clean:
	rm -f dlars.*.txt dlasso.*.txt *.diff

