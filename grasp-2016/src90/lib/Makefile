

SUBDIR = libcommons lib9290 liblapack libmcp90 librang90 librang90_inter
TARGETS	 = install

$(TARGETS):
	@for i in $(SUBDIR) ;\
	do \
		echo "Current directory: $$i" ; \
	                cd $$i ; \
			make; \
echo "....Leaving: $$i" ; echo; echo;\
		cd .. ; \
	done
clean:
	@for i in $(SUBDIR) ;\
        do \
                echo "Current directory: $$i" ; \
                        cd $$i ; \
                        $(MAKE) $@ ; \
                cd .. ; \
        done


