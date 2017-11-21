SUBDIR = rci_symb rci_symb_old rci_print rcsfgenerate_sym  #\
#         rci_mpi rwfnestimate rcsl rnucleus  \
#         rcsfgenerate rcsfmerge \
#         rmixcondens \
#         rcsfdelete \
#         rmcdhf \
#         rbiotransform rtransition \
#         rangular rci rhfs rhfs_lsj ris \
#         rangular_mpi rmcdhf_mpi rci_mpi rbiotransform_mpi rtransition_mpi rhfsd rhfsd_mpi HF rmcdhf_mpi_old rci_mpi_old

# rcsfzerofirst from src90

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

clean :
	@for i in $(SUBDIR) ;\
	do \
	echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

