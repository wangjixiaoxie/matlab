
prefix = /usr/local
libdir = ${prefix}/lib
matlabdir = ${libdir}/matlab/
matlabsitedir = ${matlabdir}/site-packages/
matlabpkgdir = ${matlabsitedir}/uisonganal/

MEX = /usr/bin/mex

SUBDIRS =

matlab_files = \
	b_spectranal.m \
	cepstrum.m \
	comp_scores.m \
	costrans.m \
	deltaprods.m \
	deltaprods2.m \
	disp_spectra.m \
	disp_spectra2.m \
	findmax1.m \
	old_distchck.m \
	old_normpdf.m \
	psd_get_xrange.m \
	psdanal.m \
	psdupdate.m \
	sinprods.m \
	spect_comp.m \
	spect_length.m \
	spectranal.m \
	split_spect.m \
	temp_plot.m \
	win_spect.m

mex_targets =

CLEAN_FILES =

all_targets = \
	${mex_targets}

all: all-subdirs $(all_targets)

all-subdirs: dummy
	@for i in $(SUBDIRS); do \
	  make -C $$i ; \
	done

install: dummy
	install -d ${matlabpkgdir}
	for i in ${matlab_files}; do \
	  install -m644 $$i ${matlabpkgdir}; \
	done
	for i in ${mex_targets}; do \
	  install -m755 $$i ${matlabpkgdir}; \
	done

read_song.mexglx: read_song.c
	$(MEX) read_song.c

clean-subdirs: dummy
	@for i in $(SUBDIRS); do \
	  make -C $$i clean; \
	done

clean: clean-subdirs
	-rm -f *.o
	rm -f $(CLEAN_FILES)

dummy:

