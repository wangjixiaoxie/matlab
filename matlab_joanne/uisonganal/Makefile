
prefix = /usr/local
libdir = ${prefix}/lib
matlabdir = ${libdir}/matlab/
matlabsitedir = ${matlabdir}/site-packages/
matlabpkgdir = ${matlabsitedir}/uisonganal/

MEX = /usr/bin/mex

SUBDIRS = matprint spect

matlab_data = \
	uifig.fig \
	uisalayout3.fig \
	uisasetparams.fig \
	uisasetup.fig

matlab_files = \
	bandpass.m \
	batch_setup.m \
	centercb.m \
	disp_amp.m \
	disp_idx_spect.m \
	disp_song.m \
	edit_note.m \
	filecb.m \
	get_print_vals.m \
	get_scale_vals.m \
	get_segvals.m \
	get_thresh_val.m \
	get_xrange.m \
	get_yrange.m \
	handseg2binary.m \
	imagesc1.m \
	isint.m \
	label.m \
	load_n_plot.m \
	load_n_plot_spec_batch.m \
	load_n_plot_spec_only.m \
	load_uisonganal_results.m \
	loadraw.m \
	make_current.m \
	make_map.m \
	makecol.m \
	makerow.m \
	move_left.m \
	move_right.m \
	myginput.m \
	mysign.m \
	paramcb.m \
	print_song.m \
	read_filt.m \
	read_labv_file.m \
	read_note.m \
	read_obsdata.m \
	read_rawfile42c.m \
	read_spect.m \
	save_data.m \
	scale_spect.m \
	segment.m \
	set_print_vals.m \
	set_scale_vals.m \
	songanalpath.m \
	soundin.m \
	tabdlg.m \
	uibatch.m \
	uifig.m \
	uigetdir_java.m \
	uimove.m \
	uiprint_sa.m \
	uipsdanal.m \
	uiresegment.m \
	uisalayout3.m \
	uisasetparams.m \
	uisasetup.m \
	uiscale.m \
	uisonganal.m \
	uisongcontrols.m \
	uispectcontrols.m \
	uizoom.m \
	uizoomcontrols.m \
	wavread16.m \
	win_spect.m \
	write_filt.m \
	write_spect.m \
	zero_crossings.m \
	zoom_x.m \
	zoom_y.m

mex_targets = \
	read_song.mexglx

CLEAN_FILES = \
	read_song.mexglx

all_targets = \
	${mex_targets}

all: all-subdirs $(all_targets)

all-subdirs: dummy
	@for i in $(SUBDIRS); do \
	  make -C $$i ; \
	done

install: dummy
	install -d ${matlabpkgdir}
	for i in ${matlab_files} ${matlab_data}; do \
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

