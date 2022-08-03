
SRCDIRS  = Lib FFTW Models occlib \
           Phantom GCM GRIB ICO ECHAM NCEP MPAS HMC Wavelib \
           Diff occnoise occstat statmerge occutil occ2nc \
           occfps occ Invert Hybrid Wave atmPrf2asc occvar

#SRCDIRS  = occlib MPAS Wavelib Wave
#SRCDIRS  = occlib MPAS NCEP Wavelib Wave
# Only rebuild MPAS libraries and then link with Wave
# SRCDIRS  = MPAS Wave

gen:
	@echo `pwd`
	@chmod u+x hdr/mak.gen
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR; pwd ;\
	  chmod u+x $$DIR.gen ; ./$$DIR.gen ;\
	  if [ -d $$DIR\_test ] ; then \
	             cd $$DIR\_test ;  ./$$DIR\_test.gen ; fi ; ) ;\
	  done       

all:
	@echo `pwd`
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) -f $$DIR.mak ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done       

test:
	@echo `pwd`
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (cd $$DIR ;\
          if [ -d $$DIR\_test ] ; then \
             cd $$DIR\_test ;  $(MAKE) -f $$DIR\_test.mak ; fi ; ) ;\
	  done       

clean:
	@echo `pwd`
	@for DIR in $(SRCDIRS) $(TESTDIRS) ;\
	  do \
	  (cd $$DIR ;\
	  rm -f *.o *.a *.mod *.x; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;\
          if [ -d $$DIR\_test ] ; then \
             cd $$DIR\_test ;  rm -f *.o *.a *.mod *.x; fi ; ) ;\
	  done

release:
	@echo `pwd`
	@chmod u+x hdr/mak.gen
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (MODE="release"; export MODE; cd $$DIR ;\
	  chmod u+x $$DIR.gen ; ./$$DIR.gen ;\
	  if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;\
	  if [ -d $$DIR\_test ] ; then \
	     cd $$DIR\_test ; chmod u+x $$DIR\_test.gen ;\
	     ./$$DIR\_test.gen; fi ; ) ;\
	  done       

debug:
	@echo `pwd`
	@chmod u+x hdr/mak.gen
	@for DIR in $(SRCDIRS) ;\
	  do \
	  (MODE="debug"; export MODE; cd $$DIR ;\
	  chmod u+x $$DIR.gen ; ./$$DIR.gen ;\
	  if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;\
	  if [ -d $$DIR\_test ] ; then \
	     cd $$DIR\_test ; chmod u+x $$DIR\_test.gen ;\
	     ./$$DIR\_test.gen; fi ; ) ;\
	  done       

