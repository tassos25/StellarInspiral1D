ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif

# copied the content from work_standard_makefile_binary for mesa/binary

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_BINARY) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = binary
endif

# BINARY_WORK_OBJS defined in makefile_header
# added more code located in WORK_SRC_DIR
IONIZATION_OBJS = ionization_def.o ion_tables_load.o mod_tables.o mod_ionization.o ion_tables_eval.o ion_table_plot.o ionization_lib.o
POSYDON_OBJS = POSYDON_run_binary_extras.o POSYDON_run_star_extras.o
OBJS = $(IONIZATION_OBJS) $(POSYDON_OBJS) $(BINARY_WORK_OBJS) binary_run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
COMMON_CODE_DIR = $(WORK_DIR)/../common_code
STAR_JOB_DIR = $(MESA_DIR)/star/job
BINARY_JOB_DIR = $(MESA_DIR)/binary/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header

# STAR_WORK_OBJS

ifeq ($(MESA_CONTRIB_DIR),)
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $<
else
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) -I $(MESA_CONTRIB_DIR)/hooks $<
endif

run_star.o: $(STAR_JOB_DIR)/run_star.f90
	$(WORK_COMPILE) $<

# NOTE: we don't use work/src/run.f


# BINARY_WORK_OBJS

binary_run.o: $(WORK_SRC_DIR)/binary_run.f90
	$(WORK_COMPILE) $<

ifeq ($(MESA_CONTRIB_DIR),)
run_binary_extras.o: $(WORK_SRC_DIR)/run_binary_extras.f90
	$(WORK_COMPILE) $<
else
run_binary_extras.o: $(WORK_SRC_DIR)/run_binary_extras.f90
	$(WORK_COMPILE) -I $(MESA_CONTRIB_DIR)/hooks $<
endif

run_binary.o: $(BINARY_JOB_DIR)/run_binary.f90
	$(WORK_COMPILE) $<

#Added rules of how to compile the ionization code
ionization_def.o: $(COMMON_CODE_DIR)/ionization/ionization_def.f90
	$(WORK_COMPILE) $<

ionization_lib.o: $(COMMON_CODE_DIR)/ionization/ionization_lib.f90
	$(WORK_COMPILE) $<

ion_table_plot.o: $(COMMON_CODE_DIR)/ionization/ion_table_plot.f90
	$(WORK_COMPILE) $<

ion_tables_eval.o: $(COMMON_CODE_DIR)/ionization/ion_tables_eval.f90
	$(WORK_COMPILE) $<

ion_tables_load.o: $(COMMON_CODE_DIR)/ionization/ion_tables_load.f90
	$(WORK_COMPILE) $<

mod_ionization.o: $(COMMON_CODE_DIR)/ionization/mod_ionization.f90
	$(WORK_COMPILE) $<

mod_tables.o: $(COMMON_CODE_DIR)/ionization/mod_tables.f90
	$(WORK_COMPILE) $<

#Added rules of how to compile the extra code from POSYDON
POSYDON_run_binary_extras.o: $(COMMON_CODE_DIR)/POSYDON/POSYDON_run_binary_extras.f90
	$(WORK_COMPILE) $<

POSYDON_run_star_extras.o: $(COMMON_CODE_DIR)/POSYDON/POSYDON_run_star_extras.f90
	$(WORK_COMPILE) $<

#A rule of how to compile additional code in WORK_SRC_DIR
%.o: $(WORK_SRC_DIR)/%.f90
	$(WORK_COMPILE) $<


clean:
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
