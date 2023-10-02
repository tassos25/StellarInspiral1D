ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif

# copied the content from work_standard_makefile

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

ifndef STAR
STAR = star
endif

# STAR_WORK_OBJS defined in makefile_header
# added more code located in WORK_SRC_DIR
IONIZATION_OBJS = ionization_def.o ion_tables_load.o mod_tables.o mod_ionization.o ion_tables_eval.o ion_table_plot.o ionization_lib.o
INSPIRAL_OBJS = CE_orbit.o CE_energy.o CE_torque.o CE_adjust_mdot.o CE_before_struct_burn_mix.o CE_after_struct_burn_mix.o CE_timestep.o CE_run_star_extras.o
OBJS = $(IONIZATION_OBJS) $(INSPIRAL_OBJS) $(RUN_STAR_EXTRA_OBJS) $(STAR_WORK_OBJS) run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
COMMON_CODE_DIR = $(WORK_DIR)/../common_code
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header

ifeq ($(SKIP_NAN_TRAPS),YES)
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<
else
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

ifeq ($(MESA_CONTRIB_DIR),)
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $<
else
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) -I $(MESA_CONTRIB_DIR)/hooks $<
endif

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

#Added rules of how to compile the inspiral code
CE_orbit.o: $(COMMON_CODE_DIR)/inspiral/CE_orbit.f90
	$(WORK_COMPILE) $<

CE_energy.o: $(COMMON_CODE_DIR)/inspiral/CE_energy.f90
	$(WORK_COMPILE) $<

CE_torque.o: $(COMMON_CODE_DIR)/inspiral/CE_torque.f90
	$(WORK_COMPILE) $<

CE_adjust_mdot.o: $(COMMON_CODE_DIR)/inspiral/CE_adjust_mdot.f90
	$(WORK_COMPILE) $<

CE_before_struct_burn_mix.o: $(COMMON_CODE_DIR)/inspiral/CE_before_struct_burn_mix.f90
	$(WORK_COMPILE) $<

CE_after_struct_burn_mix.o: $(COMMON_CODE_DIR)/inspiral/CE_after_struct_burn_mix.f90
	$(WORK_COMPILE) $<

CE_timestep.o: $(COMMON_CODE_DIR)/inspiral/CE_timestep.f90
	$(WORK_COMPILE) $<

CE_run_star_extras.o: $(COMMON_CODE_DIR)/inspiral/CE_run_star_extras.f90
	$(WORK_COMPILE) $<

#A rule of how to compile the additional code that we added in WORK_SRC_DIR
%.o: $(WORK_SRC_DIR)/%.f90
	$(WORK_COMPILE) $<

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
