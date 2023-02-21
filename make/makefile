# This is the makefile for the mesa binary library

MESA_DIR = ../..

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

MOD_DEFAULTS_DIR = $(MODULE_DIR)/defaults
MOD_OTHER_DIR = $(MODULE_DIR)/other

#################################################################
#
# SOURCE FILES

SRCS = \
   binary_def.f90 \
   binary_private_def.f90 \
   binary_jdot.f90 \
   binary_wind.f90 \
   binary_mdot.f90 \
   binary_tides.f90 \
   binary_edot.f90 \
   binary_rlobe.f90 \
   binary_ctrls_io.f90 \
   binary_irradiation.f90 \
   binary_history_specs.f90 \
   binary_history.f90 \
   binary_evolve.f90 \
   binary_timestep.f90\
	binary_job_ctrls_io.f90\
	binary_do_one_utils.f90\
   run_binary_support.f90 \
   binary_lib.f90 \
   mod_other_rlo_mdot.f90 \
   mod_other_tsync.f90 \
   mod_other_sync_spin_to_orbit.f90 \
   mod_other_mdot_edd.f90 \
   mod_other_accreted_material_j.f90 \
   mod_other_binary_jdot.f90 \
   mod_other_binary_wind_transfer.f90 \
   mod_other_binary_edot.f90 \
   mod_other_binary_extras.f90

#################################################################
#
# TARGETS

LIB = libbinary.a
LIB_OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

#################################################################
#
# COMPULATION RULES

BINARY_INCLUDES = -I$(MOD_DEFAULTS_DIR)

COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree) $(BINARY_INCLUDES)
#COMPILE = $(COMPILE_TO_TEST) $(FCfree)

COMPILE_CMD = $(COMPILE)

%.o : %.f
	$(COMPILE_CMD) $<
%.o : %.f90
	$(COMPILE_CMD) $<
%.mod : %.o
	@true
%.o : %.mod

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.a .depend

#################################################################
#
# DEPENDENCIES

# Paths to search for source & include files

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR):$(MOD_OTHER_DIR)

# Automatic dependency generation using makedepf90

nodeps : $(.DEFAULT_GOAL)

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  DEPEND_DIRS = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR):$(MOD_OTHER_DIR)

  .depend :
	makedepf90 -m %m.mod -I$(DEPEND_DIRS) $(SRCS) > .depend
  -include .depend

endif
