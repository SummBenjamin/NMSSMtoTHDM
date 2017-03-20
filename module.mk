DIR          := models/NMSSMtower
MODNAME      := NMSSMtower
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes

NMSSMtower_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NMSSMtower_MK     := \
		$(DIR)/module.mk

NMSSMtower_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

NMSSMtower_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

NMSSMtower_TWO_SCALE_MK := \
		$(NMSSMtower_TWO_SCALE_SUSY_MK) \
		$(NMSSMtower_TWO_SCALE_SOFT_MK)

NMSSMtower_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NMSSMtower_generated \
		$(DIR)/LesHouches.in.NMSSMtower \
		$(DIR)/LesHouches.in.NMSSMtower_1507.05093_TP3

NMSSMtower_GNUPLOT := \
		$(DIR)/NMSSMtower_plot_rgflow.gnuplot \
		$(DIR)/NMSSMtower_plot_spectrum.gnuplot

NMSSMtower_TARBALL := \
		$(MODNAME).tar.gz

LIBNMSSMtower_SRC :=
EXENMSSMtower_SRC :=
LLNMSSMtower_LIB  :=
LLNMSSMtower_OBJ  :=
LLNMSSMtower_SRC  :=
LLNMSSMtower_MMA  :=

LIBNMSSMtower_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBNMSSMtower_SRC += \
		$(DIR)/NMSSMtower_effective_couplings.cpp \
		$(DIR)/NMSSMtower_mass_eigenstates.cpp \
		$(DIR)/NMSSMtower_info.cpp \
		$(DIR)/NMSSMtower_input_parameters.cpp \
		$(DIR)/NMSSMtower_observables.cpp \
		$(DIR)/NMSSMtower_slha_io.cpp \
		$(DIR)/NMSSMtower_physical.cpp \
		$(DIR)/NMSSMtower_utilities.cpp \
		$(DIR)/NMSSMtower_standard_model_matching.cpp \
		$(DIR)/NMSSMtower_standard_model_two_scale_matching.cpp \
		$(DIR)/NMSSMtower_two_scale_convergence_tester.cpp \
		$(DIR)/NMSSMtower_two_scale_high_scale_constraint.cpp \
		$(DIR)/NMSSMtower_two_scale_initial_guesser.cpp \
		$(DIR)/NMSSMtower_two_scale_low_scale_constraint.cpp \
		$(DIR)/NMSSMtower_two_scale_model.cpp \
		$(DIR)/NMSSMtower_two_scale_model_slha.cpp \
		$(DIR)/NMSSMtower_two_scale_susy_parameters.cpp \
		$(DIR)/NMSSMtower_two_scale_soft_parameters.cpp \
		$(DIR)/NMSSMtower_two_scale_susy_scale_constraint.cpp
EXENMSSMtower_SRC += \
		$(DIR)/run_NMSSMtower.cpp \
		$(DIR)/run_cmd_line_NMSSMtower.cpp \
		$(DIR)/scan_NMSSMtower.cpp
LIBNMSSMtower_HDR += \
		$(DIR)/NMSSMtower_convergence_tester.hpp \
		$(DIR)/NMSSMtower_effective_couplings.hpp \
		$(DIR)/NMSSMtower_high_scale_constraint.hpp \
		$(DIR)/NMSSMtower_mass_eigenstates.hpp \
		$(DIR)/NMSSMtower_info.hpp \
		$(DIR)/NMSSMtower_initial_guesser.hpp \
		$(DIR)/NMSSMtower_input_parameters.hpp \
		$(DIR)/NMSSMtower_low_scale_constraint.hpp \
		$(DIR)/NMSSMtower_model.hpp \
		$(DIR)/NMSSMtower_model_slha.hpp \
		$(DIR)/NMSSMtower_observables.hpp \
		$(DIR)/NMSSMtower_physical.hpp \
		$(DIR)/NMSSMtower_slha_io.hpp \
		$(DIR)/NMSSMtower_spectrum_generator_interface.hpp \
		$(DIR)/NMSSMtower_spectrum_generator.hpp \
		$(DIR)/NMSSMtower_standard_model_matching.hpp \
		$(DIR)/NMSSMtower_standard_model_two_scale_matching.hpp \
		$(DIR)/NMSSMtower_susy_scale_constraint.hpp \
		$(DIR)/NMSSMtower_utilities.hpp \
		$(DIR)/NMSSMtower_two_scale_convergence_tester.hpp \
		$(DIR)/NMSSMtower_two_scale_high_scale_constraint.hpp \
		$(DIR)/NMSSMtower_two_scale_initial_guesser.hpp \
		$(DIR)/NMSSMtower_two_scale_low_scale_constraint.hpp \
		$(DIR)/NMSSMtower_two_scale_model.hpp \
		$(DIR)/NMSSMtower_two_scale_model_slha.hpp \
		$(DIR)/NMSSMtower_two_scale_soft_parameters.hpp \
		$(DIR)/NMSSMtower_two_scale_susy_parameters.hpp \
		$(DIR)/NMSSMtower_two_scale_susy_scale_constraint.hpp
LLNMSSMtower_SRC  += \
		$(DIR)/NMSSMtower_librarylink.cpp

LLNMSSMtower_MMA  += \
		$(DIR)/NMSSMtower_librarylink.m \
		$(DIR)/run_NMSSMtower.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NMSSMtower_TWO_SCALE_SUSY_MK)
-include $(NMSSMtower_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NMSSMtower_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSMtower_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

endif

# remove duplicates in case all algorithms are used
LIBNMSSMtower_SRC := $(sort $(LIBNMSSMtower_SRC))
EXENMSSMtower_SRC := $(sort $(EXENMSSMtower_SRC))

LIBNMSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNMSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNMSSMtower_SRC)))

EXENMSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENMSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENMSSMtower_SRC)))

EXENMSSMtower_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENMSSMtower_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENMSSMtower_SRC)))

LIBNMSSMtower_DEP := \
		$(LIBNMSSMtower_OBJ:.o=.d)

EXENMSSMtower_DEP := \
		$(EXENMSSMtower_OBJ:.o=.d)

LLNMSSMtower_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNMSSMtower_SRC)))

LLNMSSMtower_OBJ  := $(LLNMSSMtower_SRC:.cpp=.o)
LLNMSSMtower_LIB  := $(LLNMSSMtower_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNMSSMtower     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NMSSMtower := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NMSSMtower := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNMSSMtower) $(EXENMSSMtower_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NMSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSMtower_SRC) $(NMSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSMtower_HDR) $(NMSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENMSSMtower_SRC) $(NMSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNMSSMtower_SRC) $(NMSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNMSSMtower_MMA) $(NMSSMtower_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NMSSMtower_MK) $(NMSSMtower_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NMSSMtower_TWO_SCALE_MK) $(NMSSMtower_INSTALL_DIR)
ifneq ($(NMSSMtower_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NMSSMtower_SLHA_INPUT) $(NMSSMtower_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NMSSMtower_GNUPLOT) $(NMSSMtower_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNMSSMtower_DEP)
		-rm -f $(EXENMSSMtower_DEP)
		-rm -f $(LLNMSSMtower_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNMSSMtower)
		-rm -f $(LLNMSSMtower_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNMSSMtower_OBJ)
		-rm -f $(EXENMSSMtower_OBJ)
		-rm -f $(LLNMSSMtower_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNMSSMtower_SRC)
		-rm -f $(LIBNMSSMtower_HDR)
		-rm -f $(EXENMSSMtower_SRC)
		-rm -f $(LLNMSSMtower_SRC)
		-rm -f $(LLNMSSMtower_MMA)
		-rm -f $(METACODE_STAMP_NMSSMtower)
		-rm -f $(NMSSMtower_TWO_SCALE_MK)
		-rm -f $(NMSSMtower_SLHA_INPUT)
		-rm -f $(NMSSMtower_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENMSSMtower_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NMSSMtower_TARBALL) \
		$(LIBNMSSMtower_SRC) $(LIBNMSSMtower_HDR) \
		$(EXENMSSMtower_SRC) \
		$(LLNMSSMtower_SRC) $(LLNMSSMtower_MMA) \
		$(NMSSMtower_MK) $(NMSSMtower_TWO_SCALE_MK) \
		$(NMSSMtower_SLHA_INPUT) $(NMSSMtower_GNUPLOT)

$(LIBNMSSMtower_SRC) $(LIBNMSSMtower_HDR) $(EXENMSSMtower_SRC) $(LLNMSSMtower_SRC) $(LLNMSSMtower_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NMSSMtower)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NMSSMtower): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NMSSMtower)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_NMSSMtower)"
		@echo "Note: to regenerate NMSSMtower source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NMSSMtower)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NMSSMtower):
		@true
endif

$(LIBNMSSMtower_DEP) $(EXENMSSMtower_DEP) $(LLNMSSMtower_DEP) $(LIBNMSSMtower_OBJ) $(EXENMSSMtower_OBJ) $(LLNMSSMtower_OBJ) $(LLNMSSMtower_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNMSSMtower_DEP) $(EXENMSSMtower_DEP) $(LLNMSSMtower_DEP) $(LIBNMSSMtower_OBJ) $(EXENMSSMtower_OBJ) $(LLNMSSMtower_OBJ) $(LLNMSSMtower_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNMSSMtower_OBJ) $(LLNMSSMtower_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBNMSSMtower): $(LIBNMSSMtower_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNMSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNMSSMtower_LIB): $(LLNMSSMtower_OBJ) $(LIBNMSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBNMSSMtower_DEP) $(EXENMSSMtower_DEP)
ALLSRC += $(LIBNMSSMtower_SRC) $(EXENMSSMtower_SRC)
ALLLIB += $(LIBNMSSMtower)
ALLEXE += $(EXENMSSMtower_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNMSSMtower_DEP)
ALLSRC += $(LLNMSSMtower_SRC)
ALLLL  += $(LLNMSSMtower_LIB)
endif
