CXX = g++
MPICC = mpic++
CPP_FLAG = -O3 -std=c++11 -Wall -fopenmp -Wno-unknown-pragmas
ADDITIONAL_FLAG =
#ADDITIONAL_FLAG = -D_DEBUG 
#ADDITIONAL_FLAG = -D_PRINT_GRADIENT
#ADDITIONAL_FLAG = -DLOAD_FILE_PER_WORKER -D_DEBUG
CPP_INCLUDE = -I./src -I./hazytl/include/
CPP_LAST = -lpthread

BIN_DIR = bin
APP_DIR = app
SPLIT_DATA_TOOL = ${BIN_DIR}/tools/split_data.out
RANDOM_RESHUFFLE_TOOL = ${BIN_DIR}/tools/random_reshuffle.out
GENERATE_SPLIT_FILE_TOOL = ${BIN_DIR}/tools/generate_split_file.out
SRC_FILES = src/executor.h
SRC_FILES += app/linearmodel_exec.h app/linearmodel/linearmodel.h app/linearmodel/linearmodel_loader.h

STRATEGY := ASYNC ALLREDUCE ASYNCSPARSE ALLREDUCESPARSEBIG ALLREDUCESPARSESMALL ALLREDUCESPARSERECDBL ALLREDUCESPARSERING
#STRATEGY := HOGWILD SYNC ASYNC ALLREDUCE ASYNCSPARSE ALLREDUCEBINOMIAL ALLREDUCESPARSEBIG ALLREDUCESPARSESMALL ALLREDUCESPARSERECDBL
STEPSIZE := DECREASING_STEPSIZES EXPBACKOFF_STEPSIZES
REPRESENTATION := DENSE SPARSE

.PHONY: depend clean all
all: LINREG LOGIT SVM $(SPLIT_DATA_TOOL) ${RANDOM_RESHUFFLE_TOOL} ${GENERATE_SPLIT_FILE_TOOL}

clean: CLEAN-LINREG CLEAN-LOGIT CLEAN-SVM
	rm -f ${BIN_DIR}/tools/*.out

LINREGS := $(foreach SY,$(STRATEGY), $(foreach SE, $(STEPSIZE), $(foreach RN, $(REPRESENTATION), $(BIN_DIR)/LINREG_$(SY)_$(SE)_$(RN))))
LOGITS := $(foreach SY,$(STRATEGY), $(foreach SE, $(STEPSIZE), $(foreach RN, $(REPRESENTATION), $(BIN_DIR)/LOGIT_$(SY)_$(SE)_$(RN))))
SVMS := $(foreach SY,$(STRATEGY), $(foreach SE, $(STEPSIZE), $(foreach RN, $(REPRESENTATION), $(BIN_DIR)/SVM_$(SY)_$(SE)_$(RN))))

ALLEXECS = $(LINREGS) $(LOGITS) $(SVMS)

CLEAN-LINREG:
	rm -f $(LINREGS)

CLEAN-LOGIT:
	rm -f $(LOGITS)

CLEAN-SVM:
	rm -f $(SVMS)

$(LINREGS): ${APP_DIR}/linreg.cpp ${APP_DIR}/linreg/linreg_exec.h $(SRC_FILES)
	mkdir -p $(BIN_DIR)
	${MPICC} $(CPP_FLAG) ${ADDITIONAL_FLAG} $(PARAMS) $(CPP_INCLUDE) $< -o $@ $(CPP_LAST)

$(LOGITS): ${APP_DIR}/logit.cpp ${APP_DIR}/logit/logit_exec.h $(SRC_FILES)
	mkdir -p $(BIN_DIR)
	${MPICC} $(CPP_FLAG) ${ADDITIONAL_FLAG} $(PARAMS) $(CPP_INCLUDE) $< -o $@ $(CPP_LAST)

$(SVMS): ${APP_DIR}/svm.cpp ${APP_DIR}/svm/svm_exec.h $(SRC_FILES)
	mkdir -p $(BIN_DIR)
	${MPICC} $(CPP_FLAG) ${ADDITIONAL_FLAG} $(PARAMS) $(CPP_INCLUDE) $< -o $@ $(CPP_LAST)

LINREG: $(LINREGS)
LOGIT: $(LOGITS)
SVM: $(SVMS)

# Add strategy specific options
#$(foreach f, $(ALLEXECS), $(if $(findstring _HOGWILD_, $f),$f,)): src/strategy/hogwild.h

#$(foreach f, $(ALLEXECS), $(if $(findstring _SYNC_, $f),$f,)): PARAMS += -D_SYNC
#$(foreach f, $(ALLEXECS), $(if $(findstring _SYNC_, $f),$f,)): src/strategy/parameter_server_sync.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ASYNC_, $f),$f,)): PARAMS += -D_ASYNC
$(foreach f, $(ALLEXECS), $(if $(findstring _ASYNC_, $f),$f,)): src/strategy/parameter_server_async.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ASYNCSPARSE_, $f),$f,)): PARAMS += -D_ASYNCSPARSE
$(foreach f, $(ALLEXECS), $(if $(findstring _ASYNCSPARSE_, $f),$f,)): src/strategy/parameter_server_async_sparse.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCE_, $f),$f,)): PARAMS += -D_ALLREDUCE
$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCE_, $f),$f,)): src/strategy/all_reduce_sync.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSEBIG_, $f),$f,)): PARAMS += -D_ALLREDUCESPARSE -D_C_BIG
$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSEBIG_, $f),$f,)): src/strategy/all_reduce_sparse.h src/strategy/c_allreduce/c_common.h src/strategy/c_allreduce/c_allreduce_big.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSESMALL_, $f),$f,)): PARAMS += -D_ALLREDUCESPARSE -D_C_SMALL
$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSESMALL_, $f),$f,)): src/strategy/all_reduce_sparse.h src/strategy/c_allreduce/c_common.h src/strategy/c_allreduce/c_allreduce_small.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSERECDBL_, $f),$f,)): PARAMS += -D_ALLREDUCESPARSE
$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSERECDBL_, $f),$f,)): src/strategy/all_reduce_sparse.h src/strategy/c_allreduce/c_common.h src/strategy/c_allreduce/c_allreduce_recdoubling.h

$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSERING_, $f),$f,)): PARAMS += -D_ALLREDUCESPARSE -D_C_RING
$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCESPARSERING_, $f),$f,)): src/strategy/all_reduce_sparse.h src/strategy/c_allreduce/c_common.h src/strategy/c_allreduce/c_allreduce_ring.h

#$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCEBINOMIAL_, $f),$f,)): PARAMS += -D_ALLREDUCEBINOMIAL
#$(foreach f, $(ALLEXECS), $(if $(findstring _ALLREDUCEBINOMIAL_, $f),$f,)): src/strategy/all_reduce_binomial_sync.h

# Add Stepsize behavior

$(filter %_DECREASING_STEPSIZES_DENSE %_DECREASING_STEPSIZES_SPARSE,$(ALLEXECS)): PARAMS += -D_DECREASING_STEPSIZES
$(filter %_EXPBACKOFF_STEPSIZES_DENSE %_EXPBACKOFF_STEPSIZES_SPARSE,$(ALLEXECS)): PARAMS += -D_EXPBACKOFF_STEPSIZES

# Add representation

$(filter %_DENSE,$(ALLEXECS)): PARAMS += -D_DENSE
$(filter %_SPARSE,$(ALLEXECS)): PARAMS += -D_SPARSE

${SPLIT_DATA_TOOL}: src/tools/split_data.cpp
	mkdir -p $(BIN_DIR)/tools/
	${CXX} ${CPP_FLAG} ${CPP_INCLUDE} $< -o $@

${RANDOM_RESHUFFLE_TOOL}: src/tools/random_reshuffle.cpp
	mkdir -p $(BIN_DIR)/tools/
	${CXX} ${CPP_FLAG} ${CPP_INCLUDE} $< -o $@

${GENERATE_SPLIT_FILE_TOOL}: src/tools/generate_split_file.cpp
	mkdir -p $(BIN_DIR)/tools/
	${CXX} ${CPP_FLAG} ${CPP_INCLUDE} $< -o $@
