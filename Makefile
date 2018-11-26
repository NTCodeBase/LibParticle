################################################################################
GCC_PREFIX  ?=
GCC_SUFFIX  ?=
COMPILER    := $(GCC_PREFIX)g++$(GCC_SUFFIX)
ALL_CCFLAGS := -g -W -O3 -lstdc++fs
ALL_CCFLAGS += -std=c++17
ALL_CCFLAGS += $(FLAG_FLTO)

################################################################################
ROOT_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
INCLUDES += -I$(ROOT_PATH)/
INCLUDES += -I$(ROOT_PATH)/LibParticle/PartioBgeo
INCLUDES += -I$(ROOT_PATH)/../LibCommon
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/glm
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/json/single_include/nlohmann
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/spdlog/include
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/tbb_linux/include

################################################################################
OUTPUT_DIR := ../../Build/Linux
OBJ_DIR    := ../../Build/Linux/OBJS
LIB_NAME   := libLibParticle.a

LIB_SRC     := $(shell find $(ROOT_PATH)/LibParticle -name *.cpp)
LIB_OBJ     := $(patsubst %.cpp, %.o, $(LIB_SRC))
COMPILE_OBJ := $(patsubst $(ROOT_PATH)/%, $(OBJ_DIR)/%, $(LIB_OBJ))
COMPILE_OBJ_SUBDIR  := $(patsubst $(ROOT_PATH)/%, $(OBJ_DIR)/%, $(dir $(LIB_OBJ)))

################################################################################
all: create_out_dir $(COMPILE_OBJ)
	$(GCC_PREFIX)gcc-ar$(GCC_SUFFIX) -rsv $(OUTPUT_DIR)/$(LIB_NAME) $(COMPILE_OBJ)
	$(GCC_PREFIX)gcc-ranlib$(GCC_SUFFIX) $(OUTPUT_DIR)/$(LIB_NAME)

create_out_dir:
	mkdir -p $(OUTPUT_DIR)
	mkdir -p $(OBJ_DIR)
	mkdir -p $(COMPILE_OBJ_SUBDIR)

$(COMPILE_OBJ): $(OBJ_DIR)/%.o: $(ROOT_PATH)/%.cpp
	$(COMPILER) $(INCLUDES) $(ALL_CCFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)/LibParticle
	rm $(OUTPUT_DIR)/$(LIB_NAME)
