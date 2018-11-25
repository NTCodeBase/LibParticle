################################################################################
COMPILER := g++-8
ALL_CCFLAGS := -W -O3 -flto -lstdc++fs
ALL_CCFLAGS += -std=c++17

################################################################################
ROOT_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
INCLUDES += -I$(ROOT_PATH)/
INCLUDES += -I$(ROOT_PATH)/LibParticle/PartioBgeo
INCLUDES += -I$(ROOT_PATH)/../LibCommon
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/glm
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/json/single_include/nlohmann
INCLUDES += -I$(ROOT_PATH)/../LibCommon/Externals/spdlog/include

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
	ar -rsv $(OUTPUT_DIR)/$(LIB_NAME) $(COMPILE_OBJ)
	ranlib $(OUTPUT_DIR)/$(LIB_NAME)

create_out_dir:
	mkdir -p $(OUTPUT_DIR)
	mkdir -p $(OBJ_DIR)
	mkdir -p $(COMPILE_OBJ_SUBDIR)

$(COMPILE_OBJ): %.o: $(patsubst $(OBJ_DIR)/%, $(ROOT_PATH)/%, $(patsubst %.o, %.cpp, $(COMPILE_OBJ)))
	$(COMPILER) $(INCLUDES) $(ALL_CCFLAGS) -c $< -o $@

clean:
	rm -rf $(OUTPUT_DIR)/LibParticle/*
