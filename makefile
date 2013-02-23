# Copyright (c) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>
#
# This file is part of hVMC.
#
# hVMC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hVMC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hVMC.  If not, see <http://www.gnu.org/licenses/>.


# ----------- BUILD SETTINGS -----------

# disable default implicit rules
MAKEFLAGS += --no-builtin-rules
.SUFFIXES:

# find the hash of the git commit we are building from
GIT_HASH = $(shell git rev-parse --short HEAD 2> /dev/null || echo *unknown*)

# common compiler/linker flags
CXX      = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDFLAGS  = -lboost_program_options -lboost_filesystem -lboost_system
LDFLAGS += -lboost_serialization -lboost_mpi
DEFINES  = -DGIT_HASH=\"$(GIT_HASH)\"

# feature specific compiler/linker flags
BUILD   ?= RELEASE
DBLPREC ?= DISABLED
CBLAS   ?= ENABLED

ifeq ($(DBLPREC), ENABLED)
  DEFINES += -DUSE_FP_DBLPREC
endif

ifneq ($(CBLAS), DISABLED)
  DEFINES += -DUSE_CBLAS
  LDFLAGS += -lcblas
else
  DEFINES += -DEIGEN_DEFAULT_TO_ROW_MAJOR
endif

# build [debug/profile/release] specific compiler/linker flags
ifeq ($(BUILD), DEBUG)
  CXXFLAGS += -g
  DEFINES  += -DVERBOSE=1
else ifeq ($(BUILD), PROFILE)
  CXXFLAGS += -march=native -O2 -g
  DEFINES  += -DNDEBUG
else
  CXXFLAGS += -march=native -O2 -flto
  LDFLAGS  += -fuse-linker-plugin -s
  DEFINES  += -DNDEBUG
endif


# -------------- TARGETS ---------------

# generate lists of object files (all and existing only)
OBJECTS_ALL = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
OBJECTS_EXT = $(wildcard *.o)

# define implicit rule to compile a cpp file
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $< -o $@
	$(CXX) -MM $(CXXFLAGS) $(DEFINES) $< > $*.d

# hVMC executable
hvmc : $(OBJECTS_ALL)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# pull in dependency info for existing object files
-include $(OBJECTS_EXT:.o=.d)

clean :
	rm -f hvmc *.o *.d
