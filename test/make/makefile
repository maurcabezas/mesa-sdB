
MESA_DIR = ../../..

#################################################################

# STEP 1: set up the standard info

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lbinary
LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_STAR)

#################################################################

# STEP 3: build the test program

TEST_OBJS = test_binary.o
TEST = tester

$(TEST) : $(TEST_OBJS)
	$(FC) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) \
	$(LOAD_LOCAL) $(LOAD_OTHER)
	
#################################################################

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod $(TEST_LIB) $(TEST_DIR)/$(TEST)
