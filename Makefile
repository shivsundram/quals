# Top-level quals tree Makefile. Provides convenience targets for things
# that span the lastz fork and the bench/ scripts.

LASTZ_DIR     := lastz
LASTZ_SRC_DIR := $(LASTZ_DIR)/src

# Object files to link the test driver against. Mirrors the lastz srcFiles
# list from $(LASTZ_SRC_DIR)/Makefile, MINUS lastz.o (which has its own main()).
# Minimal closure of object files needed by gapped_extend.o + ydrop_sane.o.
# Everything else (lastz.o main(), infer_scores, seed_search, chain, output
# formats, ...) is excluded; the unresolved symbols those files reference
# are stubbed in bench/test_ydrop.c.
TEST_YDROP_OBJS := \
    gapped_extend.o ydrop_sane.o                                                \
    edit_script.o segment.o                                                     \
    identity_dist.o coverage_dist.o continuity_dist.o                           \
    utilities.o dna_utilities.o sequences.o

TEST_YDROP_OBJS_FULL := $(addprefix $(LASTZ_SRC_DIR)/,$(TEST_YDROP_OBJS))

TEST_YDROP_CFLAGS := -O2 -Wall -Wextra -Werror -g \
                     -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE \
                     -Dscore_type=I

TEST_YDROP_INCLUDES := -I$(LASTZ_SRC_DIR)

.PHONY: test_ydrop run_test_ydrop test_ydrop_clean lastz_objs

# Builds bench/test_ydrop, linking against the standard (non-instrumented)
# lastz objects. Depends on `make -C lastz` having been run.
test_ydrop: bench/test_ydrop

bench/test_ydrop: bench/test_ydrop.c $(TEST_YDROP_OBJS_FULL)
	gcc $(TEST_YDROP_CFLAGS) $(TEST_YDROP_INCLUDES) \
	    bench/test_ydrop.c $(TEST_YDROP_OBJS_FULL) -lm -o $@

# Force-rebuild the standard lastz objects via the lastz Makefile.
lastz_objs:
	$(MAKE) -C $(LASTZ_DIR) build_lastz

# Convenience: build then run.
run_test_ydrop: test_ydrop
	./bench/test_ydrop

test_ydrop_clean:
	rm -f bench/test_ydrop
