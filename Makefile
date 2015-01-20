ROOT := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

all:
	@make -C src ROOT=$(ROOT) OPT=1

clean:
	@make -C src ROOT=$(ROOT) clean
	@make -C tests ROOT=$(ROOT) clean

test:
	@make -C tests ROOT=$(ROOT) OPT=1
.PHONY: clean

distclean: clean
	@rm -rf $(AMORDAD_ROOT)/bin
.PHONY: distclean
