ROOT := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

all:
	@make -C src ROOT=$(ROOT) OPT=1

clean:

	@make -C src ROOT=$(ROOT) clean
	@make -C tests ROOT=$(ROOT) clean

install:
	@make -C src ROOT=$(ROOT) OPT=1 install

test:
	@make -C tests ROOT=$(ROOT) OPT=1
.PHONY: clean

distclean: clean
	@rm -rf $(ROOT)/bin
.PHONY: distclean
