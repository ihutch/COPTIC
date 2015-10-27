default 	:
	cd src; make

# Rule for non-predefined cases. Pass it up to src make.
%	:
	cd src; make $*

help	:
	@echo 'For other make commands $$ make <target> or $$ make -C <directory> [<target>]'

