===========================
CosmicFish: Information Gain
===========================

This is the CosmicFish package used to produce the results of the CosmicFish release paper.

To produce all results just issue::

	make all

in this folder.

To properly work the script requires an environment variable called ``$COSMICFISH_DIR``
that points to the cosmicfish directory.

If you do not have this environment variable either issue::

	export $COSMICFISH_DIR=path_to_cosmicfish

Or go to ``script/common.sh`` and point ``COSMICFISH_PATH`` to the right direction.

1. Makefile targets:
====================

The Makefile has several targets. Here's a breakdown of the main ones:

**Main**:

* ``all``: runs all examples. Creates the Fisher matrices and then analyses them;
* ``fisher_matrices``: just creates all the Fisher matrices;
* ``analysis``: analyses all the Fisher matrices;
* ``clean``: removes just raw results but leaves analysis results there;
* ``deep_clean``: removes all results;

**Secondary**:

* ``targets``: creates Makefile targets with the names of the parameters. After issuing this one can run only one choice of parameters;
* ``additional_script``: this runs all the python script in the script folder. By default does nothing but the user might use it to include his own python script in his pipeline;
