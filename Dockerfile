# To test 
# sudo docker build . -t exhaustive -f Dockerfile
# sudo docker run -it exhaustive /bin/bash

# Pulling from reskyner/ccp4 as can't get pandda update working
FROM reskyner/ccp4

# enivronement variable cannot be done by a source script
# due to how layers are produced in docker.
# These are copied from ccp4.setup-sh
# Relative calls muyst be made in seperate lines

ENV CCP4_MASTER=/ \
    HARVESTHOME=$HOME \
    GFORTRAN_UNBUFFERED_PRECONNECTED=Y \
    CCP4_SCR="/tmp/`whoami | tr ' \\\\' _`"

ENV CCP4=$CCP4_MASTER/ccp4

ENV CCP4I_TCLTK=$CCP4/bin \
    BALBES_ROOT=$CCP4_MASTER/BALBES \
    CBIN=$CCP4/bin \
    CLIB=$CCP4/lib \
    CLIBD=$CCP4/lib/data \
    CETC=$CCP4/etc \
    CINCL=$CCP4/include \
    CHTML=$CCP4/html \
    CEXAM=$CCP4/examples \
    CLIBS=$CCP4/lib/libccp4 \
    CPROG=$CCP4/src \
    CCP4I_TOP=$CCP4/share/ccp4i \
    CLIBD_MON=$CCP4/lib/data/monomers/ \
    CCP4_HELPDIR=$CCP4/help/

ENV MMCIFDIC=$CLIB/ccp4/cif_mmdic.lib \
    CRANK=$CCP4I_TOP/crank \
    MOSFLM_WISH=${CCP4I_TCLTK}/wish \
    CCP4_OPEN=UNKNOWN

# Make directory for repo files, and set to working dir
RUN mkdir ./exhaustive
WORKDIR ./exhaustive

# Note that running
#
# docker build .
#
# is important as:
#
# COPY requires:
#
# <src> must be relative to the source directory
# that is being built (the context of the build).
#
# This adds files from the repo
COPY . .

# Install the module
CMD ["/ccp4/bin/ccp4-python exhaustive/setup.py install"]





