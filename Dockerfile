# To test 
# sudo docker build . -t exhaustive -f Dockerfile
# sudo docker run -it exhaustive /bin/bash
#
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
# Pulling from reskyner/ccp4 as can't get pandda update working

FROM reskyner/ccp4

RUN mkdir ./exhaustive
WORKDIR ./exhaustive
COPY . . 

CMD ["ccp4/start"]
CMD ["/ccp4/bin/ccp4-python exhaustive/setup.py install"]





