FROM github.com/nelse003/ccp4_phenix_docker
ADD ./compile_test.py
RUN mkdir exhaustive_search
ADD * /exhaustive_search/
