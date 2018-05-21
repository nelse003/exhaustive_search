FROM nelse003/ccp4_phenix_docker
ADD ./compile_test.py
RUN mkdir XChemExplorer
ADD * /exhaustive_search/