FROM ubuntu

RUN df -h
RUN apt-get -qq update

RUN apt-get -qq install tzdata

RUN apt-get -qq -y install vim tar libgfortran3 gcc g++ m4 python2.7 git wget bzip2 tar expect
RUN wget http://devtools.fg.oisin.rc-harwell.ac.uk/nightly/ccp4-linux64-latest.tar.bz2
RUN bunzip2 ccp4-linux64-latest.tar.bz2
RUN mkdir ./ccp4
RUN tar -xf ccp4-linux64-latest.tar -C ./ccp4 --strip-components=1

ADD ccp4.setup-sh ./ccp4/bin

ADD pandda_update /
RUN ./pandda_update

CMD ["source", "/ccp4/bin/ccp4.setup-sh"]
CMD ["ccp4/start"]

RUN mkdir exhaustive_search
ADD * /exhaustive_search/

RUN mkdir exhaustive_search/test
ADD * /exhaustive_search/test
