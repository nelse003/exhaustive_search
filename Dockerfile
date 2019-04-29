FROM ubuntu

RUN df -h
RUN apt-get -qq update

RUN apt-get -qq install tzdata

RUN apt-get -qq -y install vim tar libgfortran3 gcc g++ m4 python2.7 git wget bzip2 tar expect
RUN wget http://devtools.fg.oisin.rc-harwell.ac.uk/nightly/7.0/ccp4-7.0-linux64-latest.tar.bz2
RUN bunzip2 ccp4-7.0-linux64-latest.tar.bz2
RUN mkdir ./ccp4
RUN tar -xf ccp4-7.0-linux64-latest.tar -C ./ccp4 --strip-components=1

ADD ccp4.setup-sh ./ccp4/bin

ADD pandda_update /
RUN ./pandda_update

CMD ["source", "/ccp4/bin/ccp4.setup-sh"]
CMD ["ccp4/start"]

CMD ["/ccp4/bin/ccp4-python setup.py install"]

RUN mkdir test
ADD * /test/
