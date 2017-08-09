FROM ubuntu:xenial

RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y python-pip

COPY ["dep", "/dep/"]

ENTRYPOINT ["sh", "/dep/entrypoint.sh"]