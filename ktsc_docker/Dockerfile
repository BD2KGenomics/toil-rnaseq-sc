FROM ubuntu:xenial

RUN echo "deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
RUN apt-get update ; apt-get install -y --allow-unauthenticated r-base libxml2-dev libcurl4-openssl-dev	

COPY ["install.r","/image/"]
RUN Rscript /image/install.r

COPY ["ktsc.r", "/image/"]
COPY ["entrypoint.sh", "/image/"]
ENTRYPOINT ["sh", "/image/entrypoint.sh"]