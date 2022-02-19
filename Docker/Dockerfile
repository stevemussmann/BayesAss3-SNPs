FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG USERNAME=ba3user
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG IMAGE_NAME=ba3snps
ARG IMAGE_TAG=1.0
ENV USER $USERNAME
ENV HOME /home/$USERNAME
ENV IMAGE_NAME $IMAGE_NAME
ENV IMAGE_TAG $IMAGE_TAG

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev git autoconf automake vim wget less libgsl-dev libboost-dev libboost-program-options-dev perl

## Install ba3-snps
RUN mkdir -p /app/scripts/python /app/scripts/perl /app/src
WORKDIR /app/src
RUN git clone https://github.com/stevemussmann/BayesAss3-SNPs.git
WORKDIR /app/src/BayesAss3-SNPs
RUN chmod a+x countLociImmanc.sh
RUN make
RUN make install

## Install ba3-snps-autotune
WORKDIR /app/scripts/python
RUN git clone https://github.com/stevemussmann/BA3-SNPS-autotune.git

## Install file converters
WORKDIR /app/scripts/perl
RUN git clone https://github.com/stevemussmann/file_converters.git

## link binaries and scripts
# make bin directory
RUN mkdir -p /app/bin
WORKDIR /app/bin
RUN ln -s /app/scripts/python/BA3-SNPS-autotune/BA3-SNPS-autotune.py
RUN ln -s /app/src/BayesAss3-SNPs/countLociImmanc.sh
RUN ln -s /app/scripts/perl/file_converters/stacksStr2immanc.pl
RUN ln -s /app/scripts/perl/file_converters/pyradStr2immanc.pl

## Move stuff around
RUN mkdir -p /app/data

RUN groupadd --gid $USER_GID $USERNAME \
	&& adduser --uid $USER_UID --gid $USER_GID --disabled-password $USERNAME \
	--gecos "First LAST,RoomNumber,WorkPhone,HomePhone" \
	&& apt-get update \
	&& chown -R $USERNAME:$USERNAME /home/$USERNAME \
	&& chown $USERNAME /app

#echo /app/bin to path
RUN echo "export PATH=/app/bin:${PATH}" >> /home/$USERNAME/.bashrc

USER $USER

WORKDIR /app/data

