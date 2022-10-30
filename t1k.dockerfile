# FROM alpine
# RUN apk add g++ make zlib-dev perl curl
FROM ubuntu:20.04
RUN apt update && \
    apt install -y  g++-7 make zlib1g-dev curl && \
    ln -s /usr/bin/g++-7 /usr/bin/g++
