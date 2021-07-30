FROM ubuntu:20.04 AS builder

ENV TZ=Asia/Seoul
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN  apt-get update
RUN apt-get --assume-yes install software-properties-common
RUN apt-get install -y gdal-bin libgdal-dev libomp-dev cmake make
RUN dpkg -l | grep -i gdal
RUN apt-get --assume-yes install build-essential
#RUN apt-get install -y gcc gpp

# Copy everything
COPY . /code

RUN mkdir /code/build 
WORKDIR /code/build
RUN cmake .. && make

FROM ubuntu:20.04 
ENV TZ=Asia/Seoul
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN  apt-get update
RUN apt-get --assume-yes install libgdal-dev libomp-dev
COPY --from=builder /code /code
WORKDIR /code
RUN ln -sfnv /code/build/dem2mesh /usr/bin/dem2mesh
