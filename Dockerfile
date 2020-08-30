FROM ubuntu:18.04
MAINTAINER CHEN Yuelong <yuelong.chen.btr@gmail.com>


ARG depends="wget build-essential python3 python3-dev python3-pip"

RUN apt update && apt install -y $depends


WORKDIR /opt/tmp/
ADD . ./
RUN pip3 install -r requirements.txt && \
    pip3 install . && \
    wget https://www.csie.ntu.edu.tw/~cjlin/libsvm/libsvm-3.24.tar.gz && \
    tar -zxvf libsvm-3.24.tar.gz && \
    cd libsvm-3.22 && \
    make clean && make && \
    cp svm-predict svm-train svm-scale /usr/local/bin/ && \
    chmod +x /usr/local/bin/*


RUN rm -rf /opt/tmp/

