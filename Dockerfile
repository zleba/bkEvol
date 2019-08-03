#Official ROOT Docker image
FROM rootproject/root-ubuntu16

# Run the following commands as super user (root):
USER root 
COPY requirements.txt .
RUN  mkdir /workdir  &&\
     apt-get update && \
     apt-get install -y wget \
                        libopenmpi-dev \
                        libhdf5-dev \
                        libopenblas-dev \
                        libboost-dev \
                        libopenblas-dev \
                        python-pip \
                        python-tk \
     && sudo -H pip install --upgrade pip \
     && sudo -H pip install --trusted-host pypi.python.org -r requirements.txt  \
     && rm -rf /var/lib/apt/lists/* 

#Instal LHAPDF
WORKDIR /
COPY installARMA.sh .
RUN ./installARMA.sh
ENV ARMADIR /arma


WORKDIR /bkevol
RUN mkdir -p /bkevol/obj  /bkevol/src  /bkevol/inc  /bkevol/pybind
COPY Makefile  .
COPY src src/
COPY inc inc/
COPY pybind pybind/

RUN make -j`nproc` bkevol.so
ENV PYTHONPATH /bkevol:$PYTHONPATH
ENV LD_PRELOAD  /usr/lib/libmpi_cxx.so

WORKDIR /workdir
#ENV PATH="/lhapdf/install/bin/:${PATH}"

# When starting the container and no command is started, run bash
CMD ["/bin/bash"]



ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password --gecos "Default user" \
            --uid ${NB_UID} ${NB_USER}

#RUN sudo -H pip install pyqt5

WORKDIR ${HOME}
USER ${NB_USER}
RUN  mkdir .jupyter && echo "c.NotebookApp.token = ''" > ${HOME}/.jupyter/jupyter_notebook_config.py
RUN  mkdir -p ${HOME}/examplesNb   ${HOME}/examplesPy ${HOME}/temp 
COPY --chown=jovyan examplesNb ${HOME}/examplesNb
COPY --chown=jovyan examplesPy ${HOME}/examplesPy
EXPOSE 8888

# When starting the container and no command is started, run bash
CMD ["jupyter", "notebook",  "--ip", "0.0.0.0"]

