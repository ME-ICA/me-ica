FROM ubuntu:20.04

# Make me-ica binaries work
RUN apt-get update
RUN apt-get -y install \
    gcc mono-mcs wget libxm4 libmotif-* libgsl23* libglw1-mesa* libglu1-mesa* tcsh python3-pip python3.9
RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23.1.0 /usr/lib/libgsl.so.19
RUN mkdir -p /src/bin/afni
WORKDIR /src/bin/afni
RUN wget https://afni.nimh.nih.gov/pub/dist/tgz/linux_ubuntu_16_64.tgz
RUN tar zxvf linux_ubuntu_16_64.tgz
ENV PATH=$PATH:/src/bin/afni/linux_ubuntu_16_64

# poetry # static version
RUN mkdir -p /app
WORKDIR /app

RUN python3.9 -m pip install poetry==1.8.2
COPY pyproject.toml /pyproject.toml
COPY poetry.lock /poetry.lock
COPY README.md /README.md
RUN poetry install --no-root

RUN ln -s $(which python3.9) /usr/bin/python
# meica files
COPY libmeica libmeica
COPY meica.py meica.py
COPY tedana.py tedana.py

ENTRYPOINT [ "poetry", "run" ]
# CMD will get overwritten by caller of the container
CMD [ "" ]
