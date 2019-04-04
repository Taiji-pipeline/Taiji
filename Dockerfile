FROM centos:7.5.1804

ENV PATH="/root/.local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/lib64:/usr/local/lib:${LD_LIBRARY_PATH}"
ENV LIBRARY_PATH="/usr/local/lib:${LIBRARY_PATH}"

RUN yum install -y gcc gcc-c++ make \
    glibc-static libstdc++-static zlib-static expat-static \
    gmp-static cairo-devel pango-devel libxml2-devel gmp-devel \
    git sudo

RUN mkdir -p ~/igraph && \
    curl -Lk http://igraph.org/nightly/get/c/igraph-0.7.1.tar.gz | \
    tar xz --strip-components=1 -C ~/igraph && \
    cd ~/igraph && ./configure && make && make install && cd .. && rm -rf ~/igraph

RUN mkdir -p ~/.local/bin
RUN curl -Lk https://www.stackage.org/stack/linux-x86_64 | \
    tar xz --strip-components=1 -C ~/.local/bin

RUN git clone https://github.com/Taiji-pipeline/Taiji.git && cd Taiji
RUN stack install