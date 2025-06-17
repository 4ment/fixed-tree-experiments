FROM continuumio/anaconda3:2022.10
LABEL "author"="Mathieu Fourment"
LABEL "company"="University of Technology Sydney"

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		build-essential \
		cmake \
		default-jdk \
		git \
		wget \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

RUN git clone https://github.com/4ment/fixed-tree-experiments /fixed-tree
RUN conda env create -f /fixed-tree/environment.yml
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate fixed-tree
WORKDIR /

RUN wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz && tar -xzf v4.0.1.tar.gz
RUN cmake -S /beagle-lib-4.0.1/ -B /beagle-lib-4.0.1/build -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF
RUN cmake --build /beagle-lib-4.0.1/build/ --target install

RUN wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.5pre_thorney_v0.1.2/BEASTv1.10.5pre_thorney_0.1.2.tgz
RUN tar -xzvf BEASTv1.10.5pre_thorney_0.1.2.tgz && ln -s /BEASTv1.10.5pre_thorney_0.1.2/bin/* /usr/local/bin

RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.2.6/iqtree-2.2.6-Linux.tar.gz
RUN tar -xzf iqtree-2.2.6-Linux.tar.gz && cp iqtree-2.2.6-Linux/bin/iqtree2 /usr/local/bin/

RUN wget https://github.com/tothuhien/lsd2/releases/download/v.2.3/lsd2_unix && chmod +x lsd2_unix && mv lsd2_unix /usr/local/bin/lsd2

RUN echo "source activate fixed-tree" > ~/.bashrc
ENV PATH /opt/conda/envs/fixed-tree/bin:$PATH