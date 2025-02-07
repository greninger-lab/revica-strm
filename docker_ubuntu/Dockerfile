FROM ubuntu:20.04

ARG SAMTOOLS_VERSION=1.15 
ARG BCFTOOLS_VERSION=1.15
ARG HTSLIB_VERSION=1.15
ARG TRIMMOMATIC_VERSION=0.39
ARG PICARD_VERSION=2.26.11
ARG SEQKIT_VERSION=2.1.0
ARG IVAR_VERSION=1.3.1
ARG SEQTK_VERSION=1.3
ARG PANDEPTH_VERSION=2.25

ARG BWA_MEM2_VERSION=2.2.1

ENV DEBIAN_FRONTEND=noninteractive 
ENV PATH=$PATH:/opt/bbmap
ENV TRIMMOMATIC=/usr/local/bin/trimmomatic.jar
ENV PICARD=/usr/local/bin/picard.jar
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

RUN apt-get update && apt-get install -y \
    openjdk-11-jdk bc python3 wget libncurses5-dev bzip2 gcc g++ make zlib1g-dev libbz2-dev liblzma-dev unzip autoconf

# ---------- 
# Samtools
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install

# ---------- 
# Bcftools
RUN cd /tmp \
    && wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar xjvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install

# ---------- 
# Trimmomatic
RUN cd /tmp \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && mv Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar /usr/local/bin/trimmomatic.jar

# ---------- 
# BWA-MEM2
RUN cd /tmp \
    && wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM2_VERSION}/bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 \
    && tar jxf bwa-mem2-${BWA_MEM2_VERSION}_x64-linux.tar.bz2 \
    && mv bwa-mem2-2.2.1_x64-linux/bwa-mem2 /usr/local/bin/

# ---------- 
# Picard
RUN cd /tmp \
    && wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
    && mv picard.jar /usr/local/bin

# ---------- 
# Seqkit
RUN cd /tmp \
    && wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz \
    && tar xzvf seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin

# ---------- 
# iVar and HTSLib
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar xjvf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /tmp \
    && wget https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz \
    && tar xzvf v${IVAR_VERSION}.tar.gz \
    && cd ivar-${IVAR_VERSION} \
    && ./autogen.sh \
    && ./configure \
    && make \
    && make install

# ---------- 
# Seqtk
RUN cd /tmp \
    && wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz \ 
    && tar xzvf v${SEQTK_VERSION}.tar.gz \ 
    && cd seqtk-${SEQTK_VERSION} \ 
    && make \
    && mv seqtk /usr/local/bin

#----------
#PanDepth
RUN cd /tmp \
    && wget https://github.com/HuiyangYu/PanDepth/releases/download/v${PANDEPTH_VERSION}/PanDepth-${PANDEPTH_VERSION}-Linux-x86_64.tar.gz \
    && tar zvxf PanDepth-${PANDEPTH_VERSION}-Linux-x86_64.tar.gz \
    && cd PanDepth-${PANDEPTH_VERSION}-Linux-x86_64 \
    && mv pandepth /usr/local/bin

# ----------
# Clean up
RUN apt autoremove -y wget libncurses5-dev bzip2 gcc g++ make zlib1g-dev libbz2-dev liblzma-dev unzip autoconf \
    && rm -rf /var/lib/apt/lists/* /tmp/*

CMD ["bash"]
