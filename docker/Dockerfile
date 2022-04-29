FROM openjdk:8-jdk-alpine

ARG SAMTOOLS_VERSION=1.15 
ARG TRIMMOMATIC_VERSION=0.39
ARG BBMAP_VERSION=38.96
ARG PICARD_VERSION=2.26.11
ARG SEQKIT_VERSION=2.1.0
ARG IVAR_VERSION=1.3.1
ARG HTSLIB_VERSION=1.15
ARG SEQTK_VERSION=1.3

ENV PATH=$PATH:/opt/bbmap

COPY java_jar_helper.sh /tmp

RUN apk add --update --no-cache ncurses libbz2 xz zlib bash bc python3 \
    && apk add --virtual=deps --update --no-cache ncurses-dev musl-dev g++ make zlib-dev bzip2-dev xz-dev autoconf automake \
    # ---------- 
    # samtools
    && cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    # ---------- 
    # trimmomatic
    && cd /tmp \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && cat java_jar_helper.sh Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar > /usr/local/bin/trimmomatic \
    && chmod +x /usr/local/bin/trimmomatic \ 
    # ---------- 
    # bbmap
    && cd /tmp \
    && wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download -O BBMap_${BBMAP_VERSION}.tar.gz \
    && tar xzvf BBMap_${BBMAP_VERSION}.tar.gz \
    && mv bbmap /opt \
    # ---------- 
    # picard
    && cd /tmp \
    && wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
    && cat java_jar_helper.sh picard.jar > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard \
    # ---------- 
    # seqkit
    && cd /tmp \
    && wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz \
    && tar xzvf seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin \
    # ---------- 
    # ivar
    && cd /tmp \
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
    && make install \
    # ---------- 
    # seqtk
    && cd /tmp \
    && wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz \ 
    && tar xzvf v${SEQTK_VERSION}.tar.gz \ 
    && cd seqtk-${SEQTK_VERSION} \ 
    && make \
    && mv seqtk /usr/local/bin \ 
    # ----------
    # Clean up
    && rm -rf /var/cache/apk/* /tmp/* \
    && apk del deps
   

CMD ["bash"]
