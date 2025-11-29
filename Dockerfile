# ------------------------------------------------------------
# 1. Base image with Python and C++ build tools
# ------------------------------------------------------------
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-setuptools \
    build-essential cmake git \
    && rm -rf /var/lib/apt/lists/*

# Create workspace for Python scripts
RUN mkdir -p /workspace
WORKDIR /workspace

# ------------------------------------------------------------
# 2. Install Python dependencies
# ------------------------------------------------------------
COPY python1/ /python1/
COPY python1_1/ /python1_1/
COPY python2/ /python2/

RUN pip3 install numpy \
    && pip3 install -r /python1/requirements.txt \
    && pip3 install -r /python1_1/requirements.txt \
    && pip3 install -r /python2/requirements.txt

RUN pip3 install numpy

# ------------------------------------------------------------
# 3. Copy the C++ repo (build + cinolib + src)
# ------------------------------------------------------------
COPY /CPP1 /CPP1

# ------------------------------------------------------------
# 4. Build triangulate_city
# ------------------------------------------------------------
RUN mkdir -p /CPP1/build \
    && cd /CPP1/build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release \
    && make -j4

# ------------------------------------------------------------
# 5. Entrypoint
# ------------------------------------------------------------
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
