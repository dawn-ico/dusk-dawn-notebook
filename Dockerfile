FROM gtclang/dawn:latest

# setup jupyter lab user:
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
ENV PATH="/usr/local/dawn/bin:${PATH}"

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# copy repo files into docker container
COPY . ${HOME}

# setup things as root:
# (because bash doesn't support comments in multi-line commands,
# we use this weird `: '...'` syntax)
RUN \
    apt-get update && \
    : 'ffmpeg is required for animations in exercises' && \
    : 'llvm & clang are required by dawn' && \
    apt-get install -y \
      ffmpeg \
      llvm \
      llvm-dev \
      clang \
      clang-format \
      libclang-dev \
      libclang-cpp10-dev && \
    : 'upgrade some pip stuff' && \
    pip install --upgrade pip setuptools wheel && \
    : 'install python dependencies' && \
    pip install --no-cache-dir \
        jupyterlab \
        matplotlib \
        dusk@git+https://github.com/dawn-ico/dusk.git && \
    : 'setup AtlasUtils' && \
    cd ${HOME}/AtlasUtils/utils/ && \
    chmod +x ./build_and_install.sh && \
    ./build_and_install.sh && \
    chown -R ${NB_UID} ${HOME} && \
    :

USER ${NB_USER}

WORKDIR ${HOME}/content

RUN \
    : 'Trust all jupyter notebooks by defaults' && \
    find . -name "*.ipynb" -type f | xargs -n1 jupyter trust && \
    :
