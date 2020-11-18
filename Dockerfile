FROM gtclang/dawn:latest

# setup jupyter lab user:
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

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
    chown -R ${NB_UID} ${HOME} && \
    apt-get update && \
    : 'ffmpeg is required for animations in exercises' && \
    : 'nodejs is required for the auto scroll extension' && \
    apt-get install -y \
      nodejs \
      npm \
      ffmpeg \
      clang-format && \
    : 'upgrade some pip stuff' && \
    pip install --upgrade pip setuptools wheel && \
    : 'install python dependencies' && \
    pip install --no-cache-dir \
        jupyterlab \
        matplotlib \
        dusk@git+https://github.com/dawn-ico/dusk.git && \
    : 'Enable auto scroll extension' && \
    jupyter labextension install @wallneradam/output_auto_scroll && \
    : 'setup AtlasUtils' && \
    cd ${HOME}/AtlasUtils/utils/ && \
    chmod +x ./build_and_install.sh && \
    ./build_and_install.sh && \
    :

USER ${NB_USER}

WORKDIR ${HOME}/content

RUN \
    : 'Trust all jupyter notebooks by defaults' && \
    find . -name "*.ipynb" -type f | xargs -n1 jupyter trust && \
    :
