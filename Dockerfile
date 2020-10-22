FROM gtclang/dawn:latest

# setup user:
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}

RUN apt-get install ffmpeg

# setup python packages & AtlasUtils:
RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir jupyterlab && \
    pip install matplotlib && \
    pip install dusk@git+https://github.com/dawn-ico/dusk.git && \
    cd ${HOME}/AtlasUtils/utils/ && chmod +x ./build_and_install.sh && ./build_and_install.sh

USER ${NB_USER}

WORKDIR ${HOME}/content
