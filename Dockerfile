FROM gtclang/dawn:latest

# setup python packages:
RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir jupyterlab && \
    pip install dusk@git+https://github.com/dawn-ico/dusk.git

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
USER ${NB_USER}
