# IMAGE needs be be set to one of the docker/dawn-env.dockerfile images
FROM gtclang/dawn:latest

ARG NB_USER=lambda
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

RUN  pip install --no-cache-dir notebook==5.* && pip install dusk@git+https://github.com/dawn-ico/dusk.git
CMD /bin/bash
