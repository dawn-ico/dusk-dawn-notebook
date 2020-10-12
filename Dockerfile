# IMAGE needs be be set to one of the docker/dawn-env.dockerfile images
FROM gtclang/dawn-env-ubuntu20.04
ARG BUILD_TYPE=Release
COPY . /usr/src/dawn
RUN /usr/src/dawn/scripts/build-and-test \
    --dawn-install-dir /usr/local/dawn \
    --parallel $(nproc) \
    --docker-env \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE

RUN pip install dusk@git+https://github.com/dawn-ico/dusk.git
CMD /bin/bash
