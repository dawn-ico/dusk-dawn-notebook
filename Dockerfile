# IMAGE needs be be set to one of the docker/dawn-env.dockerfile images
FROM gtclang/dawn

RUN pip install dusk@git+https://github.com/dawn-ico/dusk.git
CMD /bin/bash
