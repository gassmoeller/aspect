FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

MAINTAINER rene.gassmoeller@mailbox.org

#User and group of the 'jenkins' user
USER root
RUN adduser --uid 10000 jenkins

ENV OMPI_MCA_btl_base_warn_component_unused=0
