FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

MAINTAINER rene.gassmoeller@mailbox.org

USER root
ADD script.sh script.sh

RUN chmod u+x *.sh
CMD /home/dealii/script.sh
