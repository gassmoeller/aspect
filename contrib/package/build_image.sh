#!/bin/bash

docker run -v /home/rene/software/:/home/software --name=aspect-runtime --rm -it --device /dev/fuse ubuntu:14.04 bash build.sh 
