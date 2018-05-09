#!/bin/bash

id="gassmoeller/aspect-tester:stable"
echo "running tester: $id"
docker run -e BUILDS='gcc' -it $id bash
