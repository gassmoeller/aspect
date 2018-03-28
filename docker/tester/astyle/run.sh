#!/bin/bash

id="gassmoeller/aspect-tester:astyle"
echo "running tester: $id"
docker run -e BUILDS='gcc' -it $id bash
