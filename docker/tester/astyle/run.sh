#!/bin/bash

id="gassmoeller/aspect-tester:astyle"
echo "running tester: $id"
docker run -it $id bash
