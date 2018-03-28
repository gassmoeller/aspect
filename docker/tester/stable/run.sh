#!/bin/bash

id="gassmoeller/aspect-tester:8.5.0"
echo "running tester: $id"
docker run -it $id bash
