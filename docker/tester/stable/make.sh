#!/bin/bash
id="gassmoeller/aspect-tester:stable2"
echo "building: $id"
docker build --no-cache -t $id .
