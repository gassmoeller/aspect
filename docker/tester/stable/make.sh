#!/bin/bash
id="gassmoeller/aspect-tester:stable"
echo "building: $id"
docker build -t $id .
