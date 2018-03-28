#!/bin/bash
id="gassmoeller/aspect-tester:8.5.0"
echo "building: $id"
docker build -t $id .
