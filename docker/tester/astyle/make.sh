#!/bin/bash
id="gassmoeller/aspect-tester:astyle"
echo "building: $id"
docker build -t $id .
