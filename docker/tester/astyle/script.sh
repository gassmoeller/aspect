#!/bin/bash

# Automated script to compile&test ASPECT in different configurations

# Mounts:
# writable log dir: /home/bob/log/
# ASPECT source repo: /home/bob/source/ (writeable) or /source (can be read-only)

# settings from ENV variables:
# BUILDS:
#  (gcc|clang)[petsc]|astyle
# QUIET=1 (default, not) -- only print summary

# example usage:
# 1) mount readonly, you need to copy results out (prefered):
#      docker run -it --rm -v "$(pwd):/source:ro" tjhei/aspect-tester-8.4.1 /bin/bash
#      BUILDS="gcc" ./script.sh
#      docker cp CONTAINER:/home/bob/log/changes-BUILD.diff . # from outside
# 2) mount writeable, this will modify your files outside the container
#      docker run -it --rm -v "$(pwd):/home/bob/source" tjhei/aspect-tester-8.4.1 /bin/bash
#      BUILDS="clang" ./script.sh



mkdir -p ~/log

if [ -s ~/source/CMakeLists.txt ]; then
  touch ~/source/VERSION 2>/dev/null || { echo "~/source needs to be mounted R/W, aborting."; exit 1; }
fi

git clone --depth=1 https://github.com/geodynamics/aspect.git ~/source

if [ -z "$PULL_REQUEST" ]; then
  cd ~/source
  git fetch --depth=1 origin pull/${PULL_REQUEST}/head:branch
  git checkout branch
  cd ..
fi

submit="OFF"

run()
{
build=$1
desc=$2
submit=$3

LOGFILE="`pwd`/changes.diff"
cd ~/source
./doc/indent || { echo "indent FAILED"; return; }
git diff >$LOGFILE
echo "git diff >$LOGFILE"
git diff --exit-code --name-only || { echo "FAILED: `git diff --name-only`"; return; }
echo "ok"
return
}


summary=~/log/summary
indexhtml=~/log/index.html

main()
{
#clean contents:
> $summary

echo "BUILD $build:" |tee -a $summary
logfile=~/log/log-$build
mkdir -p build-$build
cd build-$build
eval run $build $build$name $submit 2>&1 | tee $logfile
if [ -s changes.diff ]; then
  cp changes.diff ~/log/changes-$build.diff
  echo "DIFFS: changes-$build.diff" | tee -a $logfile
fi
cd ..
rep "FAILED" $logfile | grep -v "FAILED: /" | grep -v "The following tests FAILED" | grep -v "FAILED: cd /" | tee -a $summary

grep "^ok$" $logfile | tee -a $summary
grep "tests passed" $logfile | tee -a $summary

sed -i 's/[[:space:]]*0 Compiler errors/ok/' $summary
sed -i 's/\([0-9]*\)% tests passed, 0 tests failed out of \([0-9]*\)/tests: \2 passed/' $summary 

sed -i 's/\([0-9]*\)% tests passed, \([0-9]*\) tests failed out of \([0-9]*\)/tests: \2 \/ \3 FAILED/' $summary 

cp $summary $indexhtml

grep -h "DIFFS: changes-" ~/log/log-* | tee -a $indexhtml
sed -i 's#$# <br/>#' $indexhtml
sed -i 's#^BUILD \(.*\):#<a href="log-\1">BUILD \1:</a>#' $indexhtml
sed -i 's#^DIFFS: \(.*diff\)#DIFFS: <a href="\1">\1</a>#' $indexhtml
}

if [ "$QUIET" == "1" ];
then
  main >/dev/null
  cat $summary
else
  main
fi

