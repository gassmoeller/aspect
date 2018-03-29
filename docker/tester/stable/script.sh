#!/bin/bash

# Automated script to compile&test ASPECT in different configurations

# Mounts (optional):
# writable log dir: /home/dealii/log/
# if this directory is mounted it will be filled with the test results
# if this directory is not mounted, it will be created internally and destroyed after the run

# settings from ENV variables:
# BUILDS:
#  (gcc|clang)[petsc]
# QUIET=1 (default, not) -- only print summary
# PULL_REQUEST: if set, will trigger a test of this pull request number
# USER, BRANCH: [not yet implemented] if set, will trigger a test of this branch of this github user
# ASPECT_DIR: [not yet implemented] if set, will trigger a test of this aspect directory (has to be mounted in the container
# if neither PULL_REQUEST nor USER, BRANCH is set, it will test current master

# example usage:
# 1) test pull request
#      docker run -it --rm -v "$(pwd):/home/dealii/log" -e BUILDS='gcc' gassmoeller/aspect-tester:8.5.0
#
# 2) mount writeable, this will modify your files outside the container
#      docker run -it --rm -v "$(pwd):/home/bob/source" gassmoeller/aspect-tester:8.5.0 /bin/bash
#      BUILDS="clang" ./script.sh


BUILDS=gcc

if [ -z "$BUILDS" ]; then
  echo 'Please specify list of builds to do in the ENV variable $BUILDS.'
  echo "Separate build with spaces. Valid options: gcc gccpetsc clang clangpetsc"
  exit 0
fi

mkdir -p ~/log

submit="OFF"

run()
{
build=$1
desc=$2
submit=$3

petsc="OFF"
if [[ $build =~ .*petsc.* ]];
then
  petsc="ON"
fi
compiler=""
if [[ $build =~ .*clang.* ]];
then
 compiler=" -D DEAL_II_DIR=~/deal.II/installed-v8.4.1-clang "
fi

echo "+ cmake -G "Ninja" $compiler -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=$petsc -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON ~/source"
time cmake -G "Ninja" $compiler -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=$petsc -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON ~/source || { echo "configure FAILED"; return; }

echo "+ ninja"
time ninja || { echo "build FAILED"; return; }

# build all libs first to work around parallel ninja quirks, ignore errors, be silent
echo "prebuilding tests..."
cd tests

time ninja -k 0 tests >/dev/null 2>&1

#echo "+ make -j $NPROC -k tests >/dev/null 2>&1"
#time make -j $NPROC -k tests >/dev/null 2>&1
cd .. 

echo "+ ctest --output-on-failure -j $NPROC"
time ctest --output-on-failure -j $NPROC || { echo "test FAILED"; }

echo "+ ninja generate_reference_output"
time ninja generate_reference_output
echo "ok"
}


summary=~/log/summary
indexhtml=~/log/index.html

main()
{
#clean contents:
> $summary

for build in $BUILDS;
do
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

  grep "FAILED" $logfile | grep -v "FAILED: /" | grep -v "The following tests FAILED" | grep -v "FAILED: cd /" | tee -a $summary

  grep "^ok$" $logfile | tee -a $summary
  grep "tests passed" $logfile | tee -a $summary
done

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

