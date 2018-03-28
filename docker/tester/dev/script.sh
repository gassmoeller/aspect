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



if [ -z "$BUILDS" ]; then
  echo 'Please specify list of builds to do in the ENV variable $BUILDS.'
  echo "Separate build with spaces. Valid options: gcc gccpetsc clang clangpetsc"
  exit 0
fi

mkdir -p ~/log

if [ -s ~/source/CMakeLists.txt ]; then
  touch ~/source/VERSION 2>/dev/null || { echo "~/source needs to be mounted R/W, aborting."; exit 1; }
else
  if [ -s /source/CMakeLists.txt ]; then
    cp -r /source ~/source
    rm -rf ~/source/CMakeCache.txt
  else
    echo "ERROR, no ASPECT mounted under ~/source/ or /source/"
    exit 1
  fi
fi


submit="OFF"

run()
{
build=$1
desc=$2
submit=$3

if [[ $build == "astyle" ]];
then
  LOGFILE="`pwd`/changes.diff"
  cd ~/source
  ./doc/indent || { echo "indent FAILED"; return; }
  git diff >$LOGFILE
  echo "git diff >$LOGFILE"
  git diff --exit-code --name-only || { echo "FAILED: `git diff --name-only`"; return; }
  echo "ok"
  return
fi

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

#echo "+ cmake -G "Ninja" $compiler -D ASPECT_USE_PETSC=$petsc -D ASPECT_RUN_ALL_TESTS=ON ~/source"
#time cmake -G "Ninja" $compiler -D ASPECT_USE_PETSC=$petsc -D ASPECT_RUN_ALL_TESTS=ON ~/source || { echo "configure FAILED"; return; }

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

