#!/bin/bash

# Automated script to compile&test ASPECT in different configurations

source_folder=$1

cd $source_folder

log_folder=$source_folder/log
mkdir -p $log_folder

submit="OFF"

logfile=${log_folder}/log
touch $logfile
summary=${log_folder}/summary
touch $summary
diff=${log_folder}/diff
touch $diff

indexhtml=${log_folder}/index.html

main()
{
#clean contents:
> $summary

./doc/indent || { echo "indent FAILED"; return 1; }
echo "git diff >$diff"
git diff >$diff

git diff --exit-code --name-only || { echo "FAILED:"; git diff; return 1; }
echo "ok"

if [ -s changes.diff ]; then
  echo "DIFFS: changes-astyle.diff" | tee -a $logfile
fi
cd ..

cp $summary $indexhtml
}

if [ "$QUIET" == "1" ];
then
  main >/dev/null
  cat $summary
else
  main
fi

