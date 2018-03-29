#!/bin/bash

# Automated script to compile&test ASPECT in different configurations

cd /drone/src/github.com/gassmoeller/aspect

submit="OFF"
mkdir -p /drone/log
summary=/drone/log/summary
indexhtml=/drone/log/index.html

main()
{
#clean contents:
> $summary

logfile=drone/log/log-astyle
LOGFILE="`pwd`/changes.diff"
./doc/indent || { echo "indent FAILED"; return 1; }
git diff >$LOGFILE
echo "git diff >$LOGFILE"
git diff --exit-code --name-only || { echo "FAILED:"; echo `git diff`; return 1; }
echo "ok"

if [ -s changes.diff ]; then
  cp changes.diff /drone/log/changes-astyle.diff
  echo "DIFFS: changes-astyle.diff" | tee -a $logfile
fi
cd ..
rep "FAILED" $logfile | grep -v "FAILED: /" | grep -v "The following tests FAILED" | grep -v "FAILED: cd /" | tee -a $summary

grep "^ok$" $logfile | tee -a $summary
grep "tests passed" $logfile | tee -a $summary

cp $summary $indexhtml

grep -h "DIFFS: changes-" /drone/log/log-* | tee -a $indexhtml
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

