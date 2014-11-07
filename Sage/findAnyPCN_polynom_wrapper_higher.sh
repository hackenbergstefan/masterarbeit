#!/bin/sh
N=$1
BORDER=$(echo "$1^4" | bc)
echo $N $BORDER
curR=2
if [ "$#" -ge 2 ]; then
  curR=$2
fi

while [ $(echo "2^$curR" | bc) -le $BORDER ]; do
  for (( c=1; c<=$BORDER; c+=100000 )); do
    d=$c+100000
    echo "curR = "$curR
    sage <<EOF
load("./findAnyPCN_trinom.spyx")
findAnyPCN_polynom_wrapper($N, fileoutput=True, startPrime=$c, stopPrime=$d, onlyR=$curR)
exit()
EOF
  done
  curR=$(( $curR + 1 ))
done

