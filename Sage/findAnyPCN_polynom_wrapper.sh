#!/bin/sh
N=$1
R=$2
BORDER=$(echo "$1^4" | bc)
if [ "$#" -eq 3 ]; then
  STARTPRIME=$3
else
  STARTPRIME=1
fi
echo $N $BORDER
for (( c=$STARTPRIME; c<=$BORDER; c+=100000 ))
do
  d=$c+100000
  sage <<EOF
load("./findAnyPCN_trinom.spyx")
findAnyPCN_polynom_wrapper($N, fileoutput=True, startPrime=$c, stopPrime=$d, onlyR=$R)
exit()
EOF
done

