#!/bin/bash

rm -f tmp.gp

echo 'set term png' >>tmp.gp
#echo 'set logscale y' >>tmp.gp
echo 'set output "out.png"' >>tmp.gp
echo "plot \"$1\" using 1:2 title '$1' with linespoints, \"$2\" using 1:2 title '$2' with linespoints" >>tmp.gp

gnuplot tmp.gp
rm -f tmp.gp

