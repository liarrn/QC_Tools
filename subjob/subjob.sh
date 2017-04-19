#!/bin/sh

for i in {011..032}
do
        cd disp-$i
        yhbatch -N 6 -n 144 -p nsfc3 -J ${i} RUN_VASP_5.4.1.sh
        cd ..
done
