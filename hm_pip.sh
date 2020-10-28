#!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
        continue
fi
echo $line
cd $line"_Output"
../deal2_new20190723.py --prefix ${line}"_hm"
cd ..
done<readme
