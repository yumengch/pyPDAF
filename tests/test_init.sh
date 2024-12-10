#!/bin/bash
# testing the array size in Python
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 100;
do
    for j in 0 1 2 3 4 5;
    do
        if [ $i -eq 0 ] && [ $j -eq 4 ]; then
            continue
        fi
        if [ $i -eq 2 ] && ([ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 3 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 4 ] && [ $j -eq 4 ]; then
            continue
        fi
        if [ $i -eq 5 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 6 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 7 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 8 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 9 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 10 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 11 ] && ([ $j -eq 2 ] || [ $j -eq 3 ]); then
            continue
        fi
        if [ $i -eq 12 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 100 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ] || [ $j -eq 5 ]); then
            continue
        fi
        /home/ia923171/miniconda3/envs/pypdaf-devel/bin/python \
        test_init.py $i $j
    done
done

for i in 0 1 4 5 6 7;
do
    /home/ia923171/miniconda3/envs/pypdaf-devel/bin/python \
    test_init.py 200 $i
done

# testing the array size in Fortran
for i in 1 3 6 7;
do
    for j in 0 1 2 3 4 5;
    do
        if [ $i -eq 3 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 6 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 7 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        ./test_init $i $j
    done
done

for i in 1 4 6 7;
do
    ./test_init  200 $i
done

for i in 0 4 5 9 10 11;
do
    for j in 0 1 2 3 4 5;
    do
        if [ $i -eq 0 ] && [ $j -eq 4 ]; then
            continue
        fi
        if [ $i -eq 4 ] && [ $j -eq 4 ]; then
            continue
        fi
        if [ $i -eq 5 ] && ([ $j -eq 1 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 9 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 10 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 11 ] && ([ $j -eq 2 ] || [ $j -eq 3 ]); then
            continue
        fi
        ./test_init $i $j
    done
done

for i in 2  8 12 100 200;
do
    for j in 0 1 2 3 4 5;
    do
        if [ $i -eq 2 ] && ([ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 8 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 12 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ]); then
            continue
        fi
        if [ $i -eq 100 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ] || [ $j -eq 5 ]); then
            continue
        fi
        if [ $i -eq 200 ] && ([ $j -eq 1 ] || [ $j -eq 2 ] || [ $j -eq 3 ] || [ $j -eq 4 ] || [ $j -eq 5 ]); then
            continue
        fi
        ./test_init $i $j
    done
done
