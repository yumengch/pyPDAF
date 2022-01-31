export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:lib
which python
mpiexec -n 9 python -u example/main.py