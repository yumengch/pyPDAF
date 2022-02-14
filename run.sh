export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:lib
mpiexec -n 8 python -u example/main.py