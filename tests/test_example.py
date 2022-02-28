import numpy as np

def test_ens_forecast():
    f_out_name = 'true_outputs_online/ens_{i_ens:02}_step{step:02}_for.txt'
    py_out_name = 'ens_{i_ens}_step{step}_for.txt'
    nt = 18
    n_ens = 4
    for step in range(1, nt, 2):
        for i in range(n_ens):
            f_out = np.loadtxt(f_out_name.format(i_ens=i+1, step=step+1))
            py_out = np.loadtxt(py_out_name.format(i_ens=i+1, step=step+1), delimiter=';')
            assert np.allclose(f_out, py_out), 'Python output and Fortran ensemble forecast differs!'

def test_ens_analysis():
    f_out_name = 'true_outputs_online/ens_{i_ens:02}_step{step:02}_ana.txt'
    py_out_name = 'ens_{i_ens}_step{step}_ana.txt'
    nt = 18
    n_ens = 4
    for step in range(1, nt, 2):
        for i in range(n_ens):
            f_out = np.loadtxt(f_out_name.format(i_ens=i+1, step=step+1))
            py_out = np.loadtxt(py_out_name.format(i_ens=i+1, step=step+1), delimiter=';')
            assert np.allclose(f_out, py_out), 'Python output and Fortran ensemble analysis differs!'

def test_state_analysis():
    f_out_name = 'true_outputs_online/state_step{step:02}_ana.txt'
    py_out_name = 'state_step{step}_ana.txt'
    nt = 18
    n_ens = 4
    for step in range(1, nt, 2):
        for i in range(n_ens):
            f_out = np.loadtxt(f_out_name.format(i_ens=i+1, step=step+1))
            py_out = np.loadtxt(py_out_name.format(i_ens=i+1, step=step+1), delimiter=';')
            assert np.allclose(f_out, py_out), 'Python output and Fortran ensemble analysis mean differs!'

def test_state_forecast():
    f_out_name = 'true_outputs_online/state_step{step:02}_for.txt'
    py_out_name = 'state_step{step}_for.txt'
    nt = 18
    n_ens = 4
    for step in range(1, nt, 2):
        for i in range(n_ens):
            f_out = np.loadtxt(f_out_name.format(i_ens=i+1, step=step+1))
            py_out = np.loadtxt(py_out_name.format(i_ens=i+1, step=step+1), delimiter=';')
            assert np.allclose(f_out, py_out), 'Python output and Fortran ensemble forecast mean differs!'