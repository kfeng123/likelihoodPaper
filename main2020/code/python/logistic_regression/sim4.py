from simulation import logistic_simulator

if __name__ == "__main__" :
    import time

    now_time = time.time()
    simulator = logistic_simulator(n =100, p=3, p0=2, true_nu=[1,0], true_xi = [10] * 1)
    res = simulator.do_it()
    with open("n_100_nu_10_xi_10", "w") as f:
        f.write("LRT:")
        f.write(str(res["LRT"]))
        f.write("\n")
        f.write("gFBF:")
        f.write(str(res["gFBF"]))
    print("time", time.time()-now_time)

    now_time = time.time()
    simulator = logistic_simulator(n =100, p=3, p0=2, true_nu=[1,0], true_xi = [20] * 1)
    res = simulator.do_it()
    with open("n_100_nu_10_xi_20", "w") as f:
        f.write("LRT:")
        f.write(str(res["LRT"]))
        f.write("\n")
        f.write("gFBF:")
        f.write(str(res["gFBF"]))
    print("time", time.time()-now_time)
