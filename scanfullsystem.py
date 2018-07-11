from pycftboot.bootstrap import *
import numpy as np
import matplotlib.pyplot as plt

max_pole = 25
max_der = 11


def test_spec(del1, del2):
    g_tab1 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0)
    g_tab2 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del2 - del1, del1 - del2)
    g_tab3 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del1 - del2, del1 - del2)

    f_tab1a = ConvolvedBlockTable(g_tab1)
    f_tab1s = ConvolvedBlockTable(g_tab1, symmetric = True)
    f_tab2a = ConvolvedBlockTable(g_tab2)
    f_tab2s = ConvolvedBlockTable(g_tab2, symmetric = True)
    f_tab3a = ConvolvedBlockTable(g_tab3)

    dim_list = [del1, del2]
    tab_list = [f_tab1a, f_tab1s, f_tab2a, f_tab2s, f_tab3a]

    m1 = [ [ [1,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m4 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [1,0,1,1] ] ]
    m7 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m8 = [ [ [0,0,0,0], [0.5,0,0,1] ], [ [0.5,0,0,1], [0,0,0,0] ] ]
    m9 = [ [ [0,1,0,0], [0.5,1,0,1] ], [ [0.5,1,0,1], [0,1,0,0] ] ]

    v1 = [m1, m4, m7, m8, m9]
    v2 = [ [0,0,0,0], [0,0,0,0], [1,4,1,0], [1,2,0,0], [-1,3,0,0] ]

    info = [ [v1, 0, "E"], [v2, 0, "O"] ]

    sdp = SDP(dim_list, tab_list, vector_types = info)
    sdp.set_bound([0, "E"], del1)
    sdp.set_bound([0, "O"], del2)
    sdp.set_option("dualErrorThreshold", 1e-15)

    return  sdp.iterate()

def test_spec_full(del1, del2):
    g_tab1 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0)
    g_tab2 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del2 - del1, del1 - del2)
    g_tab3 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del1 - del2, del1 - del2)

    f_tab1a = ConvolvedBlockTable(g_tab1)
    f_tab1s = ConvolvedBlockTable(g_tab1, symmetric = True)
    f_tab2a = ConvolvedBlockTable(g_tab2)
    f_tab2s = ConvolvedBlockTable(g_tab2, symmetric = True)
    f_tab3a = ConvolvedBlockTable(g_tab3)

    dim_list = [del1, del2]
    tab_list = [f_tab1a, f_tab1s, f_tab2a, f_tab2s, f_tab3a]

    m1 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m2 = [ [ [1,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m3 = [ [ [1,1,0,0], [0,1,0,0] ], [ [0,1,0,0], [0,1,0,0] ] ]
    m4 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m5 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [1,0,1,1] ] ]
    m6 = [ [ [0,1,0,0], [0,1,0,0] ], [ [0,1,0,0], [1,1,1,1] ] ]
    m7 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
    m8 = [ [ [0,0,0,0], [0.5,0,0,1] ], [ [0.5,0,0,1], [0,0,0,0] ] ]
    m9 = [ [ [0,1,0,0], [0.5,1,0,1] ], [ [0.5,1,0,1], [0,1,0,0] ] ]

    v1 = [m1, m2, m3, m4, m5, m6, m7, m8, m9]
    v2 = [ [1,0,0,0], [(1-2/6.0),0,0,0], [-(1+2/6.0),1,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0] ]
    v4 = [ [0,0,0,0], [0,0,0,0], [0,0,0,0], [1,0,1,1], [(1-2/3.0),0,1,1], [-(1+2/3.0),1,1,1], [0,0,0,0], [0,0,0,0], [0,0,0,0] ]
    v6 = [ [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [1,4,1,0], [-1,2,0,0], [1,3,0,0] ]

    info = [[v1,0,"S"], [v2,0,"T1"], [v4,0,"T2"], [v6,0,"VV"]]
    sdp = SDP(dim_list, tab_list, vector_types = info)

    sdp.add_point([0,"T1"], del1)
    sdp.add_point([0,"T2"], del2)

    sdp.set_bound([0,"S"], 1)
    sdp.set_bound([0,"T1"], del1)
    sdp.set_bound([0,"T2"], 1)
    sdp.set_bound([0,"VV"], del2)
    sdp.set_option("dualErrorThreshold", 1e-15)
    return  sdp.iterate()



def bisect(del1, lower, upper, threshold,full=False):
    x = .5
    while abs(upper - lower) > threshold:
        test = lower + x * (upper - lower)
        print("Trying del_F=" + str(test))
        if full:
            result = test_spec_full(del1, test)
        else:
            result = test_spec(del1, test)

        if result==False:
            upper = test
        else:
            lower = test
    return lower

del_phi_vec = np.linspace(0., 2.4, 13)
del_phi_vec[0] = 0.05
del_F_vec = 0*del_phi_vec
del_F_vec2 = 0*del_phi_vec

for i in range(len(del_phi_vec)):
    del_F_vec[i] = bisect(del_phi_vec[i], 0, 4., .01)
    del_F_vec2[i] = bisect(del_phi_vec[i], 0, 4., .01, full=True)

plt.plot(del_phi_vec, del_F_vec)
plt.plot(del_phi_vec, del_F_vec2)
plt.xlabel("Delta_1")
plt.ylabel("Delta_2")
plt.show()
