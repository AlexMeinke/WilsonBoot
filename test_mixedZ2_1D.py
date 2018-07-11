from pycftboot.bootstrap import *
import numpy as np
import matplotlib.pyplot as plt


cutoff = 0

max_pole = 20
max_der = 7

del_phi = 0.1
del_F = .8

g_tab1 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0)
g_tab2 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del_F - del_phi, del_phi - del_F)
g_tab3 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0, del_phi - del_F, del_phi - del_F)


f_tab1a = ConvolvedBlockTable(g_tab1)
f_tab1s = ConvolvedBlockTable(g_tab1, symmetric = True)
f_tab2a = ConvolvedBlockTable(g_tab2)
f_tab2s = ConvolvedBlockTable(g_tab2, symmetric = True)
f_tab3a = ConvolvedBlockTable(g_tab3)


dim_list = [del_phi, del_F]
tab_list = [f_tab1a, f_tab1s, f_tab2a, f_tab2s, f_tab3a]

m1 = [ [ [1,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
m4 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [1,0,1,1] ] ]
m7 = [ [ [0,0,0,0], [0,0,0,0] ], [ [0,0,0,0], [0,0,0,0] ] ]
m8 = [ [ [0,0,0,0], [0.5,0,0,1] ], [ [0.5,0,0,1], [0,0,0,0] ] ]
m9 = [ [ [0,1,0,0], [0.5,1,0,1] ], [ [0.5,1,0,1], [0,1,0,0] ] ]

# Singlets
v1 = [m1, m4, m7, m8, m9]

v2 = [ [0,0,0,0], [0,0,0,0], [1,4,1,0], [1,2,0,0], [-1,3,0,0] ]
#v2 = [ [0,0,0,0], [0,0,0,0] ]

info = [ [v1, 0, "E"], [v2, 0, "O"] ]


sdp = SDP(dim_list, tab_list, vector_types = info)

sdp.set_bound([0, "O"], del_F)
sdp.set_bound([0, "E"], del_phi)
sdp.set_option("dualErrorThreshold", 1e-15)

result = sdp.iterate()
