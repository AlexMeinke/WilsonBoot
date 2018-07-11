from pycftboot.bootstrap import *
cutoff = 0

sig = 0.518
eps = 1.412

max_spin = 15
max_pole = 15

g_tab1 = ConformalBlockTable(3, max_pole, max_spin, 2, 4)
g_tab2 = ConformalBlockTable(3, max_pole, max_spin, 2, 4, eps - sig, sig - eps, odd_spins = True)
g_tab3 = ConformalBlockTable(3, max_pole, max_spin, 2, 4, sig - eps, sig - eps, odd_spins = True)

f_tab1a = ConvolvedBlockTable(g_tab1)
f_tab1s = ConvolvedBlockTable(g_tab1, symmetric=True)
f_tab2a = ConvolvedBlockTable(g_tab2)
f_tab2s = ConvolvedBlockTable(g_tab2, symmetric=True)
f_tab3a = ConvolvedBlockTable(g_tab3)

dim_list = [sig, eps]
tab_list = [f_tab1a, f_tab1s, f_tab2a, f_tab2s, f_tab3a]

m1 = [ [[1,0,0,0], [0,0,0,0]], [[0,0,0,0],[0,0,0,0]] ]
m2 = [ [[0,0,0,0], [0,0,0,0]], [[0,0,0,0],[1,0,1,1]] ]
m3 = [ [[0,0,0,0], [0,0,0,0]], [[0,0,0,0],[0,0,0,0]] ]
m4 = [ [[0,0,0,0], [0.5,0,0,1]], [[0.5,0,0,1],[0,0,0,0]] ]
m5 = [ [[0,1,0,0], [0.5,1,0,1]], [[0.5,1,0,1],[0,1,0,0]] ]

v1 = [m1, m2, m3, m4, m5]

v2 = [ [0,0,0,0], [0,0,0,0], [1,4,1,0], [1,2,0,0], [-1,3,0,0] ]
v3 = [ [0,0,0,0], [0,0,0,0], [1,4,1,0], [-1,2,0,0], [1,3,0,0] ]

info = [[v1, 0, 0], [v2, 0, 1], [v3, 1, 2]]


sdp = SDP(dim_list, tab_list, vector_types = info)
sdp.set_bound([0, 0], dim_list[1])
sdp.set_bound([0, 1], 3.0)
sdp.add_point([0, 1], dim_list[0])
sdp.set_option("dualErrorThreshold", 1e-15)
result = sdp.iterate()

