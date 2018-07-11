from pycftboot.bootstrap import *



del_phi = .2
del_F = 1.4

cutoff = 0

max_pole = 40
max_der = 15


g_tab1 = ConformalBlockTable(1.01, max_pole, 1, max_der, 0)
g_tab2 = ConformalBlockTable(1.01, max_pole, 1, max_der, 0, del_F - del_phi, del_phi - del_F)
g_tab3 = ConformalBlockTable(1.01, max_pole, 1, max_der, 0, del_phi - del_F, del_phi - del_F)


f_tab1a = ConvolvedBlockTable(g_tab1)
f_tab1s = ConvolvedBlockTable(g_tab1, symmetric = True)
f_tab2a = ConvolvedBlockTable(g_tab2)
f_tab2s = ConvolvedBlockTable(g_tab2, symmetric = True)
f_tab3a = ConvolvedBlockTable(g_tab3)


dim_list = [del_phi, del_F]
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

# Singlets
v1 = [m1, m2, m3, m4, m5, m6, m7, m8, m9]

# O(6) symmetric traceless
v2 = [ [1,0,0,0], [(1-2/6.0),0,0,0], [-(1+2/6.0),1,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0] ]


# O(3) symmetric traceless
v4 = [ [0,0,0,0], [0,0,0,0], [0,0,0,0], [1,0,1,1], [(1-2/3.0),0,1,1], [-(1+2/3.0),1,1,1], [0,0,0,0], [0,0,0,0], [0,0,0,0] ]


# O(3)xO(6) vector
v6 = [ [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [1,4,1,0], [-1,2,0,0], [1,3,0,0] ]

info = [[v1,0,"S"], [v2,0,"T1"], [v4,0,"T2"], [v6,0,"VV"]]
sdp = SDP(dim_list, tab_list, vector_types = info)

bound = del_F + .01
sdp.add_point([0,"T1"], del_phi)
sdp.add_point([0,"T2"], del_F)

sdp.set_bound([0,"S"], 1)
sdp.set_bound([0,"T1"], del_phi)
sdp.set_bound([0,"T2"], del_F)
sdp.set_bound([0,"VV"], 1)
sdp.set_option("dualErrorThreshold", 1e-20)


result = sdp.iterate() 

