from pycftboot.bootstrap import *

max_pole = 30
max_der = 25

g_tab1 = ConformalBlockTable(1.001, max_pole, 1, max_der, 0)
f_tab1s = ConvolvedBlockTable(g_tab1, symmetric = True)
f_tab1a = ConvolvedBlockTable(g_tab1)

N = 6.0
table_list = [f_tab1s, f_tab1a]
v1 = [[0, 1], [1, 1], [1,0]]
v2 = [[1, 1], [1.0 - (2.0 / N), 1], [-(1.0 + (2.0 / N)), 0]]

info = [[v1, 0, "S"], [v2, 0, "T"]]

sdp = SDP(.3, table_list, vector_types = info)
sdp.set_bound([0, "T"], 0)
result = sdp.bisect(1, 3, 0.01, [0, "S"])
