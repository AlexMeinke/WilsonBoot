from pycftboot.bootstrap import *

max_pole = 30
max_der = 25

g_tab = ConformalBlockTable(1.001, max_pole, 1, max_der, 0)
f_tab = ConvolvedBlockTable(g_tab)

sdp = SDP(1., f_tab)
#result = sdp.bisect(1, 3, 0.01, 0)
result1 = sdp.opemax(2., 0)
