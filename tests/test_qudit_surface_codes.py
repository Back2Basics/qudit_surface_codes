from qudit_surface_codes import *


sigma = ((0,1,2),(3,4,5),(6,7))
alpha = ((0,3),(1,6),(2,4),(5,7))

def test_init():
    SCG = HybridSurfaceCodeGraph(sigma, alpha)

def test_tuple_of_tuples():
    SCG = HybridSurfaceCodeGraph(sigma, alpha)
    isinstance(SCG.phi, tuple)

def test_isdictionary():
    SCG = HybridSurfaceCodeGraph(sigma, alpha)
    isinstance(SCG.node_info, dict)
