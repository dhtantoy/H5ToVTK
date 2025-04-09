using H5ToVTK

filename = "colin27_coarse_boundaries"
grps = ["mesh", "boundaries"]
h5tovtk(filename, grps)