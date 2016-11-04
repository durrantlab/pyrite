import MDAnalysis
u = MDAnalysis.Universe("test.psf", "test.dcd")
print(u)
u.atoms.write("MDsims.pdb")