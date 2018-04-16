"""
Draw probabilistic graphical models for the luminosity inference tutorial.

Anthony G.A. Brown Nov 2017 - Nov 2017
<brown@strw.leidenuniv.nl>
"""

from matplotlib import rc
rc("font", family="serif", size=12)
rc("text", usetex=True)

import daft

pgm = daft.PGM([5,5], origin=[0, 0], node_unit=1.5, grid_unit=2.5)

pgm.add_node(daft.Node("obsm", r"$m_k$", 3, 2, observed=True))
pgm.add_node(daft.Node("obsplx", r"$\varpi_k$", 3, 1, observed=True))
pgm.add_node(daft.Node("obsm_err", r"$\sigma_{m,k}$", 4, 2, fixed=True))
pgm.add_node(daft.Node("obsplx_err", r"$\sigma_{\varpi,k}$", 4, 1, fixed=True))
pgm.add_node(daft.Node("r", r"$r_k$", 2, 1))
pgm.add_node(daft.Node("trueAbsMag", r"$M_k$", 3, 3))
pgm.add_node(daft.Node("rmin", r"$r_\mathrm{L}$", 0.5, 1, fixed=True))
pgm.add_node(daft.Node("rmax", r"$r_\mathrm{H}$", 0.5, 1.7, fixed=True))
pgm.add_node(daft.Node("alpha", r"$\alpha$", 0.5, 0.3, fixed=True))
pgm.add_node(daft.Node("meanAbsMag", r"$\mu_M$", 2, 4.5))
pgm.add_node(daft.Node("sigmaAbsMag", r"$\sigma_M$", 3, 4.5))
pgm.add_node(daft.Node("mlim", r"$m_\mathrm{lim}$", 4.5, 2.7, fixed=True))

pgm.add_edge("obsm_err", "obsm")
pgm.add_edge("obsplx_err", "obsplx")
pgm.add_edge("r", "obsm")
pgm.add_edge("r", "obsplx")
pgm.add_edge("trueAbsMag", "obsm")
pgm.add_edge("rmin", "r")
pgm.add_edge("rmax", "r")
pgm.add_edge("alpha", "r")
pgm.add_edge("meanAbsMag", "trueAbsMag")
pgm.add_edge("sigmaAbsMag", "trueAbsMag")
pgm.add_edge("mlim", "obsm")

pgm.add_plate(daft.Plate([1.2, 0.5, 3, 3], label=r"$k = 1, \ldots, N$", label_offset=[10,150], shift=0.0,
    bbox={'fill':False, 'lw':0}))

pgm.render()
pgm.figure.savefig("pgm_luminosity_inference_distprior.png", dpi=150)
