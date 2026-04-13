import postgkyl as pg
from postgkyl.commands import load_gk_distf
import matplotlib.pyplot as plt

distf = pg.commands.load_gk_distf(
                name="gk_lorentzian_mirror",
                species="ion",
                frame=0)
pg.data.select(distf, z0=0.0, overwrite=True)
pg.output.plot(distf)
plt.show()