##copy/paste into terminal/python

source /Users/nlaniewski/Desktop/Pyvenv_1/bin/activate

python

from fastTSNE import TSNE
from fastTSNE.callbacks import ErrorLogger
import numpy as np

datapath_x =  '/Volumes/nlaniewski/RPRC.PRISM.Flow/results_R/fSOMs/TPHE/input.forTSNE.TPHE_CD4p_PASS3.csv'
data_x = np.loadtxt(datapath_x, delimiter = ',')


embeddings = []

tsne = TSNE(
  # Let's use the fast approximation methods
  neighbors='approx', negative_gradient_method='fft', initialization='random', perplexity=50,
  n_iter=2000,
  # The embedding will be appended to the list we defined above, make sure we copy the
  # embedding, otherwise the same object reference will be stored for every iteration
  callbacks=lambda it, err, emb: embeddings.append(np.array(emb)),
  # This should be done on every iteration
  callbacks_every_iters=1,
  # -2 will use all but one core so I can look at cute cat pictures while this computes
  n_jobs=10
)

_ = tsne.fit(data_x)


np.savetxt("/Volumes/nlaniewski/RPRC.PRISM.Flow/results_R/fSOMs/TPHE/embeddings.TSNE.TPHE_CD4p_PASS3.csv", embeddings[len(embeddings)-1], delimiter = ',')

use_python("/Users/nlaniewski/Desktop/Pyvenv_1/", required = TRUE)

np <- import("numpy")
fastTSNE <- import("fastTSNE")

data_x = as.matrix(counts.merged[grep("Cluster", colnames(counts.merged))])

tsne = fastTSNE$TSNE(
  # Let's use the fast approximation methods
  neighbors='approx', negative_gradient_method='bh', initialization='random', perplexity= as.integer(50),
  n_iter=as.integer(1000),
  n_jobs=as.integer(10)
)

tmp = tsne$fit(data_x.np)

tmp.plot <- cbind(counts.merged, tmp)

ggplot(tmp.plot, aes(`1`, `2`, color = term.visit)) + 
  geom_point(size = 3) + 
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

use_python("/Users/nlaniewski/Desktop/Pyvenv_1/", required = TRUE)
openTSNE <- import("openTSNE")
sklearn <- import("sklearn.datasets")

iris = sklearn$
x, y = iris["data"], iris["target"]

which python
