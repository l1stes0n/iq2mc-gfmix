* MCMCtree control template; the workflow fills the @...@ placeholders.
* Edit this file or pass a custom one with --mcmctree-ctl FILE.

seed = -1
seqfile = @SEQFILE@
treefile = @TREEFILE@
mcmcfile = @PREFIX@.mcmctree.log
outfile = @PREFIX@.mcmctree.out
ckpfile = @PREFIX@.mcmctree.ckp
hessianfile = @HESSIANFILE@        * supported by the iqtree/paml MCMCtree

checkpoint = 1          * 0: nothing; 1: save; 2: resume
ndata = 1
seqtype = 0             * 0: nucleotides; 1: codons; 2: AAs
usedata = 2             * 2: approximate likelihood from hessianfile
clock = 2               * 1: global clock; 2: independent rates; 3: correlated rates
RootAge = <1.0          * used only if the tree has no root calibration

BDparas = 1 1 0.5       * birth rate, death rate, sampling fraction
rgene_gamma = 2 40      * mean rate prior (alpha, beta)
sigma2_gamma = 1 10     * rate-variance prior (alpha, beta)
finetune = 1: 0.1 0.1 0.1 0.01 .5   * auto (0/1): times, musigma2, rates, mixing, paras, FossilErr
print = 1               * 1: normal output; 2: verbose

burnin = 4000
sampfreq = 10
nsample = 5000

* Only needed for exact likelihood (usedata = 1); approximate likelihood ignores these.
model = 0
alpha = 0
ncatG = 0
cleandata = 0
kappa_gamma = 6 2
alpha_gamma = 1 1
