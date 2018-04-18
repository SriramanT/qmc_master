# huy-nguyen's description

This pair of papers by Prof. Zhang and others pioneered the constrained path Monte Carlo method. In them, he and the co-authors combined techniques from DMC (open-ended branching random walk in imaginary time, importance sampling) and AFQMC (Hubbard-Stratonovich transformation, Slater determinants) to do ground state projection.

However, the real deal here is the insight that by constraining the random walkers to always maintain a positive overlap with the trial wave function, they can remove the exponential signal-to-noise decay characteristic of the “fermion sign problem.” This insight was first observed by Fahy (see above) but the traditional AFQMC method makes the implementation of this constraint computationally difficult. In contrast, this constraint is very easily implemented in a branching random walk!

I’ve read the PRB article about 30 times and each time I find new insights that I haven’t noticed before. It really showed me how much work went into developing this method. Take it slow!

---

This paper is a nice addition to the pair of papers above by Zhang et al. It very clearly illuminates the similarities and differences between first- and second-quantized QMC algorithms (in this case DMC vs CPMC). Having read papers from DMC and CPMC, this solidifies my understanding of the two methods.

Fun trivia: Prof Zhang introduced me to Carlson when I visited W&M in Oct 2013. Of course, being an undergraduate, I didn’t have much to say to him but it was still fun to meet one of the pioneers of the method and to be able to associate a face to name on paper.

---

This paper by Zhang is a follow-up to the pioneering pair above (which describes ground-state CPMC). It extends the CPMC method to finite-temperature. It turns out the extension is fairly straightforward.


In this paper Zhang continues developing the CPMC method by extending it to treat interactions that couple auxiliary fields to complex numbers. This extension, known as phase-free or phaseless AFQMC, allows us to deal with any kind of interactions.§
