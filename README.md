# MuonTridents
This repository has the code needed to calculate muon trident events in coherent and incoherent interactions. See arxiv 2510.18943 and 2610.33153.

We implemented the full squared matrix element generated via FaynCalc package taking into account Bethe-Heitler ans Bremsstrahlung channels. Then we sample the phase space, spanned by seven invariant masses, using rejection sampling. 

We can calculate total cross sections or generate events. Events generated are unweighted.
