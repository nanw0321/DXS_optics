    # CRL_C2: drift
    ['op_C1_C2_L', 'f', 2.6, 'length'],

    # C2: crystal
    ['op_C2_hfn', 's', '', 'heightProfileFile'],
    ['op_C2_dim', 's', 'x', 'orientation'],
    ['op_C2_d_sp', 'f', 0.7838928390938714, 'dSpacing'],
    ['op_C2_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_C2_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_C2_psiHr', 'f', -3.244798554456362e-06, 'psiHr'],
    ['op_C2_psiHi', 'f', 1.5044880611361188e-07, 'psiHi'],
    ['op_C2_psiHBr', 'f', -3.244798554456362e-06, 'psiHBr'],
    ['op_C2_psiHBi', 'f', 1.5044880611361188e-07, 'psiHBi'],
    ['op_C2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_C2_uc', 'f', 1, 'useCase'],
    ['op_C2_ang_as', 'f', 0.08726646259971647, 'asymmetryAngle'],
    ['op_C2_nvx', 'f', 0.622180932039, 'nvx'],
    ['op_C2_nvy', 'f', 4.228e-09, 'nvy'],
    ['op_C2_nvz', 'f', -0.782873481354, 'nvz'],
    ['op_C2_tvx', 'f', 0.782873481354, 'tvx'],
    ['op_C2_tvy', 'f', 5.32e-09, 'tvy'],
    ['op_C2_ang', 'f', 0.0, 'grazingAngle'],
    ['op_C2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C2_energy', 'f', 9481.0, 'energy'],
    ['op_C2_diffractionAngle', 'f', -1.57079632, 'diffractionAngle'],





### HE meeting: pushed by DOE review. The next review probably at the end of this year
### but maybe it makes sense to quickly touch base with the review committee before that and show them the progress report
### Or put the abstracts to SRI meetings, due early March



### To-do:
### Then Gaussian but with an expected SASE bandwidth of ~30eV
### Sending a SASE pulse through
### thermal deformation needs to be implemented (beyond month-time scale)
### Jacek's code and Matt's code can be cross-referenced and give insights
### benchmark Matt's code with perfect optics and see our limits

### don't know how stable LCLS-II is yet
### but it's new, and the injector works much closer to a steady-state operation
### so should be more stable than LCLS-I, but will never reach the stability of a synchrotron

### Won't have LCLS-II data by the design review (end of year), so need a way to simulate beam fluctuations.

### Can have a weekly meeting about XPCS modelling.



### Next, high heat load mono simulation; the parameters: https://docs.google.com/spreadsheets/d/1dCrDyQOWnNd47cDOBF22fLEGHUD4JVOU/edit#gid=1087546977
### Note: the distances are in lab frame, need to translate to path-lengths in SRW by multiplying a bunch of cosines