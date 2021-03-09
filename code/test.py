    # incident beam
    wfr = v.w_res
    srwlpy.SetRepresElecField(wfr, 'f')
    plot_wfr_diagnostic(wfr, label='input', dir_plot=dir_plot, i=1)
    
    print('Propagating through HHLM1: ', end='')
    t0 = time()
    bl1 = set_optics_HHLM1(v)
    srwlpy.PropagElecField(wfr, bl1)
    print('done in', round(time() - t0, 3), 's')
    plot_wfr_diagnostic(wfr, label='after HHLM1', dir_plot=dir_plot, i=2)
    
    print('Propagating through HHLM2: ', end='')
    t0 = time()
    bl2 = set_optics_HHLM2(v, drift=drift_list[0])
    srwlpy.PropagElecField(wfr, bl2)
    print('done in', round(time() - t0, 3), 's')
    plot_wfr_diagnostic(wfr, label='after HHLM2', dir_plot=dir_plot, i=3)
    
    print('Propagating through HHLM3: ', end='')
    t0 = time()
    bl3 = set_optics_HHLM3(v, drift=drift_list[1])
    srwlpy.PropagElecField(wfr, bl3)
    print('done in', round(time() - t0, 3), 's')
    plot_wfr_diagnostic(wfr, label='after HHLM3', dir_plot=dir_plot, i=4)
    
    print('Propagating through HHLM4: ', end='')
    t0 = time()
    bl4 = set_optics_HHLM4(v, drift=drift_list[2])
    srwlpy.PropagElecField(wfr, bl4)
    print('done in', round(time() - t0, 3), 's')
    plot_wfr_diagnostic(wfr, label='after HHLM4', dir_plot=dir_plot, i=5)
    





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