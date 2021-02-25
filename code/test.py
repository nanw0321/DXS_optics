def main():
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    v.si = True
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)

    wfr = v.w_res       # initialize wavefront

    # propagation
    srwlpy.SetRepresElecField(wfr, 'f') # switch to frequency domain
    bl = set_optics(v)                  # initialize beamline
    srwlpy.PropagateElecField(wfr, bl)  # propagate through beamline

    # resizing
    srwlpy.ResizeElecField(wfr, 'f', [0, fRange, fRes]) # in frequency
    srwlpy.ResizeElecField(wfr, 't', [0, tRange, tRes]) # in time

    





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


### Try a single crystal reflection and see if there's tail lol
### Can have a weekly meeting about XPCS modelling.