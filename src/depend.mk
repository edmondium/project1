MOD0 = arrays.f
DEP0 = correlation.f correlation_wannier.f dump.f exchange.f getinput.f quasiparticle.f selfenergy.f selfenergy_com.f selfenergy_wannier.f spex.f xc_energy.f
MOD1 = file.f
DEP1 = Mutil.f bandinfo.f bethesalpeter.f checkinput.f checkmem.f continuation.f correlation.f correlation_wannier.f coulombmatrix.f dump.f exchange.f getinput.f hilbert.f iterate.f mixedbasis.f numerics.f quasiparticle.f readwrite_fleur.f readwrite_fleurR3.f readwrite_fleurR5.1.f readwrite_fleurR5.f readwrite_tmpl.f selfenergy.f selfenergy_wannier.f spectra_wannier.f spex.f susceptibility.f susceptibility_core.f susceptibility_wannier.f timer_util.f util.f wannier.f
MOD2 = freqintegral.f
DEP2 = selfenergy.f selfenergy_com.f selfenergy_wannier.f
MOD3 = global.f
DEP3 = Hwrapper.f Mutil.f Mwrapper.f bandinfo.f bethesalpeter.f checkinput.f checkmem.f continuation.f correlation.f correlation_wannier.f coulomb_wannier.f coulombmatrix.f dirac.f divergence.f dump.f dwavefproducts.f exchange.f fft.f freqintegral.f gauss.f getinput.f hilbert.f ibc.f irreps.f iterate.f key.f mixedbasis.f numerics.f overlap.f pade.f quasiparticle.f readwrite_fleur.f readwrite_fleurR3.f readwrite_fleurR5.1.f readwrite_fleurR5.f readwrite_tmpl.f selfenergy.f selfenergy_com.f selfenergy_wannier.f spectra_wannier.f spex.f susceptibility.f susceptibility_core.f susceptibility_wannier.f symmetry.f tetrahedron.f trafo.f util.f vector.f wannier.f wavefproducts.f xc_energy.f
MOD4 = Hwrapper.f
DEP4 = correlation.f correlation_wannier.f iterate.f readwrite_fleur.f readwrite_fleurR3.f readwrite_fleurR5.1.f readwrite_fleurR5.f readwrite_tmpl.f
MOD5 = key.f
DEP5 = bethesalpeter.f checkinput.f correlation.f correlation_wannier.f coulombmatrix.f getinput.f hilbert.f iterate.f mixedbasis.f quasiparticle.f selfenergy.f susceptibility.f susceptibility_core.f susceptibility_wannier.f wannier.f
MOD6 = fft.f
DEP6 = getinput.f spex.f wavefproducts.f
MOD7 = Mwrapper.f
DEP7 = Hwrapper.f Mutil.f bethesalpeter.f correlation.f correlation_wannier.f coulomb_wannier.f coulombmatrix.f divergence.f dwavefproducts.f exchange.f getinput.f irreps.f iterate.f mixedbasis.f quasiparticle.f readwrite_fleur.f readwrite_fleurR3.f readwrite_fleurR5.1.f readwrite_fleurR5.f readwrite_tmpl.f selfenergy.f selfenergy_com.f selfenergy_wannier.f spectra_wannier.f spex.f susceptibility.f susceptibility_wannier.f tetrahedron.f timer_util.f wannier.f wavefproducts.f wrapper.f
MOD8 = readwrite$(DFT).f
DEP8 = checkinput.f exchange.f getinput.f irreps.f iterate.f mixedbasis.f quasiparticle.f selfenergy.f selfenergy_com.f spex.f susceptibility.f symmetry.f wannier.f
MOD9 = timer_util.f
DEP9 = coulombmatrix.f selfenergy.f spex.f susceptibility.f wavefproducts.f wrapper.f
MOD10 = util.f
DEP10 = Mutil.f bandinfo.f bethesalpeter.f checkinput.f checkmem.f continuation.f correlation.f correlation_wannier.f coulomb_wannier.f coulombmatrix.f divergence.f dump.f dwavefproducts.f exchange.f fft.f getinput.f hilbert.f irreps.f iterate.f key.f mixedbasis.f numerics.f pade.f quasiparticle.f readwrite_fleur.f readwrite_fleurR3.f readwrite_fleurR5.1.f readwrite_fleurR5.f selfenergy.f selfenergy_com.f spectra_wannier.f spex.f susceptibility.f susceptibility_wannier.f symmetry.f tetrahedron.f timer_util.f trafo.f wannier.f wrapper.f
MOD11 = wrapper.f
DEP11 = Mutil.f bandinfo.f bethesalpeter.f checkinput.f continuation.f correlation.f coulomb_wannier.f coulombmatrix.f divergence.f dwavefproducts.f exchange.f freqintegral.f getinput.f ibc.f irreps.f iterate.f mixedbasis.f numerics.f overlap.f quasiparticle.f selfenergy.f selfenergy_com.f spectra_wannier.f susceptibility.f susceptibility_core.f susceptibility_wannier.f symmetry.f trafo.f wannier.f wavefproducts.f
$(DEP0:%.f=%.o): $(findstring $(MOD0:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP1:%.f=%.o): $(findstring $(MOD1:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP2:%.f=%.o): $(findstring $(MOD2:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP3:%.f=%.o): $(findstring $(MOD3:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP4:%.f=%.o): $(findstring $(MOD4:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP5:%.f=%.o): $(findstring $(MOD5:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP6:%.f=%.o): $(findstring $(MOD6:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP7:%.f=%.o): $(findstring $(MOD7:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP8:%.f=%.o): $(findstring $(MOD8:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP9:%.f=%.o): $(findstring $(MOD9:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP10:%.f=%.o): $(findstring $(MOD10:%.f=%.o),$(SRC:%.f=%.o)) $@
$(DEP11:%.f=%.o): $(findstring $(MOD11:%.f=%.o),$(SRC:%.f=%.o)) $@
