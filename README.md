# RuO2
Detailed data of [Phys. Rev. B 99, 184432 (2019)].

The [DFT] directory contains the results of the DFT runs, including all the input/output files.
- [RuO2_1_GGA_afterSCF.tar.gz]: The result of the NM GGA.
- [RuO2_2_GGA_U_afterSCF.tar.gz]: The result of the AFM GGA+U, with U = 2.8 eV and J = 0.2 eV.

The [TEX] directory contains the main script and the figures for the arXiv.

The [HYB] directory contains the codes for the Hartree-Fock (static mean-field) calculations.
- [codes]: Below is how to install the codes.

1. Install [SRC_dmft1.1].
[dmft]: To calculate the occupation matrix.

2. Install Hartree-Fock codes.
ifort -o init_sig init_sig.f90
ifort -o mk_kanamori_mf_h mk_kanamori_mf_h.f90
ifort -o mx_occ mx_occ.f90
ifort -o ck_mm ck_mm.f90

[init_sig]: To make the ``zero'' self-energy. Making the input, only once at the first stage.
[mk_kanamori]: To put the Kanamori parametrization into the H(R=0). It prints a hamilt file out, which is the input for [dmft].
[mx_occ]: To mix the occpation matrix, between the old and new iteration.
[ck_mm]: To check the spin magnetic moment. Run only once after convergence.

- [example]: The result of U = 1.7 eV and J = 0.2 eV.
