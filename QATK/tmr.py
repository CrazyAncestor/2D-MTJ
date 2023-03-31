# Calculate conductance for parallel spin
transmission_para = nlread('ParaCobaltGrCobalt.hdf5', TransmissionSpectrum)[0]
conductance_para_uu = transmission_para.conductance(spin=Spin.Up)
conductance_para_dd = transmission_para.conductance(spin=Spin.Down)
conductance_para = conductance_para_uu + conductance_para_dd

# Calculate conductance for anti-parallel spin
transmission_anti = nlread('AntiCobaltGrCobalt.hdf5', TransmissionSpectrum)[0]
conductance_anti_uu = transmission_anti.conductance(spin=Spin.Up)
conductance_anti_dd = transmission_anti.conductance(spin=Spin.Down)
conductance_anti = conductance_anti_uu + conductance_anti_dd

print('Conductance Parallel Spin (Siemens)')
print('Up=%8.2e, Down=%8.2e' % (conductance_para_uu.inUnitsOf(Siemens),
                                conductance_para_dd.inUnitsOf(Siemens)))
print('Total = %8.2e' % (conductance_para.inUnitsOf(Siemens)))
print()

print('Conductance Anti-Parallel Spin (Siemens)')
print('Up=%8.2e, Down=%8.2e' % (conductance_anti_uu.inUnitsOf(Siemens),
                                conductance_anti_dd.inUnitsOf(Siemens)))
print('Total = %8.2e' % (conductance_anti.inUnitsOf(Siemens)))
print()

print('TMR (optimistic)  = %8.2f percent' % \
      (100.*(conductance_para-conductance_anti)/conductance_anti))
print('TMR (pessimistic) = %8.2f percent' % \
      (100.*(conductance_para-conductance_anti)/(conductance_para+conductance_anti)))
