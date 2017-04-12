Drifts = {}
QuadLengths = {}
QuadGradients = {}
DipLengths = {}
DipAngles = {}
UndulConfigs = {}

cell_keys = ['LPWA-QAP1','QAP1','QAP1-QAP2','QAP2','QAP2-QAP3','QAP3','QAP3-IMG1','IMG1',\
  'IMG1-DIP1','DIP1','DIP1-DIP2','DIP2','DIP2-IMG2','IMG2','IMG2-DIP3','DIP3','DIP3-DIP4','DIP4',\
  'DIP4-QEM1','QEM1','QEM1-QEM2','QEM2','QEM2-IMG4','IMG4','IMG4-QEM3','QEM3','QEM3-QEM4','QEM4',\
  'QEM4-UNDL','UNDL1','UNDL2','UNDL-IMG5','IMG5']

Drifts['LPWA-QAP1'] = 0.04685
Drifts['QAP1-QAP2'] = 0.10295
Drifts['QAP2-QAP3'] = 0.09895
Drifts['QAP3-IMG1'] = 0.52085
Drifts['IMG1-DIP1'] = 0.263
Drifts['DIP1-DIP2'] = 0.2
Drifts['DIP2-IMG2'] = 0.275
Drifts['IMG2-DIP3'] = 0.275
Drifts['DIP3-DIP4'] = 0.2
Drifts['DIP4-QEM1'] = 0.327
Drifts['QEM1-QEM2'] = 0.4
Drifts['QEM2-IMG4'] = 0.15
Drifts['IMG4-QEM3'] = 0.15
Drifts['QEM3-QEM4'] = 0.35
Drifts['QEM4-UNDL'] = 0.298
Drifts['UNDL-IMG5'] = 0.951

QuadLengths['QAP1'] = 0.047
QuadLengths['QAP2'] = 0.0511
QuadLengths['QAP3'] = 0.0323
QuadLengths['QEM1'] = 0.2133
QuadLengths['QEM2'] = 0.2133
QuadLengths['QEM3'] = 0.2133
QuadLengths['QEM4'] = 0.2133

QuadGradients['QAP1'] = 172.9488
QuadGradients['QAP2'] = -169.3680
QuadGradients['QAP3'] = 150.9580

QuadGradients['QEM1'] = -1.817065
QuadGradients['QEM2'] = 3.090898
QuadGradients['QEM3'] = -4.477488
QuadGradients['QEM4'] = 2.760612

DipLengths['DIP1'] = 0.2
DipLengths['DIP2'] = 0.2
DipLengths['DIP3'] = 0.2
DipLengths['DIP4'] = 0.2

DipAngles['DIP1'] = 0.08
DipAngles['DIP2'] = -0.08
DipAngles['DIP3'] = -0.08
DipAngles['DIP4'] = 0.08

UndulConfigs['Period'] = 0.018
UndulConfigs['Strength'] = 1.62
UndulConfigs['NumPeriods'] = 54

BeamEnergy_ref = 0.167 #0.176

latt_par_string = """
###############################
R_11 = {0:.3g}, R_33 = {1:.3g},
R_56 = {2:.3g} mm,
R_226 = {3:.3g}, R_446 = {4:.3g},
R_126 = {5:.3g}, R_346 = {6:.3g}
################################
"""
