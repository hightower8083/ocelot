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

Drifts['QAP3-IMG1'] = 0.52085		##
Drifts['IMG1-DIP1'] = 0.5725		## correct

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
QuadLengths['QEM1'] = 0.2
QuadLengths['QEM2'] = 0.2
QuadLengths['QEM3'] = 0.2
QuadLengths['QEM4'] = 0.2

QuadGradients['QAP1'] = 175.8703
QuadGradients['QAP2'] = -173.622
QuadGradients['QAP3'] = 150.2504
QuadGradients['QEM1'] = -6.678787
QuadGradients['QEM2'] = 8.127195
QuadGradients['QEM3'] = -4.970518
QuadGradients['QEM4'] = 0.04438066

DipLengths['DIP1'] = 0.2
DipLengths['DIP2'] = 0.2
DipLengths['DIP3'] = 0.2
DipLengths['DIP4'] = 0.2

DipAngles['DIP1'] = 1
DipAngles['DIP2'] = -1.
DipAngles['DIP3'] = -1.
DipAngles['DIP4'] = 1.

UndulConfigs['Period'] = 0.02
UndulConfigs['Strength'] = 1.7291
UndulConfigs['NumPeriods'] = 50

BeamEnergy_ref = 0.18

latt_par_string = """
###############################
R_11 = {0:.3g}, R_33 = {1:.3g},
R_56 = {2:.3g} mm,
R_226 = {3:.3g}, R_446 = {4:.3g},
R_126 = {5:.3g}, R_346 = {6:.3g}
################################
"""
