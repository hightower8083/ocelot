lattice_elements = {}
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
Drifts['QAP3-IMG1'] = 0.25085
Drifts['IMG1-DIP1'] = 0.533 # or add 0.04
Drifts['DIP1-DIP2'] = 0.2
Drifts['DIP2-IMG2'] = 0.275
Drifts['IMG2-DIP3'] = 0.275
Drifts['DIP3-DIP4'] = 0.2
Drifts['DIP4-QEM1'] = 0.32035
Drifts['QEM1-QEM2'] = 0.3867
Drifts['QEM2-IMG4'] = 0.14335
Drifts['IMG4-QEM3'] = 0.14335
Drifts['QEM3-QEM4'] = 0.3367
Drifts['QEM4-UNDL'] = 0.29135
Drifts['UNDL-IMG5'] = 0.573

QuadLengths['QAP1'] = 0.047
QuadLengths['QAP2'] = 0.0511
QuadLengths['QAP3'] = 0.0323
QuadLengths['QEM1'] = 0.2133
QuadLengths['QEM2'] = 0.2133
QuadLengths['QEM3'] = 0.2133
QuadLengths['QEM4'] = 0.2133

QuadGradients['QAP1'] = 183.94514302272
QuadGradients['QAP2'] = -172.97923884219415
QuadGradients['QAP3'] = 149.78095430510771

QuadGradients['QEM1'] = -5.6687515941714768
QuadGradients['QEM2'] = 8.2921337356209825
QuadGradients['QEM3'] = -8.1471210041503745
QuadGradients['QEM4'] = 1.7746903401162113

DipLengths['DIP1'] = 0.2
DipLengths['DIP2'] = 0.2
DipLengths['DIP3'] = 0.2
DipLengths['DIP4'] = 0.2

DipAngles['DIP1'] = 0.02304
DipAngles['DIP2'] = -0.02304
DipAngles['DIP3'] = -0.02304
DipAngles['DIP4'] = 0.02304

UndulConfigs['Period'] = 0.01816
UndulConfigs['Strength'] = 1.80227173
UndulConfigs['NumPeriods'] = 54

BeamEnergy_ref = 0.176

lattice_elements['Drifts'] = Drifts
lattice_elements['QuadLengths'] = QuadLengths
lattice_elements['QuadGradients'] = QuadGradients
lattice_elements['DipLengths'] = DipLengths
lattice_elements['DipAngles'] = DipAngles
lattice_elements['UndulConfigs'] = UndulConfigs

latt_par_string = """
###############################
R_11 = {0:.3g}, R_33 = {1:.3g},
R_56 = {2:.3g} mm,
R_226 = {3:.3g}, R_446 = {4:.3g},
R_126 = {5:.3g}, R_346 = {6:.3g}
################################
"""
