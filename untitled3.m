load('ukf_wonoise.mat')
load('ekf_wonoise.mat')
load('pf_wonoise.mat')

ekf_wonoise.plotDeviationVersusError(ekf_wonoise.Deviations)
ukf_wonoise.plotDeviationVersusError(ukf_wonoise.Deviations)
pf_wonoise.plotDeviationVersusError(pf_wonoise.Deviations)