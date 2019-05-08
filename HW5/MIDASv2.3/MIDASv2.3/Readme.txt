The mixed frequency regression studies the explanatory power of high frequency variables on the low frequency outcome. The weights associated with high frequency regressors are usually assumed some functional form. This toolbox is a repack of the Mi(xed) Da(ta) S(ampling) regressions (MIDAS) programs written by Eric Ghysels. It supports ADL-MIDAS type regressions. It also includes two functions for GARCH-MIDAS and DCC-MIDAS estimation

Syntax:
[...] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate)
[...] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,name,value,...)
[...] = GarchMidas(y, name,value,...)
[...] = DccMidas(Data, name,value,...)


Version History

v2.1 Add MIDAS quantile regression
v2.0 Add GARCH-MIDAS and DCC-MIDAS
v1.1 Allow MIDAS leads and lags specification 'horizon' be negative. Add true out-of-sample forecast; results are restored in the last output argument 'Extended Forecast' struct. Report the approximated dates associated with the forecasted values in the output struct, in character instead of serial dates.
v1.0 First release of the repacked MIDAS regression
