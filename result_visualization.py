from pcraster import *
from pcraster.framework import *
import numpy as np
import matplotlib.pyplot as plt

# notes #
# out.map: location where discharge was measured
# observed.tss: observed discharge (same unit as discharge variable in script)

# helper functions to get values from a map
def getCellValue(Map, Row, Column):
  Value, Valid=cellvalue(Map, Row, Column)
  if Valid:
    return Value
  else:
    print('missing value in input of getCellValue')

def getCellValueAtBooleanLocation(location,map):
  # map can be any type, return value always float
  valueMap=mapmaximum(ifthen(location,scalar(map)))
  value=getCellValue(valueMap,1,1)
  return value

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('input/dem.map')

  def initial(self):
    # example how to provide parameter values to the script
    # the value of a is assigned here to the parameter aValue

    # the measurement location
    self.observationLocation=self.readmap('out')

    dem = self.readmap('input/dem')
    elevationMeteoStation = 208.1
    elevationAboveMeteoStation = dem - elevationMeteoStation
    #temperatureLapseRate = 0.005
    temperatureLapseRate = 0.0041
    self.temperatureCorrection = elevationAboveMeteoStation * temperatureLapseRate
    #self.report(self.temperatureCorrection,'output/tempCor')

    self.snow=0.0

    self.ldd=lddcreate(dem,1e31,1e31,1e31,1e31)
    #self.report(self.ldd,'output/ldd')

    # example how to calculate total precipitation
    self.totPrecip=scalar(0)

  def dynamic(self):
    precipitation = timeinputscalar('input/precip.tss',1)
    #self.report(precipitation,'output/pFromTss')
    temperatureObserved = timeinputscalar('input/temp.tss',1)
    #self.report(temperatureObserved,'output/tempObs')
    temperature= temperatureObserved - self.temperatureCorrection
    #self.report(temperature,'output/temp')

    freezing=temperature < 0.0
    #self.report(freezing,'output/freez')
    snowFall=ifthenelse(freezing,precipitation,0.0)
    #self.report(snowFall,'output/snowFall')
    rainFall=ifthenelse(pcrnot(freezing),precipitation,0.0)
    #self.report(rainFall,'output/rain')

    self.snow = self.snow+snowFall
    #self.report(self.snow,'snow')

    #potentialMelt = ifthenelse(pcrnot(freezing),temperature*0.01,0)
    potentialMelt = ifthenelse(pcrnot(freezing),temperature*0.0081,0)
    #self.report(potentialMelt,'output/pmelt')
    actualMelt = min(self.snow, potentialMelt)
    #self.report(actualMelt,'output/amelt')

    self.snow = self.snow - actualMelt
    #self.report(self.snow,'output/snow')

    runoffGenerated = actualMelt + rainFall
    #self.report(runoffGenerated,'output/rg')

    discharge=accuflux(self.ldd,runoffGenerated*cellarea())
    self.report(discharge,'output/q')

    # reading values from a map at the observation location
    runoffAtOutflowPoint=getCellValueAtBooleanLocation(self.observationLocation,discharge)
    modelled.append([runoffAtOutflowPoint,self.currentTimeStep()])
    observedMap = timeinputscalar('input/observed.tss',1)
    observedAtOutflowPoint = getCellValueAtBooleanLocation(self.observationLocation,observedMap)
    observed.append([observedAtOutflowPoint,self.currentTimeStep()])
    global squaredErrorStore
    squaredErrorStore = squaredErrorStore + ((runoffAtOutflowPoint - observedAtOutflowPoint) ** 2)


modelled=[]
observed=[]
squaredErrorStore=0

nrOfTimeSteps=181
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

meanSquaredError = squaredErrorStore/nrOfTimeSteps
print(meanSquaredError)

plt.xlabel("timestep (day)")
plt.ylabel("outflow (m/day)")
plt.plot([item[1] for item in observed],[item[0] for item in observed],linestyle='dashed',linewidth=2, label='observed')
plt.plot([item[1] for item in modelled],[item[0] for item in modelled],label='modelled')
plt.legend()
plt.savefig('calibrated.png', dpi=300)
plt.show()
