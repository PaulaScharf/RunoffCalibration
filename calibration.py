import builtins
import math

from pcraster import *
from pcraster.framework import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import plotly.express as px
import plotly.graph_objects as go

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
    self.meltRate = meltRate
    print('meltRate: ', self.meltRate)
    print('temperatureLapseRate: ', temperatureLapseRate)

    # the measurement location
    self.observationLocation=self.readmap('out')

    dem = self.readmap('dem')
    elevationMeteoStation = 208.1
    elevationAboveMeteoStation = dem - elevationMeteoStation
    self.temperatureCorrection = elevationAboveMeteoStation * temperatureLapseRate
    # temperature correction value
    # self.report(self.temperatureCorrection,'output/tempCor')

    self.snow=0.0

    self.ldd=lddcreate(dem,1e31,1e31,1e31,1e31)
    # map of localdrain direction
    # self.report(self.ldd,'ldd')

    # example how to calculate total precipitation
    self.totPrecip=scalar(0)


  def dynamic(self):
    precipitation = timeinputscalar('input/precip.tss',1)
    # observed precipitation
    # self.report(precipitation,'output/pFromTss')
    temperatureObserved = timeinputscalar('input/temp.tss',1)
    # observed temperature
    # self.report(temperatureObserved,'output/tempObs')
    temperature= temperatureObserved - self.temperatureCorrection
    # temperature corrected by elevation
    # self.report(temperature,'output/temp')

    freezing=temperature < 0.0
    # temperature where water freezes
    # self.report(freezing,'output/freez')
    snowFall=ifthenelse(freezing,precipitation,0.0)
    # precipiation that is snow
    # self.report(snowFall,'output/snowFall')
    rainFall=ifthenelse(pcrnot(freezing),precipitation,0.0)
    # precipitaion that is rain
    # self.report(rainFall,'output/rain')

    self.snow = self.snow+snowFall
    # snow store without melt
    # self.report(self.snow,'snow')

    potentialMelt = ifthenelse(pcrnot(freezing),temperature*self.meltRate,0)
    # potential melt
    # self.report(potentialMelt,'output/pmelt')
    actualMelt = min(self.snow, potentialMelt)
    # actual melted snow
    # self.report(actualMelt,'output/amelt')

    self.snow = self.snow - actualMelt
    # snow store after precipitation and melt
    # self.report(self.snow,'output/snow')

    runoffGenerated = actualMelt + rainFall
    # water as melting snow + rain
    # self.report(runoffGenerated,'output/rg')

    discharge=accuflux(self.ldd,runoffGenerated*cellarea())
    self.report(discharge,'output/q')

    # reading values from a map at the observation location
    runoffAtOutflowPoint=getCellValueAtBooleanLocation(self.observationLocation,discharge)
    observed = timeinputscalar('input/observed.tss',1)
    observedAtOutflowPoint = getCellValueAtBooleanLocation(self.observationLocation,observed)
    global squaredErrorStore
    squaredErrorStore = squaredErrorStore + ((runoffAtOutflowPoint - observedAtOutflowPoint) ** 2)

nrOfTimeSteps=181
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
squaredErrorStore = 0
results = []

## coarse resolution
'''
for i in range(1,20,1):
  meltRate=0.002*(i)
  for k in range(1,20,1):
      temperatureLapseRate = 0.001 * (k)
      squaredErrorStore = 0
      dynamicModel.run()
      meanSquaredError = squaredErrorStore/nrOfTimeSteps
      results.append([math.sqrt(meanSquaredError), meltRate,temperatureLapseRate])
      print(math.sqrt(meanSquaredError))

## fine resolution
for i in range(26,60,1):
  meltRate=0.0002*(i)
  for k in range(6,30,1):
      temperatureLapseRate = 0.0002 * (k)
      squaredErrorStore = 0
      dynamicModel.run()
      meanSquaredError = squaredErrorStore/nrOfTimeSteps
      results.append([math.sqrt(meanSquaredError), meltRate,temperatureLapseRate])
      print(math.sqrt(meanSquaredError))
          '''
## finest resolution
for i in range(78,98,1):
  meltRate=0.0001*(i)
  for k in range(32,52,1):
      temperatureLapseRate = 0.0001 * (k)
      squaredErrorStore = 0
      dynamicModel.run()
      meanSquaredError = squaredErrorStore/nrOfTimeSteps
      results.append([math.sqrt(meanSquaredError), meltRate,temperatureLapseRate])
      print(math.sqrt(meanSquaredError))

min_results = builtins.min(results, key=lambda x: x[0])
print(min_results)

fig = plt.figure()
ax = Axes3D(fig)
#ax.plot_trisurf([item[2] for item in results], [item[0] for item in results], np.array([item[1] for item in results]))
ax.plot_trisurf([item[1] for item in results], [item[2] for item in results], np.array([item[0] for item in results]))

fig = px.scatter_3d(results, [item[1] for item in results], [item[2] for item in results], [item[0] for item in results],
                    labels={'x':'snow melt rate', 'y':'tempreature lapse rate', 'z':'root mean square deviation'})
fig.update_layout(
    title={
        'text': "Finest resolution",
        'x': 0.5,
        'xanchor': 'center',
        'font_size': 30
    }
)
fig.write_html("results/scatter_finest.html")
fig.show()
