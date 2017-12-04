# - newDimensionalUnits (pyTree) -
import Converter.Internal as Internal

# Create a DimensionalUnits node
n = Internal.newDimensionalUnits(massUnit='Kilogram', lengthUnit='Meter', timeUnit='Second', temperatureUnit='Kelvin', angleUnit='Radian')
Internal.printTree(n)
#>> ['DimensionalUnits',array(shape=(5,),dtype='object',order='F'),[1 son],'DimensionalUnits_t']
#>>   |_['AdditionalUnits',array('NullNullNullNullNull',dtype='|S1'),[0 son],'AdditionalUnits_t']

# Attach it to a parent node
d = Internal.newGridCoordinates()
Internal.newDataClass('Dimensional', parent=d)
Internal.newDimensionalUnits('Kilogram', 'Meter', 'Second', 'Kelvin', 'Radian', parent=d)
