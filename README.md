# HuHoLa : Hummock - Hollow - Lawn Microtopography model
`Mires` `Peatlands` `Microtopography` `Hydromorphology` 
`Hydrology` `DEM` `Fill sinks`

## Description
HuHoLa is a simple model that relies on fill algorithms for identifying microtopography features using a digital elevation model (DEM).

Fill algorithms are used in catchment delineation for instance, to remove small depressions that may trap the water flow and mislead the delineation. In this application, it effectively fills the depressions, which means that we can potentially identify hollows that will be filled, vs. hummocks and lawns that will not be filled.

To identify hummocks, the DEM is flipped upside down by substracting it from its maximum value. Applying the fill algorithm to the inverted DEM identifies the hummocks, and everything that has not been filled in either of the DEM and the inverted DEM are considered to be lawns.

## Dependencies
`WBT` `rasterio` `numpy`

The fill algorithms are applied using Whitebox Tools

...

## How to use HuHoLa model

## Known issues
* Artificial barriers:
* Very high resolution digital elevation models:
* Tradeoff between the precision of hummocks & hollows vs. that of lawns:
