# GNSS Raw Measurements Processing

## Overview

This project provides a solution for processing GNSS (Global Navigation Satellite System) raw measurements. It includes functionalities to parse raw measurement log files, compute positioning using a naive algorithm, and generate visualization outputs in KML format.

## Requirements

To run this project, you need:

- Python 3.x
- pandas
- numpy
- matplotlib
- navpy
- gnssutils
- simplekml

You can install the dependencies using pip:

```bash
pip install pandas numpy matplotlib navpy gnssutils simplekml georinex unlzw3
```

## Usage
- 1. Clone this repository:
```bash
git clone https://github.com/LiorJerbi/AutoRobots_Ex1
```
- 2. Navigate to the project directory:
```bash
cd AutoRobots_Ex1
```

## How to run
```bash
put your gnss raw measurments txt files in the main folder(AutoRobots_Ex1)
```
- open the "GnssToPosition.py" file and change the input file at the beginning of the code
```bash
###################################################################################
################################ CHANGE FILE NAME #################################

input_filepath = os.path.join('gnss_log_2024_04_13_19_51_17.txt')

################################ CHANGE FILE NAME #################################
###################################################################################
```
- run cmd in the main folder and write the next code:
```bash
python GnssToPosition.py
```

## Output
The script generates the following outputs:
1. GnssToRoute.kml: KML file containing the computed path with timestamps for each position.
2. GNSStoPosition.csv: CSV file with the computed positions and additional columns (Pos.X, Pos.Y, Pos.Z, Lat, Lon, Alt).

## Contributors
- Yael Rosen - 209498211
- Tomer Klugman - 312723612
- Hadas Evers - 206398984
- Lior Jerbi - 314899493
