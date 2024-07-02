import pandas as pd
import numpy as np
from pykalman import KalmanFilter
import simplekml
from pyproj import Proj, transform


# Placeholder function for acquiring raw GNSS measurements from a log file
def read_gnss_log(file_path):
    data = pd.read_csv(file_path)
    return data


# Function to filter satellites based on power
def filter_satellites(satellites, min_power=30):
    filtered = satellites[satellites['CN0'] > min_power]
    return filtered


# Kalman Filter for real-time location determination
def kalman_filter(observations):
    initial_state = observations[0]
    kf = KalmanFilter(initial_state_mean=initial_state, n_dim_obs=3, n_dim_state=3)
    state_means, _ = kf.filter(observations)
    return state_means[-1]  # Return the last state estimate


# Function to detect disturbances
def detect_disturbances(satellites, threshold=50):
    disturbances = []
    prev_pos = None
    for idx, row in satellites.iterrows():
        if prev_pos is not None:
            distance = np.linalg.norm(np.array([row['Pos_X'], row['Pos_Y'], row['Pos_Z']]) - prev_pos)
            if distance > threshold:
                disturbances.append(idx)
        prev_pos = np.array([row['Pos_X'], row['Pos_Y'], row['Pos_Z']])
    return disturbances


# Function to handle detected disturbances
def handle_disturbances(satellites, disturbances):
    for idx in disturbances:
        if idx > 0:
            satellites.at[idx, 'Pos_X'] = satellites.at[idx - 1, 'Pos_X']
            satellites.at[idx, 'Pos_Y'] = satellites.at[idx - 1, 'Pos_Y']
            satellites.at[idx, 'Pos_Z'] = satellites.at[idx - 1, 'Pos_Z']
    return satellites


# Conversion function from ECEF to Geodetic coordinates (lat, lon, alt)
def ecef_to_geodetic(x, y, z):
    proj_ecef = Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    proj_latlon = Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lon, lat, alt = transform(proj_ecef, proj_latlon, x, y, z, radians=False)
    return lat, lon, alt


# Function to compute real-time location and generate output files
def compute_real_time_location(input_file, output_kml_file, output_csv_file):
    raw_data = read_gnss_log(input_file)

    # Filter satellites
    filtered_data = filter_satellites(raw_data)

    # Detect and handle disturbances
    disturbances = detect_disturbances(filtered_data)
    filtered_data = handle_disturbances(filtered_data, disturbances)

    # Extract observations (positions)
    observations = np.array(list(zip(filtered_data['Pos_X'], filtered_data['Pos_Y'], filtered_data['Pos_Z'])))

    # Filter out invalid observations
    valid_observations = observations[~np.isnan(observations).any(axis=1)]

    real_time_locations = []

    if len(valid_observations) > 0:
        corrected_position = kalman_filter(valid_observations)
        lat, lon, alt = ecef_to_geodetic(corrected_position[0], corrected_position[1], corrected_position[2])
        real_time_locations.append((lat, lon, alt))
    else:
        lat, lon, alt = None, None, None

    # Write output KML file
    kml = simplekml.Kml()
    for idx, row in filtered_data.iterrows():
        kml.newpoint(name=str(row['SatPRN (ID)']), coords=[(row['Lon'], row['Lat'])])
    kml.save(output_kml_file)

    # Write output CSV file
    filtered_data['Lat'] = lat
    filtered_data['Lon'] = lon
    filtered_data['Alt'] = alt
    filtered_data.to_csv(output_csv_file, index=False, columns=[
        'GPS time', 'SatPRN (ID)', 'Sat.X', 'Sat.Y', 'Sat.Z', 'Pseudo-Range',
        'CN0', 'GPS_Unique_Time', 'Pos_X', 'Pos_Y', 'Pos_Z', 'Lat', 'Lon', 'Alt'
    ])

    # Print real-time location
    if lat is not None and lon is not None and alt is not None:
        print(f"Real-time location: Latitude={lat}, Longitude={lon}, Altitude={alt}")
    else:
        print("Real-time location could not be determined due to insufficient valid data.")

def main():
    input_file = 'outcomes/walking.csv'
    output_kml_file = 'computed_path.kml'
    output_csv_file = 'computed_data.csv'

    compute_real_time_location(input_file, output_kml_file, output_csv_file)

main()

