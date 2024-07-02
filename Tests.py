import unittest
import numpy as np
from LogToCsvKml import least_squares , calculate_satellite_position , one_epoch , ephemeris , ecef_to_lla


class TestAlgorithm(unittest.TestCase):


    def calculate_distance(self, point1, point2):
        return np.linalg.norm(np.array(point1) - np.array(point2))
    
    def test_position(self):
        real_position = (4436894.27800659, 3085290.16479373, 3376331.62495113)
        sv_position = calculate_satellite_position(ephemeris, one_epoch['tTxSeconds'])

        xs = sv_position[['x_k', 'y_k', 'z_k']].to_numpy()
        pr = sv_position['Psuedo-range']
        b0 = 0
        x0 = np.array([0, 0, 0])
        position_from_algo, b, dp = least_squares(xs, pr, x0, b0)
        distance = self.calculate_distance(real_position, position_from_algo)
        margin = 30
        self.assertLessEqual(distance, margin)
    
    def testXYZConLLA(self):
        xyz_position = (4436894.27800659, 3085290.16479373, 3376331.62495113)
        lla_position = (32.16885701, 34.81365361 , 67.16203436)
        lla_Suppose_position = ecef_to_lla(xyz_position)
        distance = self.calculate_distance(lla_position, lla_Suppose_position)
        margin = 5
        self.assertLessEqual(distance, margin)

if __name__ == '__main__':
    unittest.main()