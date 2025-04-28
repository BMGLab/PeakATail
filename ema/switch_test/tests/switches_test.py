import unittest
from ema.switch_test.switchs import apa_switch

class TestSwitchs(unittest.TestCase):
    def test_apa_switch(self):
        # Test case 1
        cell_combine = [[0, 1], [2, 3]]
        apa_switch(cell_combine)
        # Add assertions to verify the expected behavior

        # Test case 2
        cell_combine = [[4, 5], [6, 7]]
        apa_switch(cell_combine)
        # Add assertions to verify the expected behavior

        # Test case 3
        cell_combine = [[8, 9], [10, 11]]
        apa_switch(cell_combine)
        # Add assertions to verify the expected behavior

if __name__ == '__main__':
    unittest.main()