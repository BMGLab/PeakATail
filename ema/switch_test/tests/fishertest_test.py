import unittest
import pandas as pd
from ema.switch_test.fishertest import fishertest, process_gene

class TestFisherTest(unittest.TestCase):
    def setUp(self):
        self.selected_cells = pd.DataFrame({
            'cluster1': [1, 2, 3, 4, 5],
            'cluster2': [5, 4, 3, 2, 1]
        }, index=pd.MultiIndex.from_tuples(
            [('gene1', 'pas1'), ('gene1', 'pas2'), ('gene2', 'pas1'), ('gene2', 'pas2'), ('gene3', 'pas1')],
            names=['Ensemble_ID', 'pas']
        ))
        self.result_dir = 'test_results.txt'

    def test_process_gene(self):
        gene = 'gene1'
        group = self.selected_cells.loc['gene1']
        cluster1 = 'cluster1'
        cluster2 = 'cluster2'
        result = process_gene(gene, group, cluster1, cluster2)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_fishertest(self):
        fishertest(self.selected_cells, self.result_dir)
        with open(self.result_dir, 'r') as f:
            results = f.readlines()
        self.assertEqual(len(results), 5)

if __name__ == '__main__':
    unittest.main()