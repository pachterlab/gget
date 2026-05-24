import unittest
import pandas as pd
from gget.gget_dataverse import dataverse
import os
import shutil

#TODO: Verify the test code, this is drafted using co-pilot!
class TestDataverse(unittest.TestCase):
    def test_dataverse_download(self):
        df = pd.DataFrame({
            'id': [6180617],
            'name': ['nodes'],
            'type': ['tab']
        })

        dataverse(df, 'temp_datasets')

        # Check if the file is downloaded
        self.assertTrue(os.path.exists('temp_datasets/nodes.tab'))

        # Clean up by removing the datasets folder
        shutil.rmtree('temp_datasets')

if __name__ == '__main__':
    unittest.main()