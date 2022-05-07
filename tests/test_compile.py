import unittest
# Used here to mock different operating systems
from unittest.mock import patch
import os
import shutil

from gget.compile import compile_muscle

# Get absolute package path
PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))

class TestCompilerWindows(unittest.TestCase):
    def test_compiler_windows(self):
        with self.assertRaises(OSError):
            with patch('platform.system', MagicMock(return_value="Windows")):
                compile_muscle()

class TestCompilerLinux(unittest.TestCase):
    def test_compiler_linux(self):
        with patch('platform.system', MagicMock(return_value="Linux")):
            compile_muscle()
        # Assert that muscle binary was created
        path = os.path.join(PACKAGE_PATH, "bins/compiled/muscle/src/Linux/muscle")
        self.assertIsFile(path)

    def tearDown(self):
        super(TestCompilerLinux, self).tearDown()
        # Delete muscle folder
        shutil.rmtree(os.path.join(PACKAGE_PATH, "bins/compiled/muscle"))

class TestCompilerDarwin(unittest.TestCase):
    def test_compiler_linux(self):
        with patch('platform.system', MagicMock(return_value="Darwin")):
            compile_muscle()
        # Assert that muscle binary was created
        path = os.path.join(PACKAGE_PATH, "bins/compiled/muscle/src/Darwin/muscle")
        self.assertIsFile(path)

    def tearDown(self):
        super(TestCompilerDarwin, self).tearDown()
        # Delete muscle folder
        shutil.rmtree(os.path.join(PACKAGE_PATH, "bins/compiled/muscle"))