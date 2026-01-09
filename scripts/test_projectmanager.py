# /* ---------------------------------------------------------------------
# *                                       _
# *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
# * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
# * | | | | | | (_| | |  | | | | | | (_) | |_
# * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
# *
# * Unit of Strength of Materials and Structural Analysis
# * University of Innsbruck,
# * 2020 - today
# *
# * festigkeitslehre@uibk.ac.at
# *
# * Matthias Neuner matthias.neuner@uibk.ac.at
# *
# * This file is part of the MAteRialMOdellingToolbox (marmot).
# *
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * The full text of the license can be found in the file LICENSE.md at
# * the top level directory of marmot.
# * ---------------------------------------------------------------------
# */

import os
import tempfile
import unittest
from pathlib import Path

from projectmanager import walk_modules


class TestWalkModules(unittest.TestCase):
    """Test cases for the walk_modules function in projectmanager.py"""

    def setUp(self):
        """Create a temporary directory structure for testing"""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root_path = Path(self.temp_dir.name)
        # Store original working directory to restore later
        self.original_cwd = os.getcwd()
        # Change to temp directory to use relative paths
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        """Clean up temporary directory"""
        # Restore original working directory
        os.chdir(self.original_cwd)
        self.temp_dir.cleanup()

    def test_walk_modules_finds_module_at_depth_1(self):
        """Test that walk_modules finds modules at depth 1"""
        # Create a module at depth 1
        module_dir = Path("module1")
        module_dir.mkdir()
        (module_dir / "module.cmake").touch()

        # levels=2 to find modules 1 level deep from root
        result = walk_modules(".", levels=2)
        self.assertEqual(len(result), 1)
        # The function returns paths like './module1'
        self.assertTrue(any('module1' in r for r in result))

    def test_walk_modules_respects_depth_limit(self):
        """Test that walk_modules respects the depth limit parameter"""
        # Create modules at different depths
        module1 = Path("module1")
        module1.mkdir()
        (module1 / "module.cmake").touch()

        module2 = Path("level1") / "module2"
        module2.mkdir(parents=True)
        (module2 / "module.cmake").touch()

        module3 = Path("level1") / "level2" / "module3"
        module3.mkdir(parents=True)
        (module3 / "module.cmake").touch()

        # Test with levels=2 - should find only module1 (1 separator)
        result = walk_modules(".", levels=2)
        self.assertEqual(len(result), 1)
        self.assertTrue(any('module1' in r for r in result))

        # Test with levels=3 - should find module1 and module2 (up to 2 separators)
        result = walk_modules(".", levels=3)
        self.assertEqual(len(result), 2)
        self.assertTrue(any('module1' in r for r in result))
        self.assertTrue(any('module2' in r for r in result))

        # Test with levels=4 - should find all three modules (up to 3 separators)
        result = walk_modules(".", levels=4)
        self.assertEqual(len(result), 3)
        self.assertTrue(any('module1' in r for r in result))
        self.assertTrue(any('module2' in r for r in result))
        self.assertTrue(any('module3' in r for r in result))

    def test_walk_modules_ignores_dirs_without_module_cmake(self):
        """Test that directories without module.cmake are ignored"""
        # Create directories without module.cmake
        dir1 = Path("dir1")
        dir1.mkdir()
        dir2 = Path("dir2")
        dir2.mkdir()

        # Create one directory with module.cmake
        module_dir = Path("module1")
        module_dir.mkdir()
        (module_dir / "module.cmake").touch()

        result = walk_modules(".", levels=2)
        self.assertEqual(len(result), 1)
        self.assertTrue(any('module1' in r for r in result))

    def test_walk_modules_empty_directory(self):
        """Test walk_modules on an empty directory"""
        result = walk_modules(".", levels=1)
        self.assertEqual(len(result), 0)

    def test_walk_modules_with_nested_structure(self):
        """Test walk_modules with a realistic nested module structure"""
        # Create a structure similar to the actual modules directory
        core_dir = Path("core")
        core_dir.mkdir()

        math_core = core_dir / "MathCore"
        math_core.mkdir()
        (math_core / "module.cmake").touch()

        mechanics_core = core_dir / "MechanicsCore"
        mechanics_core.mkdir()
        (mechanics_core / "module.cmake").touch()

        materials_dir = Path("materials")
        materials_dir.mkdir()

        linear_elastic = materials_dir / "LinearElastic"
        linear_elastic.mkdir()
        (linear_elastic / "module.cmake").touch()

        # Test with levels=3 to find modules 2 levels deep
        result = walk_modules(".", levels=3)
        self.assertEqual(len(result), 3)
        self.assertTrue(any('MathCore' in r for r in result))
        self.assertTrue(any('MechanicsCore' in r for r in result))
        self.assertTrue(any('LinearElastic' in r for r in result))

    def test_walk_modules_handles_trailing_slash(self):
        """Test that walk_modules handles root paths with trailing slashes"""
        module_dir = Path("module1")
        module_dir.mkdir()
        (module_dir / "module.cmake").touch()

        # Test with trailing slash
        result = walk_modules("." + os.path.sep, levels=2)
        self.assertEqual(len(result), 1)
        self.assertTrue(any('module1' in r for r in result))


if __name__ == "__main__":
    unittest.main()
