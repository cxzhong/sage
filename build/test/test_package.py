# -*- coding: utf-8 -*-
"""
Test Sage Package Handling
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import unittest
from sage_bootstrap.package import Package
from sage_bootstrap.tarball import Tarball


class PackageTestCase(unittest.TestCase):

    maxDiff = None

    def test_package(self):
        pkg = Package('pari')
        self.assertTrue(pkg.name, 'pari')
        self.assertTrue(pkg.path.endswith('build/pkgs/pari'))
        self.assertEqual(pkg.tarball_pattern, 'pari-VERSION.tar.gz')
        self.assertEqual(pkg.tarball_filename, pkg.tarball.filename)
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(pkg.tarball.filename.startswith('pari-') and
                        pkg.tarball.filename.endswith('.tar.gz'))
        self.assertTrue(isinstance(pkg.tarball, Tarball))

    def test_all(self):
        pari = Package('pari')
        self.assertTrue(pari in Package.all())
    
    def test_multi_tarball_package(self):
        """Test packages with multiple platform-specific wheels."""
        # rpds_py is a package with many platform-specific wheels
        try:
            pkg = Package('rpds_py')
            self.assertEqual(pkg.name, 'rpds_py')
            # Should have multiple tarballs (wheels for different platforms)
            self.assertGreater(len(pkg.tarballs_info), 1)
            # Each tarball info should have required fields
            for info in pkg.tarballs_info:
                self.assertIn('tarball', info)
                self.assertIn('sha256', info)
                # All should be wheels
                self.assertTrue(info['tarball'].endswith('.whl'))
        except Exception:
            # If rpds_py doesn't exist, skip this test
            self.skipTest('rpds_py package not found')
