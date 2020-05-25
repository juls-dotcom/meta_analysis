# -*- coding: utf-8 -*-

import pytest
from meta_analysis.skeleton import fib

__author__ = "Julien Hernandez Lallement"
__copyright__ = "Julien Hernandez Lallement"
__license__ = "mit"


def test_fib():
    assert fib(1) == 1
    assert fib(2) == 1
    assert fib(7) == 13
    with pytest.raises(AssertionError):
        fib(-10)
