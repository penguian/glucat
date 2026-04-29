import PyClical
import doctest
import pytest

def test_pyclical_doctests():
    results = doctest.testmod(PyClical, verbose=False)
    if results.failed > 0:
        pytest.fail(f"PyClical doctests failed: {results}")
