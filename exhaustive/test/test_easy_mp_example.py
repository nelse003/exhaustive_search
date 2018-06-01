import pytest
from easy_mp_example import try_easy_mp

def test_try_easy_mp():
    result = try_easy_mp()
    assert result[-1] == 10000