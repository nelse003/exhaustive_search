from exhaustive.depreceated.easy_mp_example import try_easy_mp, try_f_caller

def test_try_easy_mp():
    result = try_easy_mp()
    assert result[-1] == 10000

def test_try_f_caller():
    result = try_f_caller()
    assert result[-1] == 10000
