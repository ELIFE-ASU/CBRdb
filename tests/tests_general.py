import CBRdb


def test_side_to_dict():
    tmp = CBRdb.side_to_dict('C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('C00001 + C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('1 C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('+1 C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('-1 C00001 + 1 C00002')
    assert tmp == {'C00001': -1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('n C00001 + 1 C00002')
    assert tmp == {'C00001': 'n', 'C00002': 1}

    tmp = CBRdb.side_to_dict('n+1 C00001 + 1 C00002')
    assert tmp == {'C00001': 'n+1', 'C00002': 1}

    tmp = CBRdb.side_to_dict('n-1 C00001 + 1 C00002')
    assert tmp == {'C00001': 'n-1', 'C00002': 1}

    tmp = CBRdb.side_to_dict('0 C00001 + 1 C00002')
    assert tmp == {'C00001': 0, 'C00002': 1}

    tmp = "1 C00007 + 2 C00339 + C01438"
    assert tmp == {'C00007': 1, 'C00339': 2, 'C01438': 1}


def test_convert_formula_to_dict():
    tmp = CBRdb.convert_formula_to_dict("C2H2*BrO2")
    assert tmp == {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 1}

    tmp = CBRdb.convert_formula_to_dict("C2H2*32BrO2")
    assert tmp == {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 32}

    tmp = CBRdb.convert_formula_to_dict("Te+")
    assert tmp == {'Te': 1, '+': 1}

    tmp = CBRdb.convert_formula_to_dict("Te-")
    assert tmp == {'Te': 1, '-': 1}

    tmp = CBRdb.convert_formula_to_dict("Te+1")
    assert tmp == {'Te': 1, '+': 1}

    tmp = CBRdb.convert_formula_to_dict("Te-1")
    assert tmp == {'Te': 1, '-': 1}

    tmp = CBRdb.convert_formula_to_dict("C2H4*NO2-")
    assert tmp == {'C': 2, 'H': 4, 'N': 1, 'O': 2, '-': 1, '*': 1}




if __name__ == "__main__":
    test_side_to_dict()
    test_convert_formula_to_dict()
    print("tests_general.py passed")
