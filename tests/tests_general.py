import CBRdb


def test_side_to_dict():
    tmp = CBRdb.side_to_dict('C00001 + 1 C00002')
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


if __name__ == "__main__":
    test_side_to_dict()
    print("tests_general.py passed")
