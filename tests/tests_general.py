import CBRdb



# # Example usage
# matches = ["C00001", "C00002", "C00001"]
# coeff_out = [1, "2", "3"]
# merged_matches, merged_coeff_out = merge_duplicates(matches, coeff_out)
# print(merged_matches)  # Output: ['C00001', 'C00002']
# print(merged_coeff_out)  # Output: [4, 2]



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
    assert tmp == {'C00002': 1}

    tmp = CBRdb.side_to_dict("1 C00007 + 2 C00339 + C01438")
    assert tmp == {'C00007': 1, 'C00339': 2, 'C01438': 1}

    tmp = CBRdb.side_to_dict("C00024 + n C00083 + n C00005 + n C00004 + 2n C00080")
    assert tmp == {'C00024': 1, 'C00083': 'n', 'C00005': 'n', 'C00004': 'n', 'C00080': '2n'}

    tmp = CBRdb.side_to_dict("C00003 + C00039(n) + C02128(m)")
    assert tmp == {'C00003': 1, 'C00039': 'n', 'C02128': 'm'}

    tmp = CBRdb.side_to_dict("C00020 + C00455 + C00039(n+m)")
    assert tmp == {'C00020': 1, 'C00455': 1, 'C00039': 'n+m'}

    tmp = CBRdb.side_to_dict("C03323(m) + C03323(n)")
    assert tmp == {'C03323': 'm+n'}

    tmp = CBRdb.side_to_dict("C03323(m-1) + C03323(n+1)")
    assert tmp == {'C03323': 'm-1+n+1'}

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
