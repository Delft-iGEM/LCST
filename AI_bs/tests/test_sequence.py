from elp_lcst.sequence import parse_elp, guests_from_repeat_spec


def test_parse_simple_repeat():
    seq = "VPGVG" * 10
    p = parse_elp(seq)
    assert p.n_pentads == 10
    assert p.guests == ["V"] * 10
    assert p.is_elp
    assert p.motif_match_fraction == 1.0


def test_parse_mixed_guests():
    seq = "VPGVG" * 5 + "VPGAG" * 2 + "VPGGG" * 3
    p = parse_elp(seq)
    assert p.n_pentads == 10
    assert p.guests.count("V") == 5
    assert p.guests.count("A") == 2
    assert p.guests.count("G") == 3


def test_parse_with_flanking_tag():
    # His-tag + linker then ELP; parser should still find the repeat frame.
    seq = "MGHHHHHH" + "VPGVG" * 8
    p = parse_elp(seq)
    assert p.n_pentads >= 7
    assert set(p.guests) == {"V"}


def test_guest_spec_expansion():
    guests = guests_from_repeat_spec("V5A2G3")
    assert len(guests) == 10
    assert guests.count("V") == 5
    assert guests.count("A") == 2
    assert guests.count("G") == 3
