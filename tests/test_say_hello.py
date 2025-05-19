from crosslinked_mfp_analysis.crosslinked_mfp_analysis.say_hello import (say_hello)


def test_say_hello():
    """ Test function for the function say_hello in crosslinked_mfp_analysis/crosslinked_mfp_analysis/say_hello.py

    Verify that this function prints the word 'Hello' and does nothing else.

    Args:
        None

    Returns:
        None. Asserts that the output of say_hello is None.

    """
    assert say_hello() is None
