import imp
try:
    reload
except NameError:
    from imp import reload


def test_reload():
    import mdp
    imp.reload(mdp)
