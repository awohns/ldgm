__version__ = "undefined"
try:
    from . import _version

    __version__ = _version.version
except ImportError:  # pragma: nocover
    try:
        from setuptools_scm import get_version

        __version__ = get_version(root="..", relative_to=__file__)
    except ImportError:
        pass
