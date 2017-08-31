
import collections as _collections


# provide dnormtomo.__version__
from . import _version

_VersionInfo = _collections.namedtuple("_VersionInfo", ("major", "minor") )

__version__ = _version.version

version_info = _VersionInfo(major=_version.version_maj, minor=_version.version_min)
