"""Satellite Tracker entry point. On the first run the tracker opens in "click your location" mode."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def bat_environment():
    """what Satellite Tracker.bat sets up, for when Windows starts this file directly - a pinned taskbar button
    runs 'pythonw.exe -s start.py' (see taskbar.py), not the .bat. Only fills in what is missing, so a start from
    the .bat is unchanged."""
    runtime = os.path.dirname(sys.executable)          # the bundled python folder
    dirs = [runtime] + [os.path.join(runtime, *p) for p in (('Library', 'mingw-w64', 'bin'), ('Library', 'usr', 'bin'),
                                                            ('Library', 'bin'), ('Scripts',), ('bin',))]
    have = [os.path.normcase(os.path.normpath(p)) for p in os.environ.get('PATH', '').split(os.pathsep) if p]
    missing = [d for d in dirs if os.path.normcase(os.path.normpath(d)) not in have]
    if missing:
        os.environ['PATH'] = os.pathsep.join(missing + [os.environ.get('PATH', '')])
    plugins = os.path.join(runtime, 'Library', 'plugins')
    if os.path.isdir(plugins):
        os.environ.setdefault('QT_PLUGIN_PATH', plugins)
    os.environ.setdefault('PYTHONNOUSERSITE', '1')      # for the Pythons we start (this one has -s)
    os.chdir(HERE)                                      # the .bat runs from app\


bat_environment()
# the orbit fits are small; one maths thread avoids OpenBLAS reserving a buffer for every CPU core
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, HERE)
import iridium_tracker  # noqa: E402

if __name__ == '__main__':
    iridium_tracker.main()
