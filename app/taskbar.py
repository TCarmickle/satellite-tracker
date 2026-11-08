"""The taskbar button's identity on Windows: what a pinned button starts, and with which icon and name.

The tracker runs as pythonw.exe. SetCurrentProcessExplicitAppUserModelID (iridium_tracker.main) already gives the
running window its own button and our icon, but when that button is pinned Windows needs to know how to start
the program again; without being told, it pins pythonw.exe itself - Python's icon, and a click starts a bare
Python with no tracker. The window's "relaunch" properties tell it: our command line, icon and name.

Plain ctypes (no pywin32 / comtypes in the runtime): SHGetPropertyStoreForWindow, then IPropertyStore::SetValue
through its vtable. Never raises - a failure only means pinning falls back to Python's icon as before.
"""
import ctypes
import os
import sys

APP_ID = 'SatelliteTracker.App'      # also passed to SetCurrentProcessExplicitAppUserModelID
NAME = 'Satellite Tracker'

_VT_LPWSTR = 31
# PKEY_AppUserModel_* (propkey.h): one format id, these property ids
_FMTID = '{9F4C2855-9F79-4B39-A8D0-E1D42DE1D5F3}'
_RELAUNCH_COMMAND, _RELAUNCH_ICON, _RELAUNCH_NAME, _ID = 2, 3, 4, 5
_IID_IPROPERTYSTORE = '{886D8EEB-8CF2-4446-8D02-CDBA1DBDCF99}'


class _GUID(ctypes.Structure):
    _fields_ = [('Data1', ctypes.c_ulong), ('Data2', ctypes.c_ushort), ('Data3', ctypes.c_ushort),
                ('Data4', ctypes.c_ubyte * 8)]


class _PROPERTYKEY(ctypes.Structure):
    _fields_ = [('fmtid', _GUID), ('pid', ctypes.c_ulong)]


class _PROPVARIANT(ctypes.Structure):
    # vt + three reserved words, then a 16-byte union (on 64-bit) - only the string pointer is used here
    _fields_ = [('vt', ctypes.c_ushort), ('r1', ctypes.c_ushort), ('r2', ctypes.c_ushort), ('r3', ctypes.c_ushort),
                ('pwszVal', ctypes.c_wchar_p), ('pad', ctypes.c_void_p)]


def _guid(text):
    g = _GUID()
    ctypes.oledll.ole32.CLSIDFromString(ctypes.c_wchar_p(text), ctypes.byref(g))
    return g


def relaunch_info(app_dir):
    """(command line, icon resource) that start this folder's tracker the way Satellite Tracker.bat does"""
    pythonw = os.path.join(os.path.dirname(sys.executable), 'pythonw.exe')
    # -s: no per-user site-packages (the .bat sets PYTHONNOUSERSITE); start.py sets up the rest of the .bat's
    # environment itself
    command = '"%s" -s "%s"' % (pythonw, os.path.join(app_dir, 'start.py'))
    return command, os.path.join(app_dir, 'icon.ico') + ',0'


def set_relaunch(hwnd, app_dir):
    """give window hwnd our taskbar identity; True when Windows took every property"""
    if sys.platform != 'win32':
        return False
    store = ctypes.c_void_p()
    try:
        command, icon = relaunch_info(app_dir)
        ctypes.oledll.shell32.SHGetPropertyStoreForWindow(ctypes.c_void_p(int(hwnd)),
                                                          ctypes.byref(_guid(_IID_IPROPERTYSTORE)), ctypes.byref(store))
        vtbl = ctypes.cast(store, ctypes.POINTER(ctypes.POINTER(ctypes.c_void_p))).contents
        set_value = ctypes.WINFUNCTYPE(ctypes.HRESULT, ctypes.c_void_p, ctypes.POINTER(_PROPERTYKEY),
                                       ctypes.POINTER(_PROPVARIANT))(vtbl[6])
        commit = ctypes.WINFUNCTYPE(ctypes.HRESULT, ctypes.c_void_p)(vtbl[7])
        fmtid = _guid(_FMTID)
        # the id last: setting it is what makes the taskbar read the relaunch properties
        for pid, value in ((_RELAUNCH_COMMAND, command), (_RELAUNCH_ICON, icon), (_RELAUNCH_NAME, NAME),
                           (_ID, APP_ID)):
            set_value(store, ctypes.byref(_PROPERTYKEY(fmtid, pid)), ctypes.byref(_PROPVARIANT(vt=_VT_LPWSTR,
                                                                                              pwszVal=value)))
        commit(store)
        return True
    except (OSError, AttributeError, ValueError) as e:
        print('taskbar: %r' % e)
        return False
    finally:
        if store:
            ctypes.WINFUNCTYPE(ctypes.c_ulong, ctypes.c_void_p)(
                ctypes.cast(store, ctypes.POINTER(ctypes.POINTER(ctypes.c_void_p))).contents[2])(store)
