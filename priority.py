"""Raise this process out of Windows' background power throttling.

Measured 3 Oct 2026 on the home PC: fourteen busy worker processes were
delivering 3.96 cores of real throughput on a sixteen-thread machine, each
getting about 28% of a core, with affinity unrestricted and the power plan at
100% maximum processor state. Windows was applying EcoQoS power throttling,
which it does to processes it judges to be background work - which long detached
compute jobs look exactly like. Raising the priority class lifted throughput to
6.38 cores at ABOVE_NORMAL and 14.09 at HIGH: a factor of 3.6.

Priority alone turned out not to be enough. Runs relaunched that evening by the
S4U scheduled task (session 0) drew only ~36% of the machine with fourteen
workers, despite each worker setting HIGH. So this now also opts out of EcoQoS
explicitly, through SetProcessInformation(ProcessPowerThrottling) with the
execution-speed bit controlled and cleared - the documented switch, rather than
relying on Windows to infer it from priority.

So every wall-clock timing recorded before 3 Oct 2026 was taken at roughly a
third of this machine's capacity. CPU-time figures (core-hours) are unaffected,
since they measure processor time actually consumed.

HIGH is safe for batch compute - it sits below REALTIME and below the system's
own threads - but it will make the desktop feel slow if someone is using the
machine. Call with high=False for ABOVE_NORMAL.
"""
import ctypes
import os

ABOVE_NORMAL_PRIORITY_CLASS = 0x00008000
HIGH_PRIORITY_CLASS = 0x00000080
PROCESS_POWER_THROTTLING = 4                       # PROCESS_INFORMATION_CLASS
PROCESS_POWER_THROTTLING_CURRENT_VERSION = 1
PROCESS_POWER_THROTTLING_EXECUTION_SPEED = 0x1


class _ThrottleState(ctypes.Structure):
    _fields_ = [("Version", ctypes.c_ulong),
                ("ControlMask", ctypes.c_ulong),
                ("StateMask", ctypes.c_ulong)]


def raise_priority(high=True):
    """Opt out of EcoQoS and raise the priority class.

    Returns {"priority": bool, "ecoqos_off": bool} so callers can record whether
    each step actually took. No-op off Windows; never fatal.
    """
    out = {"priority": False, "ecoqos_off": False}
    if os.name != "nt":
        return out
    try:
        from ctypes import wintypes
        k32 = ctypes.windll.kernel32
        # Without explicit types ctypes treats the HANDLE as a 32-bit int, so the
        # pseudo-handle -1 from GetCurrentProcess reaches the next call as
        # 0x00000000FFFFFFFF, an invalid handle, and both calls fail silently.
        # That bug was live from the first version: every worker's "raise
        # priority" was a no-op, and only the PowerShell fix applied from outside
        # ever took effect - which is why task-launched runs stayed throttled.
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        k32.SetPriorityClass.restype = wintypes.BOOL
        k32.SetProcessInformation.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                              ctypes.c_void_p, wintypes.DWORD]
        k32.SetProcessInformation.restype = wintypes.BOOL
        me = k32.GetCurrentProcess()
        st = _ThrottleState(PROCESS_POWER_THROTTLING_CURRENT_VERSION,
                            PROCESS_POWER_THROTTLING_EXECUTION_SPEED, 0)
        out["ecoqos_off"] = bool(k32.SetProcessInformation(
            me, PROCESS_POWER_THROTTLING, ctypes.byref(st), ctypes.sizeof(st)))
        cls = HIGH_PRIORITY_CLASS if high else ABOVE_NORMAL_PRIORITY_CLASS
        out["priority"] = bool(k32.SetPriorityClass(me, cls))
    except Exception:
        pass
    return out
