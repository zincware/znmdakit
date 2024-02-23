import h5py
import numpy as np
import logging

log = logging.getLogger(__name__)

def convert_legacy_znh5md(filename, *, force=False):
    """Fix the issues with the legacy h5md files to comply with the current h5md standard."""
    if not force:
        log.critical("This function will modify the file in place. It will remove all observables and make some assumptions. Set 'force=True' to proceed.")
        return

    with h5py.File(filename, mode="a") as f:
        # if the first group in "particles" is not named "trajectory", rename it
        particle_group = list(f["particles"].keys())[0]
        log.critical(f"Particle groups: {particle_group}")
        # if particle_groups[0] != "trajectory":
        #     log.critical(f"Renaming {particle_groups[0]} to trajectory")
        #     f["particles"].move(particle_groups[0], "trajectory")

        # fix PBC
        # del f[f"particles/{particle_groups}/box/edges/value"]
        log.critical(f"Fixing PBC")
        f[f"particles/{particle_group}/box"].attrs['boundary'] = ["periodic", "periodic", "periodic"]
        f[f"particles/{particle_group}/box"].attrs['dimension'] = 3
        del f[f"particles/{particle_group}/box/edges/time"]
        del f[f"particles/{particle_group}/box/edges/step"]
        del f[f"particles/{particle_group}/box/boundary"]
        del f[f"particles/{particle_group}/box/dimension"]
        length = len(f[f"particles/{particle_group}/box/edges/value"])
        f[f"particles/{particle_group}/box/edges/time"] = np.linspace(0, 1, length)
        f[f"particles/{particle_group}/box/edges/step"] = np.arange(length)
        
        available_groups = list(f["particles"][particle_group].keys())
        for group in available_groups:
            # check if the group has a time / step attribute
            if "time" in f[f"particles/{particle_group}/{group}"]:
                log.critical(f"Replacing integer time with array for {group} in {particle_group}")
                length = len(f[f"particles/{particle_group}/{group}/value"])
                del f[f"particles/{particle_group}/{group}/time"]
                del f[f"particles/{particle_group}/{group}/step"]
                f[f"particles/{particle_group}/{group}/time"] = np.linspace(0, 1, length)
                f[f"particles/{particle_group}/{group}/step"] = np.arange(length)

        # now fix /observables
        log.critical(f"Removing observables")
        del f[f"observables/{particle_group}"]
        # available_groups = list(f["observables"][particle_group].keys())
        # for group in available_groups:
            # # check if the group has a time / step attribute
            # if "time" in f[f"observables/{particle_group}/{group}"]:
            #     log.critical(f"Replacing integer time with array for {group} in observables/{particle_group}")
            #     length = len(f[f"observables/{particle_group}/{group}/value"])
            #     del f[f"observables/{particle_group}/{group}/time"]
            #     del f[f"observables/{particle_group}/{group}/step"]
            #     f[f"observables/{particle_group}/{group}/time"] = np.linspace(0, 1, length)
            #     f[f"observables/{particle_group}/{group}/step"] = np.arange(length)
