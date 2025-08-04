import numpy as np
import sys

import yaml

def get_surface_area(atoms):
    return np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))

def dict_representer(dumper, data=None):
    return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data, flow_style=False)

def list_representer(dumper, data=None):
    return dumper.represent_sequence(yaml.resolver.BaseResolver.DEFAULT_SEQUENCE_TAG, data, flow_style=True)

class WDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()

def dumpYAML(data, filename, indent=4, sort_keys=False, explicit_start=True, explicit_end=True, default_flow_style=False,):
    yaml.add_representer(dict, dict_representer, Dumper=WDumper)
    yaml.add_representer(list, list_representer, Dumper=WDumper)
    with open(filename, 'w') as fp:
        yaml.dump(data, fp, Dumper=WDumper, sort_keys=sort_keys, explicit_start=explicit_start, explicit_end=explicit_end, default_flow_style=default_flow_style, indent=indent)



