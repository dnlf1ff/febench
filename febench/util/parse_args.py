import argparse

def parse_base_args(argv: list[str]| None=None): #TODO namespace
    # -----------BASE------------#
    parser = argparse.ArgumentParser(description= "cli tool")

    parser.add_argument('--calc', type=str, default='omni',help='potential name')
    parser.add_argument('--calc_type', type=str, default='sevennet',help='sevennet, esen, uma, etc')
    
    parser.add_argument('--modal', type=str, default='null',help='data modal; mpa, omat24,null')

    parser.add_argument('--potential_path', type=str, default='/data2/jinvk/MLP',help='path for potential file')

    parser.add_argument('--potential_ext', type=str, default='pth',help='extention for potential file')
    
    parser.add_argument('--config', type=str, default='./config.yaml', help='config yaml file dir')
    
    # nvm will be overrided by config.yaml
    # parser.add_argument('--carbon_config', type=str, default='./carbon_config.yaml', help='carbon config yaml file dir') 

    return parser.parse_args(argv)

